###{{{ phreg

function phreg(X::Matrix, t::Vector, status::Vector, entry::Vector, id=[]; beta=[], opt...)
    pre = EventHistory.coxprep(t,status,X,entry,id)
    p = size(X,2)
    if size(beta,1)==0 beta=vec(zeros(p,1)) end
    pl(beta,indiv=false) =
        EventHistory.coxPL(beta,pre["X"],pre["XX"],pre["sign"],pre["jumps"],indiv)
    val = EventHistory.NR(pl,beta; opt...)
    beta = vec(val[1]); grad = val[2].gradient; steps = val[3];
    U = pl(beta,true)
    I = -pinv(U[3])
    iid = U[2]*I
    V =  iid'iid
    coefmat = [beta sqrt(diag(V)) sqrt(diag(I))]
    coefmat = [coefmat 2*Distributions.cdf(Distributions.Normal(0,1),-abs(coefmat[:,1]./coefmat[:,2]))]
    ## coln = [:Estimate,:SE,:naiveSE,symbol("pvalue")]
    ## cc = convert(DataFrame,coefmat); names!(cc,coln);
    cc = CoefTable(coefmat,["Estimate","S.E","naive S.E.","P-value"],
                   repeat([""],inner=[size(coefmat,1)]))
    chaz = [pre["jumps"] pre["time"][pre["jumps"]] cumsum(1./(U[4][pre["jumps"]]))]
    EventHistoryModel("Cox",Formula(:.,:.),EventHistory.Surv,
                      vec(beta),cc,
                      iid,I,V,
                      [steps],vec(grad),X,[t status],
                      chaz
                      )
end


function phreg(formula::Formula, data::DataFrame; id=[], opt...)
    ## fm = DataFrames.Formula(formula.args[2],formula.args[3])
    M = ModelFrame(formula,data)
    S = M.df[:,1] ## convert(Vector,M.df[:,1])
    time = EventHistory.Time(S)
    status = EventHistory.Status(S)
    entry = try EventHistory.Entry(S) catch [] end
    X = ModelMatrix(M)
    X = X.m[:,collect(2:size(X.m,2))]
    res = EventHistory.phreg(X, time, status, entry; opt...)
    cnames = setdiff(coefnames(M),["(Intercept)"])
    res.coefmat.rownms=cnames
    res.formula = formula
    res.eventtype = typeof(S[1])
    res
end

###}}} phreg

###{{{ coxprep

function coxprep(exit,status,X=[],entry=[],id=[])
    Truncation = size(entry,1)==size(exit,1)
    n = size(exit,1)
    p = size(X,2)
    truncation = size(entry,1)==n
    XX = Array(Float64, n, p*p) ## Calculate XX' at each time-point
    ##  XX = Array(Float64, p, p, n)
    if size(X,1)>0
        for i=1:n
            Xi = X[i,:]
            XX[i,:] = vec(Xi'Xi)
            ## XX[:,:,i] = Xi'Xi
        end
    end
    sgn=[]
    if truncation
        sgn = Array(Int,2*n); fill!(sgn,1)
        status = [status;status]
        for i=1:n
            sgn[i] = -1
            status[i] = 0
        end
        exit = [entry;exit]
        X = [X;X]
        XX = [XX;XX]
    end
    ord0 = sortperm(status*1,rev=true)
    ord = sortperm(exit[ord0])
    ord = ord0[ord]
    if (truncation)
        sgn = sgn[ord]
    end
    if size(X,1)>0
        X = X[ord,:]
        XX = XX[ord,:]
    end
    exit = exit[ord]
    status = status[ord]
    jumps = find(status)
    Dict("X"=>X, "XX"=>XX, "jumps"=>jumps, "sign"=>sgn,
     "ord"=>ord, "time"=>exit, "id"=>[])
end

###}}} coxprep

###{{{ revcumsum
function revcumsum(A,dim=1)
    D = size(A)
    n = D[dim]
    res = similar(A)
    prev = zeros(Float64,1,size(A,2))
    for i=1:n
        idx =
        prev += A[n-i+1,:]
        res[n-i+1,:] = prev
    end
    res
end
###}}} revcumsum

###{{{ coxPL
function coxPL(beta::Vector, X::Matrix, XX::Matrix, sgn::Vector, jumps::Vector, indiv=false)
    n = size(X,1)
    p = size(X,2)
    Xb = X*beta
    eXb = map(exp,Xb)
    if size(sgn,1)==n ## Truncation
        eXb = sgn.*eXb
    end
    S0 = revcumsum(eXb)
    E = similar(X)
    D2 = similar(XX)
    S1 = Array(Float64,n,p)
    S2 = similar(XX)
    for j=1:p
        if (indiv) S1[:,j] = revcumsum(X[:,j].*eXb); end
        E[:,j] = revcumsum(X[:,j].*eXb)./S0 ## S1/S0(s)
    end
    for j=1:size(XX,2)
        if (indiv) S2[:,j] = revcumsum(XX[:,j].*eXb); end
        D2[:,j] = revcumsum(XX[:,j].* eXb)./S0 ## int S2/S0(s)
    end
    D2 = D2[jumps,:]
    E = E[jumps,:]
    grad = (X[jumps,:]-E) ## Score
    val = Xb[jumps]-log(S0[jumps]) ## Partial log-likelihood
    hess = -(reshape(sum(D2,1),p,p)-E'E)
    if indiv
        return(val,grad,hess,vec(S0),S1,S2,E)
    end
    val = sum(val); grad = sum(grad,1)
    EventHistory.D2Function(val[1],vec(grad),hess)
end

###}}} coxPL

###{{{ predict
# compute the values of a step function,
# i.e. how many of the jumps are smaller or
# equal to the eval points
function sindex(jump::Vector, eval::Vector)
    N = size(eval,1)
    n = size(jump,1)
    index = Array(Int64,N)
    j = 1::Int64
    for t=1:N
        while (j<n) & (jump[j]<eval[t])
            j += 1
        end
        index[t] = j
    end
    return(index)
end

## Returns [time surv-prob]
function predict(mm::EventHistory.EventHistoryModel; X=[]'::Matrix,time=mm.eventtime[:,1]::Vector,surv=true,order=false)
    ord = sortperm(time)
    idx = sindex(mm.chaz[:,2],time[ord,1])
    L0 = mm.chaz[:,3][idx]
    L0 = L0[sortperm(ord)]
    arraytype = size(X,1)==1
    if size(X,2)==0
        X = mm.X
        arraytype = true
    end
    if arraytype
        H = exp(X*coef(mm))
        res = L0.*H
    else
        res = Array(Union{Float64,Int64,Bool},size(time,1),size(X,1))
        for i=1:size(X,1)
            res[:,i] = L0.*exp(X[i,:]*coef(mm))
        end
    end
    if surv res = exp(-res) end
    res = [time res]
    if order
        return(res[sortperm(time),:])
    end
    return(res)
end

###}}} predict

function coef(mm::EventHistoryModel)
    mm.coef
end

function vcov(mm::EventHistoryModel)
    mm.vcov
end

function coeftable(mm::EventHistoryModel)
    mm.coefmat

end
