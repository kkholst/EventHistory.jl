###{{{ EventHistoryModel

type EventHistoryModel
    model::String
    call::Expr
    eventtype::DataType
    coef::Vector
    coefmat::DataFrame
    IC::Matrix
    invhess::Matrix
    vcov::Matrix
    opt::Vector
    grad::Vector
    X::Matrix
    eventtime::Matrix    
end

function show(io::IO, obj::EventHistoryModel)
    print(io,"\nModel: ", obj.model,",", obj.eventtype, " ", obj.call)
    n = size(obj.eventtime,1)
    events::Int = sum(obj.eventtime[:,2])
    print(io,"\nn=",n,", events=",events,"\n\n")
    print(io,obj.coefmat)
end

###}}} EventHistoryModel

###{{{ phreg

function phreg(X, t, status, entry, id=[]; beta=[], opt...)
    pre = coxprep(t,status,X,entry,id);
    p = size(X,2)
    if size(beta,1)==0 beta=vec(zeros(p,1)) end
    pl(beta,indiv=false) =
        coxPL(beta,pre["X"],pre["XX"],pre["sign"],pre["jumps"],indiv)
    val = NR(pl,beta,opt...)
    beta = vec(val[1]); grad = val[2]["grad"]'; steps = val[3]; 
    U = pl(beta,true)
    I = -pinv(U["hess"])
    iid = U["grad"]I
    V =  iid'iid
    coefmat = [beta sqrt(diag(V)) sqrt(diag(I))]
    coefmat = [coefmat 2*cdf(Normal(0,1),-abs(coefmat[:,1]./coefmat[:,2]))]
    coln = ["Estimate","S.E.","dU^-1/2","P-value"]
    EventHistoryModel("Cox",:(~1),EventHistory.Surv,
                      vec(beta),DataFrame(coefmat,coln),
                      iid,I,V,
                      [steps],vec(grad),X,[t status])
end

function phreg(formula, data; id=[], opt...)    
    M = ModelFrame(formula,data);
    S = convert(Vector,M.df[:,1])
    time = EventHistory.Time(S)
    status = EventHistory.Status(S)
    entry = try EventHistory.Entry(S) catch [] end
    X = ModelMatrix(M)
    X = X.m[:,[2:size(X.m,2)]];
    res = EventHistory.phreg(X, time, status, entry; opt...)
    return res
    res.call = formula
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
    XX = Array(Float64, n, p*p); ## Calculate XX' at each time-point
    ##  XX = Array(Float64, p, p, n);
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
    ["X"=>X, "XX"=>XX, "jumps"=>jumps, "sign"=>sgn,
     "ord"=>ord, "time"=>exit, "id"=>[]]
end
###}}} coxprep

###{{{ revcumsum
function revcumsum(A,dim=1)
    D = size(A)
    n = D[dim]
    res = similar(A)
    prev = 0;
    for i=1:n
        idx = 
        prev += A[n-i+1,:];
        res[n-i+1,:] = prev;
    end
    res
end
###}}} revcumsum

###{{{ coxPL
function coxPL(beta, X, XX, sgn, jumps, indiv=false)
    n = size(X,1);
    p = size(X,2);
    Xb = X*beta;
    eXb = map(exp,Xb);
    if size(sgn,1)==n ## Truncation
        eXb = sgn.*eXb;
    end
    S0 = revcumsum(eXb);
    E = similar(X);
    D2 = similar(XX)
    ## S1 = Array(Float64,n,p)
    ## S2 = similar(XX)
    for j=1:p
        ## S1[:,j] = revcumsum(X[:,j].*eXb);
        E[:,j] = revcumsum(X[:,j].*eXb)./S0; ## S1/S0(s)
    end;
    for j=1:size(XX,2)
        ## S2.col(j) = revcumsum(XX.col(j).*eXb);
        D2[:,j] = revcumsum(XX[:,j].* eXb)./S0 ## int S2/S0(s)
    end;
    D2 = D2[jumps,:]
    E = E[jumps,:]
    S0 = S0[jumps]
    grad = (X[jumps,:]-E); ## Score
    val = Xb[jumps]-log(S0); ## Partial log-likelihood
    hess = -(reshape(sum(D2,1),p,p)-E'E);
    if !indiv
        val = sum(val); grad = sum(grad,1)
    end
    ["val"=>val, "grad"=>grad, "hess"=>hess]    
end

###}}} coxPL

###{{{ predict


###}}}
