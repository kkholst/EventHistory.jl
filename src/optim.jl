struct D2Function
    value::Float64
    gradient::Array{Float64,1}
    hessian::Array{Float64,2}
end


###{{{ Newton-Raphson

function NR(f, x; iter=200,tol=1e-12,verbose=false,method=0,gamma=1,trace=0)
    val = f(x)::D2Function
    # S = val.grad["grad"]';  H = val["hess"]; loglik=val["val"]
    S = val.gradient::Vector{Float64};  H = val.hessian::Matrix{Float64}; loglik=val.value::Float64
    p = size(x,1)
    counter = 1
    for j=1:iter
        if verbose print("S=",S') end
        s = s0 = mean(S.*S)
        if isnan(s) error("Lack of convergence (bad starting values?)") end
        if s<tol break end
        if method==1
            W = pinv(-H + gamma*(S'S)[1]*eye(p))
        elseif method==2
            W = pinv(-H + gamma*S*S')
        else
            W = inv(-H)
        end
        if trace>0
            print(S')
        end
        delta = W*S
        x0 = x; loglik0 = loglik
        lambda=1;
        val = f(x0)
        S = val.gradient::Vector{Float64};  H = val.hessian::Matrix{Float64};
        print("x=",x0, " loglik=", val.value, " grad=", S, "\n")
        # s = mean(S.*S)
        # print(s,"\n")
        # x = x0 + lambda*delta
        while (s>=s0 || loglik<loglik0)
            x = x0+lambda*delta
            val = f(x)
            counter += 1
            ##S = val["grad"]';  H = val["hess"]; loglik=val["val"]
            S = val.gradient::Vector{Float64};  H = val.hessian::Matrix{Float64}; loglik=val.value::Float64
            s = mean(S.*S)
            if isnan(s) s=s0 end
            lambda = lambda/2.0
            if (lambda<1e-12)
                error("Lack of convergence (bad starting values?)")
            end
        end
        print("counter=", counter, "\n")
    end
    return x,val,counter
end

###}}} Newton-Raphson

