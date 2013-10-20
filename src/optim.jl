###{{{ Newton-Raphson
function NR(f, x; iter=200,tol=1e-12,verbose=false,method=0,gamma=1)    
    val = f(x)
    S = val["grad"]';  H = val["hess"]; loglik=val["val"]    
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
            W = pinv(-H)
        end
        delta = W*S
        x0 = x; loglik0 = loglik
        lambda=1;
        while s>=s0 || loglik<loglik0
            x = x0+lambda*delta
            val = f(x)
            counter += 1
            S = val["grad"]';  H = val["hess"]; loglik=val["val"]    
            s = mean(S.*S)
            if isnan(s) s=s0 end
            lambda = lambda/2.0
            if (lambda<1e-12)
                error("Lack of convergence (bad starting values?)")
            end
        end
    end
    return x,val,counter
end
###}}} Newton-Raphson
    
