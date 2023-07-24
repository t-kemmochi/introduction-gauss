function CG(A,b;x0::Vector{Float64}=zeros(length(b)),tol::Float64=1.0e-12,imax::Int64=length(b))
    r = b - A*x0 
    p = copy(r)
    x = copy(x0)
    i = 0
    res = 1.0
    while i<imax && res>tol 
        r0 = copy(r)
        q = A*p
        a = r'*p/(p'*q)
        x += a*p 
        r -= a*q 
        b = r'*r/(r0'*r0)
        p = r+b*p 
        res = sqrt(r'*r)/sqrt(b'*b)
        i += 1 
    end
    return x,i,res
end