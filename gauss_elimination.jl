function gauss(A,b)
    n = length(b)
    x = zeros(n)
    for k=1:n-1
        for i=k+1:n
            for j=k+1:n
                A[i,j] -= A[i,k]/A[k,k]*A[k,j]
            end
            b[i] -= A[i,k]/A[k,k]*b[k]
        end
    end
    x[n] = b[n]/A[n,n]
    for i=1:n-1
        for j=0:i-1
            b[n-i] -= A[n-i,n-j]*x[n-j]
        end
        x[n-i] = b[n-i]/A[n-i,n-i]
    end
    return x 
end

function gauss_rational(A::Matrix{T},b::Vector{T}) where T<:Rational
    n = length(b)
    x = zeros(Rational,n)
    for k=1:n-1
        for i=k+1:n
            for j=k+1:n
                A[i,j] -= A[i,k]/A[k,k]*A[k,j]
            end
            b[i] -= A[i,k]/A[k,k]*b[k]
        end
    end
    x[n] = b[n]/A[n,n]
    for i=1:n-1
        for j=0:i-1
            b[n-i] -= A[n-i,n-j]*x[n-j]
        end
        x[n-i] = b[n-i]/A[n-i,n-i]
    end
    return x 
end