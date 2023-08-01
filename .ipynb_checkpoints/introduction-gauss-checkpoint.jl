using LinearAlgebra, SparseArrays

# ガウスの消去法 (Gauss elimination) のアルゴリズム
function GE(A,b)
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

# ガウスの消去法を有理数のまま扱うアルゴリズム (普通はしない)
function GE_rational(A::Matrix{T},b::Vector{T}) where T<:Rational
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


# 規模の大きな問題を作る. 入力nは未知数の個数.
# 方程式は, n=4なら
# 
#     2x_1  -x_2             = 2/n^2
#     -x_1 +2x_2  -x_3       = 2/n^2
#           -x_2 +2x_3  -x_4 = 2/n^2
#                 -x_3 +2x_4 = 2/n^2
# 
# nが大きいと, ほとんどの係数が0になっている.
# 解は, x_i = i/(n+1)*(1-i/(n+1)) (i=1,2,...,n)
function make_example(n)
    u0 = 2*ones(n)
    u1 = -ones(n-1)
    A = spdiagm(0=>u0,1=>u1,-1=>u1)
    b = 2*ones(n)/(n+1)^2
    return A,b 
end


# 共役勾配法 (きょうやくこうばいほう, conjugate gradient method) のアルゴリズム
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