using LinearAlgebra

function gram_schmidt(U::Array)
    siz = size(U)
    n = siz[2]
    for i in 1:n
        for j in 1:i - 1
            a = U[:, j]'*U[:, i]
            U[:, i] += -a*U[:, j]
        end
        U[:, i] /= norm(U[:, i]) 
    end
end

function print_mat(A)
    dim_r = size(A)[1]
    dim_c = size(A)[2]

    for i in 1:dim_r
        for j in 1:dim_c
            if abs(A[i, j]) > 1e-10
                print(A[i,j]," ")
            else
                print("0.0 ")
            end
        end
        println()
    end
end

function main()
    A = [
        1.0 2.0 3.0;
        2.0 -2.5 3.4;
        3.0 3.4 5.0;
        ]

    B = [
        0.0 0.0 0.0;
        0.0 -2.5 3.4;
        0.0 3.4 5.0;
        ]

    E, U = eigen(B)

    tmp = U[:,2]
    U[:,2] = U[:,1]
    U[:,1] = tmp

    #println(E)
    #print_mat(U)
    print_mat(U'*A*U)
    #print_mat(U'*B*U)
    vec = rand(3)
    println(U*vec)
end
#main()

function  test()
    α = [3.0;-2.0;5.1]
    β = [0.5;3.0;4.0]
    He = zeros(3,3) + Diagonal(α)

    ω = 1.0
    G0 = 0.0
    for i in 1:3
        G0 += β[i]^2.0/(ω - α[i])
    end
    println(G0)

    U = rand(3,3)
    gram_schmidt(U)
    U[:,1] = β
    println(U*β)
    He = U'*He*U
    Ginv = inv(ω*I - He)
    println(Ginv[1,1])
end
test()