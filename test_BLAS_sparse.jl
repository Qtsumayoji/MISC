using LinearAlgebra
using SparseArrays
using MKLSparse

BLAS.set_num_threads(4)

function Ax(A,x)
    #return mul!(x, A, x)
    return A*x
end

function set_sparse_matrix(n)
    println("allocation")
    v = Float64[]
    r = Int64[]
    c = Int64[]

    println("set sparse matrix")
    @time for i in 1:n
        for j in 1:i-1
            x = rand()
            if x < 0.002
                push!(v, x)
                push!(r, i)
                push!(c, j)
            end 
        end
        x = rand()
        push!(v, x)
        push!(r, i)
        push!(c, i)
    end
    @time A = sparse(r, c, v)
    println("nonzeros = ",nnz(A)/n/n*100," %")

    return A
end

function test_sparse(n, nloop)

    A = set_sparse_matrix(n)
    B = set_sparse_matrix(n)

    @time vec = rand(n)

    println("matrix vector product")
    @time Ax(A, vec)

    println("matrix vector product*",nloop," loop")
    @time for i in 1:nloop
        A*vec
        A*vec
        A*vec
        A*vec
        A*vec
        A*vec
        A*vec
        A*vec
        #Ax(A, vec)
        #A*B
    end
end

function main()
    nloop = 100*10*10
    nGB = 1
    #n = Int64(nGB*2^27)

    nMB = 1
    n = Int64(nMB*2^17)

    nKB = 200
    n = Int64(nKB*2^7)

    n=4900

    println("dim=", n)
    println("")

    println("sparse")
    test_sparse(n, nloop)
end

main()