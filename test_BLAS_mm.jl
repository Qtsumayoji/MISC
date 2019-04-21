using LinearAlgebra
using SparseArrays

BLAS.set_num_threads(4)

function test_sparse(n, nloop)
    println("allocation")
    @time A = rand(n, n)
    @time B = rand(n, n)

    println("matrix matrix product")
    @time A*B

    println("matrix vector product*",nloop," loop")
    @time for i in 1:nloop
        A*B
    end
end

function main()
    nloop = 1000

    nGB = 1
    #n = Int64(nGB*2^27)

    nMB = 1
    n = Int64(nMB*2^17)

    nKB = 200
    n = Int64(nKB*2^7)

    n = 5000

    println("dim=", n)
    println("")

    println("sparse")
    test_sparse(n, nloop)
end

main()