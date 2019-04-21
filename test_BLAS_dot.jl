using LinearAlgebra

BLAS.set_num_threads(4)

function float(n, nloop)
    println("allocation")
    @time vec1 = rand(n)
    @time vec2 = rand(n)

    println("dot product")
    @time vec1'*vec2

    println("dot product*",nloop," loop")
    @time for i in 1:nloop
        vec1'*vec2
    end
end

function complex(n, nloop)
    println("allocation")
    @time vec1 = rand(n) + im*rand(n)
    @time vec2 = rand(n) + im*rand(n)

    println("dot product")
    @time vec1'*vec2

    println("dot product*",nloop," loop")
    @time for i in 1:nloop
        #vec1'*vec2
        vcvc(vec1, vec2)
    end
end

function vcvc(x1, x2)
    #return real(x1)'*real(x1) - imag(x2)'*imag(x2) + im*(real(x1)'*imag(x2) + imag(x1)'*real(x2)) 
    return x1'*x2
    
end

function main()
    nloop = 50

    nGB = 1
    #n = Int64(nGB*2^27)

    nMB = 500
    n = Int64(nMB*2^17)

    nKB = 50
    #n = Int64(nKB*2^7)

    println("dim=", n)
    println("")

    println("Float")
    float(n, nloop)
    println("")

    println("Complex")
    complex(n, nloop)
end

main()