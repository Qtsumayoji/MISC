using LinearAlgebra
using PyPlot
using StatsBase

# fermion signの計算に使用
bitflag = Array{Int64}(undef, 64)
for i in 1:64
    bitflag[i] = ~((1 << i) - 1)
    #println(bits(bitfrag[i]))
end

function make_n_basis(Ns::Int64, Ne::Int64, nup::Int64, ndown::Int64)
    nstate = 4^Ns
    dim = binomial(Ns, nup)*binomial(Ns, ndown)
    basis = zeros(Int64, dim)
    cnt = 1
    
    @inbounds @simd for i in 0:nstate-1
        if count_ones(i) == Ne
            if count_ones(i & bitflag[Ns]) == ndown
                if count_ones(i & ~bitflag[Ns]) == nup
                    basis[cnt] = i
                    cnt += 1
                    #println(string(i,base=2))
                end
            end
        end
    end

    return basis
end

function main()
    min = 2
    max = 28
    x = zeros(max - min + 1)
    y = similar(x)

    Ne = 2
    nup = 1
    ndown = 1

    make_n_basis(10, Ne, nup, ndown)
    
    for i in min:max
        Ns = i
        eptime = @elapsed make_n_basis(Ns, Ne, nup, ndown)
        println(i," ",eptime)
        x[i - min + 1] = i
        y[i - min + 1] = eptime
    end

    PyPlot.plot(x, y)
    PyPlot.show()

end

main()