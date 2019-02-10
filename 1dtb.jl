using LinearAlgebra
using SparseArrays
using Distributions
using PyCall

@pyimport pylab as plt
@pyimport seaborn as sns

function main()
    Nx = 10
    a = 1.0
    x = [a*i for i in 1:Nx]
    g1 = 2*pi/a
    kx = [n/Nx*g1 for n in 0:Nx-1]

    nlink = 2
    link_list = [[] for i in 1:Nx]
    for i in 1:Nx
        if i == 1
            push!(link_list[i], i + 1)
            push!(link_list[i], Nx)
            continue
        end
        if i == Nx
            push!(link_list[i], 1)
            push!(link_list[i], i - 1)
            continue
        end
        push!(link_list[i], i + 1)
        push!(link_list[i], i - 1)
    end

    Hk = zeros(Complex, Nx, Nx)

    t = 1.0
    for i in 1:Nx
        k = kx[i]
        link1 = link_list[1][1]
        link2 = link_list[1][2]
        x0 = x[1]
        x1 = x[link1]
        x2 = x[link2]
        Hk[i, i] += -t*(exp(-im*k*(x1 - x0)) + exp(-im*k*(x2 - x0)))
    end
    E, U = eigen(Hk)
    Eband = E
    plt.plot(kx, Eband,"x")
    plt.show()

end

main()