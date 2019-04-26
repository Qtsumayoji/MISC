using LinearAlgebra
using PyCall

@pyimport pylab as plt
@pyimport seaborn as sns

Nx = 10
Ny = 10
V = Nx*Ny

#a1 = [sqrt(3.0)/2.0; 1.0/2.0; 0.0]
#a2 = [sqrt(3.0)/2.0; -1.0/2.0; 0.0]
#a3 = [0.0; 0.0; 1.0]
a1 = [1.0; 0.0; 0.0]
a2 = [0.0; 1.0; 0.0]
a3 = [0.0; 0.0; 0.0]

v = a1'*cross(a2, a3)
g1 = 2.0*π*cross(a2, a3)/v
g2 = 2.0*π*cross(a3, a1)/v
g3 = 2.0*π*cross(a2, a3)/v

t = 1.0
U = 2.0

function calc_ϵk_sq(kx::Float64, ky::Float64, t::Float64)
    return -2.0*t*(cos(kx) + cos(ky))
end

function calc_fermi_dist(β, ϵ, μ)
    return 1.0/(exp(β*(ϵ - μ)) + 1.0)
end

function calc_Σ_HF(k, β, μ, Nx, Ny)
    Σ_HF = 0.0
    for i in 1:Nx
        k1 = (i - 1)/Nx*g1
        for j in 1:Ny
            k2 = (j - 1)/Ny*g2
            k = k1 + k2
            kx = k[1]
            ky = k[2]
            ϵk = calc_ϵk_sq(kx, ky, t)
            Σ_HF += U*calc_fermi_dist(β, ϵk, μ)
        end
    end
    
    return  Σ_HF/V
end

function main()
    
end
main()