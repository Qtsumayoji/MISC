using LinearAlgebra

using PyCall

@pyimport pylab as plt
@pyimport seaborn as sns

function JAB(J, kx, ky)
    return J*(cos(kx) + cos(ky))
end

function JAA(J1, J2, kx, ky)
    return 2.0*J1*cos(kx)*cos(ky) + J2*(cos(2.0*kx) + cos(2.0*ky))
end

function main()
    a1 = [1.0; 0.0; 0.0]
    a2 = [0.0; 1.0; 0.0]
    a3 = [0.0; 0.0; 0.0]
    
    v = a1'*cross(a2, a3)
    g1 = 2.0*π*cross(a2, a3)/v
    g2 = 2.0*π*cross(a3, a1)/v
    g3 = 2.0*π*cross(a2, a3)/v

    J1 = 1.0
    J2 = -0.5*J1
    J3 = 0.3*J1

    Nx = 200
    Ny = Nx
    Qx = range(-pi, stop=pi, length=Nx)
    Qy = range(-pi, stop=pi, length=Ny)
    ϵk = zeros(Nx, Ny)
    ϵkx = zeros(Nx)
    ΔJ = JAB(J1, 0.0, 0.0) - JAA(J2, J3, 0.0, 0.0)
    δ = 1e-10

    @time for i in 1:Nx
        kx = Qx[i]
        ky = 0.0
        ϵkx[i] = sqrt((ΔJ + JAA(J2, J3, kx, ky))^2.0 - JAB(J1, kx, ky)^2.0 + δ) - sqrt(δ)
        for j in 1:Ny
            ky = Qy[j]
            ϵk[i, j] = sqrt((ΔJ + JAA(J2, J3, kx, ky))^2.0 - JAB(J1, kx, ky)^2.0 + δ) - sqrt(δ)
        end
    end


    plt.figure(figsize=(18,8))
    plt.subplot(121)
    plt.pcolormesh(Qx, Qy, ϵk, cmap="Blues")
    plt.colorbar()

    plt.subplot(122)
    plt.plot(Qx, ϵkx)

    plt.show()
end
main()