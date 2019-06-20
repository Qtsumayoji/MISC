using LinearAlgebra
#using PyPlot
using Plots
pyplot()

t = 1.0

a1 = [-0.5; -0.5*sqrt(3.0); 0.0]
a2 = [1.0; 0.0; 0.0]
a3 = [-0.5; 0.5*sqrt(3.0); 0.0]

u1 = [2.0; 0.0; 0.0]
u2 = [1.0; sqrt(3.0); 0.0]
u3 = [0.0; 0.0; 1.0]

V = dot(u1, cross(u2, u3))

g1 = 2.0*pi/V*cross(u2, u3)
g2 = 2.0*pi/V*cross(u3, u1)
g3 = 2.0*pi/V*cross(u1, u2)

function E_12(k)
    c1 = cos(k'*a1)
    c2 = cos(k'*a2)
    c3 = cos(k'*a3)
    tmp = sqrt(4.0*(c1^2.0 + c2^2.0 + c3^2.0) - 3.0 + 1e-13)
    return -t*(1.0 + tmp), -t*(1.0 - tmp)
end

function E_3(k)
    return 2.0*t
end

function plot3d()
    N = 300
    E1 = zeros(N, N)
    E2 = zeros(N, N)
    E3 = zeros(N, N)
    kx = range(-pi, stop=pi, length=N)
    ky = range(-pi, stop=pi, length=N)

    @time for j in 1:N
        for i in 1:N
            #k = g1*(i - 1)/N + g2*(j - 1)/N - gb
            k = [kx[i]; ky[j]; 0.0]
            E1[i,j], E2[i,j] = E_12(k)
            E3[i, j] = E_3(k)
        end
    end
    @time plot(kx, ky, [E1', E2', E3'], st=:surface, dpi=300)
    savefig("kagome")
end

function plot_FS(ϵF)
    Δ = 0.01*0.5
    N = 1000
    E = Float64[]
    kx = range(-pi, stop=pi, length=N)
    ky = range(-pi, stop=pi, length=N)
    x = Float64[]
    y = Float64[]

    @time for j in 1:N
        for i in 1:N
            #k = g1*(i - 1)/N + g2*(j - 1)/N - gb
            k = [kx[i]; ky[j]; 0.0]
            Ek1, Ek2 = E_12(k)
            if ϵF - Δ <= Ek1 <= ϵF + Δ
                push!(E, Ek1)
                push!(x, k[1])
                push!(y, k[2])
            end
        end
    end
    #@time plot(kx, ky, [E1'], st=:surface, dpi=300)


    PyPlot.figure(figsize=(12, 10))
    PyPlot.scatter(x, y)
    PyPlot.xlabel("qx", size=20)
    PyPlot.ylabel("qy", size=20)
    PyPlot.savefig("E_k_kagome")
    PyPlot.show()
end

function DOS()
    η = 0.01
    Nx = 400
    Ny = 400
    N = Nx*Ny
    NΩ = 1000
    Ω = range(-5.0, stop=3.0, length=NΩ)
    G = zeros(ComplexF64, NΩ)

    Ek1 = zeros(Nx, Ny)
    Ek2 = zeros(Nx, Ny)
    Ek3 = zeros(Nx, Ny)
    
    for m in 1:Nx
        k1 = m/Nx*g1 
        for n in 1:Ny
            k2 = n/Ny*g2
            k = k1 + k2
            Ek1[m, n], Ek2[m, n] =  E_12(k)
            Ek3[m, n] = E_3(k)
        end
    end

    @time for i in 1:NΩ
        ω = Ω[i]
        for m in 1:Nx
            for n in 1:Ny
                E1 = Ek1[m, n]
                E2 = Ek2[m, n]
                E3 = Ek3[m, n]
                
                G[i] += 1.0/(ω - E1 + im*η)
                G[i] += 1.0/(ω - E2 + im*η)
                G[i] += 1.0/(ω - E3 + im*η)
            end
        end
    end
    G = G/N
    DOS = -1.0/pi*imag(G)

    Δω = Ω[2] - Ω[1]
    ϵF = -1.0
    n = 0.0
    for i in 1:NΩ
        ω = Ω[i]
        n += Δω*DOS[i]

        if n >= 0.75
            println("n=3/4, ϵF=",ω)
            ϵF = ω
            break
        end
    end
    println(ϵF)
    
    PyPlot.figure(figsize=(12, 10))
    PyPlot.plot(Ω, DOS)
    PyPlot.ylim(0.0, 2.0)
    PyPlot.xlabel("E", size=20)
    PyPlot.ylabel("DOS", size=20)
    PyPlot.savefig("DOS_kagome.png")
    PyPlot.show()
end
#plot_FS()
#DOS()
plot_FS(-2.0)