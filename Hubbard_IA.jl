using PyPlot

PyPlot.rc("font",family ="Times New Roman")

t = 1.0
U = 1.0
μ = U/2.0
n_σ = 0.5
η = 0.1*t

Nx = 1000

function Σ_atom(ω)
    return U*n_σ + U^2.0*n_σ*(1.0 - n_σ)/(ω + μ - U*(1 - n_σ))
end

function ϵ_k(k)
    return -2.0*t*cos(k) - μ
end

function exact()
    g = 2.0*π
    Q = zeros(Nx)
    for i in 1:Nx
        Q[i] = (i-1)/Nx*g - g/2.0
    end

    G1 = zeros(Complex, Nx)
    G2 = zeros(Complex, Nx)
    ω1 = zeros(Nx)
    ω2 = zeros(Nx)

    for i in 1:Nx
        q = Q[i]
        ϵ_q = ϵ_k(q)
        b = 2.0*μ - ϵ_q - U
        c = μ*(ϵ_q - μ + U*n_σ) - U*(1 - n_σ)*(ϵ_q - μ + 2.0*U*n_σ)
        x1 = (-b + sqrt(b^2.0 - 4.0*c))/2.0
        x2 = (-b - sqrt(b^2.0 - 4.0*c))/2.0
        ω1[i] = x1
        ω2[i] = x2
    end

    PyPlot.plot(Q, ω1)
    PyPlot.plot(Q, ω2)
end

function numerical()
    g = 2.0*π
    Nk = Nx + 1
    Q = zeros(Nk)
    for i in 1:Nk
        Q[i] = (i-1)/Nx*g - g/2.0
    end

    NΩ = 1000
    Ω = range(-3t, stop=3t, length = NΩ)

    G = zeros(Complex, Nk, NΩ)

    for i in 1:Nk
        q = Q[i]
        ϵ_q = ϵ_k(q)
        for j in 1:NΩ
            ω = Ω[j]
            G[i, j] = 1.0/(ω - ϵ_q + μ - Σ_atom(ω) + im*η)
        end
    end
    
    X = zeros(Nk + 1, NΩ + 1)
    Y = zeros(Nk + 1, NΩ + 1)
    ΔX = Q[2] - Q[1]
    ΔY = Ω[2] - Ω[1]
    for i in 1:NΩ + 1
        for m in 1:Nk + 1
            X[m, i] = Q[1] + (m - 1.5)*ΔX
            Y[m, i] = Ω[1] + (i - 1.5)*ΔY
        end
    end

    # color plotの確認用
    #G[1,1] = -im*1.0
    #G[end,1] = -im*2.0
    #G[end,end] = -im*3.0
    #G[1,end] = -im*4.0

    PyPlot.pcolormesh(X, Y, -1.0/pi*imag(G))
    PyPlot.colorbar()
end

function main()
    PyPlot.figure(figsize=(8,6))

    exact()
    numerical()

    PyPlot.xlabel("k", size=20)
    PyPlot.ylabel("ω", size=20)
    PyPlot.tight_layout()
    PyPlot.savefig("mott_gap.png")

    PyPlot.show()
end

main()