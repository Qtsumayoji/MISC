using PyPlot
using LinearAlgebra

function calc_Hk(k)
    return -2.0*cos(k)
end

function calc_integral(x, y)
    dim = length(x)
    dx = x[2] - x[1]
    I = 0.0

    for i in 1:dim-1
        I += y[i+1] + y[i]
    end
    I *= 0.5*dx

    return I
end

function convert_Matsubara_freq(Ω, G, ωn)
    NΩ = length(Ω)
    ImG_ωn = 0.0
    Integrands = zeros(NΩ)

    for i in 1:NΩ
        ω = Ω[i]
        Integrands[i] = ωn/(ω^2.0 + ωn^2.0)*imag(G[i])
    end
    
    ImG_ωn = 1.0/pi*calc_integral(Ω, Integrands)

    return ImG_ωn
end


function main()
    η = 0.01
    Ns = 100

    A = Float64[]
    B = ComplexF64[]
    g = 2π
    NΩ = 100000
    Ω = range(-3.0, stop=3.0, length=NΩ)
    Hk = zeros(Ns)
    K = zeros(Ns)
    Σ_imp = zeros(NΩ)
    G_local = zeros(ComplexF64, NΩ)
    G_bath = zeros(ComplexF64, NΩ)

    for m in 1:Ns
        K[m] = -g/2 + g*(m - 1)/Ns
        Hk[m] = calc_Hk(K[m])
    end

    @time for i in 1:NΩ
        ω = Ω[i]
        Σ = Σ_imp[i]
        for m in 1:Ns
            G_local[i] += 1.0/(ω - Hk[m] + Σ - im*η)
        end
        G_local[i] /= Ns
        G_bath[i] = 1.0/(Σ_imp[i] + 1.0/G_local[i])
    end

    is_up = true
    is_down = false
    for i in 1:NΩ-1
        tmp = imag(G_local[i+1] - G_local[i])

        if tmp < 0.0 && is_up
            push!(A, Ω[i])
            push!(B, sqrt(-im*G_local[i]*η))
            println("A=",Ω[i]," B=",1.0/pi*imag(G_local[i]))
            is_up = false
            is_down = true
        end

        if tmp > 0.0 && is_down
            is_up = true
            is_down = false
        end
    end

    println("len=",length(A))

    #PyPlot.plot(Ωn, ImG_ωn)
    #PyPlot.show()

    G_cf = zeros(ComplexF64, NΩ)
    for i in 1:NΩ
        ω = Ω[i]
        I = 0.0
        for j in 1:length(A)
            I += B[j]^2.0/(ω - A[j] - im*η)
        end
        G_cf[i] = 1.0/pi*I
    end

    G_imp = zeros(ComplexF64, NΩ)
    for i in 1:NΩ
        ω = Ω[i]
        I = ω - A[1]
        for j in 2:length(A)
            I -= B[j]^2.0/(ω - A[j] - im*η)
        end
        I = 1.0/I
        G_imp[i] = 1.0/pi*I
    end

    #PyPlot.plot(Ω, real(G_cf),label="re")
    PyPlot.plot(Ω, imag(G_imp),label="cf im")

    PyPlot.plot(Ω, 1.0/pi*imag(G_local), label="local Im")
    #PyPlot.plot(Ω, 1.0/pi*real(G_local), label="Re")
    PyPlot.legend()
    PyPlot.show()

end
main()

function convert_MF()
    η = 0.001
    Ns = 1000

    A = Float64[]
    B = ComplexF64[]
    g = 2π
    NΩ = 100000
    Ω = range(-3.0, stop=3.0, length=NΩ)
    Hk = zeros(Ns)
    K = zeros(Ns)
    Σ_imp = zeros(NΩ)
    G_local = zeros(ComplexF64, NΩ)
    G_bath = zeros(ComplexF64, NΩ)

    for m in 1:Ns
        K[m] = -g/2 + g*(m - 1)/Ns
        Hk[m] = calc_Hk(K[m])
    end

    @time for i in 1:NΩ
        ω = Ω[i]
        Σ = Σ_imp[i]
        for m in 1:Ns
            G_local[i] += 1.0/(ω - Hk[m] + Σ - im*η)
        end
        G_local[i] /= Ns
        G_bath[i] = 1.0/(Σ_imp[i] + 1.0/G_local[i])
    end
    
    T = 0.001
    N_MF = 1000
    Ωn = zeros(N_MF)
    ImG_ωn = zeros(N_MF)
    for i in 1:N_MF
        n = i - 1
        ωn = (2.0*n + 1)*pi*T
        Ωn[i] = ωn
        ImG_ωn[i] = convert_Matsubara_freq(Ω, G_local, ωn)
    end
end