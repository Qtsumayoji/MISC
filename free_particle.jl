using Distributions
using PyCall
@pyimport pylab as plt
@pyimport seaborn as sns

include("./model.jl")

#unit cell内のサイト数
const ns = 1
const Nx = 6
const Ny = 1
const Ns = Nx*Ny*ns
const nstate = 4^Ns
const filling = 1

#fermion factorの計算に使用
bitfrag = Array{Int64}(32)
for i in 1:32
    bitfrag[i] = ~((1 << i) - 1)
    #println(bits(bitfrag[i]))
end

function calc_basis()
    basis = [[] for i in 1:32]
    for i in 1:nstate
        push!(basis[count_ones(i)], i)
    end
    return basis
end

function calc_fermion_factor(i::Int64, state::Int64)
    return (-1)^(count_ones(state & bitfrag[i]))
end

function test_calc_fermion_factor()
    test1 = 21343154
    println(bits(test1))
    println(calc_fermion_factor(4, test1)," ",count_ones(test1 & bitfrag[4]))
    test2 =122321
    println(bits(test2))
    println(calc_fermion_factor(7, test2)," ",count_ones(test2 & bitfrag[7]))
end

#c_i*a_j
function ciaj(i::Int64, j::Int64, state::Int64)
    if state & (1 << (j - 1)) == 0
        return 0, 0
    else
        if state & (1 << (i - 1)) == 0
            #do マスクを作っていおいてxorで高速化
            sign　= calc_fermion_factor(j, state)
            state = xor(state , (1 << (j - 1)))
            sign *= calc_fermion_factor(i, state)
            state = xor(state , (1 << (i - 1)))
            return sign, state
        else
            return 0, 0
        end
    end
end

function test_ciaj()
    println(bits(21343154))
    sig,state = ciaj(3,5,21343154)
    println(bits(state))
    println(sig)
end

function calc_Hkij(i::Int64, j::Int64, state::Int64, t::Float64, H_mat::Array{Float64,2}, sisab::Dict)
    #up spin
    sig, ket = ciaj(i, j, state)
    if ket != 0
        id1 = sisab[state]
        id2 = sisab[ket]
        H_mat[id1, id2] += -t*sig
    end
    sig, ket = ciaj(j, i, state)
    if ket != 0
        id1 = sisab[state]
        id2 = sisab[ket]
        H_mat[id1, id2] += -t*sig
    end

    #down spin
    sig, ket = ciaj(i + Ns, j + Ns, state)
    if ket != 0
        id1 = sisab[state]
        id2 = sisab[ket]
        H_mat[id1, id2] += -t*sig
    end
    sig, ket = ciaj(j + Ns, i + Ns, state)
    if ket != 0    
        id1 = sisab[state]
        id2 = sisab[ket]
        H_mat[id1, id2] += -t*sig
    end
end

#n↓n↑
function coulmb_repulsion(i::Int64, state::Int64)
    if state & (1 << (i - 1)) == 0
        return 0
    else
        if state & (1 << (i - 1 + Ns)) == 0
            return 0
        else
            return 1
        end
    end
end

#n
function onsite_energy(i::Int64, state::Int64)
    if state & (1 << (i - 1)) == 0
        return 0
    else
        return 1
    end
end

function calc_Hv(state::Int64, U::Float64, μ::Float64, H_mat::Array{Float64, 2}, sisab::Dict)
    id1 = sisab[state]
    for i in 1:Ns
        #n↑n↓
        ket = coulmb_repulsion(i, state)
        if ket != 0
            H_mat[id1, id1] += U
        end
        #n↑
        ket = onsite_energy(i, state)
        if ket != 0
            H_mat[id1, id1] += -μ
        end
        #n↓
        ket = onsite_energy(i + Ns, state)
        if ket != 0
            H_mat[id1, id1] += -μ
        end
    end
end

#Hが二次形式で表せるときは局所粒子数密度ni=∑_l |Uil^2 f(El)とかける
#nagaiさんのノート　"演算子の二次形式で表現されたハミルトニアンの物理量について"より
function calc_local_particle_density(i, U, E, β)
    ni = 0.

    dim = length(E)
    for l in 1:dim
        ni += U[i, l]^2/(exp(β*E[l]) + 1.)
    end

    return ni
end

function calc_hermiltonian_dense(basis, link_list)
    L = length(basis)
    println("dim = ",L)
    H_mat = zeros(Float64, L, L)
    sisab = Dict()
    for i in 1:L
        sisab[basis[i]] = i
        #println(bits(basis[i]))
    end

    t = 1.0
    U = 0.
    μ = 0.

    println("Calculate Hermitian")
    @time for state in basis
        for i in 1:Ns
            link = link_list[i]
            for j in link
                calc_Hkij(i, j, state, t, H_mat, sisab)
            end
        end
        calc_Hv(state, U, μ, H_mat, sisab)
    end

    H_mat = Symmetric(H_mat)

    return H_mat
end

function testHubbardModelOnSquareLattice()
    link_mat, link_list, pos = Model.square_lattice(Nx, Ny)
    unit_vec = Model.get_square_lattice_unit_vec()
    Model.show_links(link_mat)

    basis = calc_basis()
    basis = basis[Integer(Ns*filling)]
    L = length(basis)
    H_mat = calc_hermiltonian_dense(basis, link_list)

    E, U = @time eig(H_mat)

    sns.set_style("white")
    plt.figure(figsize=(10, 8))
    plt.subplot(211)
    plt.hist(E, bins = 100)

    n_dens = zeros(Float64, Ns)
    β = 0.5
    for i in 1:Ns
        n_dens[i] = calc_local_particle_density(i, U, E, β)
    end
    println(sum(n_dens), " ", Ns*filling)

    #plotがうまくいかないっす
    plt.subplot(212)
    plt.plot([1:Ns, n_dens])
    plt.show()
end


function main()
    testHubbardModelOnSquareLattice()
end

main()