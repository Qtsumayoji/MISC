using StatsBase

function bit_flip(bit::Int64, i::Int64)
    # 1 << i-1 は左にiずらすのではなく i番目にbitを立てるために1引いている
    return xor(bit , (1 << (i - 1)))
end

function bit_on(bit::Int64, i::Int64)
    return bit | (1 << (i - 1)) 
end

function bit_off(bit::Int64, i::Int64)
    return bit & ~(1 << (i - 1)) 
end

function bit_shift(bit::Int64, i::Int64, Δ::Int64)
    return bit_flip(bit_flip(bit, i), i + Δ)
end

function bit_run_i_to_j()
end

function two_particle_Sz_0_basis(Ns)
    dim = binomial(Ns, 1)*binomial(Ns, 1)
    basis = zeros(Int64, dim)

    cnt = 1
    for i in 1:Ns
        bit = 0
        bit = bit_on(bit, 1)
        bit = bit_on(bit, i + Ns)
        basis[cnt] = bit
        cnt += 1

        for j in 1:Ns-1
            bit = bit_shift(bit, j, 1)
            basis[cnt] = bit
            cnt += 1
        end
    end

    return basis
end

function one_particle_basis(Ns, updown)
    dim = Ns
    basis = zeros(Int64, dim)

    cnt = 1
    for i in 1:Ns
        bit = 0
        bit = bit_on(bit, 1 + updown)
        basis[cnt] = bit
        cnt += 1
    end

    return basis
end

function tow_particle_basis(Ns, updown)
    dim = binomial(Ns, 2)
    basis = zeros(Int64, dim)

    cnt = 1
    for i in 1:Ns - 1
        bit = 0
        bit = bit_on(bit, 1 + updown)
        bit = bit_on(bit, i + 1 + updown)
        basis[cnt] = bit
        cnt += 1

        for j in 1:i - 1
            bit = bit_shift(bit, j, 1 + updown)
            basis[cnt] = bit
            cnt += 1
        end
    end

    return basis
end

function make_one_bit(nbit)
    bit = 0
    for i in 0:nbit - 1
        bit += 2^i
    end
    return bit
end

function print_bit(bit)
    println(string(bit, base=2))
end


function test1()
    Ns = 5
    dim = binomial(Ns, 1)*binomial(Ns, 1)
    println("dim=",dim)

    bit = 0
    bit = bit_flip(bit, 1)
    bit = bit_flip(bit, Ns + 1)
    bit = bit_off(bit, 1)
    #print_bit(bit)

    bit = 0
    bit = bit_on(bit, 1)
    bit = bit_on(bit, Ns + 1)
    #print_bit(bit)

    print_bit(make_one_bit(Ns))

    bit = 0
    for i in 1:Ns
        bit = make_one_bit(Ns)
        bit = bit_off(bit, 1)
        bit = bit_off(bit, i + Ns)
        print_bit(bit)

        for j in 1:Ns-1
            bit = bit_shift(bit, j, 1)
            print_bit(bit)
        end
    end

    @time println(two_particle_Sz_0_basis(Ns))
end
test1()