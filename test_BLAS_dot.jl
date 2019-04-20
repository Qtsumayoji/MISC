using LinearAlgebra

function main()
    nGB = 1
    n = Int64(nGB*2^27)

    #nKB = 5
    #n = Int64(nKB*2^7)
    vec1 = rand(n)
    vec2 = rand(n)
    @time for i in 1:100
        @time vec1'*vec2
    end
end

main()