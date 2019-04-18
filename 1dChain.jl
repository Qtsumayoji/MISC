function main()
    Ns = 10
    # 3次近接
    for i in 1:Ns
        println(i," ",(i%Ns)+1," ",(i%Ns)+2)
    end
end

main()