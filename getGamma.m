function gamma=getGamma(order)
    if (order == 1)
        gamma = [1/2];
    elseif (order == 2)
        gamma = [-1/6 , 2/3];
    elseif (order == 3)
        gamma = [-1/16 , 0 , 9/16];
    elseif (order == 4)
        gamma = [1/90 , -2/9 , 0 , 32/45];
    elseif (order == 5)
        gamma = [1/144 , -8/63 , 0 , 0 , 625/1008];
    elseif (order == 6)
        gamma = [-1/1680 , 1/15 , -27/80 , 0 , 0 , 27/35];
    elseif (order == 7)
        gamma = [-1/2304 , 32/675 , -729/3200 , 0 , 0 , 0 , 117649/172800];
    end

end