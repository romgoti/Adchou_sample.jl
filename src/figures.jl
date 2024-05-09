using Adchou_sample
using Plots

function figure7()

    T = 250
    m1 = solve_2states(T)
    m2 = solve_2states(T,.15)

    begin
        plot(range(0,8,T), [m1[:V][:,1] m1[:cf][:,1] m1[:sf][:,1] m1[:g][:,1]],line=(2,:blue), layout = 4,
        title = ["Value Functions" "Consumption Functions" "Saving Functions" "Stationary Distributions" ], 
        label = ["V1(a)" "C1(a)" "S1(a)" "g1(a)"])
        plot!(range(0,8,T), [m1[:V][:,2] m1[:cf][:,2] m1[:sf][:,2] m1[:g][:,2]],line=(2,:red), 
        label = ["V2(a)" "C2(a)" "S2(a)" "g2(a)"])
        plot!(range(0,8,T), [m2[:V][:,1] m2[:cf][:,1] m2[:sf][:,1] m2[:g][:,1]],line=(2,:blue,:dash), 
        label = ["V1(a)" "C1(a)" "S1(a)" "g1(a)"])
        plot!(range(0,8,T), [m2[:V][:,2] m2[:cf][:,2] m2[:sf][:,2] m2[:g][:,2]],line=(2,:red,:dash), 
        label = ["V2(a)" "C2(a)" "S2(a)" "g2(a)"])
        xlabel!("Net assets position")
    end    
end

function figure10()
    T= 250
    mh = solve_housing(T)

    # generate alternative distribution depengin on history (scenario 2 and 3)
    mh2 = solve_housing(T,2)
    mh3 = solve_housing(T,3)

    begin
        plot(range(0,3,T),[mh[:g][:,1] mh2[:g][:,1] mh3[:g][:,1] mh[:adot][:,1] mh[:V][:,1] ],line = (2,:blue), 
        label = ["g1(a)" "g1(a)" "g1(a)" "Δa1(a)" "V1(a)"],layout = 5,
        ylabel = ["Density" "Density" "Density" "Asset choice" ""],
        title = [["Scenario $i" for j in 1:1, i in 1:3]... [ "Saving Policy" "Value Function"]...])
        plot!(range(0,3,T),[ mh[:g][:,2] mh2[:g][:,2] mh3[:g][:,2] mh[:adot][:,2] mh[:V][:,2]],line = (2,:red), 
        label = ["g2(a)" "g2(a)" "g2(a)" "Δa2(a)" "V2(a)" ])
        vline!([mh[:astar] mh[:astar] mh[:astar] mh[:astar]],line=(1,:lightgrey,:dash), label = "a*"  )
        xlabel!("Net assets position")
    end
    
end
