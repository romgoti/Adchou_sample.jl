
using BenchmarkTools
using Adchou_sample
using Plots

# @benchmark Adchou_sample.solve_2states(1000)
# @benchmark Adchou_sample.solve_housing()

# REPLICATE FIGURES

##### Figure 7 ######

T = 500
m1 = Adchou_sample.solve_2states(T)
m2 = Adchou_sample.solve_2states(T,.15)
# m3 = solve_2states(T,.1,3)

begin
    plot(range(0,8,T),m1[:V][:,1],line=(2,:blue), label = "V1(a)",yticks =nothing)
    plot!(range(0,8,T),m1[:V][:,2],line=(2,:red), label = "V2(a)")
    plot!(range(0,8,T),m2[:V][:,1],line=(2,:blue,:dash), label = "V1(a)")
    plot!(range(0,8,T),m2[:V][:,2],line=(2,:red,:dash), label = "V2(a)")
end

begin
    plot(range(0,8,T),m1[:sf][:,1],line=(2,:blue), label = "V1(a)")
    plot!(range(0,8,T),m1[:sf][:,2],line=(2,:red), label = "V2(a)")
    plot!(range(0,8,T),m2[:sf][:,1],line=(2,:blue,:dash), label = "V1(a)")
    plot!(range(0,8,T),m2[:sf][:,2],line=(2,:red,:dash), label = "V2(a)")
end    

begin
    plot(range(0,8,T),m1[:g][:,1],line=(2,:blue), label = "g1(a)")
    plot!(range(0,8,T),m1[:g][:,2],line=(2,:red), label = "g2(a)")
    vline!([m1[:K]],line=(2,:lightgrey), label = "K*" )
    vline!([m2[:K]],line=(2,:lightgrey,:dash), label = "K*" )
    plot!(range(0,8,T),m2[:g][:,1],line=(2,:blue,:dash), label = "g1(a)")
    plot!(range(0,8,T),m2[:g][:,2],line=(2,:red,:dash), label = "g2(a)")
end

# begin
#     plot(range(0,8,T),m1[:g][:,1],line=(2,:blue), label = "g1(a)")
#     plot!(range(0,8,T),m1[:g][:,2],line=(2,:red), label = "g2(a)")
#     vline!([m1[:K]],line=(2,:lightgrey), label = "K*" )
#     vline!([m3[:K]],line=(2,:lightgrey,:dash), label = "K*" )
#     plot!(range(0,8,T),m3[:g][:,1],line=(2,:blue,:dash), label = "g1(a)")
#     plot!(range(0,8,T),m3[:g][:,2],line=(2,:red,:dash), label = "g2(a)")
# end

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


##### Figure 10 ######

# solve the model
T= 500
mh = Adchou_sample.solve_housing(T)

# generate alternative distribution depengin on history (scenario 2 and 3)
mh2 = Adchou_sample.solve_housing(T,2)
mh3 = Adchou_sample.solve_housing(T,3)

begin
    plot(range(0,3,T),mh[:V][:,1],line=(2,:blue), label = "V1(a)",yticks =nothing)
    plot!(range(0,3,T),mh[:V][:,2],line=(2,:red), label = "V2(a)")
    ylabel!("Value Function")
    xlabel!("Net assets position")
end

begin
    plot(range(0,3,T),mh[:adot][:,1],line=(2,:blue), label = "S1(a)")
    plot!(range(0,3,T),mh[:adot][:,2],line=(2,:red), label = "S2(a)")
    hline!([0.0],line=(1,:lightgrey,:dash), label = "Δa = 0" )
    vline!([mh[:astar]],line=(1,:lightgrey,:dash), label = "a*"  )
    ylabel!("Saving policy")
    xlabel!("Net assets position")
end

begin
    plot(range(0,3,T),mh[:g][:,1],line=(2,:blue), label = "g1(a)")
    plot!(range(0,3,T),mh[:g][:,2],line=(2,:red), label = "g2(a)")
    # plot!(range(0,3,T),mh2[:g][:,1],line=(1,:blue,:dash), label = "g1(a) - High hist")
    # plot!(range(0,3,T),mh2[:g][:,2],line=(1,:red,:dash), label = "g2(a) - High hist")
    # plot!(range(0,3,T),mh3[:g][:,1],line=(1,:blue,:dot), label = "g1(a) - Low hist")
    # plot!(range(0,3,T),mh3[:g][:,2],line=(1,:red,:dot), label = "g2(a) - Low hist")
    vline!([mh[:astar]],line=(1,:lightgrey,:dash), label = "a*"  )
    ylabel!("Density")
    ylims!(0., 7.)
    xlabel!("Net assets position")
end


begin
    plot(range(0,3,T),[mh[:g][:,1] mh2[:g][:,1] mh3[:g][:,1]],line = (2,:blue), label = "g1(a)",layout = 3,
    title = ["Scenario $i" for j in 1:1, i in 1:3])
    plot!(range(0,3,T),[mh[:g][:,2] mh2[:g][:,2] mh3[:g][:,2]],line = (2,:red), label = "g2(a)",layout = 3)
    vline!([mh[:astar]],line=(1,:lightgrey,:dash), label = "a*"  )
    ylabel!("Density")
    ylims!(0., 7.)
    xlabel!("Net assets position")
end


