using Plots

function main()

    #Initialize program variables
    Config = config()
    t, dt = 0., Config.dt₀
    X = [Config.S₀, Config.I₀]

    Time = [t]
    State = [X... 1-sum(X)]

    rhs = rhs_set(Config)
    simulate::Bool = true

    #Start simulation
    println("progress\tt\tdt")
    println(trunc(Int, (100*t/Config.t₀)),"%\t",t,"\t", dt)
    while simulate

        #Check simulation end
        simulate = t < Config.t₀ ? true : false

        #Evolve system
        t = step!(t, dt, X, rhs)

        #Save system state
        append!(Time, t)
        State = vcat(State, [X... 1-sum(X)])
        println(trunc(Int, (100*t/Config.t₀)),"%\t",t,"\t", dt)

    end

    #Plot
    gr()
    plot(Time, State, title="SIR model with β=0.35,γ=0.08", label=["S" "I" "R"])
    png("graph_1.png")

end

function config()
    Config = (β=0.35,γ=0.08,
              t₀= 60.,dt₀ = 0.01,
              S₀=0.999,I₀=0.001,R₀=0.)
    println("Config parameters:")
    dump(Config)
    println()
    return Config
end

function rhs_set(Config)
    rhs = let β = Config.β, γ = Config.γ
    (t::Real, X) -> [-β*X[1]*X[2], (β*X[1]-γ)*X[2]] 
    end
    return rhs
end

function step!(t::Real, dt::Real, X, rhs)

    dX₁ = dt*rhs(t, X)
    dX₂ = dt*rhs(t+dt/2, X.+dX₁./2)
    dX₃ = dt*rhs(t+dt/2, X.+dX₂./2)
    dX₄ = dt*rhs(t+dt, X.+dX₃)

    X .+= (dX₁.+2 .*(dX₂.+dX₃).+dX₄)./6

    return t + dt
end

main()
