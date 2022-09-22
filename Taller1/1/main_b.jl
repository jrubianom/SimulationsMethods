using LinearAlgebra
using Plots

function main()

    #Initialize program variables
    Config = config()
    t, β, S∞, δ = 0., 0., 0., 0.

    X = []
    Y = []

    simulate::Bool = true

    #Start simulation
    println("progress\tβ\tδ")
    while simulate

        #Find tendency of S
        S∞, δ = tendency(Config, β)
        append!(X, β/Config.γ)
        append!(Y, S∞)

        #Print state
        println(trunc(Int, (100*β/Config.β₀)),"%\t",β,"\t", δ)
        β += Config.dβ
        simulate = β < Config.β₀ ? true : false

    end

    #Plot
    gr()
    plot(X, Y, title="S∞ vs R₀", xlabel="R₀", ylabel="S∞", ylims=(0,1), xlims=(0,Config.β₀/Config.γ))
    png("graph_b.png")

end

function config()

    dt = 0.01
    S₀  = 0.999
    I₀  = 0.001
    R₀  = 1-(S₀+I₀)
    γ   = 0.08
    β₀  = 5*γ
    dβ  = γ/100
    ϵ   = 1e-10

    Config = (β₀=β₀, dβ=dβ, γ=γ,
              dt=dt,
              S₀=S₀, I₀=I₀, R₀=R₀,
              ϵ=ϵ)
    println("Config parameters:")
    dump(Config)
    println()

    return Config
end

function rhs_set(Config, β)
    rhs = let β = β, γ = Config.γ
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

function tendency(Config, β)

    t  = 0.
    X  = [Config.S₀, Config.I₀]
    X⁻ = [Config.S₀, Config.I₀]
    δ = 0.

    rhs = rhs_set(Config, β)
    simulate::Bool = true

    #Start simulation
    while simulate

        #Evolve system
        X⁻ = deepcopy(X)
        t = step!(t, Config.dt, X, rhs)

        #Check simulation end
        δ = hypot(X[1]-X⁻[1],X[2]-X⁻[2])
        simulate = δ < Config.ϵ ? false : true

    end

    return X[1], δ

end

main()
