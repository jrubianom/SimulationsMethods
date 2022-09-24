using DelimitedFiles
using Combinatorics
using LinearAlgebra
using StaticArrays
using Plots

function main()

    #Initialize system
    Config = config()

    t, tᵛ, dt = 0., 0., Config.dt
    iteration, vis_iteration = 0, 0

    r₁ = -Config.r*Config.m₂/(Config.m₁+Config.m₂)
    r₂ =  Config.r*Config.m₁/(Config.m₁+Config.m₂)

    particles = [Particle() for ii ∈ 1:3]
    particles[1] = Particle(Config.m₁, SVector(r₁,0.,0.), SVector(0.,r₁*Config.ω,0.))
    particles[2] = Particle(Config.m₂, SVector(r₂,0.,0.), SVector(0.,r₂*Config.ω,0.))
    particles[3] = Particle(0.005, SVector((1/2)r₂,(√3/2)r₂,0.), SVector((-√3/2)r₂*Config.ω,(1/2)r₂*Config.ω,0.))
    collide!(particles, t)

    print_state(particles, atan(particles[2].R[2],particles[2].R[1]), vis_iteration)

    #Start simulation
    simulate::Bool = true
    println("progress\tt(T)\tvisualization")
    println(trunc(Int, (100*t/Config.tᶠ)),"%\t",round(t/Config.T,digits=2),"\t", vis_iteration)

    x_history = [(1/2)r₂] 
    t_history = [t]

    while (simulate)
        #Evolve system
        step!(t, dt, particles)
        t += dt
        tᵛ += dt

        θ = atan(particles[2].R[2],particles[2].R[1])
        push!(x_history, particles[3].R[1]*cos(θ)+particles[3].R[2]*sin(θ))
        push!(t_history, t/Config.T)

        #Print system
        if (tᵛ >= Config.dtᵛ)
            vis_iteration += 1; tᵛ=0.
            print_state(particles, atan(particles[2].R[2],particles[2].R[1]), vis_iteration)
        end
        println(trunc(Int, (100*t/Config.tᶠ)),"%\t",round(t/Config.T,digits=2),"\t", vis_iteration)

        #Check end
        simulate = t <= Config.tᶠ ? true : false
    end

    #Animate results
    gr()
    plot(t_history, x_history, xlims=(0.,20.))
    png("data_e/graph.png")

end

function config()

    m₁ = 1047.
    m₂ = 1.

    r = 1000.

    ω = √((m₁+m₂)/r^3)
    T = 2π/ω
    tᶠ= 20T
    dt = T/1000
    dtᵛ = T/10

    Config = (m₁=m₁,m₂=m₂,r=r,
              ω=ω, T=T,
              tᶠ=tᶠ,dt=dt,dtᵛ=dtᵛ)
    println("Config parameters:")
    dump(Config)
    println()

    return Config
end

struct Particle

    m::Float64
    R::SVector{3,Float64}
    V::SVector{3,Float64}
    F::SVector{3,Float64}

    Particle() = new(1., SVector(0.,0.,0.), SVector(0.,0.,0.), SVector(0.,0.,0.))
    Particle(m::Float64, R::SVector{3,Float64}) = new(m, R, SVector(0.,0.,0.), SVector(0.,0.,0.))
    Particle(m::Float64, R::SVector{3,Float64}, V::SVector{3,Float64}) = new(m, R, V, SVector(0.,0.,0.))
    Particle(m::Float64, R::SVector{3,Float64}, V::SVector{3,Float64}, F::SVector{3,Float64}) = new(m, R, V, F)

end

move_pos(p::Particle, dt::Float64) = Particle(p.m, p.R.+dt.*p.V, p.V, p.F)

move_vel(p::Particle, dt::Float64) = Particle(p.m, p.R, p.V.+(dt/p.m).*p.F, p.F)

set_force(p::Particle, F::SVector{3,Float64}) = Particle(p.m, p.R, p.V, F)

function collide!(particles::Array{Particle}, t::Float64)

    #Global forces
    for i ∈ eachindex(particles)
        particles[i] = set_force(particles[i], SVector(0.,0.,0.))
    end

    #Interaction forces
    for (i, j) ∈ combinations(eachindex(particles), 2) 

        R12 = particles[j].R-particles[i].R; d12 = norm(R12)

        particles[i] = set_force(particles[i], particles[i].F.+( particles[i].m*particles[j].m/d12^3).*R12)
        particles[j] = set_force(particles[j], particles[j].F.+(-particles[i].m*particles[j].m/d12^3).*R12)
    end

    return

end

function step!(t::Real, dt::Float64, particles::Array{Particle})

    const1::Float64 = 0.1786178958448091;   #ζ
    const4::Float64 = -0.2123418310626054;  #λ
    const3::Float64 = -0.06626458266981849; #χ
    const2::Float64 = (1-2*const4)/2;       #(1-2*λ)/2
    const5::Float64 = 1-2*(const3+const1);  #1+2*(ζ+χ)

    particles .= move_pos.(particles, dt*const1)
    collide!(particles, t)
    particles .= move_vel.(particles, dt*const2)
    particles .= move_pos.(particles, dt*const3)
    collide!(particles, t)
    particles .= move_vel.(particles, dt*const4)
    particles .= move_pos.(particles, dt*const5)
    collide!(particles, t)
    particles .= move_vel.(particles, dt*const4)
    particles .= move_pos.(particles, dt*const3)
    collide!(particles, t)
    particles .= move_vel.(particles, dt*const2)
    particles .= move_pos.(particles, dt*const1)

    return

end

function print_state(particles::Array{Particle}, θ::Float64, vis_iteration::Int)
    open("data_d/data"*string(vis_iteration)*".txt", "w") do io
        for p in particles
            println(io, p.m," ",
                        p.R[1]*cos(θ)+p.R[2]*sin(θ)," ",p.R[2]*cos(θ)-p.R[1]*sin(θ)," ",p.R[3]," ",
                        p.V[1]*cos(θ)+p.V[2]*sin(θ)," ",p.V[2]*cos(θ)-p.V[1]*sin(θ)," ",p.V[3]," ",
                        p.F[1]*cos(θ)+p.F[2]*sin(θ)," ",p.F[2]*cos(θ)-p.F[1]*sin(θ)," ",p.F[3])
        end
    end
end

main()
