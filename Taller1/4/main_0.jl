using DelimitedFiles
using Combinatorics
using LinearAlgebra
using StaticArrays
using Plots

function main()

    #Initialize system
    Config = config(1e9)

    t, tᵛ, dt = 0., 0., Config.dt
    iteration, vis_iteration = 0, 0

    particles = [Particle() for ii ∈ 1:3]
    particles[1] = Particle(Config.m, Config.r, Config.l, SVector(-2Config.r,Config.l), -π/12)
    particles[2] = Particle(Config.m, Config.r, Config.l, SVector(0.,Config.l), 0.)
    particles[3] = Particle(Config.m, Config.r, Config.l, SVector(2Config.r,Config.l), 0.)
    collide!(Config, particles, t)

    print_state(particles, vis_iteration)

    #Start simulation
    simulate::Bool = true
    println("progress\tt(T)\tvisualization")
    println(trunc(Int, (100*t/Config.tᶠ)),"%\t",round(t,digits=2),"\t", vis_iteration)
    
    while (simulate)
        #Evolve system
        step!(Config, t, dt, particles)
        t += dt
        tᵛ += dt

        #Print system
        if (tᵛ >= Config.dtᵛ)
            vis_iteration += 1; tᵛ=0.
            print_state(particles, vis_iteration)
        end
        println(trunc(Int, (100*t/Config.tᶠ)),"%\t",round(t,digits=2),"\t", vis_iteration)

        #Check end
        simulate = t <= Config.tᶠ ? true : false
    end

    #Animate results
    gr()
    anim = @animate for i ∈ 0:vis_iteration
        Data = readdlm("data_a/data$i.txt")
        scatter(Data[:,4],Data[:,5], 
                xlims=(-10,10), ylims=(-2,12), aspect_ratio=1,
                legend=false,
                title=string(round(i*Config.dtᵛ,digits=2)))

        for j ∈ 1:size(Data)[1]
            plot!(circle_shape(Data[j,1],Data[j,4],Data[j,5]), linecolor = :black, c = :blue, fillalpha = 0.2)
            plot!(line_shape(Data[j,2],Data[j,3],Data[j,4],Data[j,5]), linecolor = :black)
        end
    end
    gif(anim, "data_a/anim.gif", fps = 10)

end

function config(K::Float64)

    m = 100.
    l = 12.
    r = 1.5

    g = 980.

    tᶠ= 10
    dt = 0.001
    dtᵛ = 0.1

    Config = (m=m,r=r,l=l,
              g=g,K=K,
              tᶠ=tᶠ,dt=dt,dtᵛ=dtᵛ)
    println("Config parameters:")
    dump(Config)
    println()

    return Config
end

struct Particle

    m::Float64
    r::Float64
    l::Float64
    R₀::SVector{2,Float64}
    θ::Float64
    ω::Float64
    τ::Float64
    R::SVector{2,Float64}

    Particle() = new(1., 1., 1., SVector(0.,1.), 0., 0., 0., SVector(0., 0.))
    Particle(m::Float64, r::Float64, l::Float64, R₀::SVector{2,Float64}, 
             θ::Float64, ω::Float64, τ::Float64) = new(m,r,l,R₀,θ,ω,τ,R₀.+l.*SVector(sin(θ),-cos(θ)))
    Particle(m::Float64, r::Float64, l::Float64, R₀::SVector{2,Float64}, 
             θ::Float64, ω::Float64) = new(m,r,l,R₀,θ,ω,0.,R₀.+l.*SVector(sin(θ),-cos(θ)))
    Particle(m::Float64, r::Float64, l::Float64, R₀::SVector{2,Float64}, 
             θ::Float64) = new(m,r,l,R₀,θ,0.,0.,R₀.+l.*SVector(sin(θ),-cos(θ)))

end

move(p::Particle, dt::Float64) = Particle(p.m, p.r, p.l, p.R₀, p.θ+dt*p.ω, p.ω, p.τ)

move(p::Particle, dt₁::Float64, dt₂::Float64) = Particle(p.m, p.r, p.l, p.R₀, p.θ+dt₂*p.ω, p.ω+(dt₁/(p.m*p.l^2))*p.τ, p.τ)

set_torque(p::Particle, τ::Float64) = Particle(p.m, p.r, p.l, p.R₀, p.θ, p.ω, τ)

function collide!(Config, particles::Array{Particle}, t::Float64)

    #Global forces
    for i ∈ eachindex(particles)
        particles[i] = set_torque(particles[i], -particles[i].m*Config.g*sin(particles[i].θ))
    end

    #Interaction forces
    for (i, j) ∈ combinations(eachindex(particles), 2) 
        ΔR = particles[i].R - particles[j].R

        d = norm(ΔR)
        s = particles[i].r + particles[j].r - d
        if (s > 0)
            particles[i] = set_torque(particles[i], particles[i].τ + (Config.K*s^(3/2)*particles[i].l/d)*
                                     (ΔR[1]*cos(particles[i].θ)+ΔR[2]*sin(particles[i].θ)))
            particles[j] = set_torque(particles[j], particles[j].τ - (Config.K*s^(3/2)*particles[j].l/d)*
                                     (ΔR[1]*cos(particles[j].θ)+ΔR[2]*sin(particles[j].θ)))
        end
    end

    return
end

function step!(Config, t::Real, dt::Float64, particles::Array{Particle})

    const1::Float64 = 0.1786178958448091;   #ζ
    const4::Float64 = -0.2123418310626054;  #λ
    const3::Float64 = -0.06626458266981849; #χ
    const2::Float64 = (1-2*const4)/2;       #(1-2*λ)/2
    const5::Float64 = 1-2*(const3+const1);  #1+2*(ζ+χ)

    particles .= move.(particles, dt*const1)
    collide!(Config, particles, t)
    particles .= move.(particles, dt*const2, dt*const3)
    collide!(Config, particles, t)
    particles .= move.(particles, dt*const4, dt*const5)
    collide!(Config, particles, t)
    particles .= move.(particles, dt*const4, dt*const3)
    collide!(Config, particles, t)
    particles .= move.(particles, dt*const2, dt*const1)

    return

end

function print_state(particles::Array{Particle}, vis_iteration::Int)
    open("data_a/data"*string(vis_iteration)*".txt", "w") do io
        for p in particles
            println(io, p.r," ",
                        p.R₀[1]," ",p.R₀[2]," ",
                        p.R[1]," ",p.R[2]," ",
                        p.θ," ",p.ω," ",p.τ)
        end
    end
end

function circle_shape(r::Float64, x::Float64, y::Float64)
    θ = LinRange(0, 2π, 500)
    return x.+r.*sin.(θ), y.+r.*cos.(θ)
end

function line_shape(x₀::Float64, y₀::Float64, x₁::Float64, y₁::Float64)
    t = LinRange(0, 1, 500)
    return x₀.+(x₁-x₀).*t, y₀.+(y₁-y₀).*t 
end

main()
