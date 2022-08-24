using Oceananigans
using GLMakie

grid = RectilinearGrid(size=(1, 256, 256), x=(-5, 5), y=(-5, 5), z=(0, 10))

@inline circle(x, y, t) = (x^2 + y^2) < 1.0

@inline ring(x, y, t) = ((x^2 + y^2) > 0.90) &
                        ((x^2 + y^2) < 1.00)

# w★³ ∼ Q
boundary_conditions = (;
    b = FieldBoundaryConditions(bottom=FluxBoundaryCondition(circle)),
    c = FieldBoundaryConditions(bottom=FluxBoundaryCondition(ring))
)

model = NonhydrostaticModel(; grid, boundary_conditions,
                            timestepper = :RungeKutta3,
                            advection = WENO(),
                            buoyancy = BuoyancyTracer(),
                            tracers = :b)
                            #tracers = (:b, :c))

simulation = Simulation(model, Δt=1e-3, stop_iteration=1000)

wizard = TimeStepWizard(cfl=0.8)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

progress(sim) = @info string("Iter: ",  iteration(sim), " time: ", time(sim))
simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))

simulation.output_writers[:jld2] = JLD2OutputWriter(model, model.tracers,
                                                    schedule = IterationInterval(10),
                                                    filename = "convecting_tea.jld2",
                                                    overwrite_existing = true)

run!(simulation)

ct = FieldTimeSeries("convecting_tea.jld2", "b")
t = ct.times
Nt = length(t)

fig = Figure(resolution=(800, 600))
ax = Axis(fig[1, 1])
slider = Slider(fig[2, 1], range=1:Nt, startvalue=1)
n = slider.value
cⁿ = @lift interior(ct[$n], 1, :, :)
heatmap!(ax, cⁿ)

display(fig)

