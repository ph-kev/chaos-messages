using DifferentialEquations
using Plots
using WAV
using Interpolations
using BenchmarkTools
using LaTeXStrings

include("../sending_secrets_chaos.jl")

u0 = [2.2,1.3,2.0] # initial conditions 
p=[16.0;4.0;45.6] # parameters 
tspan = (0,15.0) # time span 

"""
    parameterized_lorenz_transmitter!(du,u,p,t)

Define the dynamical system based on the Lorenz system that encrypt the message. 
"""
function parameterized_lorenz_transmitter!(du,u,p,t)
    x_T,y_T,z_T = u
    σ,β,r = p
    du[1] = σ*(y_T-x_T)
    du[2] = r*x_T - y_T - 20*(x_T * z_T)
    du[3] = 5*x_T*y_T - β*z_T
end

# Set up the ODE problem (for parameterized_lorenz_transmitter!) and solve it 
prob = ODEProblem(parameterized_lorenz_transmitter!,u0,tspan,p)
sol_transmitter = solve(prob, AutoTsit5(Rodas4P()), abstol = 1e-11, reltol = 1e-11)

# Get only the x-coordinates of the solution 
x_at_time_t_transmitter(t) = sol_transmitter(t,idxs=1) 

"""
    parameterized_lorenz_receiver!(du,u_1,p,t)

Define the dynamical system based on the Lorenz system that decrypt the message. 
"""
function parameterized_lorenz_receiver!(du,u_1,p,t)
  x_R, y_R, z_R = u_1
  σ,β,r = p
  du[1] = σ*(y_R-x_R)
  du[2] = r*x_at_time_t_transmitter(t) - y_R - 20*(x_at_time_t_transmitter(t) * z_R)
  du[3] = 5*(x_at_time_t_transmitter(t)*y_R) - β*z_R
end

u0 = [10.2,7.3,6.0] # Different intial conditions 

# Set up the ODE problem (for parameterized_lorenz_receiver!) and solve it 
prob_C = ODEProblem(parameterized_lorenz_receiver!,u0,tspan,p)
sol_receiver = solve(prob_C, AutoTsit5(Rodas4P()), abstol = 1e-11, reltol = 1e-11)

# Get only the x-coordinates of the solution 
x_at_time_t_receiver(t) = sol_receiver(t,idxs=1) 

# Plot the errors and the x-coordinates of the solutions 
abs_error = error_set_up(x_at_time_t_transmitter, x_at_time_t_receiver)

error_plot = plot(abs_error, tspan..., legend = false, xlabel=L"t", ylabel=L"E(t)",linecolor="red")
display(error_plot)
x_coord_plot = plot(x_at_time_t_transmitter, tspan...,label="Transmitter",xlabel=L"t",ylabel=L"x(t)")
plot!(x_at_time_t_receiver, tspan...,label="Receiver")

combined_plot = plot(x_coord_plot, error_plot, dpi=600)
display(combined_plot)
savefig("combined_plot.png")