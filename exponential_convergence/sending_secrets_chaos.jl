using DifferentialEquations
using Plots
using WAV
using Interpolations
using BenchmarkTools
using LaTeXStrings

u0 = [2.2,1.3,2.0]
p=[16.0;4.0;45.6]
tspan = (0,15.0) 

function parameterized_lorenz_transmitter!(du,u,p,t)
    x_T,y_T,z_T = u
    σ,β,r = p
    du[1] = σ*(y_T-x_T)
    du[2] = r*x_T - y_T - 20*(x_T * z_T)
    du[3] = 5*x_T*y_T - β*z_T
end
prob = ODEProblem(parameterized_lorenz_transmitter!,u0,tspan,p)
sol_transmitter = solve(prob, AutoTsit5(Rodas4P()), abstol = 1e-11, reltol = 1e-11)

x_at_time_t_transmitter(t) = sol_transmitter(t,idxs=1) 

function parameterized_lorenz_receiver!(du,u_1,p,t)
  x_R, y_R, z_R = u_1
  σ,β,r = p
  du[1] = σ*(y_R-x_R)
  du[2] = r*x_at_time_t_transmitter(t) - y_R - 20*(x_at_time_t_transmitter(t) * z_R)
  du[3] = 5*(x_at_time_t_transmitter(t)*y_R) - β*z_R
end

u0 = [10.2,7.3,6.0]

prob_C = ODEProblem(parameterized_lorenz_receiver!,u0,tspan,p)
sol_receiver = solve(prob_C, AutoTsit5(Rodas4P()), abstol = 1e-11, reltol = 1e-11)

x_at_time_t_receiver(t) = sol_receiver(t,idxs=1) 

function error_set_up(val_of_transmitter, val_of_receiver)
  function abs_error(t)
  val1 = val_of_transmitter(t)
  val2 = val_of_receiver(t)
  return abs(val1-val2)
  end 
  return abs_error
end

abs_error = error_set_up(x_at_time_t_transmitter, x_at_time_t_receiver)

error_plot = plot(abs_error, tspan..., legend = false, xlabel=L"t", ylabel=L"E(t)",linecolor="red")
display(error_plot)

x_coord_plot = plot(x_at_time_t_transmitter, tspan...,label="Transmitter",xlabel=L"t",ylabel=L"x(t)")
plot!(x_at_time_t_receiver, tspan...,label="Receiver")

combined_plot = plot(x_coord_plot, error_plot, dpi=600)
display(combined_plot)
# savefig("combined_plot.png")