using DifferentialEquations
using Plots
using WAV
using Interpolations
using BenchmarkTools
using LaTeXStrings

include("../sending_secrets_chaos.jl")

function main()
  u0 = [2.2,1.3,2.0] # initial condition 
  p=[10.0; 0.3333333; 60.0] # parameters of the dynamical system 
  tspan = (0,2000.0) # length of time span 

  ### Talker ###
  # Create secret message 
  f(x) = 0.0 
  secret_message = create_secret_message(u0, p, tspan, f)

  ### Receiver ###
  # Decrypt message 
  decrypted_message = decrypt_secret_message(u0, p, tspan, secret_message)

  # Plot error between constant function and function we get back 
  abs_error = error_set_up(f, decrypted_message)
  display(plot(abs_error, tspan..., legend = false, xlabel=L"t", ylabel=L"E(t)"))
  savefig("error.pdf")
end

main()