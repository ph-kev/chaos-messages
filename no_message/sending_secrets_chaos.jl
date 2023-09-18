using DifferentialEquations
using Plots
using WAV
using Interpolations
using BenchmarkTools
using LaTeXStrings

include("../sending_secrets_chaos_new.jl")

function main()
  u0 = [2.2,1.3,2.0]
  p=[10.0;0.3333333;60.0]
  tspan = (0,2000.0) 

  ### Talker ###
  # Create secret message 
  f(x) = 0.0
  secret_message = create_secret_message(u0, p, tspan, f)

  ### Receiver ###
  # Decrpyt message 
  decrypted_message = decrypt_secret_message(u0, p, tspan, secret_message)

  function error_set_up(message_unencrypted, decrypted_message)
    function abs_error(t)
    val1 = message_unencrypted(t)
    val2 = decrypted_message(t)
    return abs(val1-val2)
    end 
    return abs_error
  end

  abs_error = error_set_up(f, decrypted_message)
  # Plot error between constant function and function we get back 
  display(plot(abs_error, tspan..., legend = false, xlabel=L"t", ylabel=L"E(t)"))
  savefig("error.pdf")
end

main()