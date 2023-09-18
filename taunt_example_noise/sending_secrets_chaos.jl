using DifferentialEquations
using Plots
using WAV
using Interpolations
using BenchmarkTools
using LaTeXStrings
using Random, Distributions

include("../sending_secrets_chaos_new.jl")

# Make anything into a function 
function function_set_up(f)
  function func(t)
    val = f(t)
    return val
  end
  return func
end 

function error_set_up(message_unencrypted, decrypted_message)
  function abs_error(t)
  val1 = message_unencrypted(t)
  val2 = decrypted_message(t)
  return abs(val1-val2)
  end 
  return abs_error
end

function main()
  u0 = [2.2,1.3,2.0]
  p = [16.0, 4.0, 45.6]
  tspan = (0,4.0) # length of the message is 4 seconds 

  Random.seed!(42424242)

  error_plot = plot()

  for s in [2.0, 1.0, 0.5, 0.1]
    ### Talker ###
  # Create secret message 
  message_unencrypted = convert_message_to_samples("taunt.wav")
  message_unencrypted_noisy = convert_message_to_samples_noisy("taunt.wav", s)

  secret_message = create_secret_message(u0, p, tspan, message_unencrypted_noisy)

  # Convert secret message to WAV file 
  sample, sampling_rate = wavread("taunt.wav")
  sample = vec(sample)
  num_of_samples = length(sample)
  convert_samples_to_message(secret_message, sampling_rate, num_of_samples, "tauntSecret.wav")

  ### Receiver ###
  # Get secret message 
  secret_message = convert_message_to_samples("tauntSecret.wav")

  # Decrpyt message 
  decrypted_message = decrypt_secret_message(u0, p, tspan, secret_message)

  # Get the sampling rate and number of samples 
  sample, sampling_rate = wavread("tauntSecret.wav")
  sample = vec(sample)
  num_of_samples = length(sample)

  # Convert secret_message to a wav file 
  convert_samples_to_message(decrypted_message, sampling_rate, num_of_samples, "tauntDecrpyted" * string(s) * ".wav")

  abs_error = error_set_up(message_unencrypted, decrypted_message)

  plot!(error_plot, abs_error, tspan..., xlabel=L"t", ylabel=L"E(t)", labels = L"\sigma = " * string(s), palette = :seaborn_colorblind, linewidth = 0.5, dpi = 600)
  end 

  display(error_plot)
  savefig("error_plot_noise.png")
end

main()