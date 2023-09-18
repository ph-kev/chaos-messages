using DifferentialEquations
using Plots
using WAV
using Interpolations
using BenchmarkTools
using LaTeXStrings

include("../sending_secrets_chaos.jl")

function main()
  u0 = [2.2,1.3,2.0] # initial condition
  p = [16.0, 4.0, 45.6] # parameters of the dynamical system 
  tspan = (0,4.0) # length of the message is 4 seconds 

  ### Talker ###
  # Create secret message 
  message_unencrypted = convert_message_to_samples("taunt.wav")
  secret_message = create_secret_message(u0, p, tspan, message_unencrypted)

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
  convert_samples_to_message(decrypted_message, sampling_rate, num_of_samples, "tauntDecrpyted.wav")

  abs_error = error_set_up(message_unencrypted, decrypted_message)

  # Plot error between original sound file and decrypted sound file 
  error_plot = plot(abs_error, tspan..., legend = false, xlabel=L"t", ylabel=L"E(t)", color = "red", linewidth=0.5)

  # Helper function to make anything into a function 
  function function_set_up(f)
    function func(t)
      val = f(t)
      return val
    end
    return func
  end 

  # Plot unencrypted message and encrypted message 
  sound_plot = plot(function_set_up(message_unencrypted),tspan...,label="Original", xlabel=L"t", ylabel="Amplitude", palette = :seaborn_colorblind, linewidth=0.1)
  plot!(decrypted_message,tspan...,label="Recovered", linewidth=0.5)
  combined_plot = plot(sound_plot, error_plot, dpi = 900, palette = :seaborn_colorblind)

  display(combined_plot)
  savefig("combined_error_sound_plot_paper.png")
end

main()