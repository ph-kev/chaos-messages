using DifferentialEquations
using Plots
using WAV
using Interpolations
using BenchmarkTools

function convert_message_to_samples(message::String)
  sample, sampling_rate = wavread(message)
  sample = vec(sample)
  time_arr = Float64.([(1/sampling_rate)*n for n in 1:length(sample)])
  message_unencrypted = linear_interpolation(time_arr, sample, extrapolation_bc = Line())
  return message_unencrypted
end

function create_secret_message(u0, p, tspan, message_unencrypted)
  function parameterized_lorenz_talker!(du,u,p,t)
    x_T,y_T,z_T = u
    σ,β,r = p
    du[1] = σ*(y_T-x_T)
    du[2] = r*x_T - y_T - 20*(x_T * z_T)
    du[3] = 5*x_T*y_T - β*z_T
  end
  prob = ODEProblem(parameterized_lorenz_talker!,u0,tspan,p)
  sol_talker = solve(prob, AutoTsit5(Rodas4P()), abstol = 1e-13, reltol = 1e-13)

  function secret_message(t)
    hidden = 1e-7*message_unencrypted(t)
    secret_at_time_t = sol_talker(t,idxs=1) + hidden
    return secret_at_time_t
  end

  return secret_message
end 

function decrypt_secret_message(u0, p, tspan, secret_message)
  function parameterized_lorenz_receiver!(du,u_1,p,t)
    x_R, y_R, z_R = u_1
    σ,β,r = p
    du[1] = σ*(y_R-x_R)
    du[2] = r*secret_message(t) - y_R - 20*(secret_message(t) * z_R)
    du[3] = 5*(secret_message(t)*y_R) - β*z_R
  end

  prob_C = ODEProblem(parameterized_lorenz_receiver!,u0,tspan,p)
  sol_receiver = solve(prob_C, AutoTsit5(Rodas4P()), abstol = 1e-13, reltol = 1e-13)

  function decrypted_message(t)
    return (secret_message(t) - sol_receiver(t, idxs=1))*1e7
  end
    return decrypted_message
  end

  function convert_samples_to_message(decrypted_message, sampling_rate, num_of_samples, name_of_file)
    time_arr = [(1/sampling_rate)*n for n in 1:num_of_samples]
    message = decrypted_message(time_arr)
    wavwrite(message, name_of_file, Fs=sampling_rate)
  end 


function main()
u0 = [2.2,1.3,2.0]
p=[10.0;0.3333333;60.0]
tspan = (0,5.0) # length of the message is 5 seconds 

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
end

main()