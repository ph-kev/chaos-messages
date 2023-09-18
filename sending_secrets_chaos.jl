using DifferentialEquations
using Plots
using WAV
using Interpolations
using BenchmarkTools
using LaTeXStrings

"""
    convert_message_to_samples(message::String)

Convert a message (.wav file) which is identified as the string of the file path 
and convert it into samples. 
"""
function convert_message_to_samples(message::String)
    sample, sampling_rate = wavread(message)
    sample = vec(sample)
    time_arr = Float64.([(1/sampling_rate)*n for n in 1:length(sample)])
    message_unencrypted = linear_interpolation(time_arr, sample, extrapolation_bc = Line())
    return message_unencrypted
  end

"""
    convert_message_to_samples_noisy(message::String)

Convert a message (.wav file) which is identified as the string of the file path 
and convert it into samples with added noise. The noise is generated from a Gaussian
distribution with mean 0 and standard deviation std. 
"""
function convert_message_to_samples_noisy(message::String, std)
  sample, sampling_rate = wavread(message)
  sample = vec(sample)

  # add noise 
  noise = rand(Normal(0.0, std), length(sample))
  sample += noise 

  time_arr = Float64.([(1/sampling_rate)*n for n in 1:length(sample)])
  message_unencrypted = linear_interpolation(time_arr, sample, extrapolation_bc = Line())
  return message_unencrypted
end

"""
    create_secret_message(message::String)

# Arguments
- `u0`: initial condition of the variation of the Lorenz system 
- `p`: parameters of the variation of the Lorenz system 
- `tspan` : length of the message 
- `message_unencrypted` : function that produces the message that is to be encrpyted

Encrpyt the message using a variation of the Lorenz system (see definition of 
parameterized_lorenz_talker! in this function).  
"""
function create_secret_message(u0, p, tspan, message_unencrypted)
  function parameterized_lorenz_talker!(du,u,p,t)
    x_T,y_T,z_T = u
    σ,β,r = p
    du[1] = σ*(y_T-x_T)
    du[2] = r*x_T - y_T - 20*(x_T * z_T)
    du[3] = 5*x_T*y_T - β*z_T
  end
  prob = ODEProblem(parameterized_lorenz_talker!,u0,tspan,p)
  sol_talker = solve(prob, AutoTsit5(Rodas4P()), abstol = 1e-11, reltol = 1e-11)

  function secret_message(t)
    hidden = 1e-5*message_unencrypted(t)
    secret_at_time_t = sol_talker(t,idxs=1) + hidden
    return secret_at_time_t
  end

  return secret_message
end 

"""
    decrypt_secret_message(u0, p, tspan, secret_message)

# Arguments
- `u0`: initial condition of the variation of the Lorenz system 
- `p`: parameters of the variation of the Lorenz system 
- `tspan` : length of the message 
- `secret_message` : function that produces the secret message that is to be decrypted 

Decrypt the message using a variation of the Lorenz system (see definition of 
parameterized_lorenz_receiver! in this function).  
"""
function decrypt_secret_message(u0, p, tspan, secret_message)
  function parameterized_lorenz_receiver!(du,u_1,p,t)
    x_R, y_R, z_R = u_1
    σ,β,r = p
    du[1] = σ*(y_R-x_R)
    du[2] = r*secret_message(t) - y_R - 20*(secret_message(t) * z_R)
    du[3] = 5*(secret_message(t)*y_R) - β*z_R
  end

  prob_C = ODEProblem(parameterized_lorenz_receiver!,u0,tspan,p)
  sol_receiver = solve(prob_C, AutoTsit5(Rodas4P()), abstol = 1e-11, reltol = 1e-11)
  
  function decrypted_message(t)
    return (secret_message(t) - sol_receiver(t, idxs=1))*1e5
  end
    return decrypted_message
  end

"""
    convert_samples_to_message(decrypted_message, sampling_rate, num_of_samples, name_of_file)

# Arguments
- `decrypted_message`: function that produce the decrypted message 
- `sampling_rate`: sampling rate of the resulting file 
- `num_of_samples` : number of samples of the resulting file 
- `name_of_file` : name of the resulting file 

Convert samples produced by the decrypted_message into a wav file. 
"""
function convert_samples_to_message(decrypted_message, sampling_rate, num_of_samples, name_of_file)
  time_arr = [(1/sampling_rate)*n for n in 1:num_of_samples]
  message = decrypted_message(time_arr)
  wavwrite(message, name_of_file, Fs=sampling_rate)
end 

"""
    error_set_up(message_unencrypted, decrypted_message)

Find the error between message_unencrypted and decrypted_message. 
"""
function error_set_up(message_unencrypted, decrypted_message)
  function abs_error(t)
    val1 = message_unencrypted(t)
    val2 = decrypted_message(t)
    return abs(val1-val2)
  end 
  return abs_error
end