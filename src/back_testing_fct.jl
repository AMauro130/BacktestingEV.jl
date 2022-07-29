using Distributions, Plots
using CSV
using DataFrames
using Extremes
using Optim
using ForwardDiff
using Calculus
using NLsolve
using FFTW
using StatsBase, HypothesisTests
using AdaptiveRejectionSampling
using KernelDensity
using GLM
using Interact
using Fourier
using LinearAlgebra
using CurveFit
using Decimals




# import GLM







#=

SIMPLE COMPARAISON  BETWEEN THE MEAN OF PREDICTION AND THE REAL DATA FOR EACH DAY

=#


"""
Plots the average values of the model and the real data associated.
"""
function back_real_vs_meanPred_plot(data::Vector,model::Matrix)
    return plot_samples_data_mean(data,model)
end


"""
Returns the difference between the average of a prediction set and the real data.
"""
function back_real_vs_meanpred(data::Vector,model::Matrix)
    return mean_diff(data,model)
end








#=

FOURIER TRANSFORM TO DETECT A POSSIBLE FREQUENCY OF ERROR

=#



"""
Returns the Fourier Transform of the difference between the average of a prediction set and the real data.

- Option : 1 for all values; 2 for FT plot.
"""
function back_periodicity(data::Vector,model::Matrix,option::Integer)
    freqs,F = find_period_error(data,model)
    if option == 1
        return freqs,F
    elseif option == 2
        return plot(freqs,F)
    end
end





#=

FITS GEV PARAMETERS (GEV(µ,σ,ξ))

=#


# is it helpful ?
function back_fitting_ev(data::Vector,model::Matrix,determine_threshold_algorithm,param::Integer,option::Integer)
    if option == 1
        return analyse_extreme_values(data,model,determine_threshold_algorithm,param)
    elseif option == 2
        return mean(analyse_extreme_values(data,model,determine_threshold_algorithm,param))
    elseif option == 3
        return maximum(analyse_extreme_values(data,model,determine_threshold_algorithm,param))
    elseif option == 4
        return minimum(analyse_extreme_values(data,model,determine_threshold_algorithm,param))
    end
end




"""
Returns the percentage of ev out of the whole model, depending on the chosen algorithm.

- Option : 1 for all values; 2 for the mean; 3 for the maximum; 4 for the minimum.
"""
function back_perc_ev(data::Vector,model::Matrix,determine_threshold_algorithm,option::Integer)
    if option == 1
        return perc_val_evt(data,model,determine_threshold_algorithm)
    elseif option == 2
        return mean(perc_val_evt(data,model,determine_threshold_algorithm))
    elseif option == 3
        return maximum(perc_val_evt(data,model,determine_threshold_algorithm))
    elseif option == 4
        return minimum(perc_val_evt(data,model,determine_threshold_algorithm))
    end
end






"""
Returns the EV score, depending on the chosen algorithm.
"""
function back_testingEV(data::Vector,model::Matrix,determine_threshold_algorithm)
    return ev(data,model,determine_threshold_algorithm)
end
