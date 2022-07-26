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







function back_testing_fake(data,model)
    m,n = size(model)
    if n != length(data)
        return println("size problem between your vectors")
    end

    diff = mean_diff(data,model)

    m,s,skew,kurto = moments_calculus(model)

    if minimum(kurto) >= 6          # means a long tail

    end
    #count_data_out = quantiles_test(data,model)
    res = opinions_on_the_model(data,model)
    #return count_data_out
    return res
end










#=

SIMPLE COMPARAISON  BETWEEN THE MEAN OF PREDICTION AND THE REAL DATA FOR EACH DAY

=#


# visual plot to analyse
function back_real_vs_meanPred_plot(data::Vector,model::Matrix)
    return plot_samples_data_mean(data,model)
end

# values of the difference between these values
function back_real_vs_meanpred(data::Vector,model::Matrix)
    return mean_diff(data,model)
end








#=

FOURIER TRANSFORM TO DETECT A POSSIBLE FREQUENCY OF ERROR

=#


# option 1 : values ; option 2 : FT plot
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

THE PERCENTAGE OF EV OUT OF THE WHOLE MODEL

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





#=

EV SCORE ACCORDING TO THE FORMULA

=#


# !!! most powerfull tool !!!
function back_testing(data::Vector,model::Matrix,determine_threshold_algorithm)
    return ev(data,model,determine_threshold_algorithm)
end
