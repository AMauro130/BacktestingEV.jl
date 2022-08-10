#=

THIS FILES CONTAINS THE BACK-TESTING CRITERIAS

=#


"""
Returns the difference between the average of a prediction set and the real data.
"""
function mean_diff(data,model)
    meanS = mean_samples(model)
    diff = absolute_error(data,meanS)
    return diff
end


"""
Returns the Fourier Transform of the difference between the average of a prediction set and the real data.

- It helps you to find a potential error periodicity in your model.
"""
function find_period_error(data,model)
    n = length(data)
    t0 = 1
    fs = 1
    tmax = n
    t = t0:1/fs:tmax;
    signal = mean_diff(data,model)
    F = fftshift(fft(signal))
    freqs = fftshift(fftfreq(length(t), fs))
    p1 = plot(t,signal)
    p2 = plot(freqs,abs.(F))
    return freqs,abs.(F)
end


"""
Returns the first four moments.
"""
function moments_calculus(model)
    a = applier(mean,model)
    b = applier(var,model)
    c = applier(skewness,model)
    d = applier(kurtosis,model)
    return a,b,c,d
end





#=

BACKTESTING FUNCTIONS TO COMPARE VALUES FROM THE MODEL WITH THRE REAL DATA

=#


"""
Returns the values which are part of the ev.
"""
function values_into_ev(data_point,thre)
    if data_point >= thre
        return 1
    else
        return 0
    end
end


"""
Return the EV score.

Remark :
- Needs to be applied to a stationary data (like your initial data) with the selected method
- k is the number of the elements in your future prediction set
"""
function ev(data,model,determine_threshold)
    n = length(data)
    N, = size(model)
    EV = 0
    for i in 1:n
        value_quant,q = determine_threshold(model[:,i])
        # v = quantile(model[:,i],0.99)
        inside = values_into_ev(data[i],value_quant)
        tail_ratio = q/N
        # tail_ratio = 0.99
        res = inside/(tail_ratio)
        EV += res
        # EV += inside
    end
    return EV/n
    # return EV
end





#=

EXTREME VALUES THEORY AND ANALYSE

=#


"""
Returns the parameters of Generalized Extreme Value (GEV) distribution.

- Param : 1 for μ; 2 for σ and 3 for ξ
"""
function analyse_extreme_values(data,model,determine_threshold,param)
    n = length(data)
    fi = 0
    res = []
    for i in 1:n
        extreme_values = over_threshold(model[:,i],determine_threshold)
        try
            fi = fitting_gev(extreme_values,param)
        catch err
            fi = 0
        end
        append!(res,fi)
    end
    return res
end





#=

WHAT ABOUT THIS ?? THE ASSOCIATED FILES IN NON_USE.JL

=#


"""
Returns some opinions on the model
"""
function opinions_on_the_model(data,model)
    res = []
    if is_accurate(data,model)
        append!(res,["precise"])
    else
        append!(res,["not precise"])
    end

    if is_pattern_error(data,model)
        append!(res,["pattern error found"])
    else
        append!(res,["no pattern error found"])
    end

    return res
end
