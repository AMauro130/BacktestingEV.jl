



#=

THOSE FUNCTIONS ARE CALCULUS FUNCTIONS FOR CRITERIAS

=#


#this function compare the calculate mean with the original one (based on the data)
function mean_diff(data,model)
    meanS = mean_samples(model)
    diff = absolute_error(data,meanS)
    return diff
end


#here is the fourier transform for a data (the idea is to find a pattern of error)
function find_period_error(data,model)
    n = length(data)
    t0 = 1           # start time
    fs = 1          # sampling rate
    # fs = 44100
    tmax = n          # end time
    t = t0:1/fs:tmax;
    signal = mean_diff(data,model)
    F = fftshift(fft(signal))
    freqs = fftshift(fftfreq(length(t), fs))
    p1 = plot(t,signal)
    p2 = plot(freqs,abs.(F))
    # return plot(p1,p2,layout=(2,1))
    return freqs,abs.(F)
end









#calculus of the fourth firsts moments
function moments_calculus(model)
    a = applier(mean,model)
    b = applier(var,model)
    c = applier(skewness,model)
    d = applier(kurtosis,model)
    # p1,p2,p3,p4 = plot(a),plot(b),plot(c),plot(d)
    # return plot(p1,p2,p3,p4,layout=(2,2),legend=false)
    return a,b,c,d
end












#=

EXTRACT THE EXTREME VALUES FROM A DATA SERIE
(to be finished)

=#

#you need to provide something verifying is_part_ev
function extractor(model_col)
    ev_thre, = determine_threshold2(model_col)
    return over_threshold_total(model_col,ev_thre)
end













#=

THOSE FUNCTIONS ARE "IS ..." WHICH MEANS THE RETURNS ARE TRUE OR FALSE

=#

#is the model accurate ?
function is_accurate(data,model)
    diff = mean_test(data,model)
    l = maximum(data) - minimum(data)
    thres = 0.15*l
    m = mean(diff)
    return m <= thres
end

#is there a pattern of the error ?
function is_pattern_error(data,model)
    transf_f = find_period_error(data,model)
    transf_f[1] = 0
    l = filter(x -> x>=0.2*maximum(transf_f),transf_f)
    return length(l) != 1
end

#is part of the ev ? (must be stationary first)
function is_part_ev(model_col)
    return kurtosis(model_col) > 6
end













#=

THOSE FUNCTIONS RETURNS PERCENTAGE OF A CRITERIA

=#


#returns the percentage of the mean error between the mean prediction and the original data
function percent_error_model(data,model)
    diff = mean_diff(data,model)
    perc = abs.((diff./data).*100)
    #if no error then returns false
    perc = replace(perc, NaN => 0)
    return mean(perc)
end


#this function returns the percentage of values inside of the quantile range
function quantiles_test(data,model)
    count_data_out = zeros(1,4)
    percent = zeros(2,4)
    percent[1,:]=[0.05,0.025,0.01,0.005]
    for i in 1:length(data)
        col = model[:,i]
        quan = [quantile(col,[percent[1,:][j],1-percent[1,:][j]]) for j in 1:length(percent[1,:])]
        for t in 1:length(quan)
            if (data[i]<quan[t][1])||(data[i]>quan[t][2])
                count_data_out[t]+=1
            end
        end
        #count_quantile[i] = count(j->((c[j]<a)||c[j]>b)),c)
    end
    res = (count_data_out/length(data))*100
    return 100 - sum(res)/length(res)
end













#=

REAL BACKTESTING FUNCTIONS TO COMPARE VALUES FROM THE MODEL WITH THRE REAL DATA

=#



function values_into_ev(data_point,thre)
    if data_point >= thre
        return 1
    else
        return 0
    end
end


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

# ...this one is less accurate, because of the criteria 330...
function analyse_extreme_values_dist(data,model,determine_threshold)
    extreme_dist = over_threshold_total(data,model,determine_threshold)
    res = []
    for t in 1:length(data)
        if length(extreme_dist[t])>=330
            append!(res,fitting_gev(extreme_dist[t]))
        end
        # append!(res,fitting_gev(extreme_dist[t]))
    end
    return res
end


# !!! USE THIS ONE !!!
# param 1, 2, 3 respectively for μ, σ and ξ
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

OPINIONS ON THE MODEL

=#


#this is the main opinion function
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
