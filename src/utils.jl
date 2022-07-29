

"""
Convert CSV file to matrix.
"""
function convert_to_matrix(path)
    new_file = CSV.read(path,DataFrame)
    return Matrix(new_file) #it converts to a matrix form
end













"""
Returns the absolute error.

- "the" stands for Theorical and "exp" for Experimental.
"""
function absolute_error(the,exp)
    return abs.(exp-the)
end



"""
Returns the relative error.

- "the" stands for Theorical and "exp" for Experimental.
"""
function relative_error(the,exp)
    (exp-the)./abs.(the)
end







#=

PLOT DATA AND THE RELATED SAMPLES

=#

"""
Plots the real data and the whole predictions for each day.
"""
function plot_samples_data(data,model)
    pdata = plot(data)
    p = plot(data)
    for i in 1:size(model)[1]
        plot!(model[i,:])
    end
    return plot(p,pdata,layout=(2,1),legend = false)
end


"""
Plots the difference between the real data the average prediction for each day.
"""
function plot_samples_data_mean(data,model)
    p = plot(data)
    meanD = mean_samples(model)
    plot!(vec(meanD))
    return plot(p)
end







#=

EXTREME VALUES THEORY

=#



"""
Extracts ev from the whole distribution.
"""
function over_threshold(dist,determine_threshold)
    a,b = determine_threshold(dist)
    res = filter(x -> x>=a,dist)
    return res
end


"""
Returns the percentage of ev for each sample.
"""
function perc_val_evt(data,model,determine_threshold)
    n = length(data)
    res = []
    c = 0
    for i in 1:n
        extreme_values = over_threshold(model[:,i],determine_threshold)
        c = (length(extreme_values)/length(model[:,i]))*100
        append!(res,c)
    end
    return res
end



"""
Returns the ev distribution (only with the ev) for each sample.
"""
function over_threshold_total(data,model,determine_threshold)
    m,n = size(model)
    res = []
    for t in 1:n
        append!(res, [over_threshold(model[:,t],determine_threshold)])
    end
    return res
end



# ????

#test if over_threshold_total is efficient
# ?? TO USE MAYBE ???
function test_thre(a,k)
    return [(length(a[i])./k).*100 for i in 1:length(a)]
end














#=

FOR EXTREME VALUES ANALYSIS

=#


"""
Returns the parameters of the GEV (Generalized Extreme Value) distribution μ, σ and ξ.

Parameter selection :
    choose param = 1 for μ ; 2 for σ ; 3 for ξ.

- Needs to be used on an ev distribution.
"""
function fitting_gev(ev_dist,param)
    return gevfit(ev_dist).θ̂[param]
end










#=

TYPES OF WINDOWS FOR THE MODEL

=#


"""
Returns the data from the begining to the limit "lim".
"""
function incremental_window(data,lim)
    return data[1:lim]
end


"""
Returns the data with a precise size "wsize" until "lim".
"""
function moving_window(data,lim,wsize)
    if lim-wsize < 1
        return data[1:lim]
    else
        return data[lim-wsize:lim]
    end
end










#=

FIND THE STATIONARY STATE FOR THE DATA

=#


"""
Returns the difference between elelments.

- The length is length(data) - 1
- The return data is not necessarily stationary.
"""
function diff_values(data)
    return [data[i+1]-data[i] for i in 1:length(data)-1]
end


"""
Returns a stationay data.

- This method is based the Augmented Dickey-Fuller test (ADFTest).
"""
function statio_data(data)
    statio = diff_values(data)
    shift = 1
    test_statio = ADFTest(statio,:trend,1)
    while pvalue(test_statio) > 0.01
        statio = diff_values(statio)
        shift += 1
        test_statio = ADFTest(statio,:trend,1)
    end
    return statio, shift
end












#=

CREATE NEW DISTRIBUTIONS BASED ON THE DATA

=#


"""
Returns probabilities and weights based on the histogram fitting.

- Needs to be applied to a stationary data (like your initial data).
- k is the number of the elements in your future prediction set.
"""
function new_dist_historical(data,k)
    h = GLM.fit(Histogram,data,nbins = k/10)
    val = h.edges[1]
    val = [i for i in val]
    pop!(val)
    proba = h.weights
    w = Weights(proba)
    return val,w
end



"""
Returns probabilities and weights based on the Kernel density estimation.

- Needs to be applied to a stationary data (like your initial data).
- k is the number of the elements in your future prediction set.
"""
function new_dist_kerneldensity(data,k)
    kerneldist = kde(data)
    val = [i for i in kerneldist.x]
    w = Weights(kerneldist.density)

    return val,w
end



"""
Returns the distribution from a given methods of probabilities and weights selection.

- Needs to be applied to a stationary data (like your initial data) with the selected method.
- k is the number of the elements in your future prediction set.
"""
function new_dist_applier(data,k,method_dist)
    predic = zeros(k,1)
    val,w = method_dist(data,k)
    for t in 1:k
        predic[t] = sample(val,w)
    end
    return predic
end



"""
Returns the most appropriate distribution comparing with Nomal, Laplace and Rayleigh distributions.

- Needs to be applied to a stationary data (like your initial data).
- k is the number of the elements in your future prediction set.
"""
function new_dist_fitting(data,k)
    nbr = k
    n = length(data)
    step_ = Int(round(n/10))
    predic = zeros(nbr,1)
    res1 = 0
    res2 = 0
    res3 = 0
    for i in 1:8
        statiodata_test = data[1:step_*(i+1)-1]
        statiodata_valid = data[step_*(i+1):end]
        h = GLM.fit(Histogram,statiodata_valid,nbins = Int(round(length(statiodata_valid)/20)))
        #it creates the new distribution
        d1 = Distributions.fit_mle(Normal,statiodata_test)
        d2 = Distributions.fit_mle(Laplace,statiodata_test)
        d3 = Distributions.fit_mle(Rayleigh,statiodata_test)
        #d4 = fit_mle(Exponential,statiodata_test) #if z > 0 for each value
        #then
        x = [i for i in h.edges[1]]
        pop!(x)
        c1 = pdf.(d1,x)
        c2 = pdf.(d2,x)
        c3 = pdf.(d3,x)
        w = h.weights ./ length(statiodata_valid)
        res1 += norm(w.-c1)
        res2 += norm(w.-c2)
        res3 += norm(w.-c3)
    end
    if (res1 <= res2)&&(res1 <= res3)
            predic .= rand(d1,nbr)
    elseif (res2 <= res3)
        predic .= rand(d2,nbr)
    else
        predic .= rand(d3,nbr)
    end
    return predic
end
























#=

TRANSFORM THE model TO MEAN(model)

=#


"""
Returns the average value for each prediction set.
"""
function mean_samples(model)
    n = LinearAlgebra.size(model)[2]
    meanD = zeros(1,n)
    for i in 1:n
        meanD[i] = StatsBase.mean(model[:,i])
    end
    return LinearAlgebra.vec(meanD)
end


















#=

CALCULUS OF MOMENTS FOR A DISTRIBUTION

=#


"""
Applies a function f to an array m*n (like your model).

This function is mainly used for calculating the first four moments.
"""
function applier(f,model)
    m,n = LinearAlgebra.size(model)
    res = [f(list[:,i]) for i in 1:n]
    return res
end
