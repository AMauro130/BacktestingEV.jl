

#=

CONVERT FROM CSV FILE TO MATRIX

=#

#convert csv file to matrix
function convert_to_matrix(path)
    new_file = CSV.read(path,DataFrame)
    return Matrix(new_file) #it converts to a matrix form
end













#=

DIFFERENCE BETWEEN 2 COMPARABLE VECTORS. EX : DATA AND MEAN_SAMPLES

=#

function absolute_error(the,exp)
    return abs.(exp-the)
end

function relative_error(the,exp)
    (exp-the)./abs.(the)
end







#=

PLOT DATA AND THE RELATED SAMPLES

=#

#plot the data and the whole predictions for each day
function plot_samples_data(data,model)
    pdata = plot(data)
    p = plot(data)
    for i in 1:size(model)[1]
        plot!(model[i,:])
    end
    return plot(p,pdata,layout=(2,1),legend = false)
end

#plot the difference between the mean for each data and the real data
function plot_samples_data_mean(data,model)
    p = plot(data)
    meanD = mean_samples(model)
    plot!(vec(meanD))
    return plot(p)
end







#=

EXTREME VALUES THEORY

=#

#percent ev values out of total values
function perc_evt(dist,determine_threshold)
    n = length(dist)
    a,b = determine_threshold(dist)
    c = n-b
    return (n-c)/n*100
end


# !!! TO CHOOSE BETWEEN THIS ONE AND perc_val_evt !!!

#this one is longer than perc_val_evt
function perc_val_evt2(data,model,determine_threshold)
    n = length(data)
    res = []
    for i in 1:n
        append!(res,perc_evt(model[:,i],determine_threshold))
    end
    return res
end


#extract ev from the whole dist
function over_threshold(dist,determine_threshold)
    a,b = determine_threshold(dist)
    res = filter(x -> x>=a,dist)
    return res
end


# !!! TO CHOOSE BETWEEN THIS ONE AND perc_val_evt2 !!!
# returns the percentage of ev for each day
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


#IT RETURNS THE BABY DISTRIBUTION FOR EACH DAY
function over_threshold_total(data,model,determine_threshold)
    m,n = size(model)
    res = []
    for t in 1:n
        append!(res, [over_threshold(model[:,t],determine_threshold)])
    end
    return res
end


#test if over_threshold_total is efficient
function test_thre(a,k)
    return [(length(a[i])./k).*100 for i in 1:length(a)]
end














#=

FOR EXTREME VALUES ANALYSIS

=#

#for just one column of a child dist (for predictions of a data) after threshold of evt
# param 1, 2, 3 respectively for μ, σ and ξ
function fitting_gev(child_dist,param)
    return gevfit(child_dist).θ̂[param]
end









#=

TYPES OF WINDOWS FOR THE MODEL

=#

# incremental window
function incremental_window(data,tn)
    return data[1:tn]
end

# moving window
function moving_window(data,tn)
    wsize = 100
    if tn-wsize < 1
        return data[1:tn]
    else
        return data[tn-wsize:tn]
    end
end










#=

FIND THE STATIONARY STATE FOR THE DATA

=#

# data[i+1] - data[i] !!!Not the same size!!! length() -= 1
function diff_values(data)
    return [data[i+1]-data[i] for i in 1:length(data)-1]
end

# IS THE DATA STATIONARY ?
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

#you can apply these 2 functions to a stationary data

#fit problem due to cross packages...
#non parametric
function new_dist_historical(data,k)
    h = GLM.fit(Histogram,data,nbins = k/10)
    val = h.edges[1]
    val = [i for i in val]
    pop!(val)
    proba = h.weights
    w = Weights(proba)
    return val,w
end


#takes a stationary data and Kernel density estimation
function new_dist_kerneldensity(data,k)
    kerneldist = kde(data)
    val = [i for i in kerneldist.x]
    w = Weights(kerneldist.density)

    return val,w
end


#from prob and weights
function new_dist_applier(data,k,method_dist)
    predic = zeros(k,1)
    val,w = method_dist(data,k)
    for t in 1:k
        predic[t] = sample(val,w)
    end
    return predic
end


#test some distribution and find the most convenient
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
        d1 = fit_mle(Normal,statiodata_test)
        d2 = fit_mle(Laplace,statiodata_test)
        d3 = fit_mle(Rayleigh,statiodata_test)
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


function mean_samples(model)
    n = size(model)[2]
    meanD = zeros(1,n)
    for i in 1:n
        meanD[i] = mean(model[:,i])
    end
    return vec(meanD)
end


















#=

CALCULUS OF MOMENTS FOR A DISTRIBUTION

=#


#you can apply a function (like var) to a m*n size list TO COLULN !
function applier(f,list)
    m,n = size(list)
    res = [f(list[:,i]) for i in 1:n]
    return res
end
