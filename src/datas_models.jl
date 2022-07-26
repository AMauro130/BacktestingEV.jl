
#=

DATA CONSTRUCTOR :

=#


function constructor_data_randomN(n)
    x = [i for i in 1:n]
    lambda1,lambda2,lambda3 = 1,0,0
    #lambda1,lambda2,lambda3 = 1,0,0
    z1,z2,z3 = rand(Normal(0,3)),1,1
    data = zeros(n)
    for t in 3:n
        z3 = z2
        z2 = z1
        z1 = rand(Normal(0,3))
        data[t]=data[t-1] + lambda1.*z1 + lambda2.*z2 + lambda3.*z3
    end
    return data
end

function constructor_data_sin(n)
    x = [i for i in 1:n]
    lambda1,lambda2,lambda3 = 0.5,0.3,0.1
    #lambda1,lambda2,lambda3 = 1,0,0
    z1,z2,z3 = randn()./3,1,1
    data = zeros(n)
    data[1]=10
    for t in 2:n
        z3 = z2
        z2 = z1
        z1 = randn()
        data[t] = 5*cos(200*t/(n*pi)) + lambda1.*z1 + lambda2.*z2 + lambda3.*z3 + 10
    end
    return data
end

function constructor_data_random(n)
    a = rand(Normal(5,0.2),n)
    b = rand(Exponential(50),n)
    c = rand(Gumbel(5,50),n)
    d = rand(Frechet(5,50),n)
    e = rand(Pareto(2.5),n)

    u = a.+b.+c.+d.+d.+e
    u = u.-mean(u)
    data = zeros(1,n)
    for i in 2:n
        data[i] = data[i-1]+u[i-1]
    end
    return vec(data)
end








#=

CONSTRUCTOR OF THE MODEL :

=#

function constructor_model_real(data,n,k,window)
    n = length(data)
    coefD,coefP = 0.02,0.98
    predic_final = zeros(k,n)
    statiodata = statio_data(data)[1]
    predic_final[1] = data[1]
    for t in 2:n-statio_data(data)[2]
        predic_final[:,t] = coefD*data[t-1] .+ coefP.*predic_final[:,t-1] .+ new_dist_applier(window(statiodata,t),k,new_dist_historical)
    end
    return predic_final
end

function constructor_model_histo(data,n,k,window)
    n = length(data)
    predic_final = zeros(k,n)
    statiodata = statio_data(data)[1]
    for t in 2:n-statio_data(data)[2]
        predic_final[:,t] = data[t-1] .+ new_dist_applier(window(statiodata,t),k,new_dist_kerneldensity)
        #predic_final[:,t] = data[t-1] .+ new_dist_fitting(statiodata)
    end
    return predic_final
end

function constructor_model_simple(data,n,k,window)
    multiplier = 1
    x = [i for i in 1:n]
    lambda1,lambda2,lambda3 = 1,0.8,0.6
    #lambda1,lambda2,lambda3 = 1,0,0
    z1,z2,z3 = randn(),1,1
    nbr_of_samples = k
    new_samples = zeros(nbr_of_samples,n)
    new_samples[1]=10
    for t in 2:n
        z3 = z2
        z2 = z1
        new_samples[:,t]=[data[t-1]+multiplier*(lambda1.*randn()+lambda2.*z2+lambda3.*z3) for i in 1:nbr_of_samples]
    end
    return new_samples
end

function constructor_model_wrong(data,n,k,window)
    real_model = constructor_model_histo(data,n,k,window)
    step = 5
    n = length(data)
    predic_wrong = zeros(size(real_model))

    err = (abs(maximum(data))+abs(minimum(data)))/3

    T = Int(n/step)

    for t in 1:step-1
        predic_wrong[:,Int(t*T-0.02*n):Int(t*T+0.02*n)] .= err
    end

    return real_model .+ predic_wrong
end










#=

THE FOLLOWING IS CSV CONSTRUCTOR

=#


# UNNECESSARY
function gene_data(n,construc)
    data = construc(n)
    nms = ["h $i" for i in 1:n]
    data_array = DataFrame(transpose(data),nms)
    #file_name = string(rand(1:1:100))
    file_name = string(n)
    return CSV.write("D:\\Armand\\Astrolabium\\PROJECT\\data$file_name.csv",data_array)
end

# UNNECESSARY
function gene_samples(n,data,construc)
    new_samples = construc(data,n)
    nms = ["h $i" for i in 1:n]
    samples_array = DataFrame(new_samples,nms)
    #file_name = string(rand(1:1:100))
    file_name = string(n)
    return CSV.write("D:\\Armand\\Astrolabium\\PROJECT\\samples$file_name.csv",samples_array)
end

















































#=

THE EASY WAY TO DO AND TEST

=#




# n = 400
# k = 10000


# data_fake  = constructor_data_random(n)
#
# samples_data1 = constructor_model_histo(data_fake,n,k,moving_window)
# samples_data_fake = constructor_model_histo(data_fake,n,k,incremental_window)
# samples_data3 = constructor_model_real(data_fake,n,k,incremental_window)
# samples_data_wrong1 = constructor_model_wrong(data_fake,n,k,moving_window)
# samples_data_wrong2 = constructor_model_wrong(data_fake,n,k,incremental_window)


# plot_samples_data(data_fake,samples_data1)
# plot_samples_data(data_fake,samples_data_fake)
# plot_samples_data(data_fake,samples_data_wrong1)
# plot_samples_data(data_fake,samples_data_wrong2)
# plot_samples_data(data_fake,samples_data3)

#plot percentage of numbers of samples over the threshold
# plot(test_thre(over_threshold_total(data_fake,samples_data1,determine_threshold1),k))
# plot(test_thre(over_threshold_total(data_fake,samples_data_fake,determine_threshold1),k))

# plot_samples_data_mean(data_fake,samples_data1)
# plot_samples_data_mean(data_fake,samples_data_wrong1)
# plot_samples_data_mean(data_fake,samples_data_fake)
# plot_samples_data_mean(data_fake,samples_data_wrong2)

#=

THE REAL WAY TO DO WITH CSV FILES

|

#gene_data(500,constructor_data_random)
#data = convert_to_matrix("D:\\Armand\\Astrolabium\\PROJECT\\data10000.csv")

#gene_samples(500,data,constructor_model_histo)
#model = convert_to_matrix("D:\\Armand\\Astrolabium\\PROJECT\\samples10000.csv")


=#

# find_period_error(data_fake,samples_data1)
# find_period_error(data_fake,samples_data_fake)
# find_period_error(data_fake,samples_data3)
# find_period_error(data_fake,samples_data_wrong1)
# find_period_error(data_fake,samples_data_wrong2)
