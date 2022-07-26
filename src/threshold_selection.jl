
#=

HILL ESTIMATOR APPLIED TO THE 10% EV

=#


function hill_estimator(dist,k)
    n = length(dist)
    res = sum([log(dist[n-i+1])-log(dist[n-k]) for i in 1:k])./k
    return res
end

# !!! do not apply this function to the whole data !!!
function apply_hill(dist)
    dist = sort(dist)
    n = length(dist)
    k1 = Int(round(0.1*n))
    # k1 = n
    k2 = Int(round(0.9*n))
    dist = dist[k2:end]
    res = []
    for i in 1:k1-1
        append!(res,hill_estimator(dist,i))
    end
    # return scatter(res)
    return res
end


#useful for the magic number
function over_stop(deriv,thres)
    deriv = replace(deriv, NaN => 0)     #if error
    n = length(deriv)
    deriv = abs.(deriv)
    i = 1
    while (deriv[i]<thres)&&(i<n)
        i+=1
    end
    return i
end



#=

NORMALIZATION AND DENORMALIZATION

(to work between 0 and 1, no negative values)

=#


function normalization(dist)
    min = minimum(dist)
    max = maximum(dist)
    len = abs(max-min)
    res = (dist .- min) ./ len
    return res,max,min,len
end

function denormalization(dist,max,min,len)
    res = dist .* len .+ min
    return res
end
















"""

DETERMINE THE THRESHOLD OF EV 1

(based on the area under the curve, the difference between the regression and the real estipmator)
(used with the hill estimator)

"""
function reg_new_left(dist)
    y = dist
    x = 1:length(dist)
    f1 = curve_fit(Polynomial,x,y,6)
    y1 = f1.(x)
    p1 = plot(x, y)
    return y1
end

function reg_new_right(dist)
    y = dist
    x = 1:length(dist)
    f1 = curve_fit(Polynomial,x,y,1)
    y1 = f1.(x)
    p1 = plot(x, y)
    return y1
end

function err(dist1,dist2,reg1,reg2)
    e1 = relative_error(dist1,reg1)
    e2 = relative_error(dist2,reg2)
    m1 = sum(e1)
    m2 = sum(e2)
    return abs(m1),abs(m2)
end

function find_threshold_for_real(hill)
    n = length(hill)
    i = 0
    border = Int(round(n/2))
    res = [border]
    hill1 = hill[1:border-1]
    hill2 = hill[border:end]
    reg1 = reg_new_left(hill1)
    reg2 = reg_new_right(hill2)
    m1,m2 = err(hill1,hill2,reg1,reg2)
    while ((m1<0.5*m2)||(m1>1.5*m2))&&(i<=1000)
    # while (m1<0.1*m2)||(m1>1.9*m2)
        if m1<0.5*m2
            border = border + Int(round(border/(2+i)))
        else
            border = border - Int(round(border/(2+i)))
        end
        hill1 = hill[1:border-1]
        hill2 = hill[border:end]
        reg1 = reg_new_left(hill1)
        reg2 = reg_new_right(hill2)
        m1,m2 = err(hill1,hill2,reg1,reg2)
        i+=1
        append!(res,border)
    end
    return border
    # return res
end

function determine_threshold1(dist)
    n = length(dist)
    dist_r,max,min,len = normalization(dist)
    hill = apply_hill(dist_r)
    hill = convert.(Float64,hill)
    res = find_threshold_for_real(hill)
    di = sort(dist_r)
    di = reverse(di)
    d_res = di[res]
    s = denormalization(d_res,max,min,len)
    return s,res
end






#=

DETERMINE THE THRESHOLD OF EV 2

(based on the regression of the derivative of the hill estimator)

=#

function line_reg(x,y)
    data = DataFrame(X=x, Y=y)
    res = lm(@formula(Y ~ X), data)
    # lm1 = fit(LinearModel, @formula(Y ~ X), data)
    # return res.model.pp.beta0
    # return plot(res.model.rr.mu)
    return res
end


function creator_function_regression(x,y)
    res = line_reg(x,y)                        # linear regression
    meanA,meanB = coeftable(res).cols[1,:][1]
    maxA,maxB = coeftable(res).cols[6,:][1]
    minA,minB = coeftable(res).cols[5,:][1]
    f(u) = meanB .* u .+ meanA
    f1(u) = maxB .* u .+ maxA
    f2(u) = minB .* u .+ minA
    # vUpp = filter(i -> y[i]>f1(i),x)
    # vLow = filter(i -> y[i]<f2(i),x)
    valUpLow = filter(i -> (f1(i)<y[i])||(y[i]<f2(i)),x)
    scatter(x,y)
    plot!(x,f)
    plot!(x,f1)
    # return vUpp,vLow
    return valUpLow
end

function is_ev_test(ul)
    ul_bis = diff_values(ul)
    n = length(ul_bis)
    # thre = Int(round(0.01 * n))
    thre = 5
    i = 1
    s = 0
    while (s<thre)&&(i<n)
        if ul_bis[i]<=2
            s += 1
        else
            s = 0
        end
        i += 1
    end
    return ul[i]
end

function determine_threshold2(dist)
    dist_r,max,min,len = normalization(dist)
    hill = apply_hill(dist_r)              # Hill estimator on the distribution
    deriva = diff_values(hill)             # Derivative of the estimator for detect the linear moment
    deriva = reverse(deriva)               # To use over_stop with 0.004
    n = length(deriva)
    hill = convert.(Float64,hill)
    x = [i for i in 1:length(deriva)]
    # x = [i for i in 1:length(hill)]
    ul = creator_function_regression(x,deriva)
    # ul = creator_function_regression(x,hill)
    if length(ul) == 0
        res = 0
    else
        res = is_ev_test(ul)
    end
    res = n-res                            # Number of extreme values starting from the end
    di = sort(dist_r)                       # Sorts for extract the first ev
    di = reverse(di)                         # Reverse to start with the highest values
    d_res = di[res]
    s = denormalization(d_res,max,min,len)
    return s,res
end







#=

DETERMINE THE THRESHOLD OF EV 3

(based on my magic number)

=#


function determine_threshold3(dist)
    dist_r,max,min,len = normalization(dist)
    hill = apply_hill(dist_r)              # Hill estimator on the distribution
    deriva = diff_values(hill)             # Derivative of the estimator for detect the linear moment
    deriva = reverse(deriva)               # To use over_stop with 0.004
    n = length(deriva)
    hill = convert.(Float64,hill)
    res = over_stop(deriva,0.004)        # 0.004 empirical value for detection
    res = n-res                            # Number of extreme values starting from the end
    di = sort(dist_r)                       # Sorts for extract the first ev
    di = reverse(di)
    if res == 0
        res = 1
    end                         # Reverse to start with the highest values
    d_res = di[res]
    s = denormalization(d_res,max,min,len)
    return s,res
end






#=

DETERMINE THE THRESHOLD OF EV 4

(based on deciamls, counting consecutives values)
(Algorithm 5 pdf)

=#

function param_decim(samples)
    n = length(samples)
    param = decimal.(samples)
    decim = [param[i].c for i in 1:n]
    nb_decim = abs.([param[i].q for i in 1:n])
    return decim,nb_decim
end


function max_j(T)
    n = length(T)
    decim,nb_decim = param_decim(T)
    max_nb_decim = maximum(nb_decim)
    while length(unique(T)) == n
        max_nb_decim -= 1
        T = round.(T,digits = max_nb_decim)
    end
    return max_nb_decim+1
end


function akj(T,j)
    return round.(T,digits=j)
end


function consec(T,j)
    diff1 = akj(diff_values(akj(T,j)),j)
    res_int = filter(i->diff1[i]==0,1:length(diff1))

    while length(res_int) == 0
        j -= 1
        diff1 = akj(diff_values(akj(T,j)),j)
        res_int = filter(i->diff1[i]==0,1:length(diff1))
    end
    return j
end

#first 2 consecituves elements
function thre_from_consec(T,j)
    diff1 = akj(diff_values(akj(T,j-2)),j-2)
    res_int = filter(i->diff1[i]==0,1:length(diff1))
    return res_int
end

function consecutive_values(res_int)
    n = length(res_int)
    res = []
    res_mem = []
    if n==0
        res_mem = [1]
        return res_mem[1]
    end

    i = 1
    while i<n
        if (res_int[i+1]-res_int[i]) == 1
            append!(res,res_int[i])
        else
            if length(res)>length(res_mem)
                res_mem = res
            end
            res = []
        end
        i+=1
    end

    # means all the values in res_int are consecutives, so res is res_int
    # means that there are not 3 consecutives values, so we take the first 2 consecutives values
    if length(res_mem)==0
        res_mem = res_int[1]
    end
    return res_mem[1]
end

function determine_threshold4(dist)
    n = length(dist)
    dist_r,max,min,len = normalization(dist)
    dist_r = replace(dist_r, NaN => 0)
    T =  apply_hill(dist_r)
    T = replace(T,NaN => 0)
    j0 = max_j(T)
    j = consec(T,j0)
    res_int = thre_from_consec(T,j)
    res = consecutive_values(res_int)
    d = sort(dist_r)
    r = n-res
    s = denormalization(d[r],max,min,len)
    return s,res
end



#=

DETERMINE THE THRESHOLD OF EV 5

(based on the previous algorithm but applied to the regression)

=#

function hill_estimator_reg(T)
    t = convert.(Float64,T)
    x = [i for i in 1:length(t)]
    res = line_reg(x,t)
    a,b = res.model.pp.beta0
    f(u) =b.*u.+a
    return [f(i) for i in x]
end

function T2_creator(T)
    T2 = hill_estimator_reg(T)
    return abs.(T-T2)
end

function determine_threshold5(dist)
    n = length(dist)
    dist_r,max,min,len = normalization(dist)
    dist_r = replace(dist_r, NaN => 0)
    T = apply_hill(dist_r)
    T = replace(T,NaN => 0)
    T2 = T2_creator(T)
    j0 = max_j(T2)
    j = consec(T2,j0)
    res_int = thre_from_consec(T2,j)
    res = consecutive_values(res_int)
    d = sort(dist_r)
    r = n-res
    s = denormalization(d[r],max,min,len)
    return s,res
    # return res
end






#=

DETERMINE THE THRESHOLD OF EV 6

(based on hill and beta estimator)
(Algorithm 1 pdf)

=#

function M(samples,k,j)
    n = length(samples)
    res = sum([(log(samples[n-i+1])-log(samples[n-k]))^j for i in 1:k])./k
    return res
end

function W(samples,k,tau)
    W1 = ((M(samples,k,1)^tau)-(M(samples,k,2)/2)^(tau/2))/((M(samples,k,2)/2)^(tau/2)-(M(samples,k,3)/6)^(tau/3))
    W2 = (log(M(samples,k,1))-log(M(samples,k,2)/2)/2)/((log(M(samples,k,2)/2)/2)-log(M(samples,k,3)/6)/3)
    if tau != 0
        return W1
    else
        return W2
    end
end

function rho_estimator(samples,k,tau)
    samples = sort(samples)
    res = -abs( 3*(W(samples,k,tau)-1)/(W(samples,k,tau)-3) )
    return res
end

function U(samples,i)
    n = length(samples)
    res = i*(log(samples[n-i+1])-log(samples[n-i]))
    return res
end

function d(k,alpha)
    res = sum([(i/k)^(-alpha) for i in 1:k])./k
    return res
end

function D(samples,k,alpha)
    res = sum([((i/k)^(-alpha))*U(samples,i) for i in 1:k])./k
    return res
end

function beta_estimator(samples,k,tau)
    n = length(samples)
    samples = sort(samples)
    rho_est = rho_estimator(samples,k,tau)
    # !!! overflow pb binomial !!!
    bi = (k/n)^rho_est * (d(k,rho_est)*D(samples,k,0)-D(samples,k,rho_est))/(d(k,rho_est)*D(samples,k,rho_est)-D(samples,k,2*rho_est))
    return bi
end

function step2(samples)
    n = length(samples)
    k1 = Int(round(n^0.995))
    k2 = Int(round(n^0.999))
    if k2==n
        k2 = n-1
    end
    k = [i for i in k1:k2]
    rho_res_K_0 = []
    rho_res_K_1 = []
    for i in k
        append!(rho_res_K_0,rho_estimator(samples,i,0))
        append!(rho_res_K_1,rho_estimator(samples,i,1))
    end
    ksi_0 = median(rho_res_K_0)
    ksi_1 = median(rho_res_K_1)

    I_0 = sum([(rho_estimator(samples,i,0)-ksi_0)^2 for i in k])
    I_1 = sum([(rho_estimator(samples,i,1)-ksi_1)^2 for i in k])

    if I_0 <= I_1
        tau = 0
    else
        tau = 1
    end
    return tau
end


function step3(samples)
    n = length(samples)
    k01 = Int(round(n^0.999))
    if k01==n
        k01 = n-1
    end
    tau = step2(samples)
    beta = beta_estimator(samples,k01,tau)
    rho = rho_estimator(samples,k01,tau)
    return beta,rho
end


function step4(samples)
    n = length(samples)
    beta,rho = step3(samples)
    k0h = (((1-rho)^2)*(n^(-2*rho)))/(-2*rho*beta^2)
    k0h = k0h^(1/(1-2*rho))
    if (k0h > 10000)||isnan(k0h)
        return 1
    else
        k0h = Int(floor(k0h))
        return k0h
    end
end

function step5(dist)
    n = length(dist)
    dist_r,max,min,len = normalization(dist)
    k1 = Int(round(0.1*n))
    k2 = Int(round(0.9*n))
    samples = sort(dist_r)
    samples = samples[k2:end]
    res = step4(samples)
    # return scatter(res)
    di = sort(dist_r)                       # Sorts for extract the first ev
    di = reverse(di)                         # Reverse to start with the highest values
    d_res = di[res]
    s = denormalization(d_res,max,min,len)
    return s,res
end

function determine_threshold6(dist)
    return step5(dist)
end
