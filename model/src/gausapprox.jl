using Distributions, RCall


#simulated mixture
m=0.1;

tau = 1;
theta1 = 0;
thetadiff = 5;
theta2 = theta1+thetadiff;
rmax = 2;
tau = 1;
sigmaE=0;
sigmaG=1;
beta = 0.001;
Z=0.5;

N1 = 10;
N2 = 10;

tmax = 100;


rfunc = function(mu,theta)
  r = rmax*tau*(exp(-((theta-mu)^2)/(2*((sigmaE^2 + sigmaG^2)+tau^2))))/(sqrt(tau^2+(sigmaE^2 + sigmaG^2)))
end


mu1 = theta1;
mu2 = mu1+thetadiff;



d1 = Normal(mu1,1);
d2 = Normal(mu2,1);
t1 = rand(d1,N1);
t2 = rand(d2,N2);

Pop1 = Array{Int64}(tmax);
Pop2 = Array{Int64}(tmax);
tnext1 = Array{Array{Float64}}(tmax);
tnext2 = Array{Array{Float64}}(tmax);
a_mean1 = Array{Float64}(tmax);
a_mean2 = Array{Float64}(tmax);
for t = 1:tmax
    #the number of strays leaving site 2 to site 1
    N1s = Int64(floor(N2*m));
    #the number of strays leaving site 1 to site 2
    N2s = Int64(floor(N1*m));
    
    #Proportion of the site that is local
    w1 = ((1-m)*N1)/((1-m)*N1 + m*N2);
    w2 = ((1-m)*N2)/((1-m)*N2 + m*N1);
    
    #strays leave the population
    N1 = N1 - N2s;
    N2 = N2 - N1s;
    
    t1_orig = copy(t1);
    t2_orig = copy(t2);
    
    #mixed trait distribution after dispersal
    t1strayloc = sort(unique(rand(collect(1:N1),N2s)));
    t2strayloc = sort(unique(rand(collect(1:N2),N1s)));
    
    t1_minus_strays = deleteat!(t1,t1strayloc);
    t2_minus_strays = deleteat!(t2,t2strayloc);
    
    td1 = [t1_minus_strays ; t2_orig[t2strayloc]];
    
    td2 = [t2_minus_strays ; t1_orig[t1strayloc]];
    
    Nr1 = Int64(floor((N1*rfunc(mu1,theta1) + N1s*rfunc(mu2,theta1))*exp(-beta*(N1 + N1s))));
    Nr2 = Int64(floor((N2*rfunc(mu2,theta2) + N2s*rfunc(mu1,theta2))*exp(-beta*(N2 + N2s))));
    
    #pair from mixture distribution to build recruit function 
    # dmix1 = MixtureModel(Normal[Normal(mu1, 1),Normal(mu2, 1)], [w1, (1-w1)]);
    # dmix2 = MixtureModel(Normal[Normal(mu2, 1),Normal(mu1, 1)], [w2, (1-w2)]);
    trecruit1 = (rand(td1,Nr1).+rand(td1,Nr1))/2;
    trecruit2 = (rand(td2,Nr2).+rand(td2,Nr2))/2;
    
    #propotional contribution of N1, N1s, Nr1 to site 1 population
    p11 = (N1*exp(-Z))/(exp(-Z)*(N1 + N1s) + Nr1);
    p12 = (N1s*exp(-Z))/(exp(-Z)*(N1 + N1s) + Nr1);
    
    #propotional contribution of N2, N2s, Nr2 to site 2 population
    p21 = (N2*exp(-Z))/(exp(-Z)*(N2 + N2s) + Nr2);
    p22 = (N2s*exp(-Z))/(exp(-Z)*(N2 + N2s) + Nr2);
    
    #dispersal trait distribution after mortality
    localsurvive1 = Int64(floor(N1*exp(-Z)));
    straysurvive1 = Int64(floor(N1s*exp(-Z)));
    
    localsurvive2 = Int64(floor(N2*exp(-Z)));
    straysurvive2 = Int64(floor(N2s*exp(-Z)));
    
    tdm1 = rand(td1,localsurvive1+straysurvive1);
    tdm2 = rand(td2,localsurvive2+straysurvive2);
    
    #NOW ACCOUNT FOR NATURAL SELECTION
    
    
    tnext1[t] = [tdm1;trecruit1];
    tnext2[t] = [tdm2;trecruit2];
    a_mean1[t] = w1*mean(t1) + (1-w1)*mean(t2);
    a_mean2[t] = w2*mean(t2) + (1-w2)*mean(t1);
    
    #Update traits and population sizes
    t1 = copy(tnext1[t]);
    t2 = copy(tnext2[t]);
    N1 = copy(length(tnext1[t]));
    N2 = copy(length(tnext2[t]));
    Pop1[t]=N1;
    Pop2[t]=N2;
    
end
t_mean1 = map(mean,tnext1);
t_std1 = map(std,tnext1);
t_mean2 = map(mean,tnext2);
t_std2 = map(std,tnext2);

R"hist($(tnext1[tmax]),xlim=c(-5,10))"

R"""
plot($t_mean1,$a_mean1)
lines(seq(0,2),seq(0,2))
"""

R"plot($t_std1)"

R"""
plot($t_mean1,ylim=c(-5,10))
points($t_mean2)
"""

R"""
plot($Pop1)
points($Pop2)
"""
