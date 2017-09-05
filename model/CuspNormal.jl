


b1vec = collect(-2:0.01:2);
b2vec = collect(-2:0.1:2);
tmax = 1000;
xtraj = Array{Float64}(length(b1vec),length(b2vec),tmax);
xss = Array{Float64}(length(b1vec),length(b2vec));
lambda = Array{Float64}(length(b1vec),length(b2vec));

for i=1:length(b1vec)
    for j=1:length(b2vec)
        b1 = b1vec[i];
        b2 = b2vec[j];
        x = Array{Float64}(tmax);
        x[1] = 0.1;
        for t = 2:(tmax-1)
            x[t+1] = x[t] + b1 + b2*x[t] - x[t]^3;
        end
        xtraj[i,j,:] = x;
        xss[i,j] = mean(x[tmax-100:tmax]);
        lambda[i,j] = 1 + b2 - 3*xss[i,j]^2; #$1 + 
    end
end

r = 19;
R"""
par(mfrow=c(3,1))
plot($(xtraj[206,r,100:150]),type='l')
plot($b1vec,$(lambda[:,r]),type='l')
lines(seq(-50,50),rep(1,101))
"""


R"""
plot($b1vec,$(xss[:,r]),type='l')
lines(seq(-50,50),rep(1,101))
"""




using Distributions
include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolveSS.jl")
include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinJacobian.jl")


#Analysis over m
tmax=10000;
mvec = [collect(0.0:0.001:0.1);collect(0.1:-0.001:0.001)];
n1ts = zeros(Float64,length(mvec),tmax);
n2ts = zeros(Float64,length(mvec),tmax);
n1mean=zeros(Float64,length(mvec));
n2mean=zeros(Float64,length(mvec));
n1sd=zeros(Float64,length(mvec));
n2sd=zeros(Float64,length(mvec));
aggmean=zeros(Float64,length(mvec));
aggsd=zeros(Float64,length(mvec));
x1mean=zeros(Float64,length(mvec));
x2mean=zeros(Float64,length(mvec));
pe=zeros(Float64,length(mvec));
eigs = Array(Array{Complex{Float64}},length(mvec));
maxeigs = Array{Float64}(length(mvec));
maximeigs = Array{Float64}(length(mvec));

z=0.5;
rmax=2.0;
beta=0.001;
theta1=5.0;
thetadiff=8;
tau=0.5;
h=0.5;
sigmaE=0;
sigmaG=1;
perror=0.01;

burnin=0.80
@time for i=1:length(mvec)
  
  m=mvec[i];
  
  if i == 1
      n0 = [2,2];
  else
      n0 = [n1mean[i-1],n2mean[i-1]];
  end
  
  n1, n2, x1, x2, w1, w2 = 
  KevinEvolveSS(
    n0,
    tmax,
    z,
    rmax,
    beta,
    theta1,
    thetadiff,
    tau,
    h,
    sigmaE,
    sigmaG,
    m,
    perror
    );

  n1trim = n1[Int64(floor(tmax*burnin)):tmax];
  n2trim = n2[Int64(floor(tmax*burnin)):tmax];
  x1trim = x1[Int64(floor(tmax*burnin)):tmax];
  x2trim = x2[Int64(floor(tmax*burnin)):tmax]
  
  n1ts[i,:] = n1;
  n2ts[i,:] = n2;
  
  n1mean[i] = mean(n1trim);
  n2mean[i] = mean(n2trim);
  n1sd[i] = std(n1trim);
  n2sd[i] = std(n2trim);
  
  aggmean[i] = mean(n1trim+n2trim);
  aggsd[i] = std(n1trim+n2trim);

  x1mean[i] = theta1-mean(x1trim);
  x2mean[i] = (theta1+thetadiff)-mean(x2trim);
  
  # #Calculate the Jacobian
  # Jac = KevinJacobian(mean(n1trim),mean(n2trim),mean(x1trim),mean(x2trim),
  # z,rmax,beta,theta1,thetadiff,tau,sigma,m)
  # eigs[i]=eigvals(Jac)
  # 
  # re = real(eigs[i]);
  # im = imag(eigs[i]);
  # maxeigs[i] = maximum(re);
  # maximeigs[i] = maximum(im);
  
  # pe[i] = (mean([std(n1trim),std(n2trim)])/mean([mean(n1trim),mean(n2trim)])) *
  # (1/(std(n1trim+n2trim)/mean(n1trim+n2trim)))
  
  pe[i] = mean([(std(n1trim)/mean(n1trim)),(std(n2trim)/mean(n2trim))])*
  (1/(std(n1trim+n2trim)/mean(n1trim+n2trim)))
  
end
midpoint = Int64(floor(length(n1mean)/2));
namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/figs2/fig_hysteresis.pdf");
R"""
library(RColorBrewer)
cols = brewer.pal(3,'Set1')
pdf($namespace,height=5,width=6)
plot($(mvec[1:midpoint]),$(n1mean[1:midpoint]),col=cols[1],xlab="Stray rate",ylab="Steady state",ylim=c(1,1000),type='l')
lines($(mvec[1:midpoint]),$(n2mean[1:midpoint]),col=cols[2])
lines($(mvec[midpoint+1:length(mvec)]),$(n1mean[midpoint+1:length(n1mean)]),col=cols[1],lty=2)
lines($(mvec[midpoint+1:length(mvec)]),$(n2mean[midpoint+1:length(n2mean)]),col=cols[2],lty=2)
dev.off()
"""


