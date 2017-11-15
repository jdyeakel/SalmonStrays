#Figure S2
#Figure S4

using RCall


thetascale = 2;
m=collect(0:0.01:0.5);
thetadiff=zeros(length(m));
for i = 1:length(m)
    thetadiff[i] = (1-2*m[i])/(thetascale*m[i])
end

namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/FinalDraft_rev/fig_mthetarelation.pdf");
R"""
pdf($namespace,height=5,width=6)
plot($thetadiff,$m,type='l',ylab='Straying rate m',xlab=expression(paste('Habitat heterogeneity ',Delta,theta)),ylim=c(0,0.5))
dev.off()
"""

using Distributions
using RCall

include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolve.jl")
include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinJacobian.jl")



#Analysis over m
tmax=10000;
mvec = collect(0.0:0.0001:1);
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
mineigs = Array{Float64}(length(mvec));
maximeigs = Array{Float64}(length(mvec));
minimeigs = Array{Float64}(length(mvec));

z=0.5;
rmax=2.0;
beta=0.001;
theta1=5.0;
thetadiff=4;
tau=1.0;
h=0.5;
sigmaE=0;
sigmaG=1;
perror=0.01;

burnin=0.80
@time for i=1:length(mvec)
  m=mvec[i];
  
  n1, n2, x1, x2, w1, w2 = 
  KevinEvolve(
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
  
  #Calculate the Jacobian
  Jac = KevinJacobian(mean(n1trim),mean(n2trim),mean(x1trim),mean(x2trim),
  z,rmax,beta,theta1,thetadiff,tau,h,sigmaE,sigmaG,m)
  eigs[i]=eigvals(Jac)
  
  re = real(eigs[i]);
  im = imag(eigs[i]);
  maxeigs[i] = maximum(re);
  mineigs[i] = minimum(re);
  maximeigs[i] = maximum(im);
  minimeigs[i] = minimum(im);
  
  # pe[i] = (mean([std(n1trim),std(n2trim)])/mean([mean(n1trim),mean(n2trim)])) *
  # (1/(std(n1trim+n2trim)/mean(n1trim+n2trim)))
  
  pe[i] = mean([(std(n1trim)/mean(n1trim)),(std(n2trim)/mean(n2trim))])*
  (1/(std(n1trim+n2trim)/mean(n1trim+n2trim)))
  
end

R"""
library(RColorBrewer)
cols = brewer.pal(3,'Set1')
plot($mvec,$n1mean,pch='.',col=cols[1],xlab="Stray rate",ylab="Steady state",cex=0.5)
points($mvec,$n2mean,pch='.',col=cols[2],cex=0.5)
"""

reigs = Array{Float64}(length(mvec),4);
for i=1:length(mvec)
    reigs[i,1] = real(eigs[i][1]);
    reigs[i,2] = real(eigs[i][2]);
    reigs[i,3] = real(eigs[i][3]);
    reigs[i,4] = real(eigs[i][4]);
end
#Plot Jacobian Eigenvalues
namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/figs2/fig_eigs.pdf");
R"""
library(RColorBrewer)
cols = brewer.pal(3,'Set1')
pdf($namespace,height=8,width=5)
par(mfrow=c(2,1),mai = c(0.8, 0.8, 0.1, 0.1))
plot($mvec,$n1mean,pch='.',col='black',cex=0.5,xlab='Straying rate m',ylab='Steady state biomass')
points($mvec,$n2mean,pch='.',col='black',cex=0.5)

plot($(mvec),$(reigs[:,1]),ylim=c(0,1),col='black',pch='.',cex=0.5,xlab='Straying rate m',ylab='Re[Jacobian eigenvalue]')
points($(mvec),$(reigs[:,2]),col='black',pch='.',cex=0.5)
points($(mvec),$(reigs[:,3]),col='black',pch='.',cex=0.5)
points($(mvec),$(reigs[:,4]),col='black',pch='.',cex=0.5)
lines(seq(-1,1,0.1),rep(1,21))
dev.off()
"""


using RCall
using Distributions
include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolveSS.jl")
include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinJacobian.jl")


#Analysis over m
tmax=10000;
mvec1 = collect(0.0:0.0001:0.15);
mvec = [mvec1 ; reverse(mvec1)];
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
mineigs = Array{Float64}(length(mvec));
minimeigs = Array{Float64}(length(mvec));

z=0.5;
rmax=2.0;
beta=0.001;
theta1=5.0;
thetadiff=6;
tau=1;
h=0.7;
sigmaE=0;
sigmaG=1;
perror=0.01;

burnin=0.80
@time for i=1:length(mvec)
  
  m=mvec[i];
  
  if i == 1
      n0 = [2,2];
    #   x0 = [theta1 + rand(Normal(0,0.0001)),(theta1 + thetadiff)  + rand(Normal(0,0.0001))];
    x0 = [theta1 ,(theta1 + thetadiff)];
  else
      n0 = [n1mean[i-1],n2mean[i-1]];
      x0 = [x1mean[i-1],x2mean[i-1]];
  end
  
  n1, n2, x1, x2, w1, w2 = 
  KevinEvolveSS(
    n0,
    x0,
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

  x1mean[i] = mean(x1trim);
  x2mean[i] = mean(x2trim);
  
  # #Calculate the Jacobian
  #Calculate the Jacobian
  Jac = KevinJacobian(mean(n1trim),mean(n2trim),mean(x1trim),mean(x2trim),
  z,rmax,beta,theta1,thetadiff,tau,h,sigmaE,sigmaG,m)
  eigs[i]=eigvals(Jac)
  
  re = real(eigs[i]);
  im = imag(eigs[i]);
  maxeigs[i] = maximum(re);
  mineigs[i] = minimum(re);
  maximeigs[i] = maximum(im);
  minimeigs[i] = minimum(im);
  
  # pe[i] = (mean([std(n1trim),std(n2trim)])/mean([mean(n1trim),mean(n2trim)])) *
  # (1/(std(n1trim+n2trim)/mean(n1trim+n2trim)))
  
  pe[i] = mean([(std(n1trim)/mean(n1trim)),(std(n2trim)/mean(n2trim))])*
  (1/(std(n1trim+n2trim)/mean(n1trim+n2trim)))
  
end
midpoint = Int64(floor(length(n1mean)/2));
namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/FinalDraft_rev/fig_hysteresis.pdf");
R"""
library(RColorBrewer)
cols = brewer.pal(3,'Set1')
pdf($namespace,height=5,width=6)
plot($(mvec[1:midpoint]),$(n1mean[1:midpoint]),col=cols[1],xlab="Straying rate m",ylab="Steady state biomass",ylim=c(1,1500),type='l')
lines($(mvec[1:midpoint]),$(n2mean[1:midpoint]),col=cols[2])
lines($(mvec[midpoint+1:length(mvec)]),$(n1mean[midpoint+1:length(n1mean)]),col=cols[1],lty=2)
lines($(mvec[midpoint+1:length(mvec)]),$(n2mean[midpoint+1:length(n2mean)]),col=cols[2],lty=2)
types = c('Increasing m','Decreasing m')
legend(x=0.1,y=1500,legend=types,col='black',lty=c(1,2),xpd=TRUE,cex=0.9, bty="n")
dev.off()
"""

