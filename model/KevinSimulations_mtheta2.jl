using Distributions
using RCall

include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolve_mtheta.jl")

#Analysis over m
tmax=10000;
mvec = collect(0.0:0.001:0.5);
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
thetascale=2;
tau=1.0;
h=0.2;
sigmaE=0;
sigmaG=1;
perror=0.01;

burnin=0.99
@time for i=1:length(mvec)
  m=mvec[i];
  
  n1, n2, x1, x2, w1, w2, thetadiff = 
  KevinEvolve_mtheta(
    tmax,
    z,
    rmax,
    beta,
    theta1,
    thetascale,
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
  
  pe[i] = mean([(std(n1trim)/mean(n1trim)),(std(n2trim)/mean(n2trim))])*
  (1/(std(n1trim+n2trim)/mean(n1trim+n2trim)))
  
end

R"""
library(RColorBrewer)
cols = brewer.pal(3,'Set1')
plot($mvec,$n1mean,pch='.',col=cols[1],xlab="Stray rate",ylab="Steady state",cex=0.5)
points($mvec,$n2mean,pch='.',col=cols[2],cex=0.5)
"""


#Portfolio effect plot
R"""
plot($mvec,$pe)
"""








@everywhere using Distributions
@everywhere using RCall
@everywhere using HDF5
@everywhere using JLD

@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolve_mtheta.jl")

@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/qualsfunc.jl")


#Analysis over m & theta divergence
mvec=collect(0.0:0.001:0.45);
hvec = collect(0.0:0.01:1.0);


n1mean=SharedArray(Float64,length(hvec),length(mvec));
n2mean=SharedArray(Float64,length(hvec),length(mvec));
x1mean=SharedArray(Float64,length(hvec),length(mvec));
x2mean=SharedArray(Float64,length(hvec),length(mvec));
pe=SharedArray(Float64,length(hvec),length(mvec));

tmax=10000;
z=0.5;
rmax=2.0;
beta=0.001;
theta1=5.0;
thetascale=2.0;
tau=1.0;
sigmaE=0;
sigmaG=1;


perror=0.01;

@sync @parallel for k=1:length(hvec)
  h = hvec[k];
  
  for j=1:length(mvec)
    m = mvec[j];
    thetadiff = (1-2*m)/(thetascale*m);
    
    n1, n2, x1, x2, w1, w2 = 
    KevinEvolve_mtheta(
      tmax,
      z,
      rmax,
      beta,
      theta1,
      thetascale,
      tau,
      h,
      sigmaE,
      sigmaG,
      m,
      perror
      );
    burnin=0.8
    n1trim = n1[Int64(floor(tmax*burnin)):tmax];
    n2trim = n2[Int64(floor(tmax*burnin)):tmax];
    x1trim = x1[Int64(floor(tmax*burnin)):tmax];
    x2trim = x2[Int64(floor(tmax*burnin)):tmax]
    
    n1mean[k,j] = mean(n1trim);
    n2mean[k,j] = mean(n2trim);
    x1mean[k,j] = theta1-mean(x1trim);
    x2mean[k,j] = (theta1+thetadiff)-mean(x2trim);
    
    pe[k,j] = mean([(std(n1trim)/mean(n1trim)),(std(n2trim)/mean(n2trim))])*
      (1/(std(n1trim+n2trim)/mean(n1trim+n2trim)))
  end
end

save(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data2/data_sig_h_m_mtheta.jld"),"n1mean",n1mean,"n2mean",n2mean,"x1mean",x1mean,"x2mean",x2mean,"pe",pe);


d = load(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data2/data_sig_h_m_mtheta.jld"));
#This loads the dictionary
n1mean = d["n1mean"];
n2mean = d["n2mean"];
x1mean = d["x1mean"];
x2mean = d["x2mean"];
pe = d["pe"];



bifvalue = bifdet(
n1mean,
n2mean,
indmvec,
hvec
);

namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/figs2/fig_MDPE_hm_mtheta.pdf");
R"""
library(RColorBrewer)
pal = rev(brewer.pal(9,"Blues"))
pdf($namespace,height=3,width=8)
par(mfrow=c(1,3))
image(x=$mvec,y=$hvec,z=t($(n1mean[:,:]))+t($(n2mean[:,:])),zlim=c(0,3000),col=pal,xlab='m',ylab='h',main='Total biomass')
lines($bifvalue)
image(x=$mvec,y=$hvec,z=sqrt((t($(n1mean[:,:]))-t($(n2mean[:,:])))^2),zlim=c(0,1500),col=pal,xlab='m',ylab='h',main='Biomass difference')
lines($bifvalue)
image(x=$mvec,y=$hvec,z=t($(pe[:,:])),zlim=c(1,2),col=pal,xlab='m',ylab='h',main='PE')
lines($bifvalue)
dev.off()
"""







#DDM simulations



@everywhere using Distributions
@everywhere using RCall
@everywhere using HDF5
@everywhere using JLD

@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolve_mtheta_ddm.jl")

#Analysis over m & theta divergence
indmvec=collect(0.0:0.001:0.45);
hvec = collect(0.0:0.01:1.0);


n1mean=SharedArray(Float64,length(hvec),length(indmvec));
n2mean=SharedArray(Float64,length(hvec),length(indmvec));
x1mean=SharedArray(Float64,length(hvec),length(indmvec));
x2mean=SharedArray(Float64,length(hvec),length(indmvec));
pe=SharedArray(Float64,length(hvec),length(indmvec));

tmax=10000;
z=0.5;
rmax=2.0;
beta=0.001;
theta1=5.0;
thetascale=2.0;
tau=1.0;
C=1000;
sigmaE=0;
sigmaG=1;


perror=0.01;

@sync @parallel for k=1:length(hvec)
  h = hvec[k];
  
  for j=1:length(indmvec)
    a0 = 1 - indmvec[j];
    m0 = indmvec[j];
    thetadiff = (1-2*m0)/(thetascale*m0);
    
    n1, n2, x1, x2, w1, w2 = 
    KevinEvolve_mtheta_ddm(
      tmax,
      z,
      rmax,
      beta,
      theta1,
      thetascale,
      tau,
      h,
      a0,
      C,
      sigmaE,
      sigmaG,
      perror
      );
    burnin=0.80
    n1trim = n1[Int64(floor(tmax*burnin)):tmax];
    n2trim = n2[Int64(floor(tmax*burnin)):tmax];
    x1trim = x1[Int64(floor(tmax*burnin)):tmax];
    x2trim = x2[Int64(floor(tmax*burnin)):tmax]
    
    n1mean[k,j] = mean(n1trim);
    n2mean[k,j] = mean(n2trim);
    x1mean[k,j] = theta1-mean(x1trim);
    x2mean[k,j] = (theta1+thetadiff)-mean(x2trim);
    
    pe[k,j] = mean([(std(n1trim)/mean(n1trim)),(std(n2trim)/mean(n2trim))])*
      (1/(std(n1trim+n2trim)/mean(n1trim+n2trim)))
  end
end

save(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data2/data_sig_h_mtheta_ddm.jld"),"n1mean",n1mean,"n2mean",n2mean,"x1mean",x1mean,"x2mean",x2mean,"pe",pe);


d = load(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data2/data_sig_h_mtheta_ddm.jld"));
#This loads the dictionary
n1mean = d["n1mean"];
n2mean = d["n2mean"];
x1mean = d["x1mean"];
x2mean = d["x2mean"];
pe = d["pe"];


bifvalue = bifdet(
n1mean,
n2mean,
indmvec,
hvec
);

namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/figs2/fig_MDPE_hm_mtheta_ddm.pdf");
R"""
library(RColorBrewer)
pal = rev(brewer.pal(9,"Blues"))
pdf($namespace,height=3,width=8)
par(mfrow=c(1,3))
image(x=$indmvec,y=$hvec,z=t($(n1mean[:,:]))+t($(n2mean[:,:])),zlim=c(0,3000),col=pal,xlab='m0',ylab='h',main='Total biomass')
lines($(bifvalue[1:60,:]))
image(x=$indmvec,y=$hvec,z=sqrt((t($(n1mean[:,:]))-t($(n2mean[:,:])))^2),zlim=c(0,1500),col=pal,xlab='m0',ylab='h',main='Biomass difference')
lines($(bifvalue[1:60,:]))
image(x=$indmvec,y=$hvec,z=t($(pe[:,:])),zlim=c(1,2),col=pal,xlab='m0',ylab='h',main='PE')
lines($(bifvalue[1:60,:]))
dev.off()
"""
