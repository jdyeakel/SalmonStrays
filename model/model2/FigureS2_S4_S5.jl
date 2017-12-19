#Figure S2
#Figure S4

using RCall


thetascale=20;
thetadiff=collect(0:0.01:3.0);
m = zeros(length(thetadiff));
for i=1:length(thetadiff)
    m[i] = 1/(2+thetascale*thetadiff[i])
end

namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/Manuscript/FinalDraft3/fig_mthetarelation.pdf");
R"""
pdf($namespace,height=5,width=6)
plot($thetadiff,$m,type='l',ylab='Straying rate m',xlab=expression(paste('Habitat heterogeneity ',Delta,theta)),ylim=c(0,0.5),xlim=c(0,3))
dev.off()
"""



using Distributions
using RCall

include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolve.jl")
include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinJacobian.jl")



#Analysis over m
tmax=10000;
mvec = collect(0.0:0.0001:0.3);
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


z=2.0;
rmax=2.0;
beta=0.001;
theta1=5.0;
thetadiff=2.1;
tau=1.0;
h=0.2;
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

# R"""
# library(RColorBrewer)
# cols = brewer.pal(3,'Set1')
# plot($mvec,$n1mean,pch='.',col=cols[1],xlab="Stray rate",ylab="Steady state",cex=0.5)
# points($mvec,$n2mean,pch='.',col=cols[2],cex=0.5)
# """

reigs = Array{Float64}(length(mvec),4);
for i=1:length(mvec)
    reigs[i,1] = real(eigs[i][1]);
    reigs[i,2] = real(eigs[i][2]);
    reigs[i,3] = real(eigs[i][3]);
    reigs[i,4] = real(eigs[i][4]);
end
#Plot Jacobian Eigenvalues
namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/FinalDraft3/fig_eigs.pdf");
R"""
library(RColorBrewer)
cols = brewer.pal(3,'Set1')
pdf($namespace,height=8,width=5)
par(mfrow=c(2,1),mai = c(0.8, 0.8, 0.1, 0.1))
plot($mvec,$n1mean,pch='.',col='black',cex=0.5,xlab='Straying ratio m',ylab='Steady state biomass')
points($mvec,$n2mean,pch='.',col='black',cex=0.5)

plot($(mvec),$(reigs[:,1]),ylim=c(0,1),col='black',pch='.',cex=0.5,xlab='Straying ratio m',ylab='Re[Jacobian eigenvalue]')
points($(mvec),$(reigs[:,2]),col='black',pch='.',cex=0.5)
points($(mvec),$(reigs[:,3]),col='black',pch='.',cex=0.5)
points($(mvec),$(reigs[:,4]),col='black',pch='.',cex=0.5)
lines(seq(-1,1,0.1),rep(1,21))
dev.off()
"""


@everywhere using RCall, Distributions, HDF5, JLD
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolveSS.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinJacobian.jl")


#Analysis over m
tmax=10000;
mvec1 = collect(0.0:0.01:0.4);
mvec = [mvec1 ; reverse(mvec1)];

thetadiffvec=[1.8 2.0 2.2 2.4 2.6 2.8];

n1ts = zeros(Float64,length(mvec),length(thetadiffvec),tmax);
n2ts = zeros(Float64,length(mvec),length(thetadiffvec),tmax);
n1mean=zeros(Float64,length(mvec),length(thetadiffvec));
n2mean=zeros(Float64,length(mvec),length(thetadiffvec));
n1sd=zeros(Float64,length(mvec),length(thetadiffvec));
n2sd=zeros(Float64,length(mvec),length(thetadiffvec));
aggmean=zeros(Float64,length(mvec),length(thetadiffvec));
aggsd=zeros(Float64,length(mvec),length(thetadiffvec));
x1mean=zeros(Float64,length(mvec),length(thetadiffvec));
x2mean=zeros(Float64,length(mvec),length(thetadiffvec));
pe=zeros(Float64,length(mvec),length(thetadiffvec));
# eigs = Array(Array{Complex{Float64}},length(mvec));
# maxeigs = Array{Float64}(length(mvec));
# maximeigs = Array{Float64}(length(mvec));
# mineigs = Array{Float64}(length(mvec));
# minimeigs = Array{Float64}(length(mvec));

z=2;
rmax=2.0;
beta=0.001;
theta1=5.0;

tau=1;
h=0.2;
sigmaE=0;
sigmaG=1;
perror=0.00;

burnin=0.80
@time for r=1:length(thetadiffvec)

    thetadiff=thetadiffvec[r];

    for i=1:length(mvec)
      
      m=mvec[i];
      
      if i == 1
          n0 = [2,2];
        #   x0 = [theta1 + rand(Normal(0,0.0001)),(theta1 + thetadiff)  + rand(Normal(0,0.0001))];
        x0 = [theta1 ,(theta1 + thetadiff)];
      else
          n0 = [n1mean[i-1,r],n2mean[i-1,r]];
          x0 = [x1mean[i-1,r],x2mean[i-1,r]];
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
      
      n1ts[i,r,:] = n1;
      n2ts[i,r,:] = n2;
      
      n1mean[i,r] = mean(n1trim);
      n2mean[i,r] = mean(n2trim);
      n1sd[i,r] = std(n1trim);
      n2sd[i,r] = std(n2trim);
      
      aggmean[i,r] = mean(n1trim+n2trim);
      aggsd[i,r] = std(n1trim+n2trim);

      x1mean[i,r] = mean(x1trim);
      x2mean[i,r] = mean(x2trim);
      
      # #Calculate the Jacobian
      #Calculate the Jacobian
      # Jac = KevinJacobian(mean(n1trim),mean(n2trim),mean(x1trim),mean(x2trim),
      # z,rmax,beta,theta1,thetadiff,tau,h,sigmaE,sigmaG,m)
      # eigs[i]=eigvals(Jac)
      # 
      # re = real(eigs[i]);
      # im = imag(eigs[i]);
      # maxeigs[i] = maximum(re);
      # mineigs[i] = minimum(re);
      # maximeigs[i] = maximum(im);
      # minimeigs[i] = minimum(im);
      
      # pe[i] = (mean([std(n1trim),std(n2trim)])/mean([mean(n1trim),mean(n2trim)])) *
      # (1/(std(n1trim+n2trim)/mean(n1trim+n2trim)))
      
      pe[i,r] = mean([(std(n1trim)/mean(n1trim)),(std(n2trim)/mean(n2trim))])*
      (1/(std(n1trim+n2trim)/mean(n1trim+n2trim)))
      
    end
end


midpoint = Int64(floor(length(mvec)/2));
namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/FinalDraft3/fig_hysteresisfull.pdf");
R"""
#pdf($namespace,height=8,width=12)
par(mfrow=c(2,3))
library(RColorBrewer)
cols = brewer.pal(3,'Set1')
plot($(mvec[1:midpoint]),$(n1mean[1:midpoint,1]),col=cols[1],xlab="Straying rate m",ylab="Steady state biomass",ylim=c(1,1500),type='l')
lines($(mvec[1:midpoint]),$(n2mean[1:midpoint,1]),col=cols[2])
lines($(mvec[midpoint+1:length(mvec)]),$(n1mean[midpoint+1:length(mvec),1]),col=cols[1],lty=2)
lines($(mvec[midpoint+1:length(mvec)]),$(n2mean[midpoint+1:length(mvec),1]),col=cols[2],lty=2)
types = c('Increasing m','Decreasing m')
legend(x=0.13,y=1500,legend=types,col='black',lty=c(1,2),xpd=TRUE,cex=0.9, bty="n")

plot($(mvec[1:midpoint]),$(n1mean[1:midpoint,2]),col=cols[1],xlab="Straying rate m",ylab="Steady state biomass",ylim=c(1,1500),type='l')
lines($(mvec[1:midpoint]),$(n2mean[1:midpoint,2]),col=cols[2])
lines($(mvec[midpoint+1:length(mvec)]),$(n1mean[midpoint+1:length(mvec),2]),col=cols[1],lty=2)
lines($(mvec[midpoint+1:length(mvec)]),$(n2mean[midpoint+1:length(mvec),2]),col=cols[2],lty=2)


plot($(mvec[1:midpoint]),$(n1mean[1:midpoint,3]),col=cols[1],xlab="Straying rate m",ylab="Steady state biomass",ylim=c(1,1500),type='l')
lines($(mvec[1:midpoint]),$(n2mean[1:midpoint,3]),col=cols[2])
lines($(mvec[midpoint+1:length(mvec)]),$(n1mean[midpoint+1:length(mvec),3]),col=cols[1],lty=2)
lines($(mvec[midpoint+1:length(mvec)]),$(n2mean[midpoint+1:length(mvec),3]),col=cols[2],lty=2)


plot($(mvec[1:midpoint]),$(n1mean[1:midpoint,4]),col=cols[1],xlab="Straying rate m",ylab="Steady state biomass",ylim=c(1,1500),type='l')
lines($(mvec[1:midpoint]),$(n2mean[1:midpoint,4]),col=cols[2])
lines($(mvec[midpoint+1:length(mvec)]),$(n1mean[midpoint+1:length(mvec),4]),col=cols[1],lty=2)
lines($(mvec[midpoint+1:length(mvec)]),$(n2mean[midpoint+1:length(mvec),4]),col=cols[2],lty=2)

plot($(mvec[1:midpoint]),$(n1mean[1:midpoint,5]),col=cols[1],xlab="Straying rate m",ylab="Steady state biomass",ylim=c(1,1500),type='l')
lines($(mvec[1:midpoint]),$(n2mean[1:midpoint,5]),col=cols[2])
lines($(mvec[midpoint+1:length(mvec)]),$(n1mean[midpoint+1:length(mvec),5]),col=cols[1],lty=2)
lines($(mvec[midpoint+1:length(mvec)]),$(n2mean[midpoint+1:length(mvec),5]),col=cols[2],lty=2)

plot($(mvec[1:midpoint]),$(n1mean[1:midpoint,6]),col=cols[1],xlab="Straying rate m",ylab="Steady state biomass",ylim=c(1,1500),type='l')
lines($(mvec[1:midpoint]),$(n2mean[1:midpoint,6]),col=cols[2])
lines($(mvec[midpoint+1:length(mvec)]),$(n1mean[midpoint+1:length(mvec),6]),col=cols[1],lty=2)
lines($(mvec[midpoint+1:length(mvec)]),$(n2mean[midpoint+1:length(mvec),6]),col=cols[2],lty=2)



#dev.off()
"""



# 2-dimensional search


#Analysis over m
tmax=10000;
mvec1 = collect(0.0:0.001:0.5);
mvec = [mvec1 ; reverse(mvec1)];

thetadiffvec=collect(1.5:0.005:2.5);

n1ts = Array{Float64}(length(mvec),length(thetadiffvec),tmax);
n2ts = Array{Float64}(length(mvec),length(thetadiffvec),tmax);
n1mean=Array{Float64}(length(mvec),length(thetadiffvec));
n2mean=Array{Float64}(length(mvec),length(thetadiffvec));
n1sd=Array{Float64}(length(mvec),length(thetadiffvec));
n2sd=Array{Float64}(length(mvec),length(thetadiffvec));
aggmean=Array{Float64}(length(mvec),length(thetadiffvec));
aggsd=Array{Float64}(length(mvec),length(thetadiffvec));
x1mean=Array{Float64}(length(mvec),length(thetadiffvec));
x2mean=Array{Float64}(length(mvec),length(thetadiffvec));
pe=Array{Float64}(length(mvec),length(thetadiffvec));
# eigs = Array(Array{Complex{Float64}},length(mvec));
maxeigs = Array{Float64}(length(mvec),length(thetadiffvec));
maximeigs = Array{Float64}(length(mvec),length(thetadiffvec));
mineigs = Array{Float64}(length(mvec),length(thetadiffvec));
minimeigs = Array{Float64}(length(mvec),length(thetadiffvec));

z=2.0;
rmax=2.0;
beta=0.001;
theta1=5.0;

tau=1;
h=0.2;
sigmaE=0;
sigmaG=1;
perror=0.00;

burnin=0.80
for r=1:length(thetadiffvec)

    thetadiff=thetadiffvec[r];

    for i=1:length(mvec)
      
      m=mvec[i];
      
      if i == 1
          n0 = [2,2];
        #   x0 = [theta1 + rand(Normal(0,0.0001)),(theta1 + thetadiff)  + rand(Normal(0,0.0001))];
        x0 = [theta1 ,(theta1 + thetadiff)];
      else
          n0 = [n1mean[i-1,r],n2mean[i-1,r]];
          x0 = [x1mean[i-1,r],x2mean[i-1,r]];
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
      
      n1ts[i,r,:] = n1;
      n2ts[i,r,:] = n2;
      
      n1mean[i,r] = mean(n1trim);
      n2mean[i,r] = mean(n2trim);
      n1sd[i,r] = std(n1trim);
      n2sd[i,r] = std(n2trim);
      
      aggmean[i,r] = mean(n1trim+n2trim);
      aggsd[i,r] = std(n1trim+n2trim);

      x1mean[i,r] = mean(x1trim);
      x2mean[i,r] = mean(x2trim);
      
      # #Calculate the Jacobian
      #Calculate the Jacobian
      # Jac = KevinJacobian(mean(n1trim),mean(n2trim),mean(x1trim),mean(x2trim),
      # z,rmax,beta,theta1,thetadiff,tau,h,sigmaE,sigmaG,m)
      # eigs=eigvals(Jac)
      # re = real(eigs);
      # im = imag(eigs);
      # maxeigs[i,r] = maximum(re);
      # mineigs[i,r] = minimum(re);
      # maximeigs[i,r] = maximum(im);
      # minimeigs[i,r] = minimum(im);
      
      # pe[i] = (mean([std(n1trim),std(n2trim)])/mean([mean(n1trim),mean(n2trim)])) *
      # (1/(std(n1trim+n2trim)/mean(n1trim+n2trim)))
      
      pe[i,r] = mean([(std(n1trim)/mean(n1trim)),(std(n2trim)/mean(n2trim))])*
      (1/(std(n1trim+n2trim)/mean(n1trim+n2trim)))
      
    end
end

save(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data3/data_hysteresis.jld"),"n1mean",n1mean,"n2mean",n2mean,"x1mean",x1mean,"x2mean",x2mean,"pe",pe,"maxeigs",maxeigs);



#Plot the difference in n1mean, n2mean
midpoint = Int64(floor(length(mvec)/2));

fmaxeigs = maxeigs[1:midpoint,:];
bmaxeigs = flipdim(maxeigs[midpoint+1:length(mvec),:],1);
fold = 0.98 .< fmaxeigs .< 1.00;
bold = 0.98 .< bmaxeigs .< 1.00;

fn1mean = n1mean[1:midpoint,:];
bn1mean = flipdim(n1mean[midpoint+1:length(mvec),:],1);

fn2mean = n2mean[1:midpoint,:];
bn2mean = flipdim(n2mean[midpoint+1:length(mvec),:],1);

fmvec = mvec[1:midpoint];
bmvec = fmvec;

fdiffmean = fn1mean .- fn2mean;
bdiffmean = bn1mean .- bn2mean;
fdiffmeanbin = fdiffmean .> 0.000001;
bdiffmeanbin = bdiffmean .> 0.000001;

hystdiff = fdiffmeanbin .!= bdiffmeanbin;

namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/FinalDraft3/fig_hysteresis2Dhys_eig.pdf");
R"""
#pdf($namespace,height=5,width=6)
par(mfrow=c(1,1))
image(x=$fmvec,y=$thetadiffvec,z=$fdiffmeanbin,xlab='Straying ratio m',ylab='Delta theta',col=c('white','black'))
image(x=$fmvec,y=$thetadiffvec,z=$hystdiff,xlab='Straying ratio m',ylab='Delta theta',col=c('#ffffff00','gray'),add=T)
#image(x=$fmvec,y=$thetadiffvec,z=$fold,xlab='Straying ratio m',ylab='Delta theta',col=c('#ffffff00','red'),add=T)
#image(x=$fmvec,y=$thetadiffvec,z=$bold,xlab='Straying ratio m',ylab='Delta theta',col=c('#ffffff00','red'),add=T)
#dev.off()
"""



midpoint = Int64(floor(length(mvec)/2));
namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/FinalDraft3/fig_hysteresisfull.pdf");
R"""
#pdf($namespace,height=8,width=12)
par(mfrow=c(2,3))
library(RColorBrewer)
cols = brewer.pal(3,'Set1')
plot($(mvec[1:midpoint]),$(n1mean[1:midpoint,1]),col=cols[1],xlab="Straying rate m",ylab="Steady state biomass",ylim=c(1,1500),type='l')
lines($(mvec[1:midpoint]),$(n2mean[1:midpoint,1]),col=cols[2])
lines($(mvec[midpoint+1:length(mvec)]),$(n1mean[midpoint+1:length(mvec),1]),col=cols[1],lty=2)
lines($(mvec[midpoint+1:length(mvec)]),$(n2mean[midpoint+1:length(mvec),1]),col=cols[2],lty=2)
types = c('Increasing m','Decreasing m')
legend(x=0.13,y=1500,legend=types,col='black',lty=c(1,2),xpd=TRUE,cex=0.9, bty="n")

#dev.off()
"""






#Analysis over m
tmax=10000;
mvec1 = collect(0.0:0.0001:0.2);
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
# eigs = Array(Array{Complex{Float64}},length(mvec));
# maxeigs = Array{Float64}(length(mvec));
# maximeigs = Array{Float64}(length(mvec));
# mineigs = Array{Float64}(length(mvec));
# minimeigs = Array{Float64}(length(mvec));

z=2.0;
rmax=2.0;
beta=0.001;
theta1=5.0;
thetadiff=2.1;
tau=1;
h=0.2;
sigmaE=0;
sigmaG=1;
perror=0.0;

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
  # Jac = KevinJacobian(mean(n1trim),mean(n2trim),mean(x1trim),mean(x2trim),
  # z,rmax,beta,theta1,thetadiff,tau,h,sigmaE,sigmaG,m)
  # eigs[i]=eigvals(Jac)
  # 
  # re = real(eigs[i]);
  # im = imag(eigs[i]);
  # maxeigs[i] = maximum(re);
  # mineigs[i] = minimum(re);
  # maximeigs[i] = maximum(im);
  # minimeigs[i] = minimum(im);
  
  # pe[i] = (mean([std(n1trim),std(n2trim)])/mean([mean(n1trim),mean(n2trim)])) *
  # (1/(std(n1trim+n2trim)/mean(n1trim+n2trim)))
  
  pe[i] = mean([(std(n1trim)/mean(n1trim)),(std(n2trim)/mean(n2trim))])*
  (1/(std(n1trim+n2trim)/mean(n1trim+n2trim)))
  
end


midpoint = Int64(floor(length(n1mean[:,1])/2));
namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/FinalDraft3/fig_hysteresis.pdf");
R"""
library(RColorBrewer)
cols = brewer.pal(3,'Set1')
pdf($namespace,height=5,width=6)
plot($(mvec[1:midpoint]),$(n1mean[1:midpoint]),col=cols[1],xlab="Straying rate m",ylab="Steady state biomass",ylim=c(1,500),type='l')
lines($(mvec[1:midpoint]),$(n2mean[1:midpoint]),col=cols[2])
lines($(mvec[midpoint+1:length(mvec)]),$(n1mean[midpoint+1:length(n1mean)]),col=cols[1],lty=2)
lines($(mvec[midpoint+1:length(mvec)]),$(n2mean[midpoint+1:length(n2mean)]),col=cols[2],lty=2)
dev.off()
"""






# 2-dimensional search over Z


#Analysis over m
tmax=10000;
mvec1 = collect(0.0:0.001:0.2);
mvec = [mvec1 ; reverse(mvec1)];

zvec=collect(0.1:0.01:5.0);

n1ts = zeros(Float64,length(mvec),length(zvec),tmax);
n2ts = zeros(Float64,length(mvec),length(zvec),tmax);
n1mean=zeros(Float64,length(mvec),length(zvec));
n2mean=zeros(Float64,length(mvec),length(zvec));
n1sd=zeros(Float64,length(mvec),length(zvec));
n2sd=zeros(Float64,length(mvec),length(zvec));
aggmean=zeros(Float64,length(mvec),length(zvec));
aggsd=zeros(Float64,length(mvec),length(zvec));
x1mean=zeros(Float64,length(mvec),length(zvec));
x2mean=zeros(Float64,length(mvec),length(zvec));
pe=zeros(Float64,length(mvec),length(zvec));
# eigs = Array(Array{Complex{Float64}},length(mvec));
# maxeigs = Array{Float64}(length(mvec));
# maximeigs = Array{Float64}(length(mvec));
# mineigs = Array{Float64}(length(mvec));
# minimeigs = Array{Float64}(length(mvec));


rmax=2.0;
beta=0.001;
theta1=5.0;
thetadiff = 3.2;
tau=1;
h=0.2;
sigmaE=0;
sigmaG=1;
perror=0.00;

burnin=0.80
@time for r=1:length(zvec)

    z=zvec[r];

    for i=1:length(mvec)
      
      m=mvec[i];
      
      if i == 1
          n0 = [2,2];
        #   x0 = [theta1 + rand(Normal(0,0.0001)),(theta1 + thetadiff)  + rand(Normal(0,0.0001))];
        x0 = [theta1 ,(theta1 + thetadiff)];
      else
          n0 = [n1mean[i-1,r],n2mean[i-1,r]];
          x0 = [x1mean[i-1,r],x2mean[i-1,r]];
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
      
      n1ts[i,r,:] = n1;
      n2ts[i,r,:] = n2;
      
      n1mean[i,r] = mean(n1trim);
      n2mean[i,r] = mean(n2trim);
      n1sd[i,r] = std(n1trim);
      n2sd[i,r] = std(n2trim);
      
      aggmean[i,r] = mean(n1trim+n2trim);
      aggsd[i,r] = std(n1trim+n2trim);

      x1mean[i,r] = mean(x1trim);
      x2mean[i,r] = mean(x2trim);
      
      # #Calculate the Jacobian
      #Calculate the Jacobian
      # Jac = KevinJacobian(mean(n1trim),mean(n2trim),mean(x1trim),mean(x2trim),
      # z,rmax,beta,theta1,thetadiff,tau,h,sigmaE,sigmaG,m)
      # eigs[i]=eigvals(Jac)
      # 
      # re = real(eigs[i]);
      # im = imag(eigs[i]);
      # maxeigs[i] = maximum(re);
      # mineigs[i] = minimum(re);
      # maximeigs[i] = maximum(im);
      # minimeigs[i] = minimum(im);
      
      # pe[i] = (mean([std(n1trim),std(n2trim)])/mean([mean(n1trim),mean(n2trim)])) *
      # (1/(std(n1trim+n2trim)/mean(n1trim+n2trim)))
      
      pe[i,r] = mean([(std(n1trim)/mean(n1trim)),(std(n2trim)/mean(n2trim))])*
      (1/(std(n1trim+n2trim)/mean(n1trim+n2trim)))
      
    end
end

save(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data3/data_hysteresisZ.jld"),"n1mean",n1mean,"n2mean",n2mean,"x1mean",x1mean,"x2mean",x2mean,"pe",pe);


#Plot the difference in n1mean, n2mean
midpoint = Int64(floor(length(mvec)/2));

fn1mean = n1mean[1:midpoint,:];
bn1mean = flipdim(n1mean[midpoint+1:length(mvec),:],1);

fn2mean = n2mean[1:midpoint,:];
bn2mean = flipdim(n2mean[midpoint+1:length(mvec),:],1);

fmvec = mvec[1:midpoint];
bmvec = fmvec;

fdiffmean = fn1mean .- fn2mean;
bdiffmean = bn1mean .- bn2mean;
fdiffmeanbin = fdiffmean .> 0.000001;
bdiffmeanbin = bdiffmean .> 0.000001;

hystdiff = fdiffmeanbin .!= bdiffmeanbin;

namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/FinalDraft3/fig_hysteresis2DZ.pdf");
R"""
pdf($namespace,height=5,width=6)
par(mfrow=c(1,1))
image(x=$fmvec,y=$zvec,z=$fdiffmeanbin,xlab='Straying ratio m',ylab='Z',col=c('white','black'))
image(x=$fmvec,y=$zvec,z=$hystdiff,xlab='Straying ratio m',ylab='Z',col=c('#ffffff00','gray'),add=T)
dev.off()
"""
