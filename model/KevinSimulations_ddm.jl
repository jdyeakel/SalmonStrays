using Distributions
using RCall

include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolve_ddm.jl")
include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinJacobian.jl")


tmax=10000;
z=0.5;
rmax=2.0;
beta=0.001;
theta1=5.0;
thetadiff=8;
tau=1.0;
h=1.0;
indm=0.42;
C=1000;
sigmaE=0.0;
sigmaG=1.0;
perror=0.05;

a0 = 1-indm;

n1, n2, x1, x2, w1, w2 = 
KevinEvolve_ddm(
  tmax,
  z,
  rmax,
  beta,
  theta1,
  thetadiff,
  tau,
  h,
  a0,
  C,
  sigmaE,
  sigmaG,
  perror
  );
  
  R"""
  library(RColorBrewer)
  cols = brewer.pal(3,'Set1')
  plot($(n1),type="l",col=cols[1],ylim=c(0,max(c(as.vector($n1),as.vector($n2)))))
  lines($(n2),col=cols[2])
  """
  
  R"""
  library(RColorBrewer)
  cols = brewer.pal(3,'Set1')
  plot($(x1),type="l",col=cols[1],ylim=c(4,12))
  lines($(x2),col=cols[2])
  """


#CalculateJacobian
Jac = KevinJacobian(last(n1),last(n2),last(x1),last(x2),
z,rmax,beta,theta1,thetadiff,tau,sigma,m)
eigvals(Jac)






#Analysis over a0
tmax=10000;
a0vec = collect(0.0:0.002:1.0);
n1ts = zeros(Float64,length(a0vec),tmax);
n2ts = zeros(Float64,length(a0vec),tmax);
n1mean=zeros(Float64,length(a0vec));
n2mean=zeros(Float64,length(a0vec));
n1sd=zeros(Float64,length(a0vec));
n2sd=zeros(Float64,length(a0vec));
aggmean=zeros(Float64,length(a0vec));
aggsd=zeros(Float64,length(a0vec));
x1mean=zeros(Float64,length(a0vec));
x2mean=zeros(Float64,length(a0vec));
pe=zeros(Float64,length(a0vec));
eigs = Array(Array{Complex{Float64}},length(a0vec));
maxeigs = Array{Float64}(length(a0vec));
maximeigs = Array{Float64}(length(a0vec));

z=0.5;
rmax=2.0;
beta=0.001;
theta1=5.0;
thetadiff=8;
tau=1.0;
h=0.2;
C=1000;
sigmaE=0;
sigmaG=2;
perror=0.05;

burnin=0.99
@time for i=1:length(a0vec)
  a0 = a0vec[i];
  
  n1, n2, x1, x2, w1, w2 = 
  KevinEvolve_ddm(
    tmax,
    z,
    rmax,
    beta,
    theta1,
    thetadiff,
    tau,
    h,
    a0,
    C,
    sigmaE,
    sigmaG,
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
  
  pe[i] = (mean([std(n1trim),std(n2trim)])/mean([mean(n1trim),mean(n2trim)])) *
  (1/(std(n1trim+n2trim)/mean(n1trim+n2trim)))
  
end

R"""
library(RColorBrewer)
cols = brewer.pal(3,'Set1')
plot(1-$a0vec,$n1mean,pch=16,col=cols[1],xlab="Individual Stray Rate",ylab="Steady state",cex=0.5)
points(1-$a0vec,$n2mean,pch=16,col=cols[2],cex=0.5)
"""


#Steady state plot
namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/fig_SSm.pdf");
R"""
library(RColorBrewer)
cols = brewer.pal(3,'Set1')
pdf($namespace,height=5,width=6)
plot($mvec,$n1mean,pch='.',col=cols[1],xlab="Stray rate",ylab="Steady state",cex=0.5)
points($mvec,$n2mean,pch='.',col=cols[2],cex=0.5)
dev.off()
"""
#Trait offset plot
namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/fig_Traitm.pdf");
R"""
library(RColorBrewer)
cols = brewer.pal(3,'Set1')
pdf($namespace,height=5,width=6)
plot($mvec,$x1mean,pch=".",col=cols[1],ylim=c(-5,5),xlab="Stray rate",ylab="Trait offset")
points($mvec,$x2mean,pch=".",col=cols[2])
dev.off()
"""
#Portfolio effect plot
namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/fig_pem.pdf");
R"""
library(RColorBrewer)
cols = brewer.pal(3,'Set1')
pdf($namespace,height=5,width=6)
plot($mvec,$pe,pch='.',col="black",xlab="Stray rate",ylab="Portfolio effect")
dev.off()
"""

#Portfolio effect plot
R"""
plot(1-$a0vec,$pe)
"""




R"""
par(mfrow=c(2,1))
plot($mvec,$aggmean,pch=16)
points($mvec,$n1mean,col='blue')
points($mvec,$n2mean)
"""
R"""
plot($mvec,$n1sd)
points($mvec,$n2sd,col='blue')
points($mvec,$aggsd,pch=16)
"""


maxsd=find(x->x==maximum(n1sd),n1sd)[1];
R"""
plot($(n1ts[maxsd,Int64(floor(tmax*burnin)):tmax]),type='l',ylim=c(0,900))
points($(n1ts[maxsd-1,Int64(floor(tmax*burnin)):tmax]),type='l',col='blue')
"""

#Eigenvalue plots
R"""
par(mar=c(4.1,4.1,0.9,2.1))
par(mfrow=c(2,1))
library(RColorBrewer)
cols = brewer.pal(3,'Set1')
plot($(mvec),$maxeigs,type='l',col=cols[1],xlab="Stray rate",ylab="Max re(eig)")
lines(seq(0,0.5,length.out=10),rep(1,10),lty=3)
plot($mvec,$n1mean,pch='.',col=cols[1],xlab="Stray rate",ylab="Steady state")
points($mvec,$n2mean,pch='.',col=cols[2])
"""




R"""
library(RColorBrewer)
cols = brewer.pal(3,'Set1')
plot($(mvec),$maximeigs,type='l',col=cols[1],xlab="Stray rate",ylab="max im(eig)")
#points($(mvec),$unitcircle,type='l',col='black')*/*/*/
"""


R"""
library(RColorBrewer)
cols = brewer.pal(3,'Set1')
plot($(n1mean),$maxeigs,pch=16,cex=0.5,col=cols[1],xlab="n1 mean",ylab="max re(eig)")
"""


#Steady state plot
namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/fig_SSm2.pdf");
R"""
library(RColorBrewer)
cols = brewer.pal(3,'Set1')
pdf($namespace,height=5,width=6)
plot($mvec,$n1mean,pch=".",col=cols[1],xlab="Stray rate",ylab="Steady state")
points($mvec,$n2mean,pch=".",col=cols[2])
points($mvec,$n1meanNONE,pch=".",col=cols[1])
points($mvec,$n2meanNONE,pch=".",col=cols[2])
dev.off()
"""





@everywhere using Distributions
@everywhere using RCall

@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolve_ddm.jl")

#Analysis over m & theta divergence
a0vec=collect(0.0:0.001:1);
thetadiffvec = collect(5:1:15);


n1mean=SharedArray(Float64,length(thetadiffvec),length(a0vec));
n2mean=SharedArray(Float64,length(thetadiffvec),length(a0vec));
x1mean=SharedArray(Float64,length(thetadiffvec),length(a0vec));
x2mean=SharedArray(Float64,length(thetadiffvec),length(a0vec));
pe=SharedArray(Float64,length(thetadiffvec),length(a0vec));

tmax=100000;
z=0.5;
rmax=2.0;
beta=0.001;
theta1=5.0;
tau=1.0;
h=1.0;
C=1000;
sigma=1.0;
perror=0.05;

@sync @parallel for i=1:length(thetadiffvec)
  for j=1:length(a0vec)
    
    thetadiff = thetadiffvec[i];
    a0 = a0vec[j];
    
    n1, n2, x1, x2, w1, w2 = 
    KevinEvolve_ddm(
      tmax,
      z,
      rmax,
      beta,
      theta1,
      thetadiff,
      tau,
      h,
      a0,
      C,
      sigma,
      perror
      );
    burnin=0.99
    n1trim = n1[Int64(floor(tmax*0.9)):tmax];
    n2trim = n2[Int64(floor(tmax*0.9)):tmax];
    x1trim = x1[Int64(floor(tmax*0.9)):tmax];
    x2trim = x2[Int64(floor(tmax*0.9)):tmax]
    
    n1mean[i,j] = mean(n1trim);
    n2mean[i,j] = mean(n2trim);
    x1mean[i,j] = theta1-mean(x1trim);
    x2mean[i,j] = (theta1+thetadiff)-mean(x2trim);
    
    pe[i,j] = (mean([std(n1trim),std(n2trim)])/mean([mean(n1trim),mean(n2trim)])) *
    (1/(std(n1trim+n2trim)/mean(n1trim+n2trim)))
    
  end

end



nmin=zeros(Float64,length(thetadiffvec),length(a0vec));
nmax=zeros(Float64,length(thetadiffvec),length(a0vec));
for i=1:length(thetadiffvec)
  for j=1:length(a0vec)
    nmin[i,j]=minimum([n1mean[i,j],n2mean[i,j]]);
    nmax[i,j]=maximum([n1mean[i,j],n2mean[i,j]]);
  end
end



namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/fig_SSmtheta.pdf");
R"""
library(RColorBrewer)
cols = brewer.pal(11,'Spectral')
pdf($namespace,height=5,width=6)
plot($(mvec),$(nmin[1,:]),type='l',col=cols[1],xlab="Stray rate",ylab="Steady state",xlim=c(0,0.5),ylim=c(0,max(na.omit(c(as.vector($n1mean),as.vector($n2mean))))),cex=0.5)
lines($(mvec),$(nmax[1,:]),col=cols[1])
legend(x=0.45,y=max(na.omit(c(as.vector($n1mean),as.vector($n2mean)))),legend=$thetadiffvec,col=cols,pch=16,title='Theta diff',bty='n',cex=0.7);
"""
for i=2:length(thetadiffvec)
  n1ss=(nmin[i,:]);
  n2ss=(nmax[i,:]);
  R"""
  lines($(mvec),$(n1ss),col=cols[$i])
  lines($(mvec),$(n2ss),col=cols[$i])
  """
end
R"""
dev.off()
"""


namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/fig_PEmtheta.pdf");
#Portfolio effect
R"""
library(RColorBrewer)
cols = brewer.pal(11,'Spectral')
pdf($namespace,height=5,width=6)
plot($(mvec),$(pe[1,:]),pch=16,col=cols[1],xlab="Stray rate",ylab="Portfolio effect",xlim=c(0,0.5),ylim=c(1,5),cex=0.5)
lines($(mvec),$(pe[1,:]),col=cols[1])
legend(x=0.47,y=4.1,legend=$thetadiffvec,col=cols,pch=16,title='Theta diff',cex=0.5,bty='n');
"""
for j=2:length(thetadiffvec)
  n1pe=(pe[j,:]);
  R"""
  points($(mvec),$(n1pe),pch=16,col=cols[$j],cex=0.5)
  lines($(mvec),$(n1pe),col=cols[$j])
  """
end
R"""
dev.off()
"""








#Evaluating the m value where PE is maximized across
#Habitat heterogeneity (thetadiff)
#sigma (trait variability)



@everywhere using Distributions
@everywhere using RCall
@everywhere using HDF5
@everywhere using JLD

@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolve_ddm.jl")

#Analysis over m & theta divergence
indmvec=collect(0.0:0.001:0.5);
thetadiffvec = collect(5:0.2:10);
sigmavec = collect(0.1:0.1:3.0);

tmax=50000;
z=0.5;
rmax=2.0;
beta=0.001;
theta1=5.0;
tau=1.0;
h=1.0;
C=1000;
sigmaE=0;

perror=0.05;

n1mean=SharedArray(Float64,length(thetadiffvec),length(sigmavec),length(indmvec));
n2mean=SharedArray(Float64,length(thetadiffvec),length(sigmavec),length(indmvec));
x1mean=SharedArray(Float64,length(thetadiffvec),length(sigmavec),length(indmvec));
x2mean=SharedArray(Float64,length(thetadiffvec),length(sigmavec),length(indmvec));
pe=SharedArray(Float64,length(thetadiffvec),length(sigmavec),length(indmvec));



@sync @parallel for i=1:length(thetadiffvec)
  thetadiff = thetadiffvec[i];
  
  for k=1:length(sigmavec)
    sigmaG = sigmavec[k];
    
    for j=1:length(indmvec)
      a0 = 1 - indmvec[j];
      
      n1, n2, x1, x2, w1, w2 = 
      KevinEvolve_ddm(
        tmax,
        z,
        rmax,
        beta,
        theta1,
        thetadiff,
        tau,
        h,
        a0,
        C,
        sigmaE,
        sigmaG,
        perror
        );
      burnin=0.99
      n1trim = n1[Int64(floor(tmax*burnin)):tmax];
      n2trim = n2[Int64(floor(tmax*burnin)):tmax];
      x1trim = x1[Int64(floor(tmax*burnin)):tmax];
      x2trim = x2[Int64(floor(tmax*burnin)):tmax]
      
      n1mean[i,k,j] = mean(n1trim);
      n2mean[i,k,j] = mean(n2trim);
      x1mean[i,k,j] = theta1-mean(x1trim);
      x2mean[i,k,j] = (theta1+thetadiff)-mean(x2trim);
      
      pe[i,k,j] = (mean([std(n1trim),std(n2trim)])/mean([mean(n1trim),mean(n2trim)])) *
      (1/(std(n1trim+n2trim)/mean(n1trim+n2trim)))
    end
  end
end
save(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data/data_td_sig_ddm.jld"),"n1mean",n1mean,"n2mean",n2mean,"x1mean",x1mean,"x2mean",x2mean,"pe",pe);

d = load(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data/data_td_sig_ddm.jld"));
#This loads the dictionary
n1mean = d["n1mean"];
n2mean = d["n2mean"];
x1mean = d["x1mean"];
x2mean = d["x2mean"];
pe = d["pe"];

maxPE = Array{Float64}(length(thetadiffvec),length(sigmavec));
maxPEvalue = Array{Float64}(length(thetadiffvec),length(sigmavec));
for i=1:length(thetadiffvec)
  for k=1:length(sigmavec)
    pecorr = pe[i,k,:];
    pecorr[pecorr.==Inf]=0;
    #Where is the highest PE?
    maxPEvalue[i,k] = maximum(pecorr);
    maxPE[i,k] = indmvec[find(x->x==maximum(pecorr),pecorr)][1];
  end
end

namespace= string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data/data_maxPE_ddm.csv");
R"write.table($maxPE,file=$namespace,row.names=F,col.names=F)"
namespace= string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data/data_maxPEvalue_ddm.csv");
R"write.table($maxPEvalue,file=$namespace,row.names=F,col.names=F)"

R"image($(maxPE))"    
    
R"plot($(maxPE[:,10]))"


R"plot($mvec,$(pe[1,5,:]),ylim=c(0,5))"


qualss = qualsfunc(
n1mean,n2mean,
1.0, #extinct_threshold
30.0 #similarity_threshold=
);

namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/fig_qualthetasig_ddm.pdf");
R"""
library(RColorBrewer)
pal = brewer.pal(3,"Set1")
pdf($namespace,height=4,width=14)
par(mfrow=c(1,6))
image(x=$thetadiffvec,y=$sigmavec,z=$(qualss[:,:,1]),zlim=c(0,2),col=pal,xlab='thetadiff',ylab='sigma_g',main=paste(c('m=',$(indmvec[1]))))
image(x=$thetadiffvec,y=$sigmavec,z=$(qualss[:,:,100]),zlim=c(0,2),col=pal,xlab='thetadiff',ylab='sigma_g',main=paste(c('m=',$(indmvec[100]))))
image(x=$thetadiffvec,y=$sigmavec,z=$(qualss[:,:,200]),zlim=c(0,2),col=pal,xlab='thetadiff',ylab='sigma_g',main=paste(c('m=',$(indmvec[200]))))
image(x=$thetadiffvec,y=$sigmavec,z=$(qualss[:,:,300]),zlim=c(0,2),col=pal,xlab='thetadiff',ylab='sigma_g',main=paste(c('m=',$(indmvec[300]))))
image(x=$thetadiffvec,y=$sigmavec,z=$(qualss[:,:,400]),zlim=c(0,2),col=pal,xlab='thetadiff',ylab='sigma_g',main=paste(c('m=',$(indmvec[400]))))
image(x=$thetadiffvec,y=$sigmavec,z=$(qualss[:,:,500]),zlim=c(0,2),col=pal,xlab='thetadiff',ylab='sigma_g',main=paste(c('m=',$(indmvec[500]))))
legend(x=10,y=2,legend=seq(0,2),col=pal,pch=22,xpd=TRUE,pt.bg=pal, bty="n")
dev.off()
"""

namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/fig_qualthetam_ddm.pdf");
R"""
library(RColorBrewer)
pal = brewer.pal(3,"Set1")
pdf($namespace,height=4,width=14)
par(mfrow=c(1,6))
image(x=$indmvec,y=$sigmavec,z=t($(qualss[1,:,:])),zlim=c(0,2),col=pal,xlab='m',ylab='sigma_G',main=paste(c('thetadiff=',$(thetadiffvec[1]))))
image(x=$indmvec,y=$sigmavec,z=t($(qualss[5,:,:])),zlim=c(0,2),col=pal,xlab='m',ylab='sigma_G',main=paste(c('thetadiff=',$(thetadiffvec[5]))))
image(x=$indmvec,y=$sigmavec,z=t($(qualss[10,:,:])),zlim=c(0,2),col=pal,xlab='m',ylab='sigma_G',main=paste(c('thetadiff=',$(thetadiffvec[10]))))
image(x=$indmvec,y=$sigmavec,z=t($(qualss[15,:,:])),zlim=c(0,2),col=pal,xlab='m',ylab='sigma_G',main=paste(c('thetadiff=',$(thetadiffvec[15]))))
image(x=$indmvec,y=$sigmavec,z=t($(qualss[20,:,:])),zlim=c(0,2),col=pal,xlab='m',ylab='sigma_G',main=paste(c('thetadiff=',$(thetadiffvec[20]))))
image(x=$indmvec,y=$sigmavec,z=t($(qualss[25,:,:])),zlim=c(0,2),col=pal,xlab='m',ylab='sigma_G',main=paste(c('thetadiff=',$(thetadiffvec[25]))))
legend(x=max($indmvec),y=max($sigmavec),legend=seq(0,2),col=pal,pch=22,xpd=TRUE,pt.bg=pal, bty="n")
dev.off()
"""


#Portfolio muli-panel plot
namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/fig_PEpanelthetasig_ddm.pdf");

R"""
library(RColorBrewer)
pal = rev(brewer.pal(9,"Blues"))
pdf($namespace,height=4,width=14)
par(mfrow=c(1,6))
  image(x=$indmvec,y=$sigmavec,z=t($(pe[1,:,:])),zlim=c(1,2),col=pal,xlab='m',ylab='sigma_G',main=paste(c('thetadiff=',$(thetadiffvec[1]))))
image(x=$indmvec,y=$sigmavec,z=t($(pe[5,:,:])),zlim=c(1,2),col=pal,xlab='m',ylab='sigma_G',main=paste(c('thetadiff=',$(thetadiffvec[5]))))
image(x=$indmvec,y=$sigmavec,z=t($(pe[10,:,:])),zlim=c(1,2),col=pal,xlab='m',ylab='sigma_G',main=paste(c('thetadiff=',$(thetadiffvec[10]))))
image(x=$indmvec,y=$sigmavec,z=t($(pe[15,:,:])),zlim=c(1,2),col=pal,xlab='m',ylab='sigma_G',main=paste(c('thetadiff=',$(thetadiffvec[15]))))
image(x=$indmvec,y=$sigmavec,z=t($(pe[20,:,:])),zlim=c(1,2),col=pal,xlab='m',ylab='sigma_G',main=paste(c('thetadiff=',$(thetadiffvec[20]))))
image(x=$indmvec,y=$sigmavec,z=t($(pe[25,:,:])),zlim=c(1,2),col=pal,xlab='m',ylab='sigma_G',main=paste(c('thetadiff=',$(thetadiffvec[25]))))
legend(x=max($indmvec),y=max($sigmavec),legend=seq(1,2,length.out=9),col=pal,pch=22,xpd=TRUE,pt.bg=pal, bty="n")
dev.off()
"""


#Portfolio muli-panel plot
namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/fig_MEANpanelthetasig_ddm.pdf");

R"""
library(RColorBrewer)
pal = rev(brewer.pal(9,"Blues"))
pdf($namespace,height=4,width=14)
par(mfrow=c(1,6))
  image(x=$indmvec,y=$sigmavec,z=(t($(n1mean[1,:,:]))+t($(n2mean[1,:,:]))),zlim=c(1,3000),col=pal,xlab='m',ylab='sigma_G',main=paste(c('thetadiff=',$(thetadiffvec[1]))))
image(x=$indmvec,y=$sigmavec,z=(t($(n1mean[5,:,:]))+t($(n2mean[5,:,:]))),zlim=c(1,3000),col=pal,xlab='m',ylab='sigma_G',main=paste(c('thetadiff=',$(thetadiffvec[5]))))
image(x=$indmvec,y=$sigmavec,z=(t($(n1mean[10,:,:]))+t($(n2mean[10,:,:]))),zlim=c(1,3000),col=pal,xlab='m',ylab='sigma_G',main=paste(c('thetadiff=',$(thetadiffvec[10]))))
image(x=$indmvec,y=$sigmavec,z=(t($(n1mean[15,:,:]))+t($(n2mean[15,:,:]))),zlim=c(1,3000),col=pal,xlab='m',ylab='sigma_G',main=paste(c('thetadiff=',$(thetadiffvec[15]))))
image(x=$indmvec,y=$sigmavec,z=(t($(n1mean[20,:,:]))+t($(n2mean[20,:,:]))),zlim=c(1,3000),col=pal,xlab='m',ylab='sigma_G',main=paste(c('thetadiff=',$(thetadiffvec[20]))))
image(x=$indmvec,y=$sigmavec,z=(t($(n1mean[25,:,:]))+t($(n2mean[25,:,:]))),zlim=c(1,3000),col=pal,xlab='m',ylab='sigma_G',main=paste(c('thetadiff=',$(thetadiffvec[25]))))
legend(x=max($indmvec),y=max($sigmavec),legend=seq(1,3000,length.out=9),col=pal,pch=22,xpd=TRUE,pt.bg=pal, bty="n")
dev.off()
"""




#Evaluating the m value where PE is maximized across
#h (heritability)
#sigma (trait variability)



@everywhere using Distributions
@everywhere using RCall
@everywhere using HDF5
@everywhere using JLD

@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolve_ddm.jl")

#Analysis over m & theta divergence
indmvec=collect(0.0:0.001:0.5);
sigmavec = collect(0.1:0.1:3.0);
hvec = collect(0.1:0.05:1.0);


n1mean=SharedArray(Float64,length(sigmavec),length(hvec),length(indmvec));
n2mean=SharedArray(Float64,length(sigmavec),length(hvec),length(indmvec));
x1mean=SharedArray(Float64,length(sigmavec),length(hvec),length(indmvec));
x2mean=SharedArray(Float64,length(sigmavec),length(hvec),length(indmvec));
pe=SharedArray(Float64,length(sigmavec),length(hvec),length(indmvec));

tmax=50000;
z=0.5;
rmax=2.0;
beta=0.001;
theta1=5.0;
thetadiff=5.0;
tau=1.0;
C=1000;
sigmaE=0;


perror=0.05;

@sync @parallel for i=1:length(sigmavec)
  sigmaG = sigmavec[i];
  
  for k=1:length(hvec)
    h = hvec[k];
    
    for j=1:length(indmvec)
      a0 = 1 - indmvec[j];
      
      n1, n2, x1, x2, w1, w2 = 
      KevinEvolve_ddm(
        tmax,
        z,
        rmax,
        beta,
        theta1,
        thetadiff,
        tau,
        h,
        a0,
        C,
        sigmaE,
        sigmaG,
        perror
        );
      burnin=0.99
      n1trim = n1[Int64(floor(tmax*burnin)):tmax];
      n2trim = n2[Int64(floor(tmax*burnin)):tmax];
      x1trim = x1[Int64(floor(tmax*burnin)):tmax];
      x2trim = x2[Int64(floor(tmax*burnin)):tmax]
      
      n1mean[i,k,j] = mean(n1trim);
      n2mean[i,k,j] = mean(n2trim);
      x1mean[i,k,j] = theta1-mean(x1trim);
      x2mean[i,k,j] = (theta1+thetadiff)-mean(x2trim);
      
      pe[i,k,j] = (mean([std(n1trim),std(n2trim)])/mean([mean(n1trim),mean(n2trim)])) *
      (1/(std(n1trim+n2trim)/mean(n1trim+n2trim)))
    end
  end
end
save(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data/data_sig_h_ddm.jld"),"n1mean",n1mean,"n2mean",n2mean,"x1mean",x1mean,"x2mean",x2mean,"pe",pe);


d = load(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data/data_sig_h_ddm.jld"));
#This loads the dictionary
n1mean = d["n1mean"];
n2mean = d["n2mean"];
x1mean = d["x1mean"];
x2mean = d["x2mean"];
pe = d["pe"];

maxPE = Array{Float64}(length(sigmavec),length(hvec));
maxPEvalue = Array{Float64}(length(sigmavec),length(hvec));
for i=1:length(sigmavec)
  for k=1:length(hvec)
    pecorr = pe[i,k,:];
    pecorr[pecorr.==Inf]=0;
    #Where is the highest PE?
    maxPEvalue[i,k] = maximum(pecorr);
    maxPE[i,k] = indmvec[find(x->x==maximum(pecorr),pecorr)][1];
  end
end

namespace= string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data/data_maxPE_sigh_ddm.csv");
R"write.table($maxPE,file=$namespace,row.names=F,col.names=F)"
namespace= string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data/data_maxPEvalue_sigh_ddm.csv");
R"write.table($maxPEvalue,file=$namespace,row.names=F,col.names=F)"


R"image($(maxPE))"    
    
R"plot($(maxPE[:,10]))"


R"plot($mvec,$(pe[1,5,:]),ylim=c(0,5))"


qualss = qualsfunc(
n1mean,n2mean,
1.0, #extinct_threshold
30.0 #similarity_threshold=
);

namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/fig_qualhsig_ddm.pdf");
R"""
library(RColorBrewer)
pal = brewer.pal(3,"Set1")
#pdf($namespace,height=4,width=14)
par(mfrow=c(1,6))
image(x=$sigmavec,y=$hvec,z=$(qualss[:,:,1]),zlim=c(0,2),col=pal,xlab='sigma_g',ylab='heritability',main=paste(c('m=',$(mvec[1]))))
image(x=$sigmavec,y=$hvec,z=$(qualss[:,:,100]),zlim=c(0,2),col=pal,xlab='sigma_g',ylab='heritability',main=paste(c('m=',$(mvec[100]))))
image(x=$sigmavec,y=$hvec,z=$(qualss[:,:,200]),zlim=c(0,2),col=pal,xlab='sigma_g',ylab='heritability',main=paste(c('m=',$(mvec[200]))))
image(x=$sigmavec,y=$hvec,z=$(qualss[:,:,300]),zlim=c(0,2),col=pal,xlab='sigma_g',ylab='heritability',main=paste(c('m=',$(mvec[300]))))
image(x=$sigmavec,y=$hvec,z=$(qualss[:,:,400]),zlim=c(0,2),col=pal,xlab='sigma_g',ylab='heritability',main=paste(c('m=',$(mvec[400]))))
image(x=$sigmavec,y=$hvec,z=$(qualss[:,:,450]),zlim=c(0,2),col=pal,xlab='sigma_g',ylab='heritability',main=paste(c('m=',$(mvec[450]))))
legend(x=max($sigmavec),y=max($hvec),legend=seq(0,2),col=pal,pch=22,xpd=TRUE,pt.bg=pal, bty="n")
#dev.off()
"""


namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/fig_qualmsig_ddm.pdf");
R"""
library(RColorBrewer)
pal = brewer.pal(3,"Set1")
pdf($namespace,height=4,width=12)
par(mfrow=c(1,5))
image(x=$indmvec,y=$sigmavec,z=t($(qualss[:,1,:])),zlim=c(0,2),col=pal,xlab='m',ylab='sigma_G',main=paste(c('h=',$(hvec[1]))))
image(x=$indmvec,y=$sigmavec,z=t($(qualss[:,5,:])),zlim=c(0,2),col=pal,xlab='m',ylab='sigma_G',main=paste(c('h=',$(hvec[5]))))
image(x=$indmvec,y=$sigmavec,z=t($(qualss[:,10,:])),zlim=c(0,2),col=pal,xlab='m',ylab='sigma_G',main=paste(c('h=',$(hvec[10]))))
image(x=$indmvec,y=$sigmavec,z=t($(qualss[:,15,:])),zlim=c(0,2),col=pal,xlab='m',ylab='sigma_G',main=paste(c('h=',$(hvec[15]))))
image(x=$indmvec,y=$sigmavec,z=t($(qualss[:,19,:])),zlim=c(0,2),col=pal,xlab='m',ylab='sigma_G',main=paste(c('h=',$(hvec[19]))))
legend(x=max($indmvec),y=max($sigmavec),legend=seq(0,2),col=pal,pch=22,xpd=TRUE,pt.bg=pal, bty="n")
dev.off()
"""



#Evaluating the m value where PE is maximized across
#Habitat heterogeneity (thetadiff)
#h (trait heritability)



@everywhere using Distributions
@everywhere using RCall
@everywhere using HDF5
@everywhere using JLD

@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolve_ddm.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/qualsfunc.jl")

#Analysis over m & theta divergence
indmvec=collect(0.0:0.001:0.5);
thetadiffvec = collect(5:0.1:10);
hvec = collect(0.1:0.1:1.0);


n1mean=SharedArray(Float64,length(thetadiffvec),length(hvec),length(indmvec));
n2mean=SharedArray(Float64,length(thetadiffvec),length(hvec),length(indmvec));
x1mean=SharedArray(Float64,length(thetadiffvec),length(hvec),length(indmvec));
x2mean=SharedArray(Float64,length(thetadiffvec),length(hvec),length(indmvec));
pe=SharedArray(Float64,length(thetadiffvec),length(hvec),length(indmvec));

tmax=50000;
z=0.5;
rmax=2.0;
beta=0.001;
theta1=5.0;
tau=1.0;
C=1000;
sigmaE=0.0;
sigmaG=1.0;

perror=0.05;

@sync @parallel for i=1:length(thetadiffvec)
  thetadiff = thetadiffvec[i];
  
  for k=1:length(hvec)
    h = hvec[k];
    
    for j=1:length(indmvec)
      a0 = 1 - indmvec[j];
      
      n1, n2, x1, x2, w1, w2 = 
      KevinEvolve_ddm(
        tmax,
        z,
        rmax,
        beta,
        theta1,
        thetadiff,
        tau,
        h,
        a0,
        C,
        sigmaE,
        sigmaG,
        perror
        );
      burnin=0.99
      n1trim = n1[Int64(floor(tmax*burnin)):tmax];
      n2trim = n2[Int64(floor(tmax*burnin)):tmax];
      x1trim = x1[Int64(floor(tmax*burnin)):tmax];
      x2trim = x2[Int64(floor(tmax*burnin)):tmax]
      
      n1mean[i,k,j] = mean(n1trim);
      n2mean[i,k,j] = mean(n2trim);
      x1mean[i,k,j] = theta1-mean(x1trim);
      x2mean[i,k,j] = (theta1+thetadiff)-mean(x2trim);
      
      pe[i,k,j] = (mean([std(n1trim),std(n2trim)])/mean([mean(n1trim),mean(n2trim)])) *
      (1/(std(n1trim+n2trim)/mean(n1trim+n2trim)))
    end
  end
end
save(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data/data_td_h_ddm.jld"),"n1mean",n1mean,"n2mean",n2mean,"x1mean",x1mean,"x2mean",x2mean,"pe",pe);


d = load(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data/data_td_h_ddm.jld"));
#This loads the dictionary
n1mean = d["n1mean"];
n2mean = d["n2mean"];
x1mean = d["x1mean"];
x2mean = d["x2mean"];
pe = d["pe"];

maxPE = Array{Float64}(length(thetadiffvec),length(hvec));
maxPEvalue = Array{Float64}(length(thetadiffvec),length(hvec));
for i=1:length(thetadiffvec)
  for k=1:length(hvec)
    pecorr = pe[i,k,:];
    pecorr[pecorr.==Inf]=0;
    #Where is the highest PE?
    maxPEvalue[i,k] = maximum(pecorr);
    maxPE[i,k] = indmvec[find(x->x==maximum(pecorr),pecorr)][1];
  end
end

namespace= string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data/data_maxPE_h_ddm.csv");
R"write.table($maxPE,file=$namespace,row.names=F,col.names=F)"
namespace= string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data/data_maxPEvalue_h_ddm.csv");
R"write.table($maxPEvalue,file=$namespace,row.names=F,col.names=F)"


R"image($(maxPE))"    
    
R"plot($(maxPE[:,10]))"


R"plot($mvec,$(pe[1,5,:]),ylim=c(0,5))"


qualss = qualsfunc(
n1mean,n2mean,
1.0, #extinct_threshold
30.0 #similarity_threshold=
);

namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/fig_qualthetah_ddm.pdf");
R"""
library(RColorBrewer)
pal = brewer.pal(3,"Set1")
pdf($namespace,height=4,width=14)
par(mfrow=c(1,6))
image(x=$thetadiffvec,y=$hvec,z=$(qualss[:,:,1]),zlim=c(0,2),col=pal,xlab='thetadiff',ylab='heritability',main=paste(c('m=',$(indmvec[1]))))
image(x=$thetadiffvec,y=$hvec,z=$(qualss[:,:,100]),zlim=c(0,2),col=pal,xlab='thetadiff',ylab='heritability',main=paste(c('m=',$(indmvec[100]))))
image(x=$thetadiffvec,y=$hvec,z=$(qualss[:,:,200]),zlim=c(0,2),col=pal,xlab='thetadiff',ylab='heritability',main=paste(c('m=',$(indmvec[200]))))
image(x=$thetadiffvec,y=$hvec,z=$(qualss[:,:,300]),zlim=c(0,2),col=pal,xlab='thetadiff',ylab='heritability',main=paste(c('m=',$(indmvec[300]))))
image(x=$thetadiffvec,y=$hvec,z=$(qualss[:,:,400]),zlim=c(0,2),col=pal,xlab='thetadiff',ylab='heritability',main=paste(c('m=',$(indmvec[400]))))
image(x=$thetadiffvec,y=$hvec,z=$(qualss[:,:,450]),zlim=c(0,2),col=pal,xlab='thetadiff',ylab='heritability',main=paste(c('m=',$(indmvec[450]))))
legend(x=10,y=2,legend=seq(0,2),col=pal,pch=22,xpd=TRUE,pt.bg=pal, bty="n")
dev.off()
"""
