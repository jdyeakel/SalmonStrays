@everywhere using Distributions
@everywhere using RCall
@everywhere using HDF5
@everywhere using JLD

@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolveExtinct.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/qualsfunc.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/timeSS.jl")


#Analysis over m & theta divergence
mvec=collect(0.0:0.001:0.45);
hvec = collect(0.0:0.01:1.0);


n1mean=SharedArray(Float64,length(hvec),length(mvec));
n2mean=SharedArray(Float64,length(hvec),length(mvec));
x1mean=SharedArray(Float64,length(hvec),length(mvec));
x2mean=SharedArray(Float64,length(hvec),length(mvec));
pe=SharedArray(Float64,length(hvec),length(mvec));
rt=SharedArray(Float64,length(hvec),length(mvec));


tmax=10000;
z=0.5;
rmax=2.0;
beta=0.001;
theta1=5.0;
thetadiff=5.0;
tau=1.0;
sigmaE=0;
sigmaG=1;
extpop="both";
refuge=0.01;
t_ext = Int64(round(tmax/2));
perror=0.01;

@sync @parallel for k=1:length(hvec)
    
    h = hvec[k];
    
    for j=1:length(mvec)
      m = mvec[j];
      
      n1, n2, x1, x2, w1, w2 = 
      KevinEvolveExtinct(tmax, 
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
      perror,
      extpop,
      t_ext,
      refuge
      );
      
      t_ss, relaxtime = timeSS(n1,n2,t_ext);
      
      rt[k,j] = relaxtime;
      
      burnin=0.9
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
save(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data2Ext/data_sig_h_m.jld"),"n1mean",n1mean,"n2mean",n2mean,"x1mean",x1mean,"x2mean",x2mean,"pe",pe,"rt",rt);


d = load(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data2Ext/data_sig_h_m.jld"));
#This loads the dictionary
n1mean = d["n1mean"];
n2mean = d["n2mean"];
x1mean = d["x1mean"];
x2mean = d["x2mean"];
pe = d["pe"];
rt = d["rt"];


namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/figs2/fig_ExtMDPE_hm.pdf");
R"""
library(RColorBrewer)
pal = rev(brewer.pal(9,"Blues"))
pdf($namespace,height=3,width=12)
par(mfrow=c(1,4))
image(x=$mvec,y=$hvec,z=t($(n1mean[:,:]))+t($(n2mean[:,:])),zlim=c(0,3000),col=pal,xlab='m',ylab='h',main='Total biomass')
image(x=$mvec,y=$hvec,z=sqrt((t($(n1mean[:,:]))-t($(n2mean[:,:])))^2),zlim=c(0,1500),col=pal,xlab='m',ylab='h',main='Biomass difference')
image(x=$mvec,y=$hvec,z=t($(pe[:,:])),zlim=c(1,2),col=pal,xlab='m',ylab='h',main='PE')
image(x=$mvec,y=$hvec,z=log(t($(rt[:,:]))),zlim=c(log(15),log(100)),col=pal,xlab='m',ylab='h',main='Return time')
dev.off()
"""

namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/Manuscript/figs2/fig_rtvspe.pdf");
R"""
library(RColorBrewer)
pal = rev(brewer.pal(9,"Blues"))
pdf($namespace,height=4,width=5)
plot($(pe[1:50,:]),$(rt[1:50,:]),log='y',xlim=c(1,4),ylim=c(18,200),pch='.',xlab='PE',ylab='Return time')
dev.off()
"""

#THETADIFF = 8


n1mean=SharedArray(Float64,length(hvec),length(mvec));
n2mean=SharedArray(Float64,length(hvec),length(mvec));
x1mean=SharedArray(Float64,length(hvec),length(mvec));
x2mean=SharedArray(Float64,length(hvec),length(mvec));
pe=SharedArray(Float64,length(hvec),length(mvec));
rt=SharedArray(Float64,length(hvec),length(mvec));

tmax=100000;
z=0.5;
rmax=2.0;
beta=0.001;
theta1=5.0;
thetadiff=8.0;
tau=1.0;
sigmaE=0;
sigmaG=1;
popext="both";
refuge=0.01;
t_ext = Int64(round(tmax/2));
perror=0.01;

@sync @parallel for k=1:length(hvec)
  h = hvec[k];
  
  for j=1:length(mvec)
    m = mvec[j];
    
    n1, n2, x1, x2, w1, w2 = 
    KevinEvolveExtinct(tmax, 
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
    perror,
    extpop,
    t_ext,
    refuge
    );
    
    t_ss, relaxtime = timeSS(n1,n2,t_ext);
    
    rt[k,j] = relaxtime;
    
    burnin=0.9
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
save(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data2Ext/data_sig_h_m_theta8.jld"),"n1mean",n1mean,"n2mean",n2mean,"x1mean",x1mean,"x2mean",x2mean,"pe",pe,"rt",rt);


d = load(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data2Ext/data2_sig_h_m_theta8.jld"));
#This loads the dictionary
n1mean = d["n1mean"];
n2mean = d["n2mean"];
x1mean = d["x1mean"];
x2mean = d["x2mean"];
pe = d["pe"];
rt = d["rt"];



namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/figs2/fig_ExtMDPE_hm_theta8.pdf");
R"""
library(RColorBrewer)
pal = rev(brewer.pal(9,"Blues"))
pdf($namespace,height=3,width=12)
par(mfrow=c(1,4))
image(x=$mvec,y=$hvec,z=t($(n1mean[:,:]))+t($(n2mean[:,:])),zlim=c(0,3000),col=pal,xlab='m',ylab='h',main='Total biomass')
image(x=$mvec,y=$hvec,z=sqrt((t($(n1mean[:,:]))-t($(n2mean[:,:])))^2),zlim=c(0,1500),col=pal,xlab='m',ylab='h',main='Biomass difference')
image(x=$mvec,y=$hvec,z=t($(pe[:,:])),zlim=c(1,2),col=pal,xlab='m',ylab='h',main='PE')
image(x=$mvec,y=$hvec,z=log(t($(rt[:,:]))),zlim=c(log(15),log(100)),col=pal,xlab='m',ylab='h',main='Return time')
dev.off()
"""

#THETADIFF = 3


n1mean=SharedArray(Float64,length(hvec),length(mvec));
n2mean=SharedArray(Float64,length(hvec),length(mvec));
x1mean=SharedArray(Float64,length(hvec),length(mvec));
x2mean=SharedArray(Float64,length(hvec),length(mvec));
pe=SharedArray(Float64,length(hvec),length(mvec));
rt=SharedArray(Float64,length(hvec),length(mvec));

tmax=10000;
z=0.5;
rmax=2.0;
beta=0.001;
theta1=5.0;
thetadiff=3.0;
tau=1.0;
sigmaE=0;
sigmaG=1;
popext="both";
refuge=0.01;
t_ext = Int64(round(tmax/2));
perror=0.01;

@sync @parallel for k=1:length(hvec)
    
  h = hvec[k];
  
  for j=1:length(mvec)
    m = mvec[j];
    
    n1, n2, x1, x2, w1, w2 = 
    KevinEvolveExtinct(tmax, 
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
    perror,
    extpop,
    t_ext,
    refuge
    );
    
    t_ss, relaxtime = timeSS(n1,n2,t_ext);
    
    rt[k,j] = relaxtime;
    
    burnin=0.9
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
save(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data2Ext/data_sig_h_m_theta3.jld"),"n1mean",n1mean,"n2mean",n2mean,"x1mean",x1mean,"x2mean",x2mean,"pe",pe,"rt",rt);


d = load(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data2Ext/data_sig_h_m_theta3.jld"));
#This loads the dictionary
n1mean = d["n1mean"];
n2mean = d["n2mean"];
x1mean = d["x1mean"];
x2mean = d["x2mean"];
pe = d["pe"];
rt = d["rt"];


namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/figs2/fig_ExtMDPE_hm_theta3.pdf");
R"""
library(RColorBrewer)
pal = rev(brewer.pal(9,"Blues"))
pdf($namespace,height=3,width=12)
par(mfrow=c(1,4))
image(x=$mvec,y=$hvec,z=t($(n1mean[:,:]))+t($(n2mean[:,:])),zlim=c(0,3000),col=pal,xlab='m',ylab='h',main='Total biomass')
image(x=$mvec,y=$hvec,z=sqrt((t($(n1mean[:,:]))-t($(n2mean[:,:])))^2),zlim=c(0,1500),col=pal,xlab='m',ylab='h',main='Biomass difference')
image(x=$mvec,y=$hvec,z=t($(pe[:,:])),zlim=c(1,2),col=pal,xlab='m',ylab='h',main='PE')
image(x=$mvec,y=$hvec,z=log(t($(rt[:,:]))),zlim=c(log(15),log(100)),col=pal,xlab='m',ylab='h',main='Return time')
dev.off()
"""
