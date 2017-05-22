using Distributions
using RCall

include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolve.jl")
include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinJacobian.jl")


tmax=10000;
z=0.5;
rmax=2.0;
beta=0.001;
theta1=5.0;
thetadiff=5.0;
tau=1.0;
h=0.49;
sigmaE=0;
sigmaG=1.0;
m=0.35;
perror=0.01;
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
  m, #m
  perror
  );
  
  R"""
  par(mfrow=c(1,2))
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
thetadiff=5;
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

R"""
library(RColorBrewer)
cols = brewer.pal(3,'Set1')
plot($mvec,$n1mean,pch='.',col=cols[1],xlab="Stray rate",ylab="Steady state",cex=0.5)
points($mvec,$n2mean,pch='.',col=cols[2],cex=0.5)
"""


#Steady state plot
namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/figs/fig_SSm.pdf");
R"""
library(RColorBrewer)
cols = brewer.pal(3,'Set1')
pdf($namespace,height=5,width=6)
plot($mvec,$n1mean,pch='.',col=cols[1],xlab="Stray rate",ylab="Steady state",cex=0.5)
points($mvec,$n2mean,pch='.',col=cols[2],cex=0.5)
dev.off()
"""
#Trait offset plot
namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/figs/fig_Traitm.pdf");
R"""
library(RColorBrewer)
cols = brewer.pal(3,'Set1')
pdf($namespace,height=5,width=6)
plot($mvec,$x1mean,pch=".",col=cols[1],ylim=c(-5,5),xlab="Stray rate",ylab="Trait offset")
points($mvec,$x2mean,pch=".",col=cols[2])
dev.off()
"""
#Portfolio effect plot
namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/figs/fig_pem.pdf");
R"""
library(RColorBrewer)
cols = brewer.pal(3,'Set1')
pdf($namespace,height=5,width=6)
plot($mvec,$pe,pch='.',col="black",xlab="Stray rate",ylab="Portfolio effect")
dev.off()
"""

#Portfolio effect plot
R"""
plot($mvec,$pe)
"""

#Portfolio effect plot in pieces
R"""
par(mfrow=c(2,1))
plot($mvec,$n1sd/$aggsd)
points($mvec,$n2sd/$aggsd)
plot($mvec,$aggmean/mean($n1mean))
points($mvec,$aggmean/mean($n2mean))
"""

#Portfolio effect plot
R"""
plot($mvec,$n1sd/$n1mean)
points($mvec,$n2sd/$n2mean,col='blue')
points($mvec,$aggsd/$aggmean,pch=16)
"""


R"""
par(mfrow=c(2,1))
plot($mvec,$aggmean,pch=16)
points($mvec,$n1mean,col='blue')
points($mvec,$n2mean)
plot($mvec,$n1sd,ylim=c(0,7))
points($mvec,$n2sd,col='blue')
points($mvec,$aggsd,pch=16)
"""




R"""
plot($n1mean[1:500],($n1sd[1:500])^2,log='xy')
points($n2mean[1:500],($n2sd[1:500])^2)
"""


#Portfolio effect plot
R"""
plot($mvec,$n1sd/$n1mean)
points($mvec,$n2sd/$n2mean,col='blue')
points($mvec,$aggsd/$aggmean,pch=16)
"""
n1trim = n1mean[200:500];
n1sdtrim = n1sd[200:500];
lmean = n1trim[find(x->x<3.1,n1sdtrim)];
lsd = n1sdtrim[find(x->x<3.1,n1sdtrim)];
umean = n1trim[find(x->x>3.1,n1sdtrim)];
usd = n1sdtrim[find(x->x>3.1,n1sdtrim)];

R"""
plot(log($(lmean)),($lsd))
m1 = lm($lsd ~ log($lmean));
summary(m1)
"""

R"""
plot(log($(umean)),($usd))
m1 = lm($usd ~ log($umean));
summary(m1)
"""

R"""
par(mfrow=c(2,1))
plot($mvec,$aggmean,pch=16)
points($mvec,$n1mean,col='blue')
points($mvec,$n2mean)
"""
R"""
plot($mvec,$n1sd,ylim=c(0,7))
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
namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/figs/fig_SSm2.pdf");
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



####################


#Evaluating the m value where PE is maximized across
#h (heritability)
#sigma (trait variability)



@everywhere using Distributions
@everywhere using RCall
@everywhere using HDF5
@everywhere using JLD

@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolve.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/qualsfunc.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/bifdet.jl")


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
thetadiff=5.0;
tau=1.0;
sigmaE=0;
sigmaG=1;

perror=0.01;

@sync @parallel for k=1:length(hvec)
    
    h = hvec[k];
    
    for j=1:length(mvec)
      m = mvec[j];
      
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
save(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data2/data_sig_h_m.jld"),"n1mean",n1mean,"n2mean",n2mean,"x1mean",x1mean,"x2mean",x2mean,"pe",pe);


d = load(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data2/data_sig_h_m.jld"));
#This loads the dictionary
n1mean = d["n1mean"];
n2mean = d["n2mean"];
x1mean = d["x1mean"];
x2mean = d["x2mean"];
pe = d["pe"];
dExt = load(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data2Ext/data_sig_h_m.jld"));
pe_ext = dExt["pe"];
rt_ext = dExt["rt"];

@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/bifdet.jl")
bifvalue = bifdet(
n1mean,
n2mean,
mvec,
hvec
);


namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/figs2/fig_MDPE_hm.pdf");
R"""
library(RColorBrewer)
pal = rev(brewer.pal(9,"Blues"))
pdf($namespace,height=3,width=10)
par(mfrow=c(1,4))
image(x=$mvec,y=$hvec,z=t($(n1mean[:,:]))+t($(n2mean[:,:])),zlim=c(0,3000),col=pal,xlab='m',ylab='h',main='Total biomass')
points($(bifvalue),type='l',cex=1)
image(x=$mvec,y=$hvec,z=sqrt((t($(n1mean[:,:]))-t($(n2mean[:,:])))^2),zlim=c(0,1500),col=pal,xlab='m',ylab='h',main='Biomass difference')
points($(bifvalue),type='l',cex=1)
image(x=$mvec,y=$hvec,z=t($(pe[:,:])),zlim=c(1,2),col=pal,xlab='m',ylab='h',main='PE')
points($(bifvalue),type='l',cex=1)
plot($(rt_ext[1:50,:]),$(pe_ext[1:50,:]),log='x',ylim=c(1,3),xlim=c(18,300),pch='.',xlab='Recovery time',ylab='PE')
dev.off()
"""

R"""
x<- 1:10
  y<- 1:15
  z<- outer( x,y,"+") 
  image.plot(x,y,z) 

# or 
  obj<- list( x=x,y=y,z=z)
  image.plot(obj, legend.lab="Sverdrups")

# add some points on diagonal using standard plot function
#(with some clipping beyond 10 anticipated)
  points( 5:12, 5:12, pch="X", cex=3)
  """



#Average pe over m
medpe = Array{Float64}(length(mvec));
pena = pe;
pena[find(x->x==true,isnan(pe))] = 1;
for i=1:length(mvec)
  medpe[i] = median(pena[:,i])
end


namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/figs2/fig_traitdiff.pdf");
traitdiff = abs(-1*(x1mean[:,:]-theta1) - (-1*(x2mean[:,:]-(theta1+thetadiff))));
R"""
pdf($namespace,height=4,width=5)
plot($mvec,$(traitdiff[26,:]),type='l',xlim=c(0,0.5),ylab=expression(paste('Phenotypic diversity, ',Delta,mu,'*')),xlab='m')
lines($mvec,$(traitdiff[51,:]),type='l')
lines($mvec,$(traitdiff[76,:]),type='l')
lines($mvec,$(traitdiff[101,:]),type='l')
text(rep(0.49,4),$(traitdiff[[26 51 76 101],451]),c('0.25','0.50','0.75','1.00'),cex=0.8)
text(0.48,2.1,expression(paste(h^2)),cex=0.8)
dev.off()
"""




#THETADIFF = 8

@everywhere using Distributions
@everywhere using RCall
@everywhere using HDF5
@everywhere using JLD

@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolve.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/qualsfunc.jl")


#Analysis over m & theta divergence
mvec=collect(0.0:0.001:0.45);
hvec = collect(0.0:0.01:1.0);


n1mean=SharedArray(Float64,length(hvec),length(mvec));
n2mean=SharedArray(Float64,length(hvec),length(mvec));
x1mean=SharedArray(Float64,length(hvec),length(mvec));
x2mean=SharedArray(Float64,length(hvec),length(mvec));
pe=SharedArray(Float64,length(hvec),length(mvec));

tmax=100000;
z=0.5;
rmax=2.0;
beta=0.001;
theta1=5.0;
thetadiff=8.0;
tau=1.0;
sigmaE=0;
sigmaG=1;

perror=0.01;

@sync @parallel for k=1:length(hvec)
  h = hvec[k];
  
  for j=1:length(mvec)
    m = mvec[j];
    
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
save(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data2/data_sig_h_m_theta8.jld"),"n1mean",n1mean,"n2mean",n2mean,"x1mean",x1mean,"x2mean",x2mean,"pe",pe);


d = load(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data2/data_sig_h_m_theta8.jld"));
#This loads the dictionary
n1mean = d["n1mean"];
n2mean = d["n2mean"];
x1mean = d["x1mean"];
x2mean = d["x2mean"];
pe = d["pe"];


bifvalue = bifdet(
n1mean,
n2mean,
mvec,
hvec
);

namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/figs2/fig_MDPE_hm_theta8.pdf");
R"""
library(RColorBrewer)
pal = rev(brewer.pal(9,"Blues"))
pdf($namespace,height=3,width=8)
par(mfrow=c(1,3))
image(x=$mvec,y=$hvec,z=t($(n1mean[:,:]))+t($(n2mean[:,:])),zlim=c(0,3000),col=pal,xlab='m',ylab='h',main='Total biomass')
points($bifvalue,type='l')
image(x=$mvec,y=$hvec,z=sqrt((t($(n1mean[:,:]))-t($(n2mean[:,:])))^2),zlim=c(0,1500),col=pal,xlab='m',ylab='h',main='Biomass difference')
points($bifvalue,type='l')
image(x=$mvec,y=$hvec,z=t($(pe[:,:])),zlim=c(1,2),col=pal,xlab='m',ylab='h',main='PE')
points($bifvalue,type='l')
dev.off()
"""

#Average pe over m
medpe = Array{Float64}(length(mvec));
pena = pe;
pena[find(x->x==true,isnan(pe))] = 1;
for i=1:length(mvec)
  medpe[i] = median(pena[:,i])
end


namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/figs2/fig_traitdiff_theta8.pdf");
traitdiff = abs(-1*(x1mean[:,:]-theta1) - (-1*(x2mean[:,:]-(theta1+thetadiff))));
R"""
pdf($namespace,height=4,width=5)
plot($mvec,$(traitdiff[26,:]),type='l',xlim=c(0,0.5),ylab=expression(paste('Phenotypic diversity, ',Delta,mu)),xlab='m')
lines($mvec,$(traitdiff[51,:]),type='l')
lines($mvec,$(traitdiff[76,:]),type='l')
lines($mvec,$(traitdiff[101,:]),type='l')
text(rep(0.49,4),$(traitdiff[[26 51 76 101],451]),c('0.25','0.50','0.75','1.00'),cex=0.8)
text(0.48,2.1,expression(paste(h^2)),cex=0.8)
dev.off()
"""


#THETADIFF = 3

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
thetadiff=3.0;
tau=1.0;
sigmaE=0;
sigmaG=1;

perror=0.01;

@sync @parallel for k=1:length(hvec)
    
  h = hvec[k];
  
  for j=1:length(mvec)
    m = mvec[j];
    
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
save(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data2/data_sig_h_m_theta3.jld"),"n1mean",n1mean,"n2mean",n2mean,"x1mean",x1mean,"x2mean",x2mean,"pe",pe);


d = load(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data2/data_sig_h_m_theta3.jld"));
#This loads the dictionary
n1mean = d["n1mean"];
n2mean = d["n2mean"];
x1mean = d["x1mean"];
x2mean = d["x2mean"];
pe = d["pe"];



bifvalue = bifdet(
n1mean,
n2mean,
mvec,
hvec
);

namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/figs2/fig_MDPE_hm_theta3.pdf");
R"""
library(RColorBrewer)
pal = rev(brewer.pal(9,"Blues"))
pdf($namespace,height=3,width=8)
par(mfrow=c(1,3))
image(x=$mvec,y=$hvec,z=t($(n1mean[:,:]))+t($(n2mean[:,:])),zlim=c(0,3000),col=pal,xlab='m',ylab='h',main='Total biomass')
points($(bifvalue[1:23,:]),type='l')
image(x=$mvec,y=$hvec,z=sqrt((t($(n1mean[:,:]))-t($(n2mean[:,:])))^2),zlim=c(0,1500),col=pal,xlab='m',ylab='h',main='Biomass difference')
points($(bifvalue[1:23,:]),type='l')
image(x=$mvec,y=$hvec,z=t($(pe[:,:])),zlim=c(1,2),col=pal,xlab='m',ylab='h',main='PE')
points($(bifvalue[1:23,:]),type='l')
dev.off()
"""

#Average pe over m
medpe = Array{Float64}(length(mvec));
pena = pe;
pena[find(x->x==true,isnan(pe))] = 1;
for i=1:length(mvec)
  medpe[i] = median(pena[:,i])
end


namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/figs2/fig_traitdiff_theta3.pdf");
traitdiff = abs(-1*(x1mean[:,:]-theta1) - (-1*(x2mean[:,:]-(theta1+thetadiff))));
R"""
pdf($namespace,height=4,width=5)
plot($mvec,$(traitdiff[26,:]),type='l',xlim=c(0,0.5),ylab=expression(paste('Phenotypic diversity, ',Delta,mu)),xlab='m')
lines($mvec,$(traitdiff[51,:]),type='l')
lines($mvec,$(traitdiff[76,:]),type='l')
lines($mvec,$(traitdiff[101,:]),type='l')
text(rep(0.49,4),$(traitdiff[[26 51 76 101],451]),c('0.25','0.50','0.75','1.00'),cex=0.8)
text(0.48,2.1,expression(paste(h^2)),cex=0.8)
dev.off()
"""
