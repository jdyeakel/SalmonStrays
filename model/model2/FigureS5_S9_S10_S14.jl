using RCall
using Distributions
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolveSS.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinJacobian.jl")


#Analysis over m
tmax=10000;
mvec1 = collect(0.0:0.0001:0.5);
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
thetadiff=2;
tau=1;
h=0.2;
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
namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/FinalDraft3/fig_traj.pdf");
R"""
library(RColorBrewer)
cols = brewer.pal(3,'Set1')
pdf($namespace,height=5,width=6)
plot($(mvec[1:midpoint]),$(n1mean[1:midpoint]),col=cols[1],xlab="Straying (m)",ylab="Steady state biomass",xlim=c(0,0.5),ylim=c(1,400),type='l')
points($(mvec[1:midpoint]),$(n2mean[1:midpoint]),col=cols[2])
points($(mvec[midpoint+1:length(mvec)]),$(n1mean[midpoint+1:length(n1mean)]),col=cols[1],lty=2)
points($(mvec[midpoint+1:length(mvec)]),$(n2mean[midpoint+1:length(n2mean)]),col=cols[2],lty=2)
text(0.04,380,'RI', xpd=TRUE)
text(0.18,380,'RII', xpd=TRUE)
types = c('Increasing m','Decreasing m')
#legend(x=0.1,y=500,legend=types,col='black',lty=c(1,2),xpd=TRUE,cex=0.9, bty="n")
dev.off()
"""





@everywhere using Distributions
@everywhere using RCall
@everywhere using HDF5
@everywhere using JLD

@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolve.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolve_ddm.jl")


@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolveExtinct.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolveExtinct_ddm.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/timeSS.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/movingaverage.jl")


#Analysis over m & theta divergence
mvec = collect(0.0001:0.001:0.3);
thetadiffvec = collect(0.0:0.1:6.0);
hvec = collect(0.1:0.1:0.8);
reps = 50;

thetadiffarray=SharedArray(Float64,reps,length(thetadiffvec),length(hvec),length(mvec));

n1mean=SharedArray(Float64,reps,length(thetadiffvec),length(hvec),length(mvec));
n2mean=SharedArray(Float64,reps,length(thetadiffvec),length(hvec),length(mvec));
x1mean=SharedArray(Float64,reps,length(thetadiffvec),length(hvec),length(mvec));
x2mean=SharedArray(Float64,reps,length(thetadiffvec),length(hvec),length(mvec));
pe=SharedArray(Float64,reps,length(thetadiffvec),length(hvec),length(mvec));
n1mean_ddm=SharedArray(Float64,reps,length(thetadiffvec),length(hvec),length(mvec));
n2mean_ddm=SharedArray(Float64,reps,length(thetadiffvec),length(hvec),length(mvec));
x1mean_ddm=SharedArray(Float64,reps,length(thetadiffvec),length(hvec),length(mvec));
x2mean_ddm=SharedArray(Float64,reps,length(thetadiffvec),length(hvec),length(mvec));
pe_ddm=SharedArray(Float64,reps,length(thetadiffvec),length(hvec),length(mvec));

rt=SharedArray(Float64,reps,length(thetadiffvec),length(hvec),length(mvec));
rt_ddm=SharedArray(Float64,reps,length(thetadiffvec),length(hvec),length(mvec));
m1mean=SharedArray(Float64,reps,length(thetadiffvec),length(hvec),length(mvec));
m2mean=SharedArray(Float64,reps,length(thetadiffvec),length(hvec),length(mvec));

tmax=10000;
z=0.5;
rmax=2.0;
beta=0.001;
theta1=5.0;
tau=1.0;
sigmaE=0;
sigmaG=1;
perror=0.01;
refuge=0.01;
extpop="both";
burnin = 0.8;
C=1000;

# n1v = SharedArray{Float64}(reps,length(mvec),length(pvec),tmax);
# n2v = SharedArray{Float64}(reps,length(mvec),length(pvec),tmax);
# x1v = SharedArray{Float64}(reps,length(mvec),length(pvec),tmax);
# x2v = SharedArray{Float64}(reps,length(mvec),length(pvec),tmax);

t_ext = Int64(round(tmax/2));


@sync @parallel for r=1:reps

  for k=1:length(thetadiffvec)
      
      thetadiff=thetadiffvec[k];
      
      for i=1:length(hvec)
        
        h = hvec[i];
      
        for j=1:length(mvec)
          m = mvec[j];
          a0 = 1-m;
          
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
          n1trim = n1[Int64(floor(tmax*burnin)):tmax];
          n2trim = n2[Int64(floor(tmax*burnin)):tmax];
          x1trim = x1[Int64(floor(tmax*burnin)):tmax];
          x2trim = x2[Int64(floor(tmax*burnin)):tmax];
          
          n1mean[r,k,i,j] = mean(n1trim);
          n2mean[r,k,i,j] = mean(n2trim);
          x1mean[r,k,i,j] = theta1-mean(x1trim);
          x2mean[r,k,i,j] = (theta1+thetadiff)-mean(x2trim);
          
          pe[r,k,i,j] = mean([(std(n1trim)/mean(n1trim)),(std(n2trim)/mean(n2trim))])*
          (1/(std(n1trim+n2trim)/mean(n1trim+n2trim)));
          
          t_ss, relaxtime = timeSS(n1,n2,t_ext);
          
          rt[r,k,i,j] = relaxtime;
          thetadiffarray[r,k,i,j] = thetadiff;
          
          n1_ddm, n2_ddm, x1_ddm, x2_ddm, w1_ddm, w2_ddm, m1_ddm, m2_ddm = 
          KevinEvolveExtinct_ddm(
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
          perror,
          extpop,
          t_ext,
          refuge
          )
          n1trim_ddm = n1_ddm[Int64(floor(tmax*burnin)):tmax];
          n2trim_ddm = n2_ddm[Int64(floor(tmax*burnin)):tmax];
          x1trim_ddm = x1_ddm[Int64(floor(tmax*burnin)):tmax];
          x2trim_ddm = x2_ddm[Int64(floor(tmax*burnin)):tmax];
          
          n1mean_ddm[r,k,i,j] = mean(n1trim_ddm);
          n2mean_ddm[r,k,i,j] = mean(n2trim_ddm);
          x1mean_ddm[r,k,i,j] = theta1-mean(x1trim_ddm);
          x2mean_ddm[r,k,i,j] = (theta1+thetadiff)-mean(x2trim_ddm);
          
          pe_ddm[r,k,i,j] = mean([(std(n1trim)/mean(n1trim)),(std(n2trim)/mean(n2trim))])*
          (1/(std(n1trim+n2trim)/mean(n1trim+n2trim)))
        
          t_ss_ddm, relaxtime_ddm = timeSS(n1_ddm,n2_ddm,t_ext);
          
          rt_ddm[r,k,i,j] = relaxtime_ddm;
          
          #Steady state stray rate?
          m1trim = m1_ddm[Int64(floor(tmax*burnin)):tmax-1];
          m2trim = m2_ddm[Int64(floor(tmax*burnin)):tmax-1];
          
          m1mean[r,k,i,j] = mean(m1trim);
          m2mean[r,k,i,j] = mean(m2trim);
        
        end
      end
  end
end

save(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data3/data_habitathetero_ext.jld"),"n1mean",n1mean,"n2mean",n2mean,"x1mean",x1mean,"x2mean",x2mean,"pe",pe,"rt",rt,"n1mean_ddm",n1mean_ddm,"n2mean_ddm",n2mean_ddm,"x1mean_ddm",x1mean_ddm,"x2mean_ddm",x2mean_ddm,"rt_ddm",rt_ddm,"m1mean",m1mean,"m2mean",m2mean,"pe_ddm",pe_ddm,"thetadiffarray",thetadiffarray);

d = load(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data3/data_habitathetero_ext.jld"));
#This loads the dictionary
n1mean = d["n1mean"];
n2mean = d["n2mean"];
x1mean = d["x1mean"];
x2mean = d["x2mean"];
pe = d["pe"];
rt = d["rt"];
n1mean_ddm = d["n1mean_ddm"];
n2mean_ddm = d["n2mean_ddm"];
x1mean_ddm = d["x1mean_ddm"];
x2mean_ddm = d["x2mean_ddm"];
pe_ddm = d["pe_ddm"];
rt_ddm = d["rt_ddm"];
m1mean = d["m1mean"];
m2mean = d["m2mean"];
thetadiffarray = d["thetadiffarray"];



mediandiff = zeros(Float64,length(hvec),length(thetadiffvec));
for i=1:length(hvec)
  for k=1:length(thetadiffvec)
    absdiff=abs(n1mean[k,i,:] - n2mean[k,i,:])
    mediandiff[i,k] = mean(absdiff[!isnan(absdiff)]);
  end
end

mediandiff_ddm = zeros(Float64,length(hvec),length(thetadiffvec));
for i=1:length(hvec)
  for k=1:length(thetadiffvec)
    absdiff=abs(n1mean_ddm[k,i,:] - n2mean_ddm[k,i,:])
    mediandiff_ddm[i,k] = mean(absdiff[!isnan(absdiff)]);
  end
end

namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/FinalDraft3/fig_thetadiffN.pdf");
R"""
pdf($namespace,height=4,width=5)
library(RColorBrewer)
pal = brewer.pal(9,"Spectral")
plot($(thetadiffvec),$(mediandiff[1,:]),type='l',xlab=expression(paste('Habitat heterogeneity ',Delta,theta)),ylab=expression(paste('Median difference ',Delta,'N')),col=pal[1],ylim=c(0,max($mediandiff_ddm)),lwd=2)
lines($(thetadiffvec),$(mediandiff_ddm[1,:]),lty=2,col=pal[1],lwd=2)
legend(x=0,y=1000,legend=$hvec,col=pal,pch=22,xpd=TRUE,pt.bg=pal,cex=0.6, bty="n",title=expression(paste(h^2)))
"""
for i=2:length(hvec)
  R"""
  lines($(thetadiffvec),$(mediandiff[i,:]),col=pal[$i],lwd=2)
  lines($(thetadiffvec),$(mediandiff_ddm[i,:]),lty=2,col=pal[$i],lwd=2)
  """
end
R"""
dev.off()
"""

window=20;
ma_mvec = movingaverage(mvec,window);
ma_rt3 = movingaverage(mean([rt[31,i,:] for i=1:5]),window);
ma_rt5 = movingaverage(mean([rt[51,i,:] for i=1:5]),window);
ma_rt8 = movingaverage(mean([rt[81,i,:] for i=1:5]),window);



namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/FinalDraft3/fig_relaxtheta2.pdf");
R"""
library(RColorBrewer)
pal = brewer.pal(3,'Set1')
pdf($namespace,height=4,width=10)
par(mfrow=c(1,2),mai = c(0.8, 0.8, 0.4, 0.1))
plot($ma_mvec,$(ma_rt5),type='l',log='y',col=pal[2],ylim=c(20,100),xlab='m',ylab='Recovery time',lwd=2)
lines($ma_mvec,$(ma_rt8),col=pal[3],lwd=2)

legend(x=0.45,y=35,legend=c(5,8),col=pal[2:3],pch=22,xpd=TRUE,pt.bg=pal[2:3],cex=0.8, bty="n",title=expression(paste(Delta,theta)))
plot($ma_mvec,$(ma_rt5),type='l',log='y',col=pal[2],ylim=c(20,100),xlab='m',ylab='Recovery time',lwd=2)
points($(mean([m1mean[51,i,:] for i=1:5])),$(mean([rt_ddm[51,i,:] for i=1:5])),pch=16,cex=.4,col=pal[2])
points($(mean([m2mean[51,i,:] for i=1:5])),$(mean([rt_ddm[51,i,:] for i=1:5])),pch=16,cex=.4,col=pal[2])
for (i in 1:length($(mean([m1mean[51,i,:] for i=1:5])))) {
segments($(mean([m1mean[51,i,:] for i=1:5]))[i],$(mean([rt_ddm[51,i,:] for i=1:5]))[i],$(mean([m2mean[51,i,:] for i=1:5]))[i],$(mean([rt_ddm[51,i,:] for i=1:5]))[i],col=paste(pal[2],75,sep=''))
}
lines($ma_mvec,$(ma_rt8),col=pal[3],lwd=2)
points($(mean([m1mean[81,i,:] for i=1:5])),$(mean([rt_ddm[81,i,:] for i=1:5])),pch=16,cex=.4,col=pal[3])
points($(mean([m2mean[81,i,:] for i=1:5])),$(mean([rt_ddm[81,i,:] for i=1:5])),pch=16,cex=.4,col=pal[3])
for (i in 1:length($(mean([m1mean[51,i,:] for i=1:5])))) {
segments($(mean([m1mean[81,i,:] for i=1:5]))[i],$(mean([rt_ddm[81,i,:] for i=1:5]))[i],$(mean([m2mean[81,i,:] for i=1:5]))[i],$(mean([rt_ddm[81,i,:] for i=1:5]))[i],col=paste(pal[3],75,sep=''))
}
legend(x=0.45,y=35,legend=c(5,8),col=pal[2:3],pch=22,xpd=TRUE,pt.bg=pal[2:3],cex=0.8, bty="n",title=expression(paste(Delta,theta)))
dev.off()
"""





####################
# Looking at bifurcation across a range of asym  values
# asym=0 to 0.05 is up to 2% difference between rmax and beta 

@everywhere using Distributions
@everywhere using RCall
@everywhere using HDF5
@everywhere using JLD
# 
# @everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolve_asym.jl")
# @everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/qualsfunc.jl")
# @everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/bifdet.jl")

@everywhere include("$(homedir())/src/KevinEvolve_asym.jl")
@everywhere include("$(homedir())/src/qualsfunc.jl")
@everywhere include("$(homedir())/src/bifdet.jl")


asymvec = collect(0.1:0.01:0.4);

#Analysis over m
tmax=10000;
mvec = collect(0.0:0.0005:0.5);
# n1ts = zeros(Float64,length(mvec),tmax);
# n2ts = zeros(Float64,length(mvec),tmax);
n1mean=SharedArray(Float64,length(asymvec),length(mvec));
n2mean=SharedArray(Float64,length(asymvec),length(mvec));
n1sd=SharedArray(Float64,length(asymvec),length(mvec));
n2sd=SharedArray(Float64,length(asymvec),length(mvec));
aggmean=SharedArray(Float64,length(asymvec),length(mvec));
aggsd=SharedArray(Float64,length(asymvec),length(mvec));
x1mean=SharedArray(Float64,length(asymvec),length(mvec));
x2mean=SharedArray(Float64,length(asymvec),length(mvec));
pe=SharedArray(Float64,length(asymvec),length(mvec));
# eigs = SharedArray(Array{Complex{Float64}},length(asymvec),length(mvec));
# maxeigs = SharedArray{Float64}(length(asymvec),length(mvec));
# maximeigs = SharedArray{Float64}(length(asymvec),length(mvec));

@time @sync @parallel for a=1:length(asymvec)
    asym = asymvec[a];
    z=2;
    rmax=2.0;
    beta=0.001;
    theta1=5.0;
    thetadiff=2.0;
    tau=1.0;
    h=0.2;
    sigmaE=0;
    sigmaG=1;
    perror=0.00;

    burnin=0.80
    for i=1:length(mvec)
      m=mvec[i];
      
      n1, n2, x1, x2, w1, w2 = 
      KevinEvolve_asym(
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
        perror,
        asym
        );

      n1trim = n1[Int64(floor(tmax*burnin)):tmax];
      n2trim = n2[Int64(floor(tmax*burnin)):tmax];
      x1trim = x1[Int64(floor(tmax*burnin)):tmax];
      x2trim = x2[Int64(floor(tmax*burnin)):tmax]
      
    #   n1ts[i,:] = n1;
    #   n2ts[i,:] = n2;
      
      n1mean[a,i] = mean(n1trim);
      n2mean[a,i] = mean(n2trim);
      n1sd[a,i] = std(n1trim);
      n2sd[a,i] = std(n2trim);
      
      aggmean[a,i] = mean(n1trim+n2trim);
      aggsd[a,i] = std(n1trim+n2trim);

      x1mean[a,i] = theta1-mean(x1trim);
      x2mean[a,i] = (theta1+thetadiff)-mean(x2trim);
      
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
      
      pe[a,i] = mean([(std(n1trim)/mean(n1trim)),(std(n2trim)/mean(n2trim))])*
      (1/(std(n1trim+n2trim)/mean(n1trim+n2trim)))
  end
end

# namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/FinalDraft3/fig_density2.pdf");
namespace = string("$(homedir())/fig_density2.pdf");
R"""
library(RColorBrewer)
pdf($namespace,height=5,width=6)
cols = rev(colorRampPalette(brewer.pal(11, "Spectral"))(length($asymvec)))
plot($mvec,$(n1mean[length(asymvec),:]),pch='.',col=cols[1],xlab="Straying ratio (m)",ylab="Steady state",ylim=c(0,1000))
points($mvec,$(n2mean[length(asymvec),:]),pch='.',col=cols[1])
colseq = $(collect(1:10:length(asymvec)));
legend(x=0.42,y=500,legend=$(asymvec)[colseq],col=cols[colseq],pch=22,xpd=TRUE,pt.bg=cols[colseq],cex=0.8, bty="n",title=expression(paste(Asymmetry)))
"""
for a=length(asymvec)-1:-1:1
    R"""
    points($mvec,$(n1mean[a,:]),pch='.',col=cols[$a])
    points($mvec,$(n2mean[a,:]),pch='.',col=cols[$a])
    """
end
R"dev.off()"


# namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/FinalDraft3/fig_density2.pdf");
namespace = string("$(homedir())/fig_density3.pdf");
R"""
library(RColorBrewer)
pdf($namespace,height=5,width=6)
cols = rev(colorRampPalette(brewer.pal(11, "Spectral"))(length($asymvec)))
plot($mvec,$(n1mean[length(asymvec),:]),pch='.',col=cols[1],xlab="Straying ratio (m)",ylab="Steady state",ylim=c(0,1000))
points($mvec,$(n2mean[length(asymvec),:]),pch='.',col=cols[1])
colseq = $(collect(1:5:length(asymvec)));
legend(x=0.42,y=1000,legend=$(asymvec)[colseq],col=cols[colseq],pch=22,xpd=TRUE,pt.bg=cols[colseq],cex=0.8, bty="n",title=expression(paste(Asymmetry)))
"""
for a=length(asymvec)-1:-1:1
    R"""
    points($mvec,$(n1mean[a,:]),pch='.',col=cols[$a])
    points($mvec,$(n2mean[a,:]),pch='.',col=cols[$a])
    """
end
R"dev.off()"
