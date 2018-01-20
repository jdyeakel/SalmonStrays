using Distributions, RCall


@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolve.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolve_ddm.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinJacobian.jl")

@everywhere include("$(homedir())/src/KevinEvolve.jl")
@everywhere include("$(homedir())/src/KevinEvolve_ddm.jl")
@everywhere include("$(homedir())/src/KevinJacobian.jl")

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

n1ts_ddm = zeros(Float64,length(mvec),tmax);
n2ts_ddm = zeros(Float64,length(mvec),tmax);
n1mean_ddm=zeros(Float64,length(mvec));
n2mean_ddm=zeros(Float64,length(mvec));
n1sd_ddm=zeros(Float64,length(mvec));
n2sd_ddm=zeros(Float64,length(mvec));
aggmean_ddm=zeros(Float64,length(mvec));
aggsd_ddm=zeros(Float64,length(mvec));
x1mean_ddm=zeros(Float64,length(mvec));
x2mean_ddm=zeros(Float64,length(mvec));
m1mean_ddm=zeros(Float64,length(mvec));
m2mean_ddm=zeros(Float64,length(mvec));
pe_ddm=zeros(Float64,length(mvec));

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
C=300;

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

  n1trim = n1[Int64(floor(tmax*burnin)):tmax-1];
  n2trim = n2[Int64(floor(tmax*burnin)):tmax-1];
  x1trim = x1[Int64(floor(tmax*burnin)):tmax-1];
  x2trim = x2[Int64(floor(tmax*burnin)):tmax-1]
  
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
  
  a0 = 1-m;
  
  n1_ddm, n2_ddm, x1_ddm, x2_ddm, w1_ddm, w2_ddm, m1_ddm, m2_ddm = 
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
  
  n1trim_ddm = n1_ddm[Int64(floor(tmax*burnin)):tmax-1];
  n2trim_ddm = n2_ddm[Int64(floor(tmax*burnin)):tmax-1];
  x1trim_ddm = x1_ddm[Int64(floor(tmax*burnin)):tmax-1];
  x2trim_ddm = x2_ddm[Int64(floor(tmax*burnin)):tmax-1]
  
  n1ts_ddm[i,:] = n1_ddm;
  n2ts_ddm[i,:] = n2_ddm;
  
  n1mean_ddm[i] = mean(n1trim_ddm);
  n2mean_ddm[i] = mean(n2trim_ddm);
  n1sd_ddm[i] = std(n1trim_ddm);
  n2sd_ddm[i] = std(n2trim_ddm);
  
  aggmean_ddm[i] = mean(n1trim_ddm+n2trim_ddm);
  aggsd_ddm[i] = std(n1trim_ddm+n2trim_ddm);

  x1mean_ddm[i] = theta1-mean(x1trim_ddm);
  x2mean_ddm[i] = (theta1+thetadiff)-mean(x2trim_ddm);
  
  m1mean_ddm[i] = mean(m1_ddm[Int64(floor(tmax*burnin)):tmax-1]);
  m2mean_ddm[i] = mean(m2_ddm[Int64(floor(tmax*burnin)):tmax-1]);
  
  pe_ddm[i] = mean([(std(n1trim)/mean(n1trim)),(std(n2trim)/mean(n2trim))])*
  (1/(std(n1trim+n2trim)/mean(n1trim+n2trim)))
  
end


R"""
library(RColorBrewer)
pal = brewer.pal(5,'Set1')
plot($mvec,$n1mean,pch=16,col=pal[4])
points($mvec,$n2mean,pch=16,col=pal[4])
points($mvec,$n1mean_ddm,pch=16,col=pal[1])
points($mvec,$n2mean_ddm,pch=16,col=pal[1])
"""



R"""
library(RColorBrewer)
pal = brewer.pal(5,'Set1')
plot($mvec,$n1mean,pch=16,col=pal[4],ylim=c(0,350))
points($mvec,$n2mean,pch=16,col=pal[4])
points($m1mean_ddm,$n1mean_ddm,pch=16,col=pal[1])
points($m2mean_ddm,$n2mean_ddm,pch=16,col=pal[1])
"""


R"""
library(RColorBrewer)
pal = brewer.pal(5,'Set1')
plot($mvec,$n1mean+$n2mean,pch=16,col=pal[4],ylim=c(0,350))
points($(mapslices(mean,[m1mean_ddm m2mean_ddm],2)),$n1mean_ddm+$n2mean_ddm,pch=16,col=pal[1])
"""



R"""
library(RColorBrewer)
pal = brewer.pal(5,'Set1')
plot($mvec,abs($n1mean-$n2mean),pch=16,col=pal[4],ylim=c(0,400))
points($m1mean_ddm,abs($n1mean_ddm-$n2mean_ddm),pch=16,col=pal[1])
points($m2mean_ddm,abs($n1mean_ddm-$n2mean_ddm),pch=16,col=pal[1])
for (i in 1:length($n1mean_ddm)) {
  segments($(m1mean_ddm)[i],abs($n1mean_ddm-$n2mean_ddm)[i],$(m2mean_ddm)[i],abs($n1mean_ddm-$n2mean_ddm)[i],col=pal[1])
}
"""

R"""
L = list()
L[[1]] = abs($n1mean-$n2mean);
L[[2]] = abs($n1mean_ddm-$n2mean_ddm)
boxplot(L)
"""

@everywhere using Distributions, RCall, JLD, HDF5
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolveSS.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolveSS_ddm.jl")

#Analysis over m and C
tmax=10000;
mvec1 = collect(0.0:0.0001:0.5);
mvec = [mvec1 ; reverse(mvec1)];
cexpvec = collect(1:0.01:5);
# n1ts = zeros(Float64,length(mvec),tmax);
# n2ts = zeros(Float64,length(mvec),tmax);
n1mean=SharedArray{Float64}(length(mvec));
n2mean=SharedArray{Float64}(length(mvec));
n1sd=SharedArray{Float64}(length(mvec));
n2sd=SharedArray{Float64}(length(mvec));
aggmean=SharedArray{Float64}(length(mvec));
aggsd=SharedArray{Float64}(length(mvec));
x1mean=SharedArray{Float64}(length(mvec));
x2mean=SharedArray{Float64}(length(mvec));
pe=SharedArray{Float64}(length(mvec));

# n1ts_ddm = zeros(Float64,length(mvec),length(cexpvec),tmax);
# n2ts_ddm = zeros(Float64,length(mvec),length(cexpvec),tmax);
n1mean_ddm=SharedArray{Float64}(length(mvec),length(cexpvec));
n2mean_ddm=SharedArray{Float64}(length(mvec),length(cexpvec));
n1sd_ddm=SharedArray{Float64}(length(mvec),length(cexpvec));
n2sd_ddm=SharedArray{Float64}(length(mvec),length(cexpvec));
aggmean_ddm=SharedArray{Float64}(length(mvec),length(cexpvec));
aggsd_ddm=SharedArray{Float64}(length(mvec),length(cexpvec));
x1mean_ddm=SharedArray{Float64}(length(mvec),length(cexpvec));
x2mean_ddm=SharedArray{Float64}(length(mvec),length(cexpvec));
m1mean_ddm=SharedArray{Float64}(length(mvec),length(cexpvec));
m2mean_ddm=SharedArray{Float64}(length(mvec),length(cexpvec));
pe_ddm=SharedArray{Float64}(length(mvec),length(cexpvec));

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
@sync @parallel for j=1:length(cexpvec)
  
  C = 10^(Float64(cexpvec[j]));
  
  for i=1:length(mvec)
      
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

      n1trim = n1[Int64(floor(tmax*burnin)):tmax-1];
      n2trim = n2[Int64(floor(tmax*burnin)):tmax-1];
      x1trim = x1[Int64(floor(tmax*burnin)):tmax-1];
      x2trim = x2[Int64(floor(tmax*burnin)):tmax-1]
      
      # n1ts[i,:] = n1;
      # n2ts[i,:] = n2;
      
      n1mean[i] = mean(n1trim);
      n2mean[i] = mean(n2trim);
      n1sd[i] = std(n1trim);
      n2sd[i] = std(n2trim);
      
      aggmean[i] = mean(n1trim+n2trim);
      aggsd[i] = std(n1trim+n2trim);

      x1mean[i] = theta1-mean(x1trim);
      x2mean[i] = (theta1+thetadiff)-mean(x2trim);
      

        
        a0 = 1-m;
        
        if i == 1
            n0_ddm = [2,2];
          #   x0 = [theta1 + rand(Normal(0,0.0001)),(theta1 + thetadiff)  + rand(Normal(0,0.0001))];
          x0_ddm = [theta1 ,(theta1 + thetadiff)];
        else
            n0_ddm = [n1mean_ddm[i-1,j],n2mean_ddm[i-1,j]];
            x0_ddm = [x1mean_ddm[i-1,j],x2mean_ddm[i-1,j]];
        end
        
        n1_ddm, n2_ddm, x1_ddm, x2_ddm, w1_ddm, w2_ddm, m1_ddm, m2_ddm = 
        KevinEvolveSS_ddm(
        n0_ddm,
        x0_ddm,
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
        
        n1trim_ddm = n1_ddm[Int64(floor(tmax*burnin)):tmax-1];
        n2trim_ddm = n2_ddm[Int64(floor(tmax*burnin)):tmax-1];
        x1trim_ddm = x1_ddm[Int64(floor(tmax*burnin)):tmax-1];
        x2trim_ddm = x2_ddm[Int64(floor(tmax*burnin)):tmax-1]
        
        # n1ts_ddm[i,j,:] = n1_ddm;
        # n2ts_ddm[i,j,:] = n2_ddm;
        
        n1mean_ddm[i,j] = mean(n1trim_ddm);
        n2mean_ddm[i,j] = mean(n2trim_ddm);
        n1sd_ddm[i,j] = std(n1trim_ddm);
        n2sd_ddm[i,j] = std(n2trim_ddm);
        
        aggmean_ddm[i,j] = mean(n1trim_ddm+n2trim_ddm);
        aggsd_ddm[i,j] = std(n1trim_ddm+n2trim_ddm);

        x1mean_ddm[i,j] = mean(x1trim_ddm);
        x2mean_ddm[i,j] = mean(x2trim_ddm);
        
        m1mean_ddm[i,j] = mean(m1_ddm[Int64(floor(tmax*burnin)):tmax-1]);
        m2mean_ddm[i,j] = mean(m2_ddm[Int64(floor(tmax*burnin)):tmax-1]);
        
        pe_ddm[i,j] = mean([(std(n1trim)/mean(n1trim)),(std(n2trim)/mean(n2trim))])*
        (1/(std(n1trim+n2trim)/mean(n1trim+n2trim)))
    
  end
  
end

save(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data3/data_meandiff.jld"),"n1mean",n1mean,"n2mean",n2mean,"x1mean",x1mean,"x2mean",x2mean,"n1mean_ddm",n1mean_ddm,"n2mean_ddm",n2mean_ddm,"x1mean_ddm",x1mean_ddm,"x2mean_ddm",x2mean_ddm,"m1mean_ddm",m1mean_ddm,"m2mean_ddm",m2mean_ddm,"pe_ddm",pe_ddm);


d = load(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data3/data_meandiff.jld"));
n1mean = d["n1mean"];
n2mean = d["n2mean"];
x1mean = d["x1mean"];
x2mean = d["x2mean"];
n1mean_ddm = d["n1mean_ddm"];
n2mean_ddm = d["n2mean_ddm"];
x1mean_ddm = d["x1mean_ddm"];
x2mean_ddm = d["x2mean_ddm"];
m1mean_ddm = d["m1mean_ddm"];
m2mean_ddm = d["m2mean_ddm"];
pe_ddm = d["pe_ddm"];

midpoint = find(x->x==0.25,mvec)[1];
endpoint = find(x->x==0.5,mvec)[1];
diff_low = abs(n1mean_ddm[1:midpoint,:] .- n2mean_ddm[1:midpoint,:]);
diff_high = abs(n1mean_ddm[midpoint+1:endpoint,:] .- n2mean_ddm[midpoint+1:endpoint,:]);
meandiff_low = mapslices(mean,diff_low,1);
meandiff_high = mapslices(mean,diff_high,1);


less = find(x->x<=0.25,mvec);
more = find(x->x>0.25,mvec);
diff_low = abs(n1mean_ddm[less,:] .- n2mean_ddm[less,:]);
diff_high = abs(n1mean_ddm[more,:] .- n2mean_ddm[more,:]);
meandiff_low = mapslices(mean,diff_low,1);
meandiff_high = mapslices(mean,diff_high,1);



jvec=Array{Int64}(3);
jvec[1]=Int64(find(x->x==1,cexpvec)[1]);
jvec[2]=Int64(find(x->x==2,cexpvec)[1]);
jvec[3]=Int64(find(x->x==3,cexpvec)[1]);


namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/FinalDraft3/fig_meandiff2.pdf");
R"""
library(RColorBrewer)
pal = brewer.pal(5,'Set1')
pdf($namespace,height=6,width=7)
par(fig = c(0,1,0,1),mai = c(0.8, 0.8, 0.1, 0.1))
plot(10^$cexpvec,$meandiff_low,log='x',type='l',lwd=3,col=pal[2],cex=0.8,ylim=c(0,600),xlab='C',ylab='Mean difference in steady state densities')
points(10^$cexpvec,$meandiff_low,col=pal[2],cex=1)
lines(10^$cexpvec,$meandiff_high,lwd=3,col=pal[1],cex=0.8)
lines(10^$cexpvec,rep($(mean(abs(n1mean[less]-n2mean[less]))),length($cexpvec)),lty=2,col=pal[2])
text(10^$(cexpvec[length(cexpvec)])-50000,$(mean(abs(n1mean[1:midpoint]-n2mean[1:midpoint])))+20,expression(paste('Low m, ',m[0],sep='')),col=pal[2])
lines(10^$cexpvec,rep($(mean(abs(n1mean[more]-n2mean[more]))),length($cexpvec)),lty=2,col=pal[1])
text(x=10^$(cexpvec[length(cexpvec)])-50000,y=$(mean(abs(n1mean[midpoint+1:length(mvec)]-n2mean[midpoint+1:length(mvec)])))-35,expression(paste('High m, ',m[0],sep='')),col=pal[1])

cvec = 10^$(cexpvec[jvec])

#Inset figs
par(fig = c(0.04,0.44, 0.7, 0.995), new = T)
plot($mvec,$n1mean,pch='.',cex=1,col=pal[4],xaxt='n',yaxt='n',ann=F,ylim=c(0,400))
points($mvec,$n2mean,pch='.',cex=1,col=pal[4])
points($(mvec[1:endpoint]),$(n1mean_ddm[1:endpoint,jvec[1]]),pch='.',cex=1,col=pal[3])
points($(mvec[1:endpoint]),$(n2mean_ddm[1:endpoint,jvec[1]]),pch='.',cex=1,col=pal[3])
text(0.25,380,paste('C=',cvec[1],sep=''),xpd=T,cex=0.9)
text(-0.06,200,'N*',xpd=T)
text(0.25,-60,expression(paste('m, ',m[0])),xpd=T)
lines(seq(0,$(mvec[midpoint]),length.out=10),rep(0,10),col=pal[2],lwd=3)
lines(seq($(mvec[midpoint]),0.5,length.out=10),rep(0,10),col=pal[1],lwd=3)
segments($(mvec[collect(1:5:length(mvec))]),$(n1mean_ddm[:,jvec[1]][collect(1:5:length(mvec))]),$(mvec[collect(1:5:length(mvec))]),$(n2mean_ddm[:,jvec[1]][collect(1:5:length(mvec))]),col=paste(pal[3],99,sep=''),lwd=0.5)

par(fig = c(0.32,0.72, 0.7, 0.995), new = T)
plot($mvec,$n1mean,pch='.',cex=2,col=pal[4],xaxt='n',yaxt='n',ann=F,ylim=c(0,400))
points($mvec,$n2mean,pch='.',cex=2,col=pal[4])
points($mvec,$(n1mean_ddm[:,jvec[2]]),pch='.',cex=1,col=pal[3])
points($mvec,$(n2mean_ddm[:,jvec[2]]),pch='.',cex=1,col=pal[3])
text(0.25,380,paste('C=',cvec[2],sep=''),xpd=T,cex=0.9)
text(0.25,-60,expression(paste('m, ',m[0])),xpd=T)
lines(seq(0,$(mvec[midpoint]),length.out=10),rep(0,10),col=pal[2],lwd=3)
lines(seq($(mvec[midpoint]),0.5,length.out=10),rep(0,10),col=pal[1],lwd=3)
segments($(mvec[collect(1:50:length(mvec))]),$(n1mean_ddm[:,jvec[2]][collect(1:50:length(mvec))]),$(mvec[collect(1:50:length(mvec))]),$(n2mean_ddm[:,jvec[2]][collect(1:50:length(mvec))]),col=paste(pal[3],99,sep=''),lwd=0.5)

par(fig = c(0.6,0.99, 0.7, 0.995), new = T)
plot($mvec,$n1mean,pch='.',cex=2,col=pal[4],xaxt='n',yaxt='n',ann=F,ylim=c(0,400))
points($mvec,$n2mean,pch='.',cex=2,col=pal[4])
points($mvec,$(n1mean_ddm[:,jvec[3]]),pch='.',cex=1,col=pal[3])
points($mvec,$(n2mean_ddm[:,jvec[3]]),pch='.',cex=1,col=pal[3])
text(0.25,380,paste('C=',cvec[3],sep=''),xpd=T,cex=0.9)
text(0.25,-60,expression(paste('m, ',m[0])),xpd=T)
lines(seq(0,$(mvec[midpoint]),length.out=10),rep(0,10),col=pal[2],lwd=3)
lines(seq($(mvec[midpoint]),0.5,length.out=10),rep(0,10),col=pal[1],lwd=3)
segments($(mvec[collect(1:50:length(mvec))]),$(n1mean_ddm[:,jvec[3]][collect(1:50:length(mvec))]),$(mvec[collect(1:50:length(mvec))]),$(n2mean_ddm[:,jvec[3]][collect(1:50:length(mvec))]),col=paste(pal[3],99,sep=''),lwd=0.5)

#Second row

#Inset figs
par(fig = c(0.04,0.44, 0.5, 0.8), new = T)
plot($mvec,$n1mean,pch='.',cex=2,col=pal[4],xaxt='n',yaxt='n',ann=F,ylim=c(0,400))
points($mvec,$n2mean,pch='.',cex=2,col=pal[4])
points($(m1mean_ddm[:,jvec[1]]),$(n1mean_ddm[:,jvec[1]]),pch='.',cex=1,col=pal[3])
points($(m2mean_ddm[:,jvec[1]]),$(n2mean_ddm[:,jvec[1]]),pch='.',cex=1,col=pal[3])
text(-0.06,200,'N*',xpd=T)
text(0.25,-60,'m, m*',xpd=T)
segments($(m1mean_ddm[:,jvec[1]][collect(1:50:length(mvec))]),$(n1mean_ddm[:,jvec[1]][collect(1:50:length(mvec))]),$(m2mean_ddm[:,jvec[1]][collect(1:50:length(mvec))]),$(n2mean_ddm[:,jvec[1]][collect(1:50:length(mvec))]),col=paste(pal[3],99,sep=''),lwd=0.5)

par(fig = c(0.32,0.72, 0.5, 0.8), new = T)
plot($mvec,$n1mean,pch='.',cex=2,col=pal[4],xaxt='n',yaxt='n',ann=F,ylim=c(0,400))
points($mvec,$n2mean,pch='.',cex=2,col=pal[4])
points($(m1mean_ddm[:,jvec[2]]),$(n1mean_ddm[:,jvec[2]]),pch='.',cex=1,col=pal[3])
points($(m2mean_ddm[:,jvec[2]]),$(n2mean_ddm[:,jvec[2]]),pch='.',cex=1,col=pal[3])
text(0.25,-60,'m, m*',xpd=T)
segments($(m1mean_ddm[:,jvec[2]][collect(1:50:length(mvec))]),$(n1mean_ddm[:,jvec[2]][collect(1:50:length(mvec))]),$(m2mean_ddm[:,jvec[2]][collect(1:50:length(mvec))]),$(n2mean_ddm[:,jvec[2]][collect(1:50:length(mvec))]),col=paste(pal[3],99,sep=''),lwd=0.5)

par(fig = c(0.6,0.99, 0.5, 0.8), new = T)
plot($mvec,$n1mean,pch='.',cex=2,col=pal[4],xaxt='n',yaxt='n',ann=F,ylim=c(0,400))
points($mvec,$n2mean,pch='.',cex=2,col=pal[4])
points($(m1mean_ddm[:,jvec[3]]),$(n1mean_ddm[:,jvec[3]]),pch='.',cex=1,col=pal[3])
points($(m2mean_ddm[:,jvec[3]]),$(n2mean_ddm[:,jvec[3]]),pch='.',cex=1,col=pal[3])
text(0.25,-60,'m, m*',xpd=T)
segments($(m1mean_ddm[:,jvec[3]][collect(1:50:length(mvec))]),$(n1mean_ddm[:,jvec[3]][collect(1:50:length(mvec))]),$(m2mean_ddm[:,jvec[3]][collect(1:50:length(mvec))]),$(n2mean_ddm[:,jvec[3]][collect(1:50:length(mvec))]),col=paste(pal[3],99,sep=''),lwd=0.5)
dev.off()
"""


lines(10^$cexpvec,rep($(mean(abs(n1mean[1:midpoint]-n2mean[1:midpoint]))),length($cexpvec)),lty=2,col=pal[2])
text(10^$(cexpvec[length(cexpvec)])-50000,$(mean(abs(n1mean[1:midpoint]-n2mean[1:midpoint])))+20,expression(paste('Low m, ',m[0],sep='')),col=pal[2])
lines(10^$cexpvec,rep($(mean(abs(n1mean[midpoint+1:endpoint]-n2mean[midpoint+1:endpoint]))),length($cexpvec)),lty=2,col=pal[1])
text(x=10^$(cexpvec[length(cexpvec)])-50000,y=$(mean(abs(n1mean[midpoint+1:length(mvec)]-n2mean[midpoint+1:length(mvec)])))+20,expression(paste('High m, ',m[0],sep='')),col=pal[1])


R"""
library(RColorBrewer)
pal = brewer.pal(5,'Set1')
par(fig = c(0,1,0,1))
plot(10^$cexpvec,$meandiff_low,log='x',pch=16,col=pal[2],ylim=c(0,300),xlab='C',ylab='Mean difference in steady state densities')
points(10^$cexpvec,$meandiff_high,pch=16,col=pal[1])
lines(10^$cexpvec,rep($(mean(abs(n1mean[1:midpoint]-n2mean[1:midpoint]))),length($cexpvec)),lty=2,col=pal[2])
text(10^$(cexpvec[length(cexpvec)])-20000,$(mean(abs(n1mean[1:midpoint]-n2mean[1:midpoint])))+10,'Low m',col=pal[2])
lines(10^$cexpvec,rep($(mean(abs(n1mean[midpoint+1:length(mvec)]-n2mean[midpoint+1:length(mvec)]))),length($cexpvec)),lty=2,col=pal[1])
text(x=10^$(cexpvec[length(cexpvec)])-20000,y=$(mean(abs(n1mean[midpoint+1:length(mvec)]-n2mean[midpoint+1:length(mvec)])))+10,'High m',col=pal[1])

cvec = 10^$(cexpvec[jvec])

#Inset figs
par(fig = c(0.01,0.4, 0.6, 0.99), new = T)
plot($mvec,$n1mean,pch='.',cex=2,col=pal[4],xaxt='n',yaxt='n',ann=F,ylim=c(0,400))
points($mvec,$n2mean,pch='.',cex=2,col=pal[4])
points($(m1mean_ddm[:,jvec[1]]),$(n1mean_ddm[:,jvec[1]]),pch='.',cex=2,col=pal[3])
points($(m2mean_ddm[:,jvec[1]]),$(n2mean_ddm[:,jvec[1]]),pch='.',cex=2,col=pal[3])
text(0.25,370,paste('C=',cvec[1],sep=''),xpd=T)
text(0.25,-60,'m, m*',xpd=T)
lines(seq(0,0.25,length.out=10),rep(0,10),col=pal[2],lwd=3)
lines(seq(0.25,0.5,length.out=10),rep(0,10),col=pal[1],lwd=3)

par(fig = c(0.3,0.7, 0.6, 0.99), new = T)
plot($mvec,$n1mean,pch='.',cex=2,col=pal[4],xaxt='n',yaxt='n',ann=F,ylim=c(0,400))
points($mvec,$n2mean,pch='.',cex=2,col=pal[4])
points($(m1mean_ddm[:,jvec[2]]),$(n1mean_ddm[:,jvec[2]]),pch='.',cex=2,col=pal[3])
points($(m2mean_ddm[:,jvec[2]]),$(n2mean_ddm[:,jvec[2]]),pch='.',cex=2,col=pal[3])
text(0.25,370,paste('C=',cvec[2],sep=''),xpd=T)
text(0.25,-60,'m, m*',xpd=T)
lines(seq(0,0.25,length.out=10),rep(0,10),col=pal[2],lwd=3)
lines(seq(0.25,0.5,length.out=10),rep(0,10),col=pal[1],lwd=3)


par(fig = c(0.6,0.99, 0.6, 0.99), new = T)
plot($mvec,$n1mean,pch='.',cex=2,col=pal[4],xaxt='n',yaxt='n',ann=F,ylim=c(0,400))
points($mvec,$n2mean,pch='.',cex=2,col=pal[4])
points($(m1mean_ddm[:,jvec[3]]),$(n1mean_ddm[:,jvec[3]]),pch='.',cex=2,col=pal[3])
points($(m2mean_ddm[:,jvec[3]]),$(n2mean_ddm[:,jvec[3]]),pch='.',cex=2,col=pal[3])
text(0.25,370,paste('C=',cvec[3],sep=''),xpd=T)
text(0.25,-60,'m, m*',xpd=T)
lines(seq(0,0.25,length.out=10),rep(0,10),col=pal[2],lwd=3)
lines(seq(0.25,0.5,length.out=10),rep(0,10),col=pal[1],lwd=3)

"""



########################
# PE and Recovery as a function of m0 and C
########################

@everywhere using Distributions
@everywhere using RCall
@everywhere using HDF5
@everywhere using JLD

@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolveExtinct_ddm.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/timeSS.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/qualsfunc.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/bifdet.jl")


#Analysis over m and C
tmax=10000;
reps = 20;
mvec = collect(0.0:0.005:0.5);
cexpvec = collect(1:0.01:5);
pvec = ["small","large","both"];

# n1ts_ddm = SharedArray{Float64}(reps,length(mvec),length(cexpvec),tmax);
# n2ts_ddm = SharedArray{Float64}(reps,length(mvec),length(cexpvec),tmax);
n1mean_ddm=SharedArray{Float64}(reps,length(mvec),length(cexpvec),length(pvec));
n2mean_ddm=SharedArray{Float64}(reps,length(mvec),length(cexpvec),length(pvec));
n1sd_ddm=SharedArray{Float64}(reps,length(mvec),length(cexpvec),length(pvec));
n2sd_ddm=SharedArray{Float64}(reps,length(mvec),length(cexpvec),length(pvec));
aggmean_ddm=SharedArray{Float64}(reps,length(mvec),length(cexpvec),length(pvec));
aggsd_ddm=SharedArray{Float64}(reps,length(mvec),length(cexpvec),length(pvec));
x1mean_ddm=SharedArray{Float64}(reps,length(mvec),length(cexpvec),length(pvec));
x2mean_ddm=SharedArray{Float64}(reps,length(mvec),length(cexpvec),length(pvec));
m1mean_ddm=SharedArray{Float64}(reps,length(mvec),length(cexpvec),length(pvec));
m2mean_ddm=SharedArray{Float64}(reps,length(mvec),length(cexpvec),length(pvec));
pe_ddm=SharedArray{Float64}(reps,length(mvec),length(cexpvec),length(pvec));
rt_ddm=SharedArray{Float64}(reps,length(mvec),length(cexpvec),length(pvec));

z=2;
rmax=2.0;
beta=0.001;
theta1=5.0;
thetadiff=2.0;
tau=1.0;
h=0.2;
sigmaE=0;
sigmaG=1;
perror=0.01;
refuge=0.01;
t_ext = Int64(round(tmax/2));

burnin=0.80
@sync @parallel for r=1:reps

    for i=1:length(mvec)
      
      m=mvec[i];
      a0 = 1-m;
      
      for j=1:length(cexpvec)
          
          C = 10^(Float64(cexpvec[j]));
          
          for k=1:length(pvec)
              
              extpop=pvec[k];
              
              n1_ddm_pre, n2_ddm_pre, x1_ddm_pre, x2_ddm_pre, w1_ddm_pre, w2_ddm_pre, m1_ddm_pre, m2_ddm_pre = 
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
              );
              
              #Set n1 to be high, n2 to be low
              #This is to align dominant and subordinate populations when we take the means across reps. This doesn't matter when there are no reps.
              if n1_ddm_pre[tmax-1] < n2_ddm_pre[tmax-1]
                  n1_ddm = copy(n2_ddm_pre);
                  n2_ddm = copy(n1_ddm_pre);
                  x1_ddm = copy(x2_ddm_pre);
                  x2_ddm = copy(x1_ddm_pre);
                  w1_ddm = copy(w2_ddm_pre);
                  w2_ddm = copy(w1_ddm_pre);
                  m1_ddm = copy(m2_ddm_pre);
                  m2_ddm = copy(m1_ddm_pre);
              else
                  n1_ddm = copy(n1_ddm_pre);
                  n2_ddm = copy(n2_ddm_pre);
                  x1_ddm = copy(x1_ddm_pre);
                  x2_ddm = copy(x2_ddm_pre);
                  w1_ddm = copy(w1_ddm_pre);
                  w2_ddm = copy(w2_ddm_pre);
                  m1_ddm = copy(m1_ddm_pre);
                  m2_ddm = copy(m2_ddm_pre);
              end
              
              
              t_ss, relaxtime = timeSS(n1_ddm,n2_ddm,t_ext);
              
              rt_ddm[r,i,j,k] = relaxtime;
              
              n1trim_ddm = n1_ddm[Int64(floor(tmax*burnin)):tmax-1];
              n2trim_ddm = n2_ddm[Int64(floor(tmax*burnin)):tmax-1];
              x1trim_ddm = x1_ddm[Int64(floor(tmax*burnin)):tmax-1];
              x2trim_ddm = x2_ddm[Int64(floor(tmax*burnin)):tmax-1]
              
              # n1ts_ddm[r,i,j,:] = n1_ddm;
              # n2ts_ddm[r,i,j,:] = n2_ddm;
              
              n1mean_ddm[r,i,j,k] = mean(n1trim_ddm);
              n2mean_ddm[r,i,j,k] = mean(n2trim_ddm);
              n1sd_ddm[r,i,j,k] = std(n1trim_ddm);
              n2sd_ddm[r,i,j] = std(n2trim_ddm);
              
              aggmean_ddm[r,i,j,k] = mean(n1trim_ddm+n2trim_ddm);
              aggsd_ddm[r,i,j,k] = std(n1trim_ddm+n2trim_ddm);
              
              x1mean_ddm[r,i,j,k] = theta1-mean(x1trim_ddm);
              x2mean_ddm[r,i,j,k] = (theta1+thetadiff)-mean(x2trim_ddm);
              
              m1mean_ddm[r,i,j,k] = mean(m1_ddm[Int64(floor(tmax*burnin)):tmax-1]);
              m2mean_ddm[r,i,j,k] = mean(m2_ddm[Int64(floor(tmax*burnin)):tmax-1]);
              
              pe_ddm[r,i,j,k] = mean([(std(n1trim_ddm)/mean(n1trim_ddm)),(std(n2trim_ddm)/mean(n2trim_ddm))])*
              (1/(std(n1trim_ddm+n2trim_ddm)/mean(n1trim_ddm+n2trim_ddm)))
              
          end
          
      end
      
    end

end

save(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data3/data_4regimes_rt_pe.jld"),"n1mean_ddm",n1mean_ddm,"n2mean_ddm",n2mean_ddm,"x1mean_ddm",x1mean_ddm,"x2mean_ddm",x2mean_ddm,"m1mean_ddm",m1mean_ddm,"m2mean_ddm",m2mean_ddm,"pe_ddm",pe_ddm,"rt_ddm",rt_ddm);


d = load(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data3/data_4regimes_rt_pe.jld"));
n1mean_ddm = d["n1mean_ddm"];
n2mean_ddm = d["n2mean_ddm"];
x1mean_ddm = d["x1mean_ddm"];
x2mean_ddm = d["x2mean_ddm"];
m1mean_ddm = d["m1mean_ddm"];
m2mean_ddm = d["m2mean_ddm"];
pe_ddm = d["pe_ddm"];
rt_ddm = d["rt_ddm"];


#Take means over reps
extpop = "both"; #small large both
indext = find(x->x==extpop,pvec)[1];

mpe_ddm = mapslices(mean,pe_ddm,1)[1,:,:,indext];
mrt_ddm = mapslices(mean,rt_ddm,1)[1,:,:,indext];

namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/FinalDraft3/fig_rtpe_ddm.pdf");
R"""
library(RColorBrewer)
library(fields)
pal = rev(brewer.pal(9,"Blues"))
pal2 = rev(brewer.pal(11,"Spectral"))
pdf($namespace,height=4,width=10)
par(mfrow=c(1,2),mai = c(0.8, 0.8, 0.3, 0.9))
image.plot(x=$mvec,y=$cexpvec,z=$mpe_ddm,zlim=c(min($(mpe_ddm[!isnan(mpe_ddm)])),2),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
text(x=0.25,y=5.2,'PE',xpd=T)
image.plot(x=$mvec,y=$cexpvec,z=log($mrt_ddm,10),zlim=c(min(log($mrt_ddm,10)),2.2),col=pal2,xlab=expression(paste(m[0])),ylab='')
text(x=0.25,y=5.2,expression(paste(log[10],' Recovery time')),xpd=T)
dev.off()
"""

R"""
par(mfrow=c(2,2))
plot($mvec,$(mrt_ddm[:,1]),type='l',log='y')
plot($mvec,$(mrt_ddm[:,101]),type='l',log='y')
plot($mvec,$(mrt_ddm[:,201]),type='l',log='y')
plot($mvec,$(mrt_ddm[:,301]),type='l',log='y')
"""

R"""
plot($(mrt_ddm[:,401]),$(mpe_ddm[:,401]),log='xy')
"""


#Small/large collapse

#Take means over reps
extpop = "small"; #small large both
indext = find(x->x==extpop,pvec)[1];

mpe_ddm_s = mapslices(mean,pe_ddm,1)[1,:,:,indext];
mrt_ddm_s = mapslices(mean,rt_ddm,1)[1,:,:,indext];


#Take means over reps
extpop = "large"; #small large both
indext_l = find(x->x==extpop,pvec)[1];

mpe_ddm_l = mapslices(mean,pe_ddm,1)[1,:,:,indext_l];
mrt_ddm_l = mapslices(mean,rt_ddm,1)[1,:,:,indext_l];


namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/FinalDraft3/fig_rtpe_ddmsl.pdf");
R"""
library(RColorBrewer)
library(fields)
pal = rev(brewer.pal(9,"Blues"))
pal2 = rev(brewer.pal(11,"Spectral"))
pdf($namespace,height=10,width=10)
par(mfrow=c(2,2),mai = c(0.8, 0.8, 0.3, 0.9))
image.plot(x=$mvec,y=$cexpvec,z=$mpe_ddm_s,zlim=c(min($(mpe_ddm_s[!isnan(mpe_ddm_s)])),2),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
text(x=0.25,y=5.2,'Portfolio effect',xpd=T)
text(x=0.4,y=4.9,'Subordinate extinct',xpd=T,col='white')
image.plot(x=$mvec,y=$cexpvec,z=log($mrt_ddm_s,10),zlim=c(min(log($mrt_ddm_s,10)),2.2),col=pal2,xlab=expression(paste(m[0])),ylab='')
text(x=0.25,y=5.2,expression(paste(log[10],' Recovery time')),xpd=T)
text(x=0.4,y=4.9,'Subordinate extinct',xpd=T,col='black')
image.plot(x=$mvec,y=$cexpvec,z=$mpe_ddm_l,zlim=c(min($(mpe_ddm_l[!isnan(mpe_ddm_l)])),2),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
text(x=0.25,y=5.2,'Portfolio effect',xpd=T)
text(x=0.4,y=4.9,'Dominant extinct',xpd=T,col='white')
image.plot(x=$mvec,y=$cexpvec,z=log($mrt_ddm_l,10),zlim=c(min(log($mrt_ddm_l,10)),2.2),col=pal2,xlab=expression(paste(m[0])),ylab='')
text(x=0.25,y=5.2,expression(paste(log[10],' Recovery time')),xpd=T)
text(x=0.4,y=4.9,'Dominant extinct',xpd=T,col='black')
dev.off()
"""

namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/FinalDraft3/fig_rtpe_ddmsl.pdf");
R"""
library(RColorBrewer)
library(fields)
pal = rev(brewer.pal(9,"Blues"))
pal2 = brewer.pal(11,"Spectral")
pdf($namespace,height=5,width=10)
par(mfrow=c(1,2),mai = c(0.8, 0.8, 0.3, 0.9))
image.plot(x=$mvec,y=$cexpvec,z=log($mrt_ddm_s,10),zlim=c(min(log($mrt_ddm_s,10)),2.2),col=pal2,xlab=expression(paste(m[0])),ylab='')
text(x=0.25,y=5.2,expression(paste(log[10],' Recovery time')),xpd=T)
image.plot(x=$mvec,y=$cexpvec,z=log($mrt_ddm_l,10),zlim=c(min(log($mrt_ddm_l,10)),2.2),col=pal2,xlab=expression(paste(m[0])),ylab='')
text(x=0.25,y=5.2,expression(paste(log[10],' Recovery time')),xpd=T)
dev.off()
"""


##############
#Figure 4
##############

namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/FinalDraft3/fig_rtpe_ddm.pdf");
R"""
library(RColorBrewer)
library(fields)
pal = rev(brewer.pal(9,"Blues"))
pal2 = rev(brewer.pal(11,"Spectral"))
pdf($namespace,height=8,width=10)
par(mfrow=c(2,2),mai = c(0.8, 0.8, 0.3, 0.9))
image.plot(x=$mvec,y=$cexpvec,z=$mpe_ddm,zlim=c(min($(mpe_ddm[!isnan(mpe_ddm)])),2),col=pal,xlab=expression(paste('Individual straying ',m[0])),ylab=expression(paste(log[10],' C')))
text(x=0.25,y=5.2,'Portfolio effect',xpd=T)
image.plot(x=$mvec,y=$cexpvec,z=log($mrt_ddm,10),zlim=c(min(log($mrt_ddm,10)),2.2),col=pal2,xlab=expression(paste('Individual straying ',m[0])),ylab=expression(paste(log[10],' C')))
text(x=0.43,y=4.8,'Near-collapse',xpd=T,col='black')
text(x=0.25,y=5.2,expression(paste(log[10],'Recovery time')),xpd=T)
image.plot(x=$mvec,y=$cexpvec,z=log($mrt_ddm_s,10),zlim=c(min(log($mrt_ddm_s,10)),2.2),col=pal2,xlab=expression(paste('Individual straying ',m[0])),ylab=expression(paste(log[10],' C')))
text(x=0.25,y=5.2,expression(paste(log[10],'Recovery time')),xpd=T)
text(x=0.4,y=4.8,'Subordinate extinct',xpd=T,col='black')
image.plot(x=$mvec,y=$cexpvec,z=log($mrt_ddm_l,10),zlim=c(min(log($mrt_ddm_l,10)),2.2),col=pal2,xlab=expression(paste('Individual straying ',m[0])),ylab=expression(paste(log[10],' C')))
text(x=0.25,y=5.2,expression(paste(log[10],'Recovery time')),xpd=T)
text(x=0.42,y=4.8,'Dominant extinct',xpd=T,col='black')
dev.off()
"""


#Take means over reps
extpop = "both"; #small large both
indext = find(x->x==extpop,pvec)[1];
n1mean_b = mapslices(mean,n1mean_ddm,1)[1,:,:,indext];
n2mean_b = mapslices(mean,n2mean_ddm,1)[1,:,:,indext];
#Take means over reps
extpop = "small"; #small large both
indext = find(x->x==extpop,pvec)[1];
n1mean_s = mapslices(mean,n1mean_ddm,1)[1,:,:,indext];
n2mean_s = mapslices(mean,n2mean_ddm,1)[1,:,:,indext];
#Take means over reps
extpop = "large"; #small large both
indext = find(x->x==extpop,pvec)[1];
n1mean_l = mapslices(mean,n1mean_ddm,1)[1,:,:,indext];
n2mean_l = mapslices(mean,n2mean_ddm,1)[1,:,:,indext];

#Compare final steady state between 3 disturbance regimes
R"""
library(RColorBrewer)
library(fields)
pal = rev(brewer.pal(9,"Blues"))
par(mfrow=c(3,2),mai = c(0.8, 0.8, 0.3, 0.9))
image.plot(x=$mvec,y=$cexpvec,z=log($n1mean_b,10),zlim=c(1,2.5),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
image.plot(x=$mvec,y=$cexpvec,z=log($n2mean_b,10),zlim=c(1,2.5),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
image.plot(x=$mvec,y=$cexpvec,z=log($n1mean_s,10),zlim=c(1,2.5),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
image.plot(x=$mvec,y=$cexpvec,z=log($n2mean_s,10),zlim=c(1,2.5),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
image.plot(x=$mvec,y=$cexpvec,z=log($n1mean_l,10),zlim=c(1,2.5),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
image.plot(x=$mvec,y=$cexpvec,z=log($n2mean_l,10),zlim=c(1,2.5),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
"""


R"""
par(mfrow = c(1,2))
image.plot(x=$mvec,y=$cexpvec,z=abs($n1mean_b-$n2mean_b),zlim=c(0,500),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
image.plot(x=$mvec,y=$cexpvec,z=abs($n1mean_l-$n2mean_l),zlim=c(0,500),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
"""


R"""
library(fields)
library(RColorBrewer)
pal = rev(brewer.pal(9,"Blues"))
par(mfrow = c(3,2))
image.plot(x=$mvec,y=$cexpvec,z=abs($n1mean_b),zlim=c(0,400),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
image.plot(x=$mvec,y=$cexpvec,z=abs($n2mean_b),zlim=c(0,400),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
image.plot(x=$mvec,y=$cexpvec,z=abs($n1mean_l),zlim=c(0,400),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
image.plot(x=$mvec,y=$cexpvec,z=abs($n2mean_l),zlim=c(0,400),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
image.plot(x=$mvec,y=$cexpvec,z=abs($n1mean_s),zlim=c(0,400),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
image.plot(x=$mvec,y=$cexpvec,z=abs($n2mean_s),zlim=c(0,400),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
"""


R"""
par(mfrow = c(3,1))
image.plot(x=$mvec,y=$cexpvec,z=abs($n1mean_b+$n2mean_b),zlim=c(0,500),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
image.plot(x=$mvec,y=$cexpvec,z=abs($n1mean_s+$n2mean_s),zlim=c(0,500),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
image.plot(x=$mvec,y=$cexpvec,z=abs($n1mean_l+$n2mean_l),zlim=c(0,500),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
"""




#Bifurcation Plot as a function of m and C


# 2-dimensional search

@everywhere using RCall, Distributions, HDF5, JLD
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolveSS_ddm.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinJacobian.jl")



#Analysis over m
tmax=10000;
mvec1 = collect(0.0:0.001:0.5);
mvec = [mvec1 ; reverse(mvec1)];

cexpvec = collect(1:0.01:5);

# n1ts = Array{Float64}(length(mvec),length(cexpvec),tmax);
# n2ts = Array{Float64}(length(mvec),length(cexpvec),tmax);
n1mean=Array{Float64}(length(mvec),length(cexpvec));
n2mean=Array{Float64}(length(mvec),length(cexpvec));
n1sd=Array{Float64}(length(mvec),length(cexpvec));
n2sd=Array{Float64}(length(mvec),length(cexpvec));
aggmean=Array{Float64}(length(mvec),length(cexpvec));
aggsd=Array{Float64}(length(mvec),length(cexpvec));
x1mean=Array{Float64}(length(mvec),length(cexpvec));
x2mean=Array{Float64}(length(mvec),length(cexpvec));
pe=Array{Float64}(length(mvec),length(cexpvec));
# eigs = Array(Array{Complex{Float64}},length(mvec));
maxeigs = Array{Float64}(length(mvec),length(cexpvec));
maximeigs = Array{Float64}(length(mvec),length(cexpvec));
mineigs = Array{Float64}(length(mvec),length(cexpvec));
minimeigs = Array{Float64}(length(mvec),length(cexpvec));

z=2.0;
rmax=2.0;
beta=0.001;
theta1=5.0;
thetadiff=2.0;

tau=1;
h=0.2;
sigmaE=0;
sigmaG=1;
perror=0.00;

burnin=0.80
for r=1:length(cexpvec)

    C = 10^(Float64(cexpvec[r]));

    for i=1:length(mvec)
      
      m=mvec[i];
      a0 = 1-m;
      
      if i == 1
          n0 = [2,2];
        #   x0 = [theta1 + rand(Normal(0,0.0001)),(theta1 + thetadiff)  + rand(Normal(0,0.0001))];
        x0 = [theta1 ,(theta1 + thetadiff)];
      else
          n0 = [n1mean[i-1,r],n2mean[i-1,r]];
          x0 = [x1mean[i-1,r],x2mean[i-1,r]];
      end
      
      n1, n2, x1, x2, w1, w2, m1, m2 = 
      KevinEvolveSS_ddm(
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
      
      # n1ts[i,r,:] = n1;
      # n2ts[i,r,:] = n2;
      
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

save(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data3/data_hysteresis_ddm.jld"),"n1mean",n1mean,"n2mean",n2mean,"x1mean",x1mean,"x2mean",x2mean,"pe",pe);



#Plot the difference in n1mean, n2mean
midpoint = Int64(floor(length(mvec)/2));

# fmaxeigs = maxeigs[1:midpoint,:];
# bmaxeigs = flipdim(maxeigs[midpoint+1:length(mvec),:],1);
# fold = 0.98 .< fmaxeigs .< 1.00;
# bold = 0.98 .< bmaxeigs .< 1.00;

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

namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/FinalDraft3/fig_hysteresis_ddm.pdf");
R"""
pdf($namespace,height=5,width=6)
par(mfrow=c(1,1))
image(x=$fmvec,y=$cexpvec,z=$fdiffmeanbin,xlab='Straying ratio m',ylab=expression(paste(log[10],' C')),col=c('white','black'))
image(x=$fmvec,y=$cexpvec,z=$hystdiff,xlab='Straying ratio m',ylab=expression(paste(log[10],' C')),col=c('#ffffff00','gray'),add=T)
dev.off()
"""






#####################
# SECTION 3 - Habitat heterogeneity
#####################

# bifurcations


@everywhere using RCall, Distributions, HDF5, JLD
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolveSS_ddm.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinJacobian.jl")



#Analysis over m
tmax=10000;
mvec1 = collect(0.0:0.01:1.0);
mvec = [mvec1 ; reverse(mvec1)];
cexpvec = collect(1:0.01:5);

thetadiffvec = [1.8 2.0 2.2];

# n1ts = Array{Float64}(length(mvec),length(cexpvec),tmax);
# n2ts = Array{Float64}(length(mvec),length(cexpvec),tmax);
n1mean=SharedArray{Float64}(length(mvec),length(cexpvec),length(thetadiffvec));
n2mean=SharedArray{Float64}(length(mvec),length(cexpvec),length(thetadiffvec));
n1sd=SharedArray{Float64}(length(mvec),length(cexpvec),length(thetadiffvec));
n2sd=SharedArray{Float64}(length(mvec),length(cexpvec),length(thetadiffvec));
aggmean=SharedArray{Float64}(length(mvec),length(cexpvec),length(thetadiffvec));
aggsd=SharedArray{Float64}(length(mvec),length(cexpvec),length(thetadiffvec));
x1mean=SharedArray{Float64}(length(mvec),length(cexpvec),length(thetadiffvec));
x2mean=SharedArray{Float64}(length(mvec),length(cexpvec),length(thetadiffvec));
pe=SharedArray{Float64}(length(mvec),length(cexpvec),length(thetadiffvec));
# eigs = Array(Array{Complex{Float64}},length(mvec));

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
@sync @parallel for r=1:length(cexpvec)

    C = 10^(Float64(cexpvec[r]));
        
    for th=1:length(thetadiffvec)
                  
        thetadiff=thetadiffvec[th];

        for i=1:length(mvec)
          
          m=mvec[i];
          a0 = 1-m;
          
              if i == 1
                  n0 = [2,2];
                #   x0 = [theta1 + rand(Normal(0,0.0001)),(theta1 + thetadiff)  + rand(Normal(0,0.0001))];
                x0 = [theta1 ,(theta1 + thetadiff)];
              else
                  n0 = [n1mean[i-1,r,th],n2mean[i-1,r,th]];
                  x0 = [x1mean[i-1,r,th],x2mean[i-1,r,th]];
              end
              
              n1, n2, x1, x2, w1, w2, m1, m2 = 
              KevinEvolveSS_ddm(
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
              
              # n1ts[i,r,:] = n1;
              # n2ts[i,r,:] = n2;
              
              n1mean[i,r,th] = mean(n1trim);
              n2mean[i,r,th] = mean(n2trim);
              n1sd[i,r,th] = std(n1trim);
              n2sd[i,r,th] = std(n2trim);
              
              aggmean[i,r,th] = mean(n1trim+n2trim);
              aggsd[i,r,th] = std(n1trim+n2trim);

              x1mean[i,r,th] = mean(x1trim);
              x2mean[i,r,th] = mean(x2trim);
              
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
              
              pe[i,r,th] = mean([(std(n1trim)/mean(n1trim)),(std(n2trim)/mean(n2trim))])*
              (1/(std(n1trim+n2trim)/mean(n1trim+n2trim)))
          
         end
          
    end
end

save(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data3/data_hysteresis_ddm_theta.jld"),"n1mean",n1mean,"n2mean",n2mean,"x1mean",x1mean,"x2mean",x2mean,"pe",pe);

d=load(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data3/data_hysteresis_ddm_theta.jld"));
n1mean = d["n1mean"];
n2mean = d["n2mean"];
x1mean = d["x1mean"];
x2mean = d["x2mean"];
pe = d["pe"];

#Plot the difference in n1mean, n2mean
midpoint = Int64(floor(length(mvec)/2));

#Low heterogeneity
fn1mean1 = copy(n1mean[1:midpoint,:,1]);
bn1mean1 = copy(flipdim(n1mean[midpoint+1:length(mvec),:,1],1));
fn2mean1 = copy(n2mean[1:midpoint,:,1]);
bn2mean1 = copy(flipdim(n2mean[midpoint+1:length(mvec),:,1],1));
fmvec = mvec[1:midpoint];
bmvec = fmvec;
fdiffmean1 = fn1mean1 .- fn2mean1;
bdiffmean1 = bn1mean1 .- bn2mean1;
fdiffmeanbin1 = fdiffmean1 .> 0.000001;
bdiffmeanbin1 = bdiffmean1 .> 0.000001;
hystdiff1 = fdiffmeanbin1 .!= bdiffmeanbin1;

#Intermediate heterogeneity
fn1mean2 = copy(n1mean[1:midpoint,:,2]);
bn1mean2 = copy(flipdim(n1mean[midpoint+1:length(mvec),:,2],1));
fn2mean2 = copy(n2mean[1:midpoint,:,2]);
bn2mean2 = copy(flipdim(n2mean[midpoint+1:length(mvec),:,2],1));
fdiffmean2 = fn1mean2 .- fn2mean2;
bdiffmean2 = bn1mean2 .- bn2mean2;
fdiffmeanbin2 = fdiffmean2 .> 0.000001;
bdiffmeanbin2 = bdiffmean2 .> 0.000001;
hystdiff2 = fdiffmeanbin2 .!= bdiffmeanbin2;

#High heterogeneity
fn1mean3 = copy(n1mean[1:midpoint,:,3]);
bn1mean3 = copy(flipdim(n1mean[midpoint+1:length(mvec),:,3],1));
fn2mean3 = copy(n2mean[1:midpoint,:,3]);
bn2mean3 = copy(flipdim(n2mean[midpoint+1:length(mvec),:,3],1));
fdiffmean3 = fn1mean3 .- fn2mean3;
bdiffmean3 = bn1mean3 .- bn2mean3;
fdiffmeanbin3 = fdiffmean3 .> 0.000001;
bdiffmeanbin3 = bdiffmean3 .> 0.000001;
hystdiff3 = fdiffmeanbin3 .!= bdiffmeanbin3;

namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/FinalDraft3/fig_hysteresis_ddm_theta.pdf");
R"""
pdf($namespace,height=5,width=12)
par(mfrow=c(1,3))
image(x=$fmvec,y=$cexpvec,z=$fdiffmeanbin1,xlim=c(0,0.5),xlab='Straying ratio m',ylab=expression(paste(log[10],' C')),col=c('white','black'),main=expression(paste('Low heterogeneity (',Delta,theta,'=1.8)',sep='')))
image(x=$fmvec,y=$cexpvec,z=$hystdiff1,xlab='Straying ratio m',ylab=expression(paste(log[10],' C')),col=c('#ffffff00','gray'),add=T)

image(x=$fmvec,y=$cexpvec,z=$fdiffmeanbin2,xlim=c(0,0.5),xlab='Straying ratio m',ylab=expression(paste(log[10],' C')),col=c('white','black'),main=expression(paste('Intermediate heterogeneity (',Delta,theta,'=2.0)',sep='')))
image(x=$fmvec,y=$cexpvec,z=$hystdiff2,xlab='Straying ratio m',ylab=expression(paste(log[10],' C')),col=c('#ffffff00','gray'),add=T)

image(x=$fmvec,y=$cexpvec,z=$fdiffmeanbin3,xlim=c(0,0.5),xlab='Straying ratio m',ylab=expression(paste(log[10],' C')),col=c('white','black'),main=expression(paste('High heterogeneity (',Delta,theta,'=2.2)',sep='')))
image(x=$fmvec,y=$cexpvec,z=$hystdiff3,xlab='Straying ratio m',ylab=expression(paste(log[10],' C')),col=c('#ffffff00','gray'),add=T)

dev.off()
"""

#Main, sub version

namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/FinalDraft3/fig_hysteresis_ddm.pdf");
R"""
pdf($namespace,height=4,width=10)
#layout(matrix(c(1,2,1,3), 2, 2, byrow = TRUE))
par(mfrow=c(1,3),mai = c(0.8, 0.7, 0.3, 0.1))

image(x=$fmvec,y=$cexpvec,z=$fdiffmeanbin1,xlim=c(0,0.5),xlab=expression(paste('Individual straying ',m[0],sep='')),ylab=expression(paste(log[10],' C')),col=c('white','black'),main=expression(paste('Low heterogeneity (',Delta,theta,'=1.8)',sep='')))
image(x=$fmvec,y=$cexpvec,z=$hystdiff1,xlab='Straying ratio m',ylab=expression(paste(log[10],' C')),col=c('#ffffff00','gray'),add=T)
text(par('usr')[1]-0.1,5.15,'(a)', xpd=TRUE,cex=1.5)

image(x=$fmvec,y=$cexpvec,z=$fdiffmeanbin2,xlim=c(0,0.5),xlab=expression(paste('Individual straying ',m[0],sep='')),ylab=expression(paste(log[10],' C')),col=c('white','black'),main=expression(paste('Intermediate heterogeneity (',Delta,theta,'=2.0)',sep='')))
image(x=$fmvec,y=$cexpvec,z=$hystdiff2,xlab='Straying ratio m',ylab=expression(paste(log[10],' C')),col=c('#ffffff00','gray'),add=T)
text(par('usr')[1]-0.1,5.15,'(b)', xpd=TRUE,cex=1.5)


image(x=$fmvec,y=$cexpvec,z=$fdiffmeanbin3,xlim=c(0,0.5),xlab=expression(paste('Individual straying ',m[0],sep='')),ylab=expression(paste(log[10],' C')),col=c('white','black'),main=expression(paste('High heterogeneity (',Delta,theta,'=2.2)',sep='')))
image(x=$fmvec,y=$cexpvec,z=$hystdiff3,xlab='Straying ratio m',ylab=expression(paste(log[10],' C')),col=c('#ffffff00','gray'),add=T)
text(par('usr')[1]-0.1,5.15,'(c)', xpd=TRUE,cex=1.5)
dev.off()
"""

#m0 linked to delta theta


@everywhere using Distributions
@everywhere using RCall
@everywhere using HDF5
@everywhere using JLD

@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolveExtinct_ddm.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/timeSS.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/qualsfunc.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/bifdet.jl")


#Analysis over m and C
tmax=10000;
reps = 20;
mvec = collect(0.0:0.01:0.5);
cexpvec = collect(1:0.01:5);
pvec = ["small","large","both"];

# n1ts_ddm = SharedArray{Float64}(reps,length(mvec),length(cexpvec),tmax);
# n2ts_ddm = SharedArray{Float64}(reps,length(mvec),length(cexpvec),tmax);
n1mean_ddm=SharedArray{Float64}(reps,length(mvec),length(cexpvec),length(pvec));
n2mean_ddm=SharedArray{Float64}(reps,length(mvec),length(cexpvec),length(pvec));
n1sd_ddm=SharedArray{Float64}(reps,length(mvec),length(cexpvec),length(pvec));
n2sd_ddm=SharedArray{Float64}(reps,length(mvec),length(cexpvec),length(pvec));
aggmean_ddm=SharedArray{Float64}(reps,length(mvec),length(cexpvec),length(pvec));
aggsd_ddm=SharedArray{Float64}(reps,length(mvec),length(cexpvec),length(pvec));
x1mean_ddm=SharedArray{Float64}(reps,length(mvec),length(cexpvec),length(pvec));
x2mean_ddm=SharedArray{Float64}(reps,length(mvec),length(cexpvec),length(pvec));
m1mean_ddm=SharedArray{Float64}(reps,length(mvec),length(cexpvec),length(pvec));
m2mean_ddm=SharedArray{Float64}(reps,length(mvec),length(cexpvec),length(pvec));
pe_ddm=SharedArray{Float64}(reps,length(mvec),length(cexpvec),length(pvec));
rt_ddm=SharedArray{Float64}(reps,length(mvec),length(cexpvec),length(pvec));

z=2;
rmax=2.0;
beta=0.001;
theta1=5.0;
thetascale=5.0;
tau=1.0;
h=0.2;
sigmaE=0;
sigmaG=1;
perror=0.01;
refuge=0.01;
t_ext = Int64(round(tmax/2));

burnin=0.80
@sync @parallel for r=1:reps

    for i=1:length(mvec)
      
      m=mvec[i];
      a0 = 1-m;
      
      thetadiff = (1-2*m)/(thetascale*m);
      
      for j=1:length(cexpvec)
          
          C = 10^(Float64(cexpvec[j]));
          
          for k=1:length(pvec)
              
              extpop=pvec[k];
              
              n1_ddm_pre, n2_ddm_pre, x1_ddm_pre, x2_ddm_pre, w1_ddm_pre, w2_ddm_pre, m1_ddm_pre, m2_ddm_pre = 
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
              );
              
              #Set n1 to be high, n2 to be low
              #This is to align dominant and subordinate populations when we take the means across reps. This doesn't matter when there are no reps.
              if n1_ddm_pre[tmax-1] < n2_ddm_pre[tmax-1]
                  n1_ddm = copy(n2_ddm_pre);
                  n2_ddm = copy(n1_ddm_pre);
                  x1_ddm = copy(x2_ddm_pre);
                  x2_ddm = copy(x1_ddm_pre);
                  w1_ddm = copy(w2_ddm_pre);
                  w2_ddm = copy(w1_ddm_pre);
                  m1_ddm = copy(m2_ddm_pre);
                  m2_ddm = copy(m1_ddm_pre);
              else
                  n1_ddm = copy(n1_ddm_pre);
                  n2_ddm = copy(n2_ddm_pre);
                  x1_ddm = copy(x1_ddm_pre);
                  x2_ddm = copy(x2_ddm_pre);
                  w1_ddm = copy(w1_ddm_pre);
                  w2_ddm = copy(w2_ddm_pre);
                  m1_ddm = copy(m1_ddm_pre);
                  m2_ddm = copy(m2_ddm_pre);
              end
              
              
              t_ss, relaxtime = timeSS(n1_ddm,n2_ddm,t_ext);
              
              rt_ddm[r,i,j,k] = relaxtime;
              
              n1trim_ddm = n1_ddm[Int64(floor(tmax*burnin)):tmax-1];
              n2trim_ddm = n2_ddm[Int64(floor(tmax*burnin)):tmax-1];
              x1trim_ddm = x1_ddm[Int64(floor(tmax*burnin)):tmax-1];
              x2trim_ddm = x2_ddm[Int64(floor(tmax*burnin)):tmax-1]
              
              # n1ts_ddm[r,i,j,:] = n1_ddm;
              # n2ts_ddm[r,i,j,:] = n2_ddm;
              
              n1mean_ddm[r,i,j,k] = mean(n1trim_ddm);
              n2mean_ddm[r,i,j,k] = mean(n2trim_ddm);
              n1sd_ddm[r,i,j,k] = std(n1trim_ddm);
              n2sd_ddm[r,i,j] = std(n2trim_ddm);
              
              aggmean_ddm[r,i,j,k] = mean(n1trim_ddm+n2trim_ddm);
              aggsd_ddm[r,i,j,k] = std(n1trim_ddm+n2trim_ddm);
              
              x1mean_ddm[r,i,j,k] = theta1-mean(x1trim_ddm);
              x2mean_ddm[r,i,j,k] = (theta1+thetadiff)-mean(x2trim_ddm);
              
              m1mean_ddm[r,i,j,k] = mean(m1_ddm[Int64(floor(tmax*burnin)):tmax-1]);
              m2mean_ddm[r,i,j,k] = mean(m2_ddm[Int64(floor(tmax*burnin)):tmax-1]);
              
              pe_ddm[r,i,j,k] = mean([(std(n1trim_ddm)/mean(n1trim_ddm)),(std(n2trim_ddm)/mean(n2trim_ddm))])*
              (1/(std(n1trim_ddm+n2trim_ddm)/mean(n1trim_ddm+n2trim_ddm)))
              
          end
          
      end
      
    end

end

save(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data3/data_4regimes_rt_pe_thetam.jld"),"n1mean_ddm",n1mean_ddm,"n2mean_ddm",n2mean_ddm,"x1mean_ddm",x1mean_ddm,"x2mean_ddm",x2mean_ddm,"m1mean_ddm",m1mean_ddm,"m2mean_ddm",m2mean_ddm,"pe_ddm",pe_ddm,"rt_ddm",rt_ddm);



d = load(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data3/data_4regimes_rt_pe_thetam.jld"));
n1mean_ddm = d["n1mean_ddm"];
n2mean_ddm = d["n2mean_ddm"];
x1mean_ddm = d["x1mean_ddm"];
x2mean_ddm = d["x2mean_ddm"];
m1mean_ddm = d["m1mean_ddm"];
m2mean_ddm = d["m2mean_ddm"];
pe_ddm = d["pe_ddm"];
rt_ddm = d["rt_ddm"];


#Take means over reps
extpop = "both"; #small large both
indext = find(x->x==extpop,pvec)[1];

mpe_ddm = mapslices(mean,pe_ddm,1)[1,:,:,indext];
mrt_ddm = mapslices(mean,rt_ddm,1)[1,:,:,indext];

namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/FinalDraft3/fig_rtpe_ddm_thetam.pdf");
R"""
library(RColorBrewer)
library(fields)
pal = rev(brewer.pal(9,"Blues"))
pal2 = brewer.pal(11,"Spectral")
pdf($namespace,height=4,width=6)
par(mfrow=c(1,1),mai = c(0.8, 0.8, 0.3, 0.9))
#image.plot(x=$mvec,y=$cexpvec,z=$mpe_ddm,zlim=c(0,8),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
#text(x=0.25,y=5.2,'PE',xpd=T)
image.plot(x=$mvec,y=$cexpvec,z=log($mrt_ddm,10),xlim=c(0,0.3),zlim=c(min(log($mrt_ddm,10)),2.0),col=pal2,xlab=expression(paste(m[0])),ylab='')
text(x=0.15,y=5.2,expression(paste(log[10],' Recovery time')),xpd=T)
dev.off()
"""

R"""
par(mfrow=c(2,2))
plot($mvec,$(mrt_ddm[:,1]),type='l',log='y')
plot($mvec,$(mrt_ddm[:,101]),type='l',log='y')
plot($mvec,$(mrt_ddm[:,201]),type='l',log='y')
plot($mvec,$(mrt_ddm[:,301]),type='l',log='y')
"""

R"""
plot($(mrt_ddm[:,401]),$(mpe_ddm[:,401]),log='xy')
"""


#Small/large collapse

#Take means over reps
extpop = "small"; #small large both
indext = find(x->x==extpop,pvec)[1];

mpe_ddm_s = mapslices(mean,pe_ddm,1)[1,:,:,indext];
mrt_ddm_s = mapslices(mean,rt_ddm,1)[1,:,:,indext];


#Take means over reps
extpop = "large"; #small large both
indext_l = find(x->x==extpop,pvec)[1];

mpe_ddm_l = mapslices(mean,pe_ddm,1)[1,:,:,indext_l];
mrt_ddm_l = mapslices(mean,rt_ddm,1)[1,:,:,indext_l];


namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/FinalDraft3/fig_rtpe_ddmsl.pdf");
R"""
library(RColorBrewer)
library(fields)
pal = rev(brewer.pal(9,"Blues"))
pal2 = brewer.pal(11,"Spectral")
pdf($namespace,height=10,width=10)
par(mfrow=c(2,2),mai = c(0.8, 0.8, 0.3, 0.9))
image.plot(x=$mvec,y=$cexpvec,z=$mpe_ddm_s,zlim=c(min($(mpe_ddm_s[!isnan(mpe_ddm_s)])),2),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
text(x=0.25,y=5.2,'Portfolio effect',xpd=T)
text(x=0.4,y=4.9,'Subordinate extinct',xpd=T,col='white')
image.plot(x=$mvec,y=$cexpvec,z=log($mrt_ddm_s,10),zlim=c(min(log($mrt_ddm_s,10)),2.2),col=pal2,xlab=expression(paste(m[0])),ylab='')
text(x=0.25,y=5.2,expression(paste(log[10],' Recovery time')),xpd=T)
text(x=0.4,y=4.9,'Subordinate extinct',xpd=T,col='black')
image.plot(x=$mvec,y=$cexpvec,z=$mpe_ddm_l,zlim=c(min($(mpe_ddm_l[!isnan(mpe_ddm_l)])),2),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
text(x=0.25,y=5.2,'Portfolio effect',xpd=T)
text(x=0.4,y=4.9,'Dominant extinct',xpd=T,col='white')
image.plot(x=$mvec,y=$cexpvec,z=log($mrt_ddm_l,10),zlim=c(min(log($mrt_ddm_l,10)),2.2),col=pal2,xlab=expression(paste(m[0])),ylab='')
text(x=0.25,y=5.2,expression(paste(log[10],' Recovery time')),xpd=T)
text(x=0.4,y=4.9,'Dominant extinct',xpd=T,col='black')
dev.off()
"""

namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/FinalDraft3/fig_rtpe_ddmsl.pdf");
R"""
library(RColorBrewer)
library(fields)
pal = rev(brewer.pal(9,"Blues"))
pal2 = brewer.pal(11,"Spectral")
pdf($namespace,height=5,width=10)
par(mfrow=c(1,2),mai = c(0.8, 0.8, 0.3, 0.9))
image.plot(x=$mvec,y=$cexpvec,z=log($mrt_ddm_s,10),zlim=c(min(log($mrt_ddm_s,10)),2.2),col=pal2,xlab=expression(paste(m[0])),ylab='')
text(x=0.25,y=5.2,expression(paste(log[10],' Recovery time')),xpd=T)
image.plot(x=$mvec,y=$cexpvec,z=log($mrt_ddm_l,10),zlim=c(min(log($mrt_ddm_l,10)),2.2),col=pal2,xlab=expression(paste(m[0])),ylab='')
text(x=0.25,y=5.2,expression(paste(log[10],' Recovery time')),xpd=T)
dev.off()
"""


#Take means over reps
extpop = "both"; #small large both
indext = find(x->x==extpop,pvec)[1];
n1mean_b = mapslices(mean,n1mean_ddm,1)[1,:,:,indext];
n2mean_b = mapslices(mean,n2mean_ddm,1)[1,:,:,indext];
#Take means over reps
extpop = "small"; #small large both
indext = find(x->x==extpop,pvec)[1];
n1mean_s = mapslices(mean,n1mean_ddm,1)[1,:,:,indext];
n2mean_s = mapslices(mean,n2mean_ddm,1)[1,:,:,indext];
#Take means over reps
extpop = "large"; #small large both
indext = find(x->x==extpop,pvec)[1];
n1mean_l = mapslices(mean,n1mean_ddm,1)[1,:,:,indext];
n2mean_l = mapslices(mean,n2mean_ddm,1)[1,:,:,indext];

#Compare final steady state between 3 disturbance regimes
R"""
library(RColorBrewer)
library(fields)
pal = rev(brewer.pal(9,"Blues"))
par(mfrow=c(3,2),mai = c(0.8, 0.8, 0.3, 0.9))
image.plot(x=$mvec,y=$cexpvec,z=log($n1mean_b,10),zlim=c(1,2.5),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
image.plot(x=$mvec,y=$cexpvec,z=log($n2mean_b,10),zlim=c(1,2.5),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
image.plot(x=$mvec,y=$cexpvec,z=log($n1mean_s,10),zlim=c(1,2.5),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
image.plot(x=$mvec,y=$cexpvec,z=log($n2mean_s,10),zlim=c(1,2.5),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
image.plot(x=$mvec,y=$cexpvec,z=log($n1mean_l,10),zlim=c(1,2.5),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
image.plot(x=$mvec,y=$cexpvec,z=log($n2mean_l,10),zlim=c(1,2.5),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
"""

R"""
par(mfrow = c(1,1))
image.plot(x=$mvec,y=$cexpvec,z=abs($n1mean_b-$n2mean_b),xlim=c(0,0.3),zlim=c(0,500),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
"""


R"""
par(mfrow = c(1,2))
image.plot(x=$mvec,y=$cexpvec,z=$n1mean_b,xlim=c(0,0.3),zlim=c(0,500),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
image.plot(x=$mvec,y=$cexpvec,z=$n2mean_b,xlim=c(0,0.3),zlim=c(0,500),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
"""


R"""
library(fields)
library(RColorBrewer)
pal = rev(brewer.pal(9,"Blues"))
par(mfrow = c(3,2))
image.plot(x=$mvec,y=$cexpvec,z=abs($n1mean_b),zlim=c(0,400),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
image.plot(x=$mvec,y=$cexpvec,z=abs($n2mean_b),zlim=c(0,400),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
image.plot(x=$mvec,y=$cexpvec,z=abs($n1mean_l),zlim=c(0,400),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
image.plot(x=$mvec,y=$cexpvec,z=abs($n2mean_l),zlim=c(0,400),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
image.plot(x=$mvec,y=$cexpvec,z=abs($n1mean_s),zlim=c(0,400),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
image.plot(x=$mvec,y=$cexpvec,z=abs($n2mean_s),zlim=c(0,400),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
"""


R"""
par(mfrow = c(3,1))
image.plot(x=$mvec,y=$cexpvec,z=abs($n1mean_b+$n2mean_b),zlim=c(0,500),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
image.plot(x=$mvec,y=$cexpvec,z=abs($n1mean_s+$n2mean_s),zlim=c(0,500),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
image.plot(x=$mvec,y=$cexpvec,z=abs($n1mean_l+$n2mean_l),zlim=c(0,500),col=pal,xlab=expression(paste(m[0])),ylab=expression(paste(log[10],' C')))
"""



#Bifurcation Plot as a function of m and C when m0 ~ thetadiff


# 2-dimensional search

@everywhere using RCall, Distributions, HDF5, JLD
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolveSS_ddm.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinJacobian.jl")

#Analysis over m
tmax=10000;
thetadiffvec1 = collect(0.0:0.01:3.0);
thetadiffvec = [thetadiffvec1 ; reverse(thetadiffvec1)];
cexpvec = collect(1:0.01:5);

# n1ts = Array{Float64}(length(mvec),length(cexpvec),tmax);
# n2ts = Array{Float64}(length(mvec),length(cexpvec),tmax);
n1mean=SharedArray{Float64}(length(thetadiffvec),length(cexpvec));
n2mean=SharedArray{Float64}(length(thetadiffvec),length(cexpvec));
n1sd=SharedArray{Float64}(length(thetadiffvec),length(cexpvec));
n2sd=SharedArray{Float64}(length(thetadiffvec),length(cexpvec));
aggmean=SharedArray{Float64}(length(thetadiffvec),length(cexpvec));
aggsd=SharedArray{Float64}(length(thetadiffvec),length(cexpvec));
x1mean=SharedArray{Float64}(length(thetadiffvec),length(cexpvec));
x2mean=SharedArray{Float64}(length(thetadiffvec),length(cexpvec));
pe=SharedArray{Float64}(length(thetadiffvec),length(cexpvec));

z=2.0;
rmax=2.0;
beta=0.001;
theta1=5.0;
thetascale=20.0;

tau=1;
h=0.2;
sigmaE=0;
sigmaG=1;
perror=0.00;

burnin=0.80
@sync @parallel for r=1:length(cexpvec)

    C = 10^(Float64(cexpvec[r]));

    for i=1:length(thetadiffvec)
      
      thetadiff = thetadiffvec[i];
      m=1/(2+thetascale*thetadiff);
      a0 = 1-m;
      
      if i == 1
          n0 = [2,2];
        #   x0 = [theta1 + rand(Normal(0,0.0001)),(theta1 + thetadiff)  + rand(Normal(0,0.0001))];
        x0 = [theta1 ,(theta1 + thetadiff)];
      else
          n0 = [n1mean[i-1,r],n2mean[i-1,r]];
          x0 = [x1mean[i-1,r],x2mean[i-1,r]];
      end
      
      n1, n2, x1, x2, w1, w2, m1, m2 = 
      KevinEvolveSS_ddm(
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
      
      # n1ts[i,r,:] = n1;
      # n2ts[i,r,:] = n2;
      
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

save(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data3/data_hysteresis_ddm_mtheta.jld"),"n1mean",n1mean,"n2mean",n2mean,"x1mean",x1mean,"x2mean",x2mean,"pe",pe);

d = load(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data3/data_hysteresis_ddm_mtheta.jld"));
n1mean = d["n1mean"];
n2mean = d["n2mean"];

#Plot the difference in n1mean, n2mean
midpoint = Int64(floor(length(thetadiffvec)/2));

# fmaxeigs = maxeigs[1:midpoint,:];
# bmaxeigs = flipdim(maxeigs[midpoint+1:length(mvec),:],1);
# fold = 0.98 .< fmaxeigs .< 1.00;
# bold = 0.98 .< bmaxeigs .< 1.00;

fn1mean = n1mean[1:midpoint,:];
bn1mean = flipdim(n1mean[midpoint+1:length(thetadiffvec),:],1);

fn2mean = n2mean[1:midpoint,:];
bn2mean = flipdim(n2mean[midpoint+1:length(thetadiffvec),:],1);

fthetadiffvec = thetadiffvec[1:midpoint];
bthetadiffvec = fthetadiffvec;

fdiffmean = fn1mean .- fn2mean;
bdiffmean = bn1mean .- bn2mean;
fdiffmeanbin = fdiffmean .> 0.000001;
bdiffmeanbin = bdiffmean .> 0.000001;

hystdiff = fdiffmeanbin .!= bdiffmeanbin;


thetascale=20;
thetadiff=collect(0:0.01:3.0);
mt = zeros(length(thetadiff));
for i=1:length(thetadiff)
    mt[i] = 1/(2+thetascale*thetadiff[i])
end


namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/Manuscript/FinalDraft3/fig_mthetarelation.pdf");
R"""
pdf($namespace,height=5,width=6)
plot($thetadiff,$m,type='l',ylab='Straying rate m',xlab=expression(paste('Habitat heterogeneity ',Delta,theta)),ylim=c(0,0.5),xlim=c(0,3))
dev.off()
"""


namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/FinalDraft3/fig_hysteresis_ddm_mtheta.pdf");
R"""
pdf($namespace,height=5,width=6)

par(fig = c(0,1,0,1),mai = c(0.9, 0.9, 0.2, 0.4))
image(x=$fthetadiffvec,y=$cexpvec,z=$fdiffmeanbin,xlim=c(3,1.0),xlab=expression(paste('Habitat heterogeneity ',Delta,theta,sep='')),ylab=expression(paste(log[10],' C')),col=c('white','black'))
image(x=$fthetadiffvec,y=$cexpvec,z=$hystdiff,xlim=c(3,1.0),xlab='Straying ratio m',ylab=expression(paste(log[10],' C')),col=c('#ffffff00','gray'),add=T)
text(x=3,y=0.3,expression(paste('low ',m[0])),xpd=T)
text(x=1,y=0.3,expression(paste('high ',m[0])),xpd=T)

par(fig = c(0.5,0.98, 0.53, 0.96), new = T)
plot($thetadiff,$mt,type='l',ylab='',xlab='',ylim=c(0.01,0.05),xlim=c(3,1),)
text(x=2,y=-0.014,expression(paste(Delta,theta)),xpd=T)
text(x=3.8,y=0.03,expression(paste(m[0])),xpd=T)
dev.off()
"""



# Just habitat heterogeneity
@everywhere using Distributions
@everywhere using RCall
@everywhere using HDF5
@everywhere using JLD

@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolveExtinct_ddm.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/timeSS.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/qualsfunc.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/bifdet.jl")


#Analysis over m and C
tmax=10000;
reps = 50;
pvec = ["small","large","both"];
thetavec = collect(1.5:0.01:2.5);

# n1ts_ddm = SharedArray{Float64}(reps,length(mvec),length(cexpvec),tmax);
# n2ts_ddm = SharedArray{Float64}(reps,length(mvec),length(cexpvec),tmax);
n1mean_ddm=SharedArray{Float64}(reps,length(thetavec),length(pvec));
n2mean_ddm=SharedArray{Float64}(reps,length(thetavec),length(pvec));
n1sd_ddm=SharedArray{Float64}(reps,length(thetavec),length(pvec));
n2sd_ddm=SharedArray{Float64}(reps,length(thetavec),length(pvec));
aggmean_ddm=SharedArray{Float64}(reps,length(thetavec),length(pvec));
aggsd_ddm=SharedArray{Float64}(reps,length(thetavec),length(pvec));
x1mean_ddm=SharedArray{Float64}(reps,length(thetavec),length(pvec));
x2mean_ddm=SharedArray{Float64}(reps,length(thetavec),length(pvec));
m1mean_ddm=SharedArray{Float64}(reps,length(thetavec),length(pvec));
m2mean_ddm=SharedArray{Float64}(reps,length(thetavec),length(pvec));
pe_ddm=SharedArray{Float64}(reps,length(thetavec),length(pvec));
rt_ddm=SharedArray{Float64}(reps,length(thetavec),length(pvec));

z=2;
rmax=2.0;
beta=0.001;
theta1=5.0;
# thetadiff=2.0;
tau=1.0;
h=0.2;
sigmaE=0;
sigmaG=1;
perror=0.01;
refuge=0.01;
t_ext = Int64(round(tmax/2));
C = 10^2.5;

burnin=0.80
@sync @parallel for r=1:reps

    m=rand(Normal(0.1,0.05));
    a0 = 1-m;

    for i=1:length(thetavec)
      
        thetadiff = thetavec[i];
      
          for k=1:length(pvec)
              
              extpop=pvec[k];
              
              n1_ddm_pre, n2_ddm_pre, x1_ddm_pre, x2_ddm_pre, w1_ddm_pre, w2_ddm_pre, m1_ddm_pre, m2_ddm_pre = 
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
              );
              
              #Set n1 to be high, n2 to be low
              #This is to align dominant and subordinate populations when we take the means across reps. This doesn't matter when there are no reps.
              if n1_ddm_pre[tmax-1] < n2_ddm_pre[tmax-1]
                  n1_ddm = copy(n2_ddm_pre);
                  n2_ddm = copy(n1_ddm_pre);
                  x1_ddm = copy(x2_ddm_pre);
                  x2_ddm = copy(x1_ddm_pre);
                  w1_ddm = copy(w2_ddm_pre);
                  w2_ddm = copy(w1_ddm_pre);
                  m1_ddm = copy(m2_ddm_pre);
                  m2_ddm = copy(m1_ddm_pre);
              else
                  n1_ddm = copy(n1_ddm_pre);
                  n2_ddm = copy(n2_ddm_pre);
                  x1_ddm = copy(x1_ddm_pre);
                  x2_ddm = copy(x2_ddm_pre);
                  w1_ddm = copy(w1_ddm_pre);
                  w2_ddm = copy(w2_ddm_pre);
                  m1_ddm = copy(m1_ddm_pre);
                  m2_ddm = copy(m2_ddm_pre);
              end
              
              
              t_ss, relaxtime = timeSS(n1_ddm,n2_ddm,t_ext);
              
              rt_ddm[r,i,k] = relaxtime;
              
              n1trim_ddm = n1_ddm[Int64(floor(tmax*burnin)):tmax-1];
              n2trim_ddm = n2_ddm[Int64(floor(tmax*burnin)):tmax-1];
              x1trim_ddm = x1_ddm[Int64(floor(tmax*burnin)):tmax-1];
              x2trim_ddm = x2_ddm[Int64(floor(tmax*burnin)):tmax-1]
              
              # n1ts_ddm[r,i,j,:] = n1_ddm;
              # n2ts_ddm[r,i,j,:] = n2_ddm;
              
              n1mean_ddm[r,i,k] = mean(n1trim_ddm);
              n2mean_ddm[r,i,k] = mean(n2trim_ddm);
              n1sd_ddm[r,i,k] = std(n1trim_ddm);
              n2sd_ddm[r,i,k] = std(n2trim_ddm);
              
              aggmean_ddm[r,i,k] = mean(n1trim_ddm+n2trim_ddm);
              aggsd_ddm[r,i,k] = std(n1trim_ddm+n2trim_ddm);
              
              x1mean_ddm[r,i,k] = theta1-mean(x1trim_ddm);
              x2mean_ddm[r,i,k] = (theta1+thetadiff)-mean(x2trim_ddm);
              
              m1mean_ddm[r,i,k] = mean(m1_ddm[Int64(floor(tmax*burnin)):tmax-1]);
              m2mean_ddm[r,i,k] = mean(m2_ddm[Int64(floor(tmax*burnin)):tmax-1]);
              
              pe_ddm[r,i,k] = mean([(std(n1trim_ddm)/mean(n1trim_ddm)),(std(n2trim_ddm)/mean(n2trim_ddm))])*
              (1/(std(n1trim_ddm+n2trim_ddm)/mean(n1trim_ddm+n2trim_ddm)))
              
          end
      
    end

end

save(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data3/data_4regimes_theta.jld"),"n1mean_ddm",n1mean_ddm,"n2mean_ddm",n2mean_ddm,"x1mean_ddm",x1mean_ddm,"x2mean_ddm",x2mean_ddm,"m1mean_ddm",m1mean_ddm,"m2mean_ddm",m2mean_ddm,"pe_ddm",pe_ddm,"rt_ddm",rt_ddm);


d = load(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data3/data_4regimes_theta.jld"));
n1mean_ddm = d["n1mean_ddm"];
n2mean_ddm = d["n2mean_ddm"];
x1mean_ddm = d["x1mean_ddm"];
x2mean_ddm = d["x2mean_ddm"];
m1mean_ddm = d["m1mean_ddm"];
m2mean_ddm = d["m2mean_ddm"];
pe_ddm = d["pe_ddm"];
rt_ddm = d["rt_ddm"];


#Take means over reps
extpop = "small"; #small large both
indext = find(x->x==extpop,pvec)[1];
mn1mean_ddm_s = mapslices(mean,n1mean_ddm,1)[1,:,indext];
mpe_ddm_s = mapslices(mean,pe_ddm,1)[1,:,indext];
mrt_ddm_s = mapslices(mean,rt_ddm,1)[1,:,indext];

extpop = "large"; #small large both
indext = find(x->x==extpop,pvec)[1];
mn1mean_ddm_l = mapslices(mean,n1mean_ddm,1)[1,:,indext];
mpe_ddm_l = mapslices(mean,pe_ddm,1)[1,:,indext];
mrt_ddm_l = mapslices(mean,rt_ddm,1)[1,:,indext];

extpop = "both"; #small large both
indext = find(x->x==extpop,pvec)[1];
mn1mean_ddm = mapslices(mean,n1mean_ddm,1)[1,:,indext];
mpe_ddm = mapslices(mean,pe_ddm,1)[1,:,indext];
mrt_ddm = mapslices(mean,rt_ddm,1)[1,:,indext];

R"""
plot($thetavec,$mrt_ddm,log='y',pch=16,col='black',ylim=c(10,250))
points($thetavec,$mrt_ddm_s,pch=16,col='lightgray')
points($thetavec,$mrt_ddm_l,pch=16,col='darkgray')
"""

R"""
plot($thetavec,$mpe_ddm,pch=16,col='black',ylim=c(1,3))
points($thetavec,$mpe_ddm_s,pch=16,col='lightgray')
points($thetavec,$mpe_ddm_l,pch=16,col='darkgray')
"""
