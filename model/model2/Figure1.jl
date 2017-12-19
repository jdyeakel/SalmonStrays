#Figure 1
@everywhere using Distributions
@everywhere using RCall

@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolve.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinJacobian.jl")

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
maximeigs = Array{Float64}(length(mvec));

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
@time for i=1:length(mvec)
  m=mvec[i];
  
  n1_pre, n2_pre, x1_pre, x2_pre, w1, w2 = 
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
    
    # if n1_pre[tmax-1] < n2_pre[tmax-1]
        n1 = copy(n2_pre);
        n2 = copy(n1_pre);
        x1 = copy(x2_pre);
        x2 = copy(x1_pre);
    # else
    #     n1 = copy(n1_pre);
    #     n2 = copy(n2_pre);
    #     x1 = copy(x1_pre);
    #     x2 = copy(x2_pre);
    # end
    
  n1trim = n1[Int64(floor(tmax*burnin)):tmax];
  n2trim = n2[Int64(floor(tmax*burnin)):tmax];
  x1trim = x1[Int64(floor(tmax*burnin)):tmax];
  x2trim = x2[Int64(floor(tmax*burnin)):tmax];
  
  
  
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
diffn = n1mean .- n2mean;
#Steady state plot
namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/Manuscript/FinalDraft3/fig_traj.pdf");
R"""
library(RColorBrewer)
cols = brewer.pal(4,'Set1')
pdf($namespace,height=8,width=5)
par(mfrow=c(2,1),mai = c(0.8, 0.8, 0.1, 0.1))
plot($mvec,$n1mean,pch='.',col=cols[4],xlab="Straying rate m",ylab="Steady state biomass",cex=0.5,ylim=c(0,max($n1mean)+200),las=1)
points($mvec,$n2mean,pch='.',col=cols[4],cex=0.5)
#arrows($(mvec[indmax(diffn)]),350,$(mvec[indmax(diffn)]),300,length=0.05,angle=40,lwd=3)
text($(mvec[indmax(diffn)])+0.01,200,'PB')
text(par('usr')[1]-0.09,1300,'(a)', xpd=TRUE)
lines(rep($(mvec[indmax(diffn)]),1001),seq(0,1000))
lines(rep($(mvec[indmax(diffn)])-0.05,1001),seq(0,1000))
text($(mvec[indmax(diffn)])-0.025,400,'regime I',cex=0.75)
text($(mvec[indmax(diffn)])+0.025,400,'regime II',cex=0.75)

plot($mvec,$x1mean,pch=".",col=cols[1],ylim=c(-5,5),xlab="Straying rate m",ylab="Trait offset",las=1)
points($mvec,$x2mean,pch=".",col=cols[2])
arrows($(mvec[indmax(diffn)]),3.4,$(mvec[indmax(diffn)]),2.5,length=0.05,angle=40,lwd=3)
text($(mvec[indmax(diffn)]),4,'PB')
text(par('usr')[1]-0.09,5,'(b)', xpd=TRUE)
dev.off()
"""
