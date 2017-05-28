using Distributions
using RCall

include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolve.jl")
include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolve_ddm.jl")

#Analysis over m
tmax=500000;
mvec = collect(0.0:0.01:0.5);

m1mean=zeros(Float64,length(mvec));
m2mean=zeros(Float64,length(mvec));
pe=zeros(Float64,length(mvec));
pe_ddm=zeros(Float64,length(mvec));

z=0.5;
rmax=2.0;
beta=0.001;
theta1=5.0;
thetadiff=5;
tau=1.0;
h=0.2;
sigmaE=0;
sigmaG=1;
perror=0.01;
C=1000;

burnin=0.80
@time for i=1:length(mvec)
  m=mvec[i];
  a0 = 1-m;
  
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
  
  pe[i] = mean([(std(n1trim)/mean(n1trim)),(std(n2trim)/mean(n2trim))])*
  (1/(std(n1trim+n2trim)/mean(n1trim+n2trim)))
  
  
  
  n1_ddm, n2_ddm, x1_ddm, x2_ddm, w1_ddm, w2_ddm, m1, m2 = 
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
    
    n1trim_ddm = n1_ddm[Int64(floor(tmax*burnin)):tmax];
    n2trim_ddm = n2_ddm[Int64(floor(tmax*burnin)):tmax];
    
    m1trim = m1[Int64(floor(tmax*burnin)):tmax-1];
    m2trim = m2[Int64(floor(tmax*burnin)):tmax-1];
    m1mean[i] = mean(m1trim);
    m2mean[i] = mean(m2trim);
    
    pe_ddm[i] = mean([(std(n1trim_ddm)/mean(n1trim_ddm)),(std(n2trim_ddm)/mean(n2trim_ddm))])*
  (1/(std(n1trim_ddm+n2trim_ddm)/mean(n1trim_ddm+n2trim_ddm)))
  
end

namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/figs2/fig_PEmvm.pdf");
R"""
pdf($namespace,height=5,width=6)
plot($mvec,$pe,type='l',lwd=2,xlab='Straying rate m, m*',ylab='PE')
points($m1mean,$pe_ddm,pch=16,cex=0.8)
points($m2mean,$pe_ddm,pch=16,cex=0.8)
for (i in 1:length($mvec)) {
  segments($m1mean[i],$pe_ddm[i],$m2mean[i],$pe_ddm[i],col=paste('#00000050'))
}
dev.off()
"""

