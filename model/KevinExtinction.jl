using Distributions
using RCall
using JLD
using HDF5

include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolve.jl")
include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolve_ddm.jl")


include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolveExtinct.jl")
include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/timeSS.jl")
include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/movingaverage.jl")



tmax=10000;
z=0.5;
rmax=2.0;
beta=0.001;
theta1=5.0;
thetadiff=5.0;
tau=1.0;
h=0.5;
sigmaE=0;
sigmaG=1;
perror=0.01;

m=0.45

extpop = "small";
t_ext = Int64(round(tmax/2));

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


R"""
library(RColorBrewer)
pal = brewer.pal(3,'Set1')
par(mfrow=c(2,1))
plot($n1,type='l',col=pal[1])
lines($n2,col=pal[2])
plot($x1,type='l',col=pal[1],ylim=c(4,12))
lines($x2,col=pal[2])
"""


t_ss, relaxtime = timeSS(n1,n2,t_ext);


R"""
library(RColorBrewer)
pal = brewer.pal(3,'Set1')
plot($n1,type='l',col=pal[1],xlim=c($t_ext,$t_ss+50))
lines($n2,col=pal[2])
lines(c($t_ss,$t_ss),c(0,10000))
"""



#Calculate relaxation time as a function of straying rate and large v small


#Analysis over m
tmax=10000;

mvec = collect(0.0001:0.001:0.5);
pvec = ["small","large","both"];


z=0.5;
rmax=2.0;
beta=0.001;
theta1=5.0;
thetadiff=5;
tau=1.0;
h=0.8;
sigmaE=0;
sigmaG=1;
perror=0.01;
refuge=0.01;

# n1v = Array{Float64}(length(mvec),tmax);
# n2v = Array{Float64}(length(mvec),tmax);

t_ext = Int64(round(tmax/2));

reps = 100;
RT = Array(Array{Float64},reps);

@time for r=1:reps
  
  rt = zeros(Float64,length(mvec),length(pvec));
  
  for i=1:length(mvec)
    
    m=mvec[i];
    
    for j=1:length(pvec);
      
      extpop = pvec[j];
      
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
      
      rt[i,j] = relaxtime;
      
      # if j==1
      #   n1v[i,:] = n1;
      #   n2v[i,:] = n2;
      # end
      
    end
  end
  
  RT[r] = rt;

end


ma_rt = mean(RT[1:reps]);

namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/figs2/fig_relax.pdf");
R"""
library(RColorBrewer)
pdf($namespace,height=4,width=5)
pal = brewer.pal(3,'Set1')
plot($mvec,$(ma_rt[:,1]),col=pal[1],type='l',log='y',cex=0.5,lwd=2,xlab='m',ylab='Return time')
lines($mvec,$(ma_rt[:,2]),col=pal[2],cex=0.5,lwd=2)
lines($mvec,$(ma_rt[:,3]),col=pal[3],cex=0.5,lwd=2)
dev.off()
"""


R"""
par(mfrow=c(2,1))
plot($(n1v[426,:]),type='l',col=pal[1],ylim=c(0,500),xlim=c(5000,5100))
lines($(n2v[426,:]),col=pal[2])
lines(c($(t_ext+rt[426,1]),$(t_ext+rt[426,1])),c(0,10000))
plot($(n1v[427,:]),type='l',col=pal[1],ylim=c(0,500),xlim=c(5000,5100))
lines($(n2v[427,:]),col=pal[2])
lines(c($(t_ext+rt[427,1]),$(t_ext+rt[427,1])),c(0,10000))
"""
