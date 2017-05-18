using Distributions
using RCall
using JLD
using HDF5

include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolve.jl")
include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolve_ddm.jl")


include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolveExtinct.jl")
include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/timeSS.jl")




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

m=0.2

extpop = "large";
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
t_ext
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



#Calculate relaxation time as a function of 

