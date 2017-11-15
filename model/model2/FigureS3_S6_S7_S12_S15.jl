#Figure S12

#Note Figures S6 and S7 have to be recreated and renamed!


using Distributions
using RCall

include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolve_mtheta.jl")


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
thetascale=2.0;
tau=1.0;
h=0.2;
sigmaE=0;
sigmaG=1;
perror=0.01;
refuge=0.01;
m=0.02;

thetadiff = (1-2*m)/(thetascale*m);

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
t_ext,
refuge
);

t_ss, relaxtime = timeSS(n1,n2,t_ext);


R"""
library(RColorBrewer)
pal = brewer.pal(3,'Set1')
par(mfrow=c(2,1),mai = c(0.8, 0.8, 0.4, 0.1))
plot($n1,type='l',col=pal[1],xlim=c($t_ext,$t_ss+50),ylim=c(0,max($n1,$n2)),ylab='Population density',xlab='Time',main=paste(c('RT=',$relaxtime)))
lines($n2,col=pal[2])
lines(c($t_ss,$t_ss),c(-10,10000))
plot($x1,type='l',col=pal[1],ylim=c(2,100),xlim=c($t_ext,$t_ss+50),ylab='Trait means',xlab='Time')
lines($x2,col=pal[2])
lines(c($t_ss,$t_ss),c(0,10000))
"""

namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/Manuscript/figs2/fig_relax_inertia.pdf");
R"""
library(RColorBrewer)
pal = brewer.pal(3,'Set1')
pdf($namespace,height=8,width=6)
par(mfrow=c(2,1),mai = c(0.8, 0.8, 0.4, 0.1))
plot($n1,type='l',col=pal[1],xlim=c($t_ext,$t_ss+50),ylim=c(0,max($n1,$n2)),ylab='Population density',xlab='Time',main=paste(c('High-density extinction; recovery time=',$relaxtime)))
lines($n2,col=pal[2])
lines(c($t_ss,$t_ss),c(-10,10000))
plot($x1,type='l',col=pal[1],ylim=c(0,40),xlim=c($t_ext,$t_ss+50),ylab='Trait means',xlab='Time')
lines($x2,col=pal[2])
lines(c($t_ss,$t_ss),c(0,10000))
dev.off()
"""

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
h=0.8;
sigmaE=0;
sigmaG=1;
perror=0.01;
refuge=0.01;
m=0.4;

extpop = "both";
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

t_ss, relaxtime = timeSS(n1,n2,t_ext);

namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/Manuscript/figs2/fig_recovery.pdf");
R"""
library(RColorBrewer)
pal = brewer.pal(3,'Set1')
pdf($namespace,height=5,width=6)
par(mfrow=c(1,1),mai = c(0.8, 0.8, 0.4, 0.1))
plot($n1,type='l',col=pal[1],xlim=c($t_ext,$t_ss+50),ylim=c(0,max($n1+$n2)),ylab='Population density',xlab='Time')
lines($n2,col=pal[2])
lines($n1+$n2,col='black')
lines(c($t_ss,$t_ss),c(-10,10000))
points($t_ss,$(n1[t_ss])+$(n2[t_ss]),pch=16)
text(x=$t_ss+18,y=1,paste(c('Recovery time=',$relaxtime),collapse=''))
#plot($x1,type='l',col=pal[1],ylim=c(5,10),xlim=c($t_ext,$t_ss+50),ylab='Trait means',xlab='Time')
#lines($x2,col=pal[2])
#lines(c($t_ss,$t_ss),c(0,10000))
dev.off()
"""
