


#############
# Calculate return time as a function of 'm' and 'Z'
#############



@everywhere using Distributions
@everywhere using RCall
@everywhere using JLD
@everywhere using HDF5

@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolve.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolve_ddm.jl")


@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolveExtinct.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolveExtinct_ddm.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/timeSS.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/movingaverage.jl")


#Analysis over m
tmax=10000;

mvec = collect(0.0001:0.001:0.3);
pvec = ["both"];
h = 0.2;
zvec = logspace(-1,log10(5),200);
reps = 50;

# z=0.5;
rmax=2.0;
beta=0.001;
theta1=5.0;
thetadiff=5;
tau=1.0;
C=1000;
sigmaE=0;
sigmaG=1;
perror=0.01;
refuge=0.01;

burnin = 0.8;

# n1v = SharedArray{Float64}(reps,length(mvec),length(pvec),tmax);
# n2v = SharedArray{Float64}(reps,length(mvec),length(pvec),tmax);
# x1v = SharedArray{Float64}(reps,length(mvec),length(pvec),tmax);
# x2v = SharedArray{Float64}(reps,length(mvec),length(pvec),tmax);

t_ext = Int64(round(tmax/2));


rt = SharedArray{Float64}(reps,length(mvec),length(pvec),length(zvec));
rt_ddm = SharedArray{Float64}(reps,length(mvec),length(pvec),length(zvec));
m1mean = SharedArray{Float64}(reps,length(mvec),length(pvec),length(zvec));
m2mean = SharedArray{Float64}(reps,length(mvec),length(pvec),length(zvec));

@sync @parallel for r=1:reps
  
  for i=1:length(mvec)
    
    m=mvec[i]; #constant straying rate
    a0 = 1-m; #individual homing rate
    
    for j=1:length(pvec);
      
      extpop = pvec[j];
      
      for k=1:length(zvec)
        
        z = zvec[k];
        
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
        
        rt[r,i,j,k] = relaxtime;
        

        # #Do the same thing for density-dependent m
        # 
        # 
        # n1_ddm, n2_ddm, x1_ddm, x2_ddm, w1_ddm, w2_ddm, m1_ddm, m2_ddm = 
        # KevinEvolveExtinct_ddm(
        # tmax, 
        # z, 
        # rmax,
        # beta,
        # theta1,
        # thetadiff,
        # tau,
        # h,
        # a0,
        # C,
        # sigmaE,
        # sigmaG,
        # perror,
        # extpop,
        # t_ext,
        # refuge
        # )
        # 
        # t_ss_ddm, relaxtime_ddm = timeSS(n1_ddm,n2_ddm,t_ext);
        # 
        # rt_ddm[r,i,j,k] = relaxtime_ddm;
        # 
        # #Steady state stray rate?
        # m1trim = m1_ddm[Int64(floor(tmax*burnin)):tmax-1];
        # m2trim = m2_ddm[Int64(floor(tmax*burnin)):tmax-1];
        # 
        # m1mean[r,i,j,k] = mean(m1trim);
        # m2mean[r,i,j,k] = mean(m2trim);
        
        
      end
    end
  end
end

save(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data3/data_relaxZ.jld"),"rt",rt,"rt_ddm",rt_ddm,"m1mean",m1mean,"m2mean",m2mean);

d = load(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data3/data_relaxZ.jld"));
#This loads the dictionary
rt = d["rt"];
rt_ddm = d["rt_ddm"];
m1mean = d["m1mean"];
m2mean = d["m2mean"];



namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/FinalDraft3/fig_relax_Z.pdf");
R"""
library(RColorBrewer)
library(fields)
pal = brewer.pal(11,"Spectral")
pdf($namespace,height=5,width=6)
par(mfrow=c(1,1))
image.plot(x=$mvec,y=$zvec,z=log($(mapslices(mean,rt[:,:,1,:],1)[1,:,:])),zlim=c(log(10),log(5000)),col=pal,xlab='m',ylab='Z',legend.lab='ln Recovery time')
#lines(seq(-1,1,length.out=10),rep(0.5,10),lty=3)
dev.off()
"""

namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/FinalDraft3/fig_pe_Z.pdf");
R"""
library(RColorBrewer)
library(fields)
pal = brewer.pal(11,"Spectral")
pdf($namespace,height=5,width=6)
par(mfrow=c(1,1))
image.plot(x=$mvec,y=$zvec,z=$(mapslices(mean,pe[:,:,1,:],1)[1,:,:]),zlim=c(1,5),col=pal,xlab='m',ylab='Z',main='Return time')
#lines(seq(-1,1,length.out=10),rep(0.5,10),lty=3)
dev.off()
"""


avgm = (mean([mapslices(mean,m2mean[:,:,3,:],1)[1,:,:],mapslices(mean,m1mean[:,:,3,:],1)[1,:,:]]));
avgm[find(x->x==true,isnan.(avgm))]=0.5;

ma1_m_ddm = mean([m1mean[i,:,3,:] for i=1:reps]);
ma2_m_ddm = mean([m2mean[i,:,3,:] for i=1:reps]);

namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/figs2/fig_relax_Z.pdf");
R"""
library(RColorBrewer)
pal = brewer.pal(9,"YlOrRd")
#pdf($namespace,height=3,width=12)*/
par(mfrow=c(1,1))
image(x=$(avgm),y=$zvec,z=log($(mapslices(mean,rt_ddm[:,:,3,:],1)[1,:,:])),zlim=c(log(15),log(5000)),col=pal,xlab='m',ylab='Z',main='Return time')
lines(seq(-1,1,length.out=10),rep(0.5,10),lty=3)
#dev.off()
"""
