# using Distributions
# using RCall
# using JLD
# using HDF5
# 
# include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolve.jl")
# include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolve_ddm.jl")
# 
# 
# include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolveExtinct.jl")
# include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/timeSS.jl")
# include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/movingaverage.jl")
# 
# 
# 
# tmax=10000;
# z=0.5;
# rmax=2.0;
# beta=0.001;
# theta1=5.0;
# thetadiff=5.0;
# tau=1.0;
# h=0.5;
# sigmaE=0;
# sigmaG=1;
# perror=0.01;
# refuge=0.01;
# m=0.25
# 
# extpop = "large";
# t_ext = Int64(round(tmax/2));
# 
# n1, n2, x1, x2, w1, w2 = 
# KevinEvolveExtinct(tmax, 
# z, 
# rmax,
# beta,
# theta1,
# thetadiff,
# tau,
# h,
# sigmaE,
# sigmaG,
# m,
# perror,
# extpop,
# t_ext,
# refuge
# );
# 
# 
# R"""
# library(RColorBrewer)
# pal = brewer.pal(3,'Set1')
# par(mfrow=c(2,1))
# plot($n1,type='l',col=pal[1],xlim=c())
# lines($n2,col=pal[2])
# plot($x1,type='l',col=pal[1],ylim=c(4,12))
# lines($x2,col=pal[2])
# """
# 
# 
# t_ss, relaxtime = timeSS(n1,n2,t_ext);
# 
# 
# R"""
# library(RColorBrewer)
# pal = brewer.pal(3,'Set1')
# par(mfrow=c(2,1),mai = c(0.8, 0.8, 0.1, 0.1))
# plot($n1,type='l',col=pal[1],xlim=c($t_ext,$t_ss+50),ylim=c(0,max($n1,$n2)),ylab='Population density',xlab='Time')
# lines($n2,col=pal[2])
# lines(c($t_ss,$t_ss),c(-10,10000))
# plot($x1,type='l',col=pal[1],ylim=c(5,10),xlim=c($t_ext,$t_ss+50),ylab='Trait means',xlab='Time')
# lines($x2,col=pal[2])
# lines(c($t_ss,$t_ss),c(0,10000))
# """


#Calculate relaxation time as a function of straying rate and large v small

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

mvec = collect(0.0001:0.001:0.5);
pvec = ["small","large","both"];
hvec = [0.2, 0.8]
reps = 100;

z=0.5;
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


rt = SharedArray{Float64}(reps,length(mvec),length(pvec),length(hvec));
rt_ddm = SharedArray{Float64}(reps,length(mvec),length(pvec),length(hvec));
m1mean = SharedArray{Float64}(reps,length(mvec),length(pvec),length(hvec));
m2mean = SharedArray{Float64}(reps,length(mvec),length(pvec),length(hvec));

@sync @parallel for r=1:reps
  
  for i=1:length(mvec)
    
    m=mvec[i]; #constant straying rate
    a0 = 1-m; #individual homing rate
    
    for j=1:length(pvec);
      
      extpop = pvec[j];
      
      for k=1:2
        
        h = hvec[k];
        
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
        

        #Do the same thing for density-dependent m
        
        
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
      
        t_ss_ddm, relaxtime_ddm = timeSS(n1_ddm,n2_ddm,t_ext);
        
        rt_ddm[r,i,j,k] = relaxtime_ddm;
        
        #Steady state stray rate?
        m1trim = m1_ddm[Int64(floor(tmax*burnin)):tmax-1];
        m2trim = m2_ddm[Int64(floor(tmax*burnin)):tmax-1];
        
        m1mean[r,i,j,k] = mean(m1trim);
        m2mean[r,i,j,k] = mean(m2trim);
        
        
      end
    end
  end
end

save(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data2/data_relax.jld"),"rt",rt,"rt_ddm",rt_ddm,"m1mean",m1mean,"m2mean",m2mean);

d = load(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data2/data_relax.jld"));
#This loads the dictionary
rt = d["rt"];
rt_ddm = d["rt_ddm"];
m1mean = d["m1mean"];
m2mean = d["m2mean"];


# save(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data2/data_relax.jld"),"n1v",n1v,"n2v",n2v,"x1v",n1v,"x2v",n2v,"rt",rt);
# 
# 
# ma_rth2 = mean([rt[i,:,:,1] for i=1:reps]);
# ma_rth8 = mean([rt[i,:,:,2] for i=1:reps]);
# 
# namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/figs2/fig_relax.pdf");
# R"""
# library(RColorBrewer)
# pdf($namespace,height=8,width=6)
# pal = brewer.pal(9,'Greys')
# palsub = pal[c(4,6,8)];
# par(mfrow=c(2,1),mai = c(0.8, 0.8, 0.1, 0.1))
# plot($mvec,$(ma_rth2[:,1]),col=palsub[1],type='l',log='y',cex=0.5,lwd=2,xlab='m',ylab='Return time',ylim=c(min($rt),max($(ma_rth2))))
# lines($mvec,$(ma_rth2[:,2]),col=palsub[2],cex=0.5,lwd=2)
# lines($mvec,$(ma_rth2[:,3]),col=palsub[3],cex=0.5,lwd=2)
# text(0.01,max($(ma_rth2)),expression(paste(h^2,'=0.2')))
# legend(x=0.44,y=20,legend=$pvec,col=palsub,pch=22,xpd=TRUE,pt.bg=palsub,cex=0.8, bty="n") #,title=expression(paste(Delta,theta))
# 
# plot($mvec,$(ma_rth8[:,1]),col=palsub[1],type='l',log='y',cex=0.5,lwd=2,xlab='m',ylab='Return time',ylim=c(min($rt),max($(ma_rth8))))
# lines($mvec,$(ma_rth8[:,2]),col=palsub[2],cex=0.5,lwd=2)
# lines($mvec,$(ma_rth8[:,3]),col=palsub[3],cex=0.5,lwd=2)
# text(0.01,max($(ma_rth8)),expression(paste(h^2,'=0.8')))
# 
# dev.off()
# """
# 
# 
# ma_rth2_ddm = mean([rt_ddm[i,:,:,1] for i=1:reps]);
# ma_rth8_ddm = mean([rt_ddm[i,:,:,2] for i=1:reps]);
# 
# #mean stray rates
# ma1_m2_ddm = mean([m1mean[i,:,:,1] for i=1:reps]);
# ma2_m2_ddm = mean([m2mean[i,:,:,1] for i=1:reps]);
# ma1_m8_ddm = mean([m1mean[i,:,:,2] for i=1:reps]);
# ma2_m8_ddm = mean([m2mean[i,:,:,2] for i=1:reps]);
# 
# 
# 
# 
# namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/figs2/fig_relax_ddm.pdf");
# R"""
# library(RColorBrewer)
# pdf($namespace,height=8,width=6)
# pal = brewer.pal(9,'Greys')
# palsub = pal[c(4,6,8)];
# par(mfrow=c(2,1),mai = c(0.8, 0.8, 0.1, 0.1))
# plot($(ma1_m2_ddm[:,1]),$(ma_rth2_ddm[:,1]),col=palsub[1],log='y',cex=0.5,pch=16,xlab='m*',ylab='Return time',ylim=c(min($rt_ddm),max($(ma_rth2_ddm))))
# points($(ma2_m2_ddm[:,1]),$(ma_rth2_ddm[:,1]),col=palsub[1],cex=0.5,pch=16)
# points($(ma1_m2_ddm[:,2]),$(ma_rth2_ddm[:,2]),col=palsub[2],cex=0.5,pch=16)
# points($(ma2_m2_ddm[:,2]),$(ma_rth2_ddm[:,2]),col=palsub[2],cex=0.5,pch=16)
# points($(ma1_m2_ddm[:,3]),$(ma_rth2_ddm[:,3]),col=palsub[3],cex=0.5,pch=16)
# points($(ma2_m2_ddm[:,3]),$(ma_rth2_ddm[:,3]),col=palsub[3],cex=0.5,pch=16)
# text(0.01,max($(ma_rth2_ddm)),expression(paste(h^2,'=0.2')))
# legend(x=0.3,y=130,legend=$pvec,col=palsub,pch=22,xpd=TRUE,pt.bg=palsub,cex=0.8, bty="n") #,title=expression(paste(Delta,theta))
# for (i in 1:length($mvec)) {
#   segments($(ma1_m2_ddm[:,1])[i],$(ma_rth2_ddm[:,1])[i],$(ma2_m2_ddm[:,1])[i],$(ma_rth2_ddm[:,1])[i],col=palsub[1])
#   segments($(ma1_m2_ddm[:,2])[i],$(ma_rth2_ddm[:,2])[i],$(ma2_m2_ddm[:,2])[i],$(ma_rth2_ddm[:,2])[i],col=palsub[2])
#   segments($(ma1_m2_ddm[:,3])[i],$(ma_rth2_ddm[:,3])[i],$(ma2_m2_ddm[:,3])[i],$(ma_rth2_ddm[:,3])[i],col=palsub[3])
# }
# plot($(ma1_m8_ddm[:,1]),$(ma_rth8_ddm[:,1]),col=palsub[1],log='y',cex=0.5,pch=16,xlab='m*',ylab='Return time',ylim=c(min($rt_ddm),max($(ma_rth8_ddm))))
# points($(ma2_m8_ddm[:,1]),$(ma_rth8_ddm[:,1]),col=palsub[1],cex=0.5,pch=16)
# points($(ma1_m8_ddm[:,2]),$(ma_rth8_ddm[:,2]),col=palsub[2],cex=0.5,pch=16)
# points($(ma2_m8_ddm[:,2]),$(ma_rth8_ddm[:,2]),col=palsub[2],cex=0.5,pch=16)
# points($(ma1_m8_ddm[:,3]),$(ma_rth8_ddm[:,3]),col=palsub[3],cex=0.5,pch=16)
# points($(ma2_m8_ddm[:,3]),$(ma_rth8_ddm[:,3]),col=palsub[3],cex=0.5,pch=16)
# text(0.01,max($(ma_rth8_ddm)),expression(paste(h^2,'=0.8')))
# for (i in 1:length($mvec)) {
#   segments($(ma1_m8_ddm[:,1])[i],$(ma_rth8_ddm[:,1])[i],$(ma2_m8_ddm[:,1])[i],$(ma_rth8_ddm[:,1])[i],col=palsub[1])
#   segments($(ma1_m8_ddm[:,2])[i],$(ma_rth8_ddm[:,2])[i],$(ma2_m8_ddm[:,2])[i],$(ma_rth8_ddm[:,2])[i],col=palsub[2])
#   segments($(ma1_m8_ddm[:,3])[i],$(ma_rth8_ddm[:,3])[i],$(ma2_m8_ddm[:,3])[i],$(ma_rth8_ddm[:,3])[i],col=palsub[3])
# }
# dev.off()
# """
# 


#COMBINED 4-PANEL plot





ma_rth2 = mean([rt[i,:,:,1] for i=1:reps]);
ma_rth8 = mean([rt[i,:,:,2] for i=1:reps]);




ma_rth2_ddm = mean([rt_ddm[i,:,:,1] for i=1:reps]);
ma_rth8_ddm = mean([rt_ddm[i,:,:,2] for i=1:reps]);

#mean stray rates
ma1_m2_ddm = mean([m1mean[i,:,:,1] for i=1:reps]);
ma2_m2_ddm = mean([m2mean[i,:,:,1] for i=1:reps]);
ma1_m8_ddm = mean([m1mean[i,:,:,2] for i=1:reps]);
ma2_m8_ddm = mean([m2mean[i,:,:,2] for i=1:reps]);


namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/figs2/fig_relax_comb.pdf");
R"""
library(RColorBrewer)
pdf($namespace,height=8,width=10)
pal = brewer.pal(9,'Greys')
palsub = pal[c(4,6,8)];
par(mfrow=c(2,2),mai = c(0.8, 0.9, 0.2, 0.1))
plot($mvec,$(ma_rth2[:,1]),col=palsub[1],type='l',log='y',cex=0.5,lwd=2,xlab='m',ylab='Recovery time',ylim=c(min($rt),max($(ma_rth2))))
lines($mvec,$(ma_rth2[:,2]),col=palsub[2],cex=0.5,lwd=2)
lines($mvec,$(ma_rth2[:,3]),col=palsub[3],cex=0.5,lwd=2)
text(0.25,max($(ma_rth2))+23,expression(paste(h^2,'=0.2')), xpd=TRUE)
text(-0.11, ## x position
     41, ## position of the low axis
     srt=90, ## angle
     labels='Constant m', ##labels
     xpd=TRUE, ## allows plotting outside the region 
     pos=2)


plot($mvec,$(ma_rth8[:,1]),col=palsub[1],type='l',log='y',cex=0.5,lwd=2,xlab='m',ylab='Recovery time',ylim=c(min($rt),max($(ma_rth8))))
lines($mvec,$(ma_rth8[:,2]),col=palsub[2],cex=0.5,lwd=2)
lines($mvec,$(ma_rth8[:,3]),col=palsub[3],cex=0.5,lwd=2)
text(0.25,max($(ma_rth8))+1600,expression(paste(h^2,'=0.8')), xpd=TRUE)
text(par('usr')[1]-0.1,max($(ma_rth2))+23,'(a)', xpd=TRUE)
legend(x=0.43,y=3500,legend=$pvec,col=palsub,pch=22,xpd=TRUE,pt.bg=palsub,cex=1, bty="n") #,title=expression(paste(Delta,theta))


plot($(ma1_m2_ddm[:,1]),$(ma_rth2_ddm[:,1]),col=palsub[1],log='y',cex=0.5,pch=16,xlab='m*',ylab='Recovery time',ylim=c(min($rt_ddm),max($(ma_rth2_ddm))))
points($(ma2_m2_ddm[:,1]),$(ma_rth2_ddm[:,1]),col=palsub[1],cex=0.5,pch=16)
points($(ma1_m2_ddm[:,2]),$(ma_rth2_ddm[:,2]),col=palsub[2],cex=0.5,pch=16)
points($(ma2_m2_ddm[:,2]),$(ma_rth2_ddm[:,2]),col=palsub[2],cex=0.5,pch=16)
points($(ma1_m2_ddm[:,3]),$(ma_rth2_ddm[:,3]),col=palsub[3],cex=0.5,pch=16)
points($(ma2_m2_ddm[:,3]),$(ma_rth2_ddm[:,3]),col=palsub[3],cex=0.5,pch=16)
#text(0.01,max($(ma_rth2_ddm)),expression(paste(h^2,'=0.2')))
for (i in 1:length($mvec)) {
  segments($(ma1_m2_ddm[:,1])[i],$(ma_rth2_ddm[:,1])[i],$(ma2_m2_ddm[:,1])[i],$(ma_rth2_ddm[:,1])[i],col=palsub[1])
  segments($(ma1_m2_ddm[:,2])[i],$(ma_rth2_ddm[:,2])[i],$(ma2_m2_ddm[:,2])[i],$(ma_rth2_ddm[:,2])[i],col=palsub[2])
  segments($(ma1_m2_ddm[:,3])[i],$(ma_rth2_ddm[:,3])[i],$(ma2_m2_ddm[:,3])[i],$(ma_rth2_ddm[:,3])[i],col=palsub[3])
}
text(-0.072, ## x position
     55, ## position of the low axis
     srt=90, ## angle
     labels='Density dependent m', ##labels
     xpd=TRUE, ## allows plotting outside the region 
     pos=2)
     
plot($(ma1_m8_ddm[:,1]),$(ma_rth8_ddm[:,1]),col=palsub[1],log='y',cex=0.5,pch=16,xlab='m*',ylab='Recovery time',ylim=c(min($rt_ddm),max($(ma_rth8_ddm))))
points($(ma2_m8_ddm[:,1]),$(ma_rth8_ddm[:,1]),col=palsub[1],cex=0.5,pch=16)
points($(ma1_m8_ddm[:,2]),$(ma_rth8_ddm[:,2]),col=palsub[2],cex=0.5,pch=16)
points($(ma2_m8_ddm[:,2]),$(ma_rth8_ddm[:,2]),col=palsub[2],cex=0.5,pch=16)
points($(ma1_m8_ddm[:,3]),$(ma_rth8_ddm[:,3]),col=palsub[3],cex=0.5,pch=16)
points($(ma2_m8_ddm[:,3]),$(ma_rth8_ddm[:,3]),col=palsub[3],cex=0.5,pch=16)
#text(0.01,max($(ma_rth8_ddm)),expression(paste(h^2,'=0.8')))
for (i in 1:length($mvec)) {
  segments($(ma1_m8_ddm[:,1])[i],$(ma_rth8_ddm[:,1])[i],$(ma2_m8_ddm[:,1])[i],$(ma_rth8_ddm[:,1])[i],col=palsub[1])
  segments($(ma1_m8_ddm[:,2])[i],$(ma_rth8_ddm[:,2])[i],$(ma2_m8_ddm[:,2])[i],$(ma_rth8_ddm[:,2])[i],col=palsub[2])
  segments($(ma1_m8_ddm[:,3])[i],$(ma_rth8_ddm[:,3])[i],$(ma2_m8_ddm[:,3])[i],$(ma_rth8_ddm[:,3])[i],col=palsub[3])
}
dev.off()
"""
