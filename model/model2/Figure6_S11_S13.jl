#Figure 6
#Figure S11
#Figure S13


@everywhere using Distributions, RCall, HDF5, JLD

@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolve_mtheta.jl")

@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/qualsfunc.jl")

@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolveExtinct.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolveExtinct_ddm.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/timeSS.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/movingaverage.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/bifdet.jl")




#Recovery time sims
#Analysis over m
tmax=10000;

mvec2 = collect(0.01:0.001:0.25);
pvec = ["small","large","both"];
reps = 100;

rt = SharedArray{Float64}(reps,length(mvec2),length(pvec));
rt_ddm = SharedArray{Float64}(reps,length(mvec2),length(pvec));
m1mean = SharedArray{Float64}(reps,length(mvec2),length(pvec));
m2mean = SharedArray{Float64}(reps,length(mvec2),length(pvec));
pe = SharedArray{Float64}(reps,length(mvec2),length(pvec));


z=0.5;
rmax=2.0;
beta=0.001;
theta1=2.0;
thetascale=5.0;
tau=1.0;
C=300;
h=0.2;
sigmaE=0;
sigmaG=1;
perror=0.01;
refuge=0.01;
burnin = 0.8;
t_ext = Int64(round(tmax/2));

@sync @parallel for r=1:reps
  for j=1:length(pvec)
    extpop = pvec[j]
    for i=1:length(mvec2)
      
      m = mvec2[i];
      a0 = 1-m;
      
      thetadiff = (1-2*m)/(thetascale*m);
      
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
      
      rt[r,i,j] = relaxtime;
      
      burnin=0.80
      n1trim = n1[Int64(floor(tmax*burnin)):tmax];
      n2trim = n2[Int64(floor(tmax*burnin)):tmax];
      x1trim = x1[Int64(floor(tmax*burnin)):tmax];
      x2trim = x2[Int64(floor(tmax*burnin)):tmax]
      
      pe[r,i,j] = mean([(std(n1trim)/mean(n1trim)),(std(n2trim)/mean(n2trim))])*
      (1/(std(n1trim+n2trim)/mean(n1trim+n2trim)))
      
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
      
      rt_ddm[r,i,j] = relaxtime_ddm;
      
      #Steady state stray rate?
      m1trim = m1_ddm[Int64(floor(tmax*burnin)):tmax-1];
      m2trim = m2_ddm[Int64(floor(tmax*burnin)):tmax-1];
      
      m1mean[r,i,j] = mean(m1trim);
      m2mean[r,i,j] = mean(m2trim);
      
    end
  end
end

save(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data3/data_relax_thetam_4.jld"),"rt",rt,"rt_ddm",rt_ddm,"m1mean",m1mean,"m2mean",m2mean,"pe",pe);



d = load(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data3/data_relax_thetam_3.jld"));
#This loads the dictionary
rt = d["rt"];
rt_ddm = d["rt_ddm"];
m1mean = d["m1mean"];
m2mean = d["m2mean"];


ma_rth = mean([rt[i,:,:] for i=1:reps]);
ma_rth_ddm = mean([rt_ddm[i,:,:] for i=1:reps]);

ma_peh = mean([pe[i,:,:] for i=1:reps]);

#mean stray rates
ma1_m_ddm = mean([m1mean[i,:,:] for i=1:reps]);
ma2_m_ddm = mean([m2mean[i,:,:] for i=1:reps]);

namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/Manuscript/FinalDraft3/fig_mtheta_rt.pdf");
R"""
library(RColorBrewer)
#library(shape)
pal = brewer.pal(9,'Greys')
pdf($namespace,height=4,width=5)
palsub = pal[c(4,6,8)];
plot($mvec2,$(ma_rth[:,1]),col=palsub[1],type='l',log='y',cex=0.5,lwd=2,xlab='m',ylab='Recovery time',xlim=c(0,0.25),ylim=c(min($rt),max($(ma_rth))),las=1)
lines($mvec2,$(ma_rth[:,2]),col=palsub[2],cex=0.5,lwd=2)
lines($mvec2,$(ma_rth[:,3]),col=palsub[3],cex=0.5,lwd=2)
#Arrows(0.119,5000,0.119,3000,code=2,arr.type ='circle')
arrows($(mvec2[indmax(ma_rth[:,1])]),$(maximum(ma_rth)+500),$(mvec2[indmax(ma_rth[:,1])]),$(maximum(ma_rth)+100),length=0.05,angle=40,lwd=3)
types = c('subordinate extinct','dominant extinct','near-collapse')
legend(x=0.14,y=5000,legend=types,col=palsub,pch=22,xpd=TRUE,pt.bg=palsub,cex=0.9, bty="n") #,title=expression(paste(Delta,theta))
text($(mvec2[indmax(ma_rth[:,1])]),$(maximum(ma_rth)+1000),'DCB')
text(0.027,90,'*',cex=2)
dev.off()
"""



namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/Manuscript/FinalDraft3/fig_relax_mtheta.pdf");
R"""
library(RColorBrewer)
pal = brewer.pal(9,'Greys')
palsub = pal[c(4,6,8)];
pdf($namespace,height=4,width=10)
par(mfrow=c(1,2),mai = c(0.8, 0.9, 0.3, 0.1))
plot($mvec2,$(ma_rth[:,1]),col=palsub[1],type='l',log='y',cex=0.5,lwd=2,xlab='m',ylab='Recovery time',xlim=c(0,0.25),ylim=c(min($rt),max($(ma_rth))))
lines($mvec2,$(ma_rth[:,2]),col=palsub[2],cex=0.5,lwd=2)
lines($mvec2,$(ma_rth[:,3]),col=palsub[3],cex=0.5,lwd=2)
arrows($(mvec2[indmax(ma_rth[:,1])]),$(maximum(ma_rth)+500),$(mvec2[indmax(ma_rth[:,1])]),$(maximum(ma_rth)+100),length=0.05,angle=40,lwd=3)
types = c('subordinate extinct','dominant extinct','near-collapse')
legend(x=0.14,y=5000,legend=types,col=palsub,pch=22,xpd=TRUE,pt.bg=palsub,cex=0.9, bty="n") #,title=expression(paste(Delta,theta))
text($(mvec2[indmax(ma_rth[:,1])]),$(maximum(ma_rth)+1000),'DCB')
text(0.027,90,'*',cex=2)

plot($(ma1_m_ddm[:,1]),$(ma_rth_ddm[:,1]),col=palsub[1],log='y',cex=0.5,pch=16,xlab='m*',ylab='Recovery time',xlim=c(0,max($(ma1_m_ddm[!isnan(ma1_m_ddm)]))),ylim=c(min($rt_ddm),200))
points($(ma2_m_ddm[:,1]),$(ma_rth_ddm[:,1]),col=palsub[1],cex=0.5,pch=16)
points($(ma1_m_ddm[:,2]),$(ma_rth_ddm[:,2]),col=palsub[2],cex=0.5,pch=16)
points($(ma2_m_ddm[:,2]),$(ma_rth_ddm[:,2]),col=palsub[2],cex=0.5,pch=16)
points($(ma1_m_ddm[:,3]),$(ma_rth_ddm[:,3]),col=palsub[3],cex=0.5,pch=16)
points($(ma2_m_ddm[:,3]),$(ma_rth_ddm[:,3]),col=palsub[3],cex=0.5,pch=16)
for (i in 1:length($mvec2)) {
 segments($(ma1_m_ddm[:,1])[i],$(ma_rth_ddm[:,1])[i],$(ma2_m_ddm[:,1])[i],$(ma_rth_ddm[:,1])[i],col=palsub[1])
 segments($(ma1_m_ddm[:,2])[i],$(ma_rth_ddm[:,2])[i],$(ma2_m_ddm[:,2])[i],$(ma_rth_ddm[:,2])[i],col=palsub[2])
 segments($(ma1_m_ddm[:,3])[i],$(ma_rth_ddm[:,3])[i],$(ma2_m_ddm[:,3])[i],$(ma_rth_ddm[:,3])[i],col=palsub[3])
}
dev.off()
"""


namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/Manuscript/FinalDraft3/fig_relax_mtheta.pdf");
R"""
library(RColorBrewer)
pal = brewer.pal(9,'Greys')
palsub = pal[c(4,6,8)];
pdf($namespace,height=4,width=10)
par(mfrow=c(1,2),mai = c(0.8, 0.9, 0.3, 0.1))
plot($mvec2,$(ma_rth[:,1]),col=palsub[1],type='l',log='y',cex=0.5,lwd=2,xlab='m',ylab='Recovery time',xlim=c(0,0.25),ylim=c(min($rt),max($(ma_rth))))
lines($mvec2,$(ma_rth[:,2]),col=palsub[2],cex=0.5,lwd=2)
lines($mvec2,$(ma_rth[:,3]),col=palsub[3],cex=0.5,lwd=2)
arrows($(mvec2[indmax(ma_rth[:,1])]),$(maximum(ma_rth)+500),$(mvec2[indmax(ma_rth[:,1])]),$(maximum(ma_rth)+100),length=0.05,angle=40,lwd=3)
types = c('subordinate extinct','dominant extinct','near-collapse')
legend(x=0.14,y=5000,legend=types,col=palsub,pch=22,xpd=TRUE,pt.bg=palsub,cex=0.9, bty="n") #,title=expression(paste(Delta,theta))
text($(mvec2[indmax(ma_rth[:,1])]),$(maximum(ma_rth)+1000),'DCB')
text(0.027,90,'*',cex=2)

plot($(mvec2),$(ma_rth_ddm[:,1]),col=palsub[1],log='y',cex=0.5,pch=16,xlab='m*',ylab='Recovery time',xlim=c(0,max($(ma1_m_ddm[!isnan(ma1_m_ddm)]))),ylim=c(min($rt_ddm),max($ma_rth_ddm)))
points($(mvec2),$(ma_rth_ddm[:,2]),col=palsub[2],cex=0.5,pch=16)
points($(mvec2),$(ma_rth_ddm[:,3]),col=palsub[3],cex=0.5,pch=16)
dev.off()
"""




#Analysis over m & theta divergence
mvec=collect(0.01:0.0001:0.1);
hvec = collect(0.01:0.001:0.3);
td = SharedArray{Float64}(length(mvec));

n1mean=SharedArray{Float64}(length(hvec),length(mvec));
n2mean=SharedArray{Float64}(length(hvec),length(mvec));
x1mean=SharedArray{Float64}(length(hvec),length(mvec));
x2mean=SharedArray{Float64}(length(hvec),length(mvec));
pe=SharedArray{Float64}(length(hvec),length(mvec));

tmax=10000;
z=0.5;
rmax=2.0;
beta=0.001;
theta1=5.0;
thetascale=10.0;
tau=1.0;
sigmaE=0;
sigmaG=1;


perror=0.01;

@sync @parallel for k=1:length(hvec)
  h = hvec[k];
  
  for j=1:length(mvec)
    m = mvec[j];
    thetadiff = (1-2*m)/(thetascale*m);
    td[j] = copy(thetadiff);
    
    n1, n2, x1, x2, w1, w2 = 
    KevinEvolve_mtheta(
      tmax,
      z,
      rmax,
      beta,
      theta1,
      thetascale,
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

# namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/Manuscript/FinalDraft3/fig_MDPE_hm_mtheta_rt.pdf");
# R"""
# library(RColorBrewer)
# pal = rev(brewer.pal(9,"Blues"))
# pdf($namespace,height=4,width=4)
# par(mfrow=c(1,1), mai = c(1.2, 0.8, 0.1, 0.1))
# image(x=$mvec,y=$hvec,z=t($(pe[:,:])),zlim=c(1,2),col=pal,xlab='m',ylab='h')
# lines($bifvalue)
# dev.off()
# """

namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/Manuscript/FinalDraft3/fig_mtheta_ss.pdf");
R"""
library(RColorBrewer)
pal = rev(brewer.pal(9,"Blues"))
pdf($namespace,height=4,width=10)
par(mfrow=c(1,2),mai = c(0.8, 0.8, 0.2, 0.2))
cols=brewer.pal(5,'Set1')
plot($mvec,$(n1mean[2,:]),pch='.',col=cols[4],xlab='m',ylab='Steady state')
points($mvec,$(n2mean[2,:]),pch='.',col=cols[4])
text(par('usr')[1]-0.016,max($(n1mean))+10,'(a)', xpd=TRUE)
image(x=$mvec,y=$hvec,z=t($(abs(n1mean-n2mean))),zlim=c(0,400),col=pal,xlab='m',ylab=expression(paste(h^2)))
text(par('usr')[1]-0.016,max($(hvec)),'(b)', xpd=TRUE)
dev.off()
"""
