

@everywhere using Distributions
@everywhere using RCall
@everywhere using HDF5
@everywhere using JLD

@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolve.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/qualsfunc.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/bifdet.jl")


#Analysis over m & theta divergence
mvec=collect(0.0:0.001:0.3);
thetadiffvec = collect(0.0:0.1:10.0);

thetadiffarray=SharedArray(Float64,length(thetadiffvec),length(mvec));

n1mean=SharedArray(Float64,length(thetadiffvec),length(mvec));
n2mean=SharedArray(Float64,length(thetadiffvec),length(mvec));
x1mean=SharedArray(Float64,length(thetadiffvec),length(mvec));
x2mean=SharedArray(Float64,length(thetadiffvec),length(mvec));
pe=SharedArray(Float64,length(thetadiffvec),length(mvec));

tmax=10000;
z=0.5;
rmax=2.0;
beta=0.001;
theta1=5.0;

tau=1.0;
sigmaE=0;
sigmaG=1;
h=0.2
perror=0.01;

@sync @parallel for k=1:length(thetadiffvec)
    
    thetadiff=thetadiffvec[k];
    
    for j=1:length(mvec)
      m = mvec[j];
      
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
      burnin=0.8
      n1trim = n1[Int64(floor(tmax*burnin)):tmax];
      n2trim = n2[Int64(floor(tmax*burnin)):tmax];
      x1trim = x1[Int64(floor(tmax*burnin)):tmax];
      x2trim = x2[Int64(floor(tmax*burnin)):tmax];
      
      n1mean[k,j] = mean(n1trim);
      n2mean[k,j] = mean(n2trim);
      x1mean[k,j] = theta1-mean(x1trim);
      x2mean[k,j] = (theta1+thetadiff)-mean(x2trim);
      thetadiffarray[k,j] = thetadiff;
      
      pe[k,j] = mean([(std(n1trim)/mean(n1trim)),(std(n2trim)/mean(n2trim))])*
      (1/(std(n1trim+n2trim)/mean(n1trim+n2trim)))
    end
end


save(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data2/data_habitathetero.jld"),"n1mean",n1mean,"n2mean",n2mean,"x1mean",x1mean,"x2mean",x2mean,"pe",pe);

d = load(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data2/data_habitathetero_ext.jld"));
#This loads the dictionary
n1mean = d["n1mean"];
n2mean = d["n2mean"];
x1mean = d["x1mean"];
x2mean = d["x2mean"];
pe = d["pe"];

mediandiff = zeros(Float64,length(thetadiffvec));
for k=1:length(thetadiffvec)
  mediandiff[k] = median(abs(n1mean[k,:] - n2mean[k,:]));
end

R"""
plot($thetadiffvec,$mediandiff,log='y',type='l')
"""

diffarray = abs(n1mean-n2mean);
R"""
plot($thetadiffarray,$diffarray,pch='.')
"""

aggarray = n1mean+n2mean;
R"""
plot($thetadiffarray,$aggarray,pch='.')
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
mvec = collect(0.0001:0.001:0.5);
thetadiffvec = collect(0.0:0.1:10.0);
hvec = collect(0.1:0.1:0.9);

thetadiffarray=SharedArray(Float64,length(thetadiffvec),length(hvec),length(mvec));

n1mean=SharedArray(Float64,length(thetadiffvec),length(hvec),length(mvec));
n2mean=SharedArray(Float64,length(thetadiffvec),length(hvec),length(mvec));
x1mean=SharedArray(Float64,length(thetadiffvec),length(hvec),length(mvec));
x2mean=SharedArray(Float64,length(thetadiffvec),length(hvec),length(mvec));
pe=SharedArray(Float64,length(thetadiffvec),length(hvec),length(mvec));
n1mean_ddm=SharedArray(Float64,length(thetadiffvec),length(hvec),length(mvec));
n2mean_ddm=SharedArray(Float64,length(thetadiffvec),length(hvec),length(mvec));
x1mean_ddm=SharedArray(Float64,length(thetadiffvec),length(hvec),length(mvec));
x2mean_ddm=SharedArray(Float64,length(thetadiffvec),length(hvec),length(mvec));
pe_ddm=SharedArray(Float64,length(thetadiffvec),length(hvec),length(mvec));

rt=SharedArray(Float64,length(thetadiffvec),length(hvec),length(mvec));
rt_ddm=SharedArray(Float64,length(thetadiffvec),length(hvec),length(mvec));
m1mean=SharedArray(Float64,length(thetadiffvec),length(hvec),length(mvec));
m2mean=SharedArray(Float64,length(thetadiffvec),length(hvec),length(mvec));

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


@sync @parallel for k=1:length(thetadiffvec)
    
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
        
        n1mean[k,i,j] = mean(n1trim);
        n2mean[k,i,j] = mean(n2trim);
        x1mean[k,i,j] = theta1-mean(x1trim);
        x2mean[k,i,j] = (theta1+thetadiff)-mean(x2trim);
        
        pe[k,i,j] = mean([(std(n1trim)/mean(n1trim)),(std(n2trim)/mean(n2trim))])*
        (1/(std(n1trim+n2trim)/mean(n1trim+n2trim)));
        
        t_ss, relaxtime = timeSS(n1,n2,t_ext);
        
        rt[k,i,j] = relaxtime;
        thetadiffarray[k,i,j] = thetadiff;
        
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
        
        n1mean_ddm[k,i,j] = mean(n1trim_ddm);
        n2mean_ddm[k,i,j] = mean(n2trim_ddm);
        x1mean_ddm[k,i,j] = theta1-mean(x1trim_ddm);
        x2mean_ddm[k,i,j] = (theta1+thetadiff)-mean(x2trim_ddm);
        
        pe_ddm[k,i,j] = mean([(std(n1trim)/mean(n1trim)),(std(n2trim)/mean(n2trim))])*
        (1/(std(n1trim+n2trim)/mean(n1trim+n2trim)))
      
        t_ss_ddm, relaxtime_ddm = timeSS(n1_ddm,n2_ddm,t_ext);
        
        rt_ddm[k,i,j] = relaxtime_ddm;
        
        #Steady state stray rate?
        m1trim = m1_ddm[Int64(floor(tmax*burnin)):tmax-1];
        m2trim = m2_ddm[Int64(floor(tmax*burnin)):tmax-1];
        
        m1mean[k,i,j] = mean(m1trim);
        m2mean[k,i,j] = mean(m2trim);
        
      end
    end
end


save(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data2/data_habitathetero_ext.jld"),"n1mean",n1mean,"n2mean",n2mean,"x1mean",x1mean,"x2mean",x2mean,"pe",pe,"rt",rt,"n1mean_ddm",n1mean_ddm,"n2mean_ddm",n2mean_ddm,"x1mean_ddm",x1mean_ddm,"x2mean_ddm",x2mean_ddm,"rt_ddm",rt_ddm,"m1mean",m1mean,"m2mean",m2mean,"pe_ddm",pe_ddm,"thetadiffarray",thetadiffarray);

d = load(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data2/data_habitathetero_ext.jld"));
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

namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/figs2/fig_thetadiffN.pdf");
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



R"""
plot($thetadiffarray,$rt,pch='.',log='y')
"""

@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/movingaverage.jl")
window=20;
ma_mvec = movingaverage(mvec,window);
ma_rt3 = movingaverage(mean([rt[31,i,:] for i=1:5]),window);
ma_rt5 = movingaverage(mean([rt[51,i,:] for i=1:5]),window);
ma_rt8 = movingaverage(mean([rt[81,i,:] for i=1:5]),window);

namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/figs2/fig_relaxtheta.pdf");
R"""
library(RColorBrewer)
pal = brewer.pal(3,'Set1')
pdf($namespace,height=8,width=5)
par(mfrow=c(2,1),mai = c(0.8, 0.8, 0.4, 0.1))
plot($ma_mvec,$(ma_rt5),type='l',log='y',col=pal[2],ylim=c(20,100),xlab='m',ylab='Recovery time',lwd=2)
lines($ma_mvec,$(ma_rt8),col=pal[3],lwd=2)
points($(mean([m1mean[51,i,:] for i=1:5])),$(mean([rt_ddm[51,i,:] for i=1:5])),pch=16,cex=.4,col=pal[2])
points($(mean([m2mean[51,i,:] for i=1:5])),$(mean([rt_ddm[51,i,:] for i=1:5])),pch=16,cex=.4,col=pal[2])
for (i in 1:length($(mean([m1mean[51,i,:] for i=1:5])))) {
segments($(mean([m1mean[51,i,:] for i=1:5]))[i],$(mean([rt_ddm[51,i,:] for i=1:5]))[i],$(mean([m2mean[51,i,:] for i=1:5]))[i],$(mean([rt_ddm[51,i,:] for i=1:5]))[i],col=pal[2])
}
legend(x=0.45,y=35,legend=c(5,8),col=pal[2:3],pch=22,xpd=TRUE,pt.bg=pal[2:3],cex=0.8, bty="n",title=expression(paste(Delta,theta)))
plot($ma_mvec,$(ma_rt5),type='l',log='y',col=pal[2],ylim=c(20,100),xlab='m',ylab='Recovery time',lwd=2)
lines($ma_mvec,$(ma_rt8),col=pal[3],lwd=2)
points($(mean([m1mean[81,i,:] for i=1:5])),$(mean([rt_ddm[81,i,:] for i=1:5])),pch=16,cex=.4,col=pal[3])
points($(mean([m2mean[81,i,:] for i=1:5])),$(mean([rt_ddm[81,i,:] for i=1:5])),pch=16,cex=.4,col=pal[3])
for (i in 1:length($(mean([m1mean[51,i,:] for i=1:5])))) {
segments($(mean([m1mean[81,i,:] for i=1:5]))[i],$(mean([rt_ddm[81,i,:] for i=1:5]))[i],$(mean([m2mean[81,i,:] for i=1:5]))[i],$(mean([rt_ddm[81,i,:] for i=1:5]))[i],col=pal[3])
}
legend(x=0.45,y=35,legend=c(5,8),col=pal[2:3],pch=22,xpd=TRUE,pt.bg=pal[2:3],cex=0.8, bty="n",title=expression(paste(Delta,theta)))
dev.off()
"""


R"""
plot($rt,$pe,pch='.',log='x',xlim=c(10,500))
"""
