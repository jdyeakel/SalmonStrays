#Figure 2
#Figure 4

@everywhere using Distributions
@everywhere using RCall
@everywhere using HDF5
@everywhere using JLD

@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolve.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolveExtinct.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/timeSS.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/qualsfunc.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/bifdet.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/movingaverage.jl")


#Analysis over m & theta divergence
mvec=collect(0.0:0.0005:0.3);
hvec = collect(0.0:0.01:0.5);


n1mean=SharedArray{Float64}(length(hvec),length(mvec));
n2mean=SharedArray{Float64}(length(hvec),length(mvec));
x1mean=SharedArray{Float64}(length(hvec),length(mvec));
x2mean=SharedArray{Float64}(length(hvec),length(mvec));
pe=SharedArray{Float64}(length(hvec),length(mvec));

n1meanE=SharedArray{Float64}(length(hvec),length(mvec));
n2meanE=SharedArray{Float64}(length(hvec),length(mvec));
x1meanE=SharedArray{Float64}(length(hvec),length(mvec));
x2meanE=SharedArray{Float64}(length(hvec),length(mvec));
peE=SharedArray{Float64}(length(hvec),length(mvec));
rtE=SharedArray{Float64}(length(hvec),length(mvec));

tmax=10000;
z=2;
rmax=2.0;
beta=0.001;
theta1=5.0;
thetadiff=2.4;
tau=1.0;
sigmaE=0;
sigmaG=1;

extpop="both";
refuge=0.01;
t_ext = Int64(round(tmax/2));

perror=0.001;

@sync @parallel for k=1:length(hvec)
    
    h = hvec[k];
    
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
      x2trim = x2[Int64(floor(tmax*burnin)):tmax]
      
      n1mean[k,j] = mean(n1trim);
      n2mean[k,j] = mean(n2trim);
      x1mean[k,j] = theta1-mean(x1trim);
      x2mean[k,j] = (theta1+thetadiff)-mean(x2trim);
      
      pe[k,j] = mean([(std(n1trim)/mean(n1trim)),(std(n2trim)/mean(n2trim))])*
      (1/(std(n1trim+n2trim)/mean(n1trim+n2trim)))
      
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
      
      rtE[k,j] = relaxtime;
      
      burnin=0.9
      n1trim = n1[Int64(floor(tmax*burnin)):tmax];
      n2trim = n2[Int64(floor(tmax*burnin)):tmax];
      x1trim = x1[Int64(floor(tmax*burnin)):tmax];
      x2trim = x2[Int64(floor(tmax*burnin)):tmax]
      
      n1meanE[k,j] = mean(n1trim);
      n2meanE[k,j] = mean(n2trim);
      x1meanE[k,j] = theta1-mean(x1trim);
      x2meanE[k,j] = (theta1+thetadiff)-mean(x2trim);
      
      peE[k,j] = mean([(std(n1trim)/mean(n1trim)),(std(n2trim)/mean(n2trim))])*
      (1/(std(n1trim+n2trim)/mean(n1trim+n2trim)))
      
    end
end
save(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data3/data_sig_h_m.jld"),"n1mean",n1mean,"n2mean",n2mean,"x1mean",x1mean,"x2mean",x2mean,"pe",pe,"n1meanE",n1meanE,"n2meanE",n2meanE,"x1meanE",x1meanE,"x2meanE",x2meanE,"peE",peE,"rtE",rtE);

d = load(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data3/data_sig_h_m.jld"));
#This loads the dictionary
n1mean = d["n1mean"];
n2mean = d["n2mean"];
x1mean = d["x1mean"];
x2mean = d["x2mean"];
pe = d["pe"];
peE = d["peE"];
rtE = d["rtE"];

rt_ext = rtE;
pe_ext = peE;

bifvalue = bifdet(
n1mean,
n2mean,
mvec,
hvec
);


namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/FinalDraft3/fig_pevsrt.pdf");
R"""
library(RColorBrewer)
pal = rev(brewer.pal(9,"Blues"))
pdf($namespace,height=5,width=5)
par(mfrow=c(1,1))
plot($(rt_ext),$(pe_ext),log='xy',ylim=c(1,3),xlim=c(18,300),pch='.',xlab='Recovery time',ylab='Portfolio effect',las=1,col='#00000025')
text(7,3.18,'(d)', xpd=TRUE)
dev.off()
"""



namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/FinalDraft3/fig_MDPE_hm.pdf");
R"""
library(RColorBrewer)
pal = rev(brewer.pal(9,"Blues"))
pdf($namespace,height=3,width=10)
par(mfrow=c(1,4))
image(x=$mvec,y=$hvec,z=t($(n1mean[:,:]))+t($(n2mean[:,:])),zlim=c(0,1000),col=pal,xlab='m',ylab=expression(paste(h^2)),las=1)
#points($(bifvalue),type='l',cex=1)
text(par('usr')[1]-0.12,1.05,'(a)', xpd=TRUE)
image(x=$mvec,y=$hvec,z=sqrt((t($(n1mean[:,:]))-t($(n2mean[:,:])))^2),zlim=c(0,500),col=pal,xlab='m',ylab=expression(paste(h^2)),las=1)
#points($(bifvalue),type='l',cex=1)
text(par('usr')[1]-0.12,1.05,'(b)', xpd=TRUE)
image(x=$mvec,y=$hvec,z=t($(pe[:,:])),zlim=c(1,2),col=pal,xlab='m',ylab=expression(paste(h^2)),las=1)
#points($(bifvalue),type='l',cex=1)
text(par('usr')[1]-0.12,1.05,'(c)', xpd=TRUE)
plot($(rt_ext[1:50,:]),$(pe_ext[1:50,:]),log='xy',ylim=c(1,3),xlim=c(18,300),pch='.',xlab='Recovery time',ylab='PE',las=1,col='#00000025')
text(7,3.18,'(d)', xpd=TRUE)
dev.off()
"""



namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/FinalDraft3/fig_traitdiff.pdf");
traitdiff = abs(-1*(x1mean[:,:]-theta1) - (-1*(x2mean[:,:]-(theta1+thetadiff))));
R"""
pdf($namespace,height=4,width=5)
plot($mvec,$(traitdiff[26,:]),type='l',xlim=c(0,0.33),ylab=expression(paste('Phenotypic diversity, ',Delta,mu,'*')),xlab='m',las=1)
lines($mvec,$(traitdiff[51,:]),type='l')
lines($mvec,$(traitdiff[76,:]),type='l')
lines($mvec,$(traitdiff[101,:]),type='l')
text(rep(0.32,4),$(traitdiff[[26 51 76 101],301]),c('0.25','0.50','0.75','1.00'),cex=0.8)
text(0.315,1.0,expression(paste(h^2)),cex=0.8)
dev.off()
"""


#Which h2 = 0.2
indh = find(x->x==0.2,hvec)[1];
pe_ma = movingaverage([mvec peE[indh,:]],10);
rt_ma = movingaverage([mvec rtE[indh,:]],5);

R"""
par(mfrow=c(1,3))
plot($mvec,$(peE[indh,:]),type='l')
plot($mvec,$(rtE[indh,:]),log='y',type='l')
plot($(peE[indh,:]),$(rtE[indh,:]),log='y')
"""

namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/FinalDraft3/fig_mpert.pdf");
R"""
pdf($namespace,height=4,width=10)
par(mfrow=c(1,2),mai = c(0.9, 0.9, 0.1, 0.1))
plot($pe_ma,type='l',xlab='Straying (m)',ylab='Portfolio effect',xlim=c(0,0.3))
plot($rt_ma,log='y',xlab='Straying (m)',ylab='Recovery time',type='l',xlim=c(0,0.3))
dev.off()
"""

