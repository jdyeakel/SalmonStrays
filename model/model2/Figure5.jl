#Figure 5


@everywhere using Distributions, RCall, JLD, HDF5

@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/movingaverage.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolve.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/KevinEvolve_ddm.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/src/qualsfunc.jl")



# Figure 5

#Analysis over m & theta divergence
mvec=collect(0.0:0.001:0.3);
hvec = collect(0.0:0.01:1.0);


n1mean=SharedArray(Float64,length(hvec),length(mvec));
n2mean=SharedArray(Float64,length(hvec),length(mvec));
x1mean=SharedArray(Float64,length(hvec),length(mvec));
x2mean=SharedArray(Float64,length(hvec),length(mvec));
pe=SharedArray(Float64,length(hvec),length(mvec));

tmax=10000;
z=0.5;
rmax=2.0;
beta=0.001;
theta1=5.0;
thetadiff=5.0;
tau=1.0;
sigmaE=0;
sigmaG=1;

perror=0.01;

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
    end
end
save(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data3/data_sig_h_m.jld"),"n1mean",n1mean,"n2mean",n2mean,"x1mean",x1mean,"x2mean",x2mean,"pe",pe);


#Analysis over m & theta divergence
mvec=collect(0.0:0.001:0.3);
hvec = collect(0.0:0.01:1.0);


n1mean=SharedArray(Float64,length(hvec),length(mvec));
n2mean=SharedArray(Float64,length(hvec),length(mvec));
x1mean=SharedArray(Float64,length(hvec),length(mvec));
x2mean=SharedArray(Float64,length(hvec),length(mvec));
pe=SharedArray(Float64,length(hvec),length(mvec));

tmax=10000;
z=0.5;
rmax=2.0;
beta=0.001;
theta1=5.0;
thetadiff=3.0;
tau=1.0;
sigmaE=0;
sigmaG=1;

perror=0.01;

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
  end
end
save(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data3/data_sig_h_m_theta3.jld"),"n1mean",n1mean,"n2mean",n2mean,"x1mean",x1mean,"x2mean",x2mean,"pe",pe);



#Analysis over m & theta divergence
mvec=collect(0.0:0.001:0.3);
hvec = collect(0.0:0.01:1.0);


n1mean=SharedArray(Float64,length(hvec),length(mvec));
n2mean=SharedArray(Float64,length(hvec),length(mvec));
x1mean=SharedArray(Float64,length(hvec),length(mvec));
x2mean=SharedArray(Float64,length(hvec),length(mvec));
pe=SharedArray(Float64,length(hvec),length(mvec));

tmax=100000;
z=0.5;
rmax=2.0;
beta=0.001;
theta1=5.0;
thetadiff=8.0;
tau=1.0;
sigmaE=0;
sigmaG=1;

perror=0.01;

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
  end
end
save(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data3/data_sig_h_m_theta8.jld"),"n1mean",n1mean,"n2mean",n2mean,"x1mean",x1mean,"x2mean",x2mean,"pe",pe);


#Analysis over m & theta divergence
indmvec=collect(0.0:0.001:0.3);
hvec = collect(0.0:0.01:1.0);


n1mean=SharedArray(Float64,length(hvec),length(indmvec));
n2mean=SharedArray(Float64,length(hvec),length(indmvec));
x1mean=SharedArray(Float64,length(hvec),length(indmvec));
x2mean=SharedArray(Float64,length(hvec),length(indmvec));
pe=SharedArray(Float64,length(hvec),length(indmvec));

tmax=10000;
z=0.5;
rmax=2.0;
beta=0.001;
theta1=5.0;
thetadiff=5.0;
tau=1.0;
C=1000;
sigmaE=0;
sigmaG=1;

perror=0.01;

@sync @parallel for k=1:length(hvec)
  h = hvec[k];
  
  for j=1:length(indmvec)
    a0 = 1 - indmvec[j];
    
    n1, n2, x1, x2, w1, w2 = 
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
    burnin=0.80
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
save(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data3/data_sig_h_m_ddm.jld"),"n1mean",n1mean,"n2mean",n2mean,"x1mean",x1mean,"x2mean",x2mean,"pe",pe);

#THETADIFF = 3

#Analysis over m & theta divergence
indmvec=collect(0.0:0.001:0.3);
hvec = collect(0.0:0.01:1.0);


n1mean=SharedArray(Float64,length(hvec),length(indmvec));
n2mean=SharedArray(Float64,length(hvec),length(indmvec));
x1mean=SharedArray(Float64,length(hvec),length(indmvec));
x2mean=SharedArray(Float64,length(hvec),length(indmvec));
pe=SharedArray(Float64,length(hvec),length(indmvec));


tmax=10000;
z=0.5;
rmax=2.0;
beta=0.001;
theta1=5.0;
thetadiff=3.0;
tau=1.0;
C=1000;
sigmaE=0;
sigmaG=1;

perror=0.01;

@sync @parallel for k=1:length(hvec)
  h = hvec[k];
  
  for j=1:length(indmvec)
    a0 = 1 - indmvec[j];
    
    n1, n2, x1, x2, w1, w2 = 
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
    burnin=0.80
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
save(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data3/data_sig_h_m_theta3_ddm.jld"),"n1mean",n1mean,"n2mean",n2mean,"x1mean",x1mean,"x2mean",x2mean,"pe",pe);


#Analysis over m & theta divergence
indmvec=collect(0.0:0.001:0.3);
hvec = collect(0.0:0.01:1.0);


n1mean=SharedArray(Float64,length(hvec),length(indmvec));
n2mean=SharedArray(Float64,length(hvec),length(indmvec));
x1mean=SharedArray(Float64,length(hvec),length(indmvec));
x2mean=SharedArray(Float64,length(hvec),length(indmvec));
pe=SharedArray(Float64,length(hvec),length(indmvec));


tmax=100000;
z=0.5;
rmax=2.0;
beta=0.001;
theta1=5.0;
thetadiff=8.0;
tau=1.0;
C=1000;
sigmaE=0;
sigmaG=1;

perror=0.01;

@sync @parallel for k=1:length(hvec)
  h = hvec[k];
  
  for j=1:length(indmvec)
    a0 = 1 - indmvec[j];
    
    n1, n2, x1, x2, w1, w2 = 
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
    burnin=0.80
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
save(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data3/data_sig_h_m_theta8_ddm.jld"),"n1mean",n1mean,"n2mean",n2mean,"x1mean",x1mean,"x2mean",x2mean,"pe",pe);


######################
#LOAD AND BUILD FIGURE
######################

mvec=collect(0.0:0.001:0.3);
# sigmavec = collect(0.1:0.1:3.0);
hvec = collect(0.0:0.01:1.0);
#Import constant m version
d_m = load(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data3/data_sig_h_m.jld"));
#This loads the dictionary
n1mean_m = d_m["n1mean"];
n2mean_m = d_m["n2mean"];
x1mean_m = d_m["x1mean"];
x2mean_m = d_m["x2mean"];
pe_m = d_m["pe"];
pena_m = pe_m;
pena_m[find(x->x==true,isnan(pe_m))] = 1;

d_m3 = load(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data3/data_sig_h_m_theta3.jld"));
#This loads the dictionary
n1mean_m3 = d_m3["n1mean"];
n2mean_m3 = d_m3["n2mean"];
x1mean_m3 = d_m3["x1mean"];
x2mean_m3 = d_m3["x2mean"];
pe_m3 = d_m3["pe"];
pena_m3 = pe_m3;
pena_m3[find(x->x==true,isnan(pe_m3))] = 1;


d_m8 = load(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data3/data_sig_h_m_theta8.jld"));
#This loads the dictionary
n1mean_m8 = d_m8["n1mean"];
n2mean_m8 = d_m8["n2mean"];
x1mean_m8 = d_m8["x1mean"];
x2mean_m8 = d_m8["x2mean"];
pe_m8 = d_m8["pe"];
pena_m8 = pe_m8;
pena_m8[find(x->x==true,isnan(pe_m8))] = 1;



#Analysis over m & theta divergence
indmvec=collect(0.0:0.001:0.3);
C=1000;
d_ddm = load(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data3/data_sig_h_m_ddm.jld"));
#This loads the dictionary
n1mean_ddm = d_ddm["n1mean"];
n2mean_ddm = d_ddm["n2mean"];
x1mean_ddm = d_ddm["x1mean"];
x2mean_ddm = d_ddm["x2mean"];
pe_ddm = d_ddm["pe"];
pena_ddm = pe_ddm;
pena_ddm[find(x->x==true,isnan(pe_ddm))] = 1;

d_ddm3 = load(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data3/data_sig_h_m_theta3_ddm.jld"));
#This loads the dictionary
n1mean_ddm3 = d_ddm3["n1mean"];
n2mean_ddm3 = d_ddm3["n2mean"];
x1mean_ddm3 = d_ddm3["x1mean"];
x2mean_ddm3 = d_ddm3["x2mean"];
pe_ddm3 = d_ddm3["pe"];
pena_ddm3 = pe_ddm3;
pena_ddm3[find(x->x==true,isnan(pe_ddm3))] = 1;


d_ddm8 = load(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data3/data_sig_h_m_theta8_ddm.jld"));
#This loads the dictionary
n1mean_ddm8 = d_ddm8["n1mean"];
n2mean_ddm8 = d_ddm8["n2mean"];
x1mean_ddm8 = d_ddm8["x1mean"];
x2mean_ddm8 = d_ddm8["x2mean"];
pe_ddm8 = d_ddm8["pe"];
pena_ddm8 = pe_ddm8;
pena_ddm8[find(x->x==true,isnan(pe_ddm8))] = 1;




#Quantitatively compare pe for constant m and pe for ddm
mstar_constant = Array(Float64,length(hvec),length(indmvec));

mstar1 = Array(Float64,length(hvec),length(indmvec));
mstar2 = Array(Float64,length(hvec),length(indmvec));
mstarmax = Array(Float64,length(hvec),length(indmvec));
mstarmin = Array(Float64,length(hvec),length(indmvec));

mstar1_3 = Array(Float64,length(hvec),length(indmvec));
mstar2_3 = Array(Float64,length(hvec),length(indmvec));
mstarmax3 = Array(Float64,length(hvec),length(indmvec));
mstarmin3 = Array(Float64,length(hvec),length(indmvec));


mstar1_8 = Array(Float64,length(hvec),length(indmvec));
mstar2_8 = Array(Float64,length(hvec),length(indmvec));
mstarmax8 = Array(Float64,length(hvec),length(indmvec));
mstarmin8 = Array(Float64,length(hvec),length(indmvec));


for j=1:length(hvec)
  for k=1:length(indmvec)
    mstar_constant[j,k] = indmvec[k];
    
    mstar1[j,k] = indmvec[k]*(1 - (n1mean_ddm[j,k]/(C+n1mean_ddm[j,k])));
    mstar2[j,k] = indmvec[k]*(1 - (n2mean_ddm[j,k]/(C+n2mean_ddm[j,k])));
    
    mstarmax[j,k] = maximum([mstar1[j,k] mstar2[j,k]]);
    mstarmin[j,k] = minimum([mstar1[j,k] mstar2[j,k]]);
    
    mstar1_3[j,k] = indmvec[k]*(1 - (n1mean_ddm3[j,k]/(C+n1mean_ddm3[j,k])));
    mstar2_3[j,k] = indmvec[k]*(1 - (n2mean_ddm3[j,k]/(C+n2mean_ddm3[j,k])));
    
    mstarmax3[j,k] = maximum([mstar1_3[j,k] mstar2_3[j,k]]);
    mstarmin3[j,k] = minimum([mstar1_3[j,k] mstar2_3[j,k]]);
    
    mstar1_8[j,k] = indmvec[k]*(1 - (n1mean_ddm8[j,k]/(C+n1mean_ddm8[j,k])));
    mstar2_8[j,k] = indmvec[k]*(1 - (n2mean_ddm8[j,k]/(C+n2mean_ddm8[j,k])));
    
    mstarmax8[j,k] = maximum([mstar1_8[j,k] mstar2_8[j,k]]);
    mstarmin8[j,k] = minimum([mstar1_8[j,k] mstar2_8[j,k]]);
  end
end


#Median values 1:51
minh=1;
maxh=51;
medpe3 = Array{Float64}(length(mvec));
medpe5 = Array{Float64}(length(mvec));
medpe8 = Array{Float64}(length(mvec));

medpe3_ddm = Array{Float64}(length(mvec));
medpe5_ddm = Array{Float64}(length(mvec));
medpe8_ddm = Array{Float64}(length(mvec));

medindm3 = Array{Float64}(length(mvec));
medindm5 = Array{Float64}(length(mvec));
medindm8 = Array{Float64}(length(mvec));

medindm3max = Array{Float64}(length(mvec));
medindm3min = Array{Float64}(length(mvec));
medindm5max = Array{Float64}(length(mvec));
medindm5min = Array{Float64}(length(mvec));
medindm8max = Array{Float64}(length(mvec));
medindm8min = Array{Float64}(length(mvec));

for i=1:length(mvec)
  
  medpe3[i] = median(pena_m3[minh:maxh,i])
  medpe5[i] = median(pena_m[minh:maxh,i])
  medpe8[i] = median(pena_m8[minh:maxh,i])
  
  medpe3_ddm[i] = median(pena_ddm3[minh:maxh,i])
  medpe5_ddm[i] = median(pena_ddm[minh:maxh,i])
  medpe8_ddm[i] = median(pena_ddm8[minh:maxh,i])
  
  medindm3[i] = mean([mstar1_3[minh:maxh,i] mstar2_3[minh:maxh,i]])
  medindm5[i] = mean([mstar1[minh:maxh,i] mstar2[minh:maxh,i]])
  medindm8[i] = mean([mstar1_8[minh:maxh,i] mstar2_8[minh:maxh,i]])
  
  
  medindm3max[i] = mean(mstarmax3[minh:maxh,i])
  medindm3min[i] = mean(mstarmin3[minh:maxh,i])
  
  medindm5max[i] = mean(mstarmax[minh:maxh,i])
  medindm5min[i] = mean(mstarmin[minh:maxh,i])
  
  medindm8max[i] = mean(mstarmax8[minh:maxh,i])
  medindm8min[i] = mean(mstarmin8[minh:maxh,i])

  
end

bifsite5 = find(x->x==maximum(abs(diff(medpe5))),abs(diff(medpe5)));
bif5 = mvec[bifsite5];
bifsite8 = find(x->x==maximum(abs(diff(medpe8))),abs(diff(medpe8)));
bif8 = mvec[bifsite8];

bifsite5ddm = find(x->x==maximum(abs(diff(medpe5_ddm))),abs(diff(medpe5_ddm)));
bif5ddm = medindm5[bifsite5];
bifsite8ddm = find(x->x==maximum(abs(diff(medpe8_ddm))),abs(diff(medpe8_ddm)));
bif8ddm = medindm8[bifsite8ddm];



window = 5;
ma_mvec = movingaverage(mvec,window);
ma_medpe3 = movingaverage(medpe3,window);
ma_medpe5 = movingaverage(medpe5,window);
ma_medpe8 = movingaverage(medpe8,window);

window = 15;
ma_mvecsmooth = movingaverage(mvec,window);
ma_medpe3smooth = movingaverage(medpe3,window);
ma_medpe5smooth = movingaverage(medpe5,window);
ma_medpe8smooth = movingaverage(medpe8,window);

window = 5;
ma_medpe3_ddm = movingaverage(medpe3_ddm,window);
ma_medpe5_ddm = movingaverage(medpe5_ddm,window);
ma_medpe8_ddm = movingaverage(medpe8_ddm,window);
ma_medindm3min = movingaverage(medindm3min,window);
ma_medindm3max = movingaverage(medindm3max,window);
ma_medindm5min = movingaverage(medindm5min,window);
ma_medindm5max = movingaverage(medindm5max,window);
ma_medindm8min = movingaverage(medindm8min,window);
ma_medindm8max = movingaverage(medindm8max,window);

mlist = collect(1:2:length(ma_mvec));
mddmlist = collect(1:10:length(ma_medpe5_ddm));

#lines($(ma_mvecsmooth),$(ma_medpe8smooth),col=pal[3],lwd=3)
#lines($(ma_mvecsmooth),$(ma_medpe3smooth),col=pal[1],lwd=3)

namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/Manuscript/FinalDraft3/fig_thetaPEmvm.pdf");
R"""
library(RColorBrewer)
pal = brewer.pal(3,'Set1')
pdf($namespace,height=8,width=5)
par(mfrow=c(2,1),mai = c(0.8, 0.8, 0.1, 0.1))
plot($(ma_mvecsmooth),$(ma_medpe5smooth),col=pal[2],type='l',ylim=c(1,max($([ma_medpe3[mlist] ma_medpe5[mlist] ma_medpe8[mlist]]))),xlab='m, m*',ylab='PE',lwd=3,las=1)
arrows(0.053,3.18,0.053,3.00,length=0.05,angle=40,lwd=3)
text(0.053,3.32,'DCB')

points($(ma_medindm5max[mddmlist]),$(ma_medpe5_ddm[mddmlist]),col=pal[2],pch=16,cex=0.8)
points($(ma_medindm5min[mddmlist]),$(ma_medpe5_ddm[mddmlist]),col=pal[2],pch=16,cex=0.8)
l = length($(ma_medindm5max[mddmlist]))
for (i in 1:l) {
  segments($(ma_medindm5max[mddmlist])[i],$(ma_medpe5_ddm[mddmlist])[i],$(ma_medindm5min[mddmlist])[i],$(ma_medpe5_ddm[mddmlist])[i],col=paste(pal[2],'60',sep=''))
}
xleft<-0.0;xright<-0.45;ybottom<-1;ytop<-3;
text(-0.05,3.4,'(a)', xpd=TRUE)

plot($(ma_mvecsmooth),$(ma_medpe5smooth),col=pal[2],type='l',ylim=c(1,max($([ma_medpe3[mlist] ma_medpe5[mlist] ma_medpe8[mlist]]))),xlab='m, m*',ylab='PE',lwd=3,las=1)
lines($(ma_mvecsmooth),$(ma_medpe8smooth),col=pal[3],lwd=3)
lines($(ma_mvecsmooth),$(ma_medpe3smooth),col=pal[1],lwd=3)
#Or moving average
points($(ma_medindm3max[mddmlist]),$(ma_medpe3_ddm[mddmlist]),col=pal[1],pch=16,cex=0.8)
points($(ma_medindm3min[mddmlist]),$(ma_medpe3_ddm[mddmlist]),col=pal[1],pch=16,cex=0.8)
points($(ma_medindm5max[mddmlist]),$(ma_medpe5_ddm[mddmlist]),col=pal[2],pch=16,cex=0.8)
points($(ma_medindm5min[mddmlist]),$(ma_medpe5_ddm[mddmlist]),col=pal[2],pch=16,cex=0.8)
points($(ma_medindm8max[mddmlist]),$(ma_medpe8_ddm[mddmlist]),col=pal[3],pch=16,cex=0.8)
points($(ma_medindm8min[mddmlist]),$(ma_medpe8_ddm[mddmlist]),col=pal[3],pch=16,cex=0.8)
l = length($(ma_medindm5max[mddmlist]))
for (i in 1:l) {
  segments($(ma_medindm3max[mddmlist])[i],$(ma_medpe3_ddm[mddmlist])[i],$(ma_medindm3min[mddmlist])[i],$(ma_medpe3_ddm[mddmlist])[i],col=paste(pal[1],'60',sep=''))
  segments($(ma_medindm5max[mddmlist])[i],$(ma_medpe5_ddm[mddmlist])[i],$(ma_medindm5min[mddmlist])[i],$(ma_medpe5_ddm[mddmlist])[i],col=paste(pal[2],'60',sep=''))
  segments($(ma_medindm8max[mddmlist])[i],$(ma_medpe8_ddm[mddmlist])[i],$(ma_medindm8min[mddmlist])[i],$(ma_medpe8_ddm[mddmlist])[i],col=paste(pal[3],'60',sep=''))
}
text(-0.05,3.4,'(b)', xpd=TRUE)
legend(x=0.27,y=3.35,legend=c(3,5,8),col=pal,pch=22,xpd=TRUE,pt.bg=pal,cex=0.8, bty="n",title=expression(paste(Delta,theta)))
dev.off()
"""

