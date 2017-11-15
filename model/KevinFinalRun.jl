


#Figure 1

#Analysis over m
tmax=10000;
mvec = collect(0.0:0.001:0.5);
n1ts = zeros(Float64,length(mvec),tmax);
n2ts = zeros(Float64,length(mvec),tmax);
n1mean=zeros(Float64,length(mvec));
n2mean=zeros(Float64,length(mvec));
n1sd=zeros(Float64,length(mvec));
n2sd=zeros(Float64,length(mvec));
aggmean=zeros(Float64,length(mvec));
aggsd=zeros(Float64,length(mvec));
x1mean=zeros(Float64,length(mvec));
x2mean=zeros(Float64,length(mvec));
pe=zeros(Float64,length(mvec));
eigs = Array(Array{Complex{Float64}},length(mvec));
maxeigs = Array{Float64}(length(mvec));
maximeigs = Array{Float64}(length(mvec));

z=0.5;
rmax=2.0;
beta=0.001;
theta1=5.0;
thetadiff=3;
tau=1.0;
h=0.2;
sigmaE=0;
sigmaG=1;
perror=0.01;

burnin=0.80
@time for i=1:length(mvec)
  m=mvec[i];
  
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
  x1trim = x1[Int64(floor(tmax*burnin)):tmax];
  x2trim = x2[Int64(floor(tmax*burnin)):tmax]
  
  n1ts[i,:] = n1;
  n2ts[i,:] = n2;
  
  n1mean[i] = mean(n1trim);
  n2mean[i] = mean(n2trim);
  n1sd[i] = std(n1trim);
  n2sd[i] = std(n2trim);
  
  aggmean[i] = mean(n1trim+n2trim);
  aggsd[i] = std(n1trim+n2trim);

  x1mean[i] = theta1-mean(x1trim);
  x2mean[i] = (theta1+thetadiff)-mean(x2trim);
  
  # #Calculate the Jacobian
  # Jac = KevinJacobian(mean(n1trim),mean(n2trim),mean(x1trim),mean(x2trim),
  # z,rmax,beta,theta1,thetadiff,tau,sigma,m)
  # eigs[i]=eigvals(Jac)
  # 
  # re = real(eigs[i]);
  # im = imag(eigs[i]);
  # maxeigs[i] = maximum(re);
  # maximeigs[i] = maximum(im);
  
  # pe[i] = (mean([std(n1trim),std(n2trim)])/mean([mean(n1trim),mean(n2trim)])) *
  # (1/(std(n1trim+n2trim)/mean(n1trim+n2trim)))
  
  pe[i] = mean([(std(n1trim)/mean(n1trim)),(std(n2trim)/mean(n2trim))])*
  (1/(std(n1trim+n2trim)/mean(n1trim+n2trim)))
  
end
#Steady state plot
namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/FinalDraft_rev/fig_traj.pdf");
R"""
library(RColorBrewer)
cols = brewer.pal(4,'Set1')
#pdf($namespace,height=8,width=5)
par(mfrow=c(2,1),mai = c(0.8, 0.8, 0.1, 0.1))
plot($mvec,$n1mean,pch='.',col=cols[4],xlab="Straying rate m",ylab="Steady state biomass",cex=0.5,ylim=c(0,max($n1mean)),las=1)
points($mvec,$n2mean,pch='.',col=cols[4],cex=0.5)
arrows($(mvec[indmax(pe)]),1200,$(mvec[indmax(pe)]),1100,length=0.05,angle=40,lwd=3)
text($(mvec[indmax(pe)]),1275,'DCB')
text(par('usr')[1]-0.09,1300,'(a)', xpd=TRUE)
plot($mvec,$x1mean,pch=".",col=cols[1],ylim=c(-5,5),xlab="Straying rate m",ylab="Trait offset",las=1)
points($mvec,$x2mean,pch=".",col=cols[2])
arrows($(mvec[indmax(pe)]),3.4,$(mvec[indmax(pe)]),2.5,length=0.05,angle=40,lwd=3)
text($(mvec[indmax(pe)]),4,'DCB')
text(par('usr')[1]-0.09,5,'(b)', xpd=TRUE)
#dev.off()
"""



# Figure 5

#Analysis over m & theta divergence
mvec=collect(0.0:0.001:0.45);
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
save(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data2/data_sig_h_m.jld"),"n1mean",n1mean,"n2mean",n2mean,"x1mean",x1mean,"x2mean",x2mean,"pe",pe);


#Analysis over m & theta divergence
mvec=collect(0.0:0.001:0.45);
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
save(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data2/data_sig_h_m_theta3.jld"),"n1mean",n1mean,"n2mean",n2mean,"x1mean",x1mean,"x2mean",x2mean,"pe",pe);


#Analysis over m & theta divergence
indmvec=collect(0.0:0.001:0.45);
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
save(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data2/data_sig_h_m_ddm.jld"),"n1mean",n1mean,"n2mean",n2mean,"x1mean",x1mean,"x2mean",x2mean,"pe",pe);

#Analysis over m & theta divergence
indmvec=collect(0.0:0.001:0.45);
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
save(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data2/data_sig_h_m_theta8_ddm.jld"),"n1mean",n1mean,"n2mean",n2mean,"x1mean",x1mean,"x2mean",x2mean,"pe",pe);

#THETADIFF = 3

#Analysis over m & theta divergence
indmvec=collect(0.0:0.001:0.45);
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
save(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data2/data_sig_h_m_theta3_ddm.jld"),"n1mean",n1mean,"n2mean",n2mean,"x1mean",x1mean,"x2mean",x2mean,"pe",pe);

######################
#LOAD AND BUILD FIGURE
######################

mvec=collect(0.0:0.001:0.45);
sigmavec = collect(0.1:0.1:3.0);
hvec = collect(0.0:0.01:1.0);
#Import constant m version
d_m = load(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data2/data_sig_h_m.jld"));
#This loads the dictionary
n1mean_m = d_m["n1mean"];
n2mean_m = d_m["n2mean"];
x1mean_m = d_m["x1mean"];
x2mean_m = d_m["x2mean"];
pe_m = d_m["pe"];
pena_m = pe_m;
pena_m[find(x->x==true,isnan(pe_m))] = 1;

d_m3 = load(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data2/data_sig_h_m_theta3.jld"));
#This loads the dictionary
n1mean_m3 = d_m3["n1mean"];
n2mean_m3 = d_m3["n2mean"];
x1mean_m3 = d_m3["x1mean"];
x2mean_m3 = d_m3["x2mean"];
pe_m3 = d_m3["pe"];
pena_m3 = pe_m3;
pena_m3[find(x->x==true,isnan(pe_m3))] = 1;


d_m8 = load(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data2/data_sig_h_m_theta8.jld"));
#This loads the dictionary
n1mean_m8 = d_m8["n1mean"];
n2mean_m8 = d_m8["n2mean"];
x1mean_m8 = d_m8["x1mean"];
x2mean_m8 = d_m8["x2mean"];
pe_m8 = d_m8["pe"];
pena_m8 = pe_m8;
pena_m8[find(x->x==true,isnan(pe_m8))] = 1;



#Analysis over m & theta divergence
indmvec=collect(0.0:0.001:0.45);
C=1000;
d_ddm = load(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data2/data_sig_h_m_ddm.jld"));
#This loads the dictionary
n1mean_ddm = d_ddm["n1mean"];
n2mean_ddm = d_ddm["n2mean"];
x1mean_ddm = d_ddm["x1mean"];
x2mean_ddm = d_ddm["x2mean"];
pe_ddm = d_ddm["pe"];
pena_ddm = pe_ddm;
pena_ddm[find(x->x==true,isnan(pe_ddm))] = 1;

d_ddm3 = load(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data2/data_sig_h_m_theta3_ddm.jld"));
#This loads the dictionary
n1mean_ddm3 = d_ddm3["n1mean"];
n2mean_ddm3 = d_ddm3["n2mean"];
x1mean_ddm3 = d_ddm3["x1mean"];
x2mean_ddm3 = d_ddm3["x2mean"];
pe_ddm3 = d_ddm3["pe"];
pena_ddm3 = pe_ddm3;
pena_ddm3[find(x->x==true,isnan(pe_ddm3))] = 1;


d_ddm8 = load(string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/model/data2/data_sig_h_m_theta8_ddm.jld"));
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

