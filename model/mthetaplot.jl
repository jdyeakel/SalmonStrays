using RCall


thetascale = 2;
m=collect(0:0.01:0.5);
thetadiff=zeros(length(m));
for i = 1:length(m)
    thetadiff[i] = (1-2*m[i])/(thetascale*m[i])
end

namespace = string("$(homedir())/Dropbox/PostDoc/2017_SalmonStrays/manuscript/FinalDraft_rev/fig_mthetarelation.pdf");
R"""
pdf($namespace,height=5,width=6)
plot($thetadiff,$m,type='l',ylab='Straying rate m',xlab=expression(paste('Habitat heterogeneity ',Delta,theta)),ylim=c(0,0.5))
dev.off()
"""

