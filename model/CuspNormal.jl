


b1vec = collect(-1:0.001:1);
b2vec = collect(-1:0.01:1);
tmax = 1000;
xtraj = Array{Float64}(length(b1vec),length(b2vec),tmax);
xss = Array{Float64}(length(b1vec),length(b2vec));
lambda = Array{Float64}(length(b1vec),length(b2vec));

for i=1:length(b1vec)
    for j=1:length(b2vec)
        b1 = b1vec[i];
        b2 = b2vec[j];
        x = Array{Float64}(tmax);
        x[1] = 0.1;
        for t = 2:(tmax-1)
            x[t+1] = x[t] + b1 + b2*x[t] - x[t]^3;
        end
        xtraj[i,j,:] = x;
        xss[i,j] = mean(x[tmax-100:tmax]);
        lambda[i,j] = 1 + b2 - 3*xss[i,j]^2; #$1 + 
    end
end

b1 = 0.12;
rb1 = find(x->x==b1,b1vec)[1];
b2 = 0.1;
r = find(x->x==b2,b2vec)[1];
R"""
par(mfrow=c(3,1))
plot($(xtraj[rb1,r,1:1000]),type='l')
plot($b1vec,$(lambda[:,r]),type='l',ylim=c(-1,1.2))
lines(seq(-50,50),rep(1,101))
plot($b1vec,$(xss[:,r]),pch=16,cex=0.5)
points(seq(-50,50),rep(1,101),pch=16,cex=0.5)
"""


b1 = 0.0;
r = find(x->x==b1,b1vec)[1];
R"""
par(mfrow=c(2,1))
plot($b2vec,$(lambda[r,:]),type='l',ylim=c(-1,2))
lines(seq(-50,50),rep(1,101))
plot($b2vec,$(xss[r,:]),pch=16,cex=0.5)
"""

