using Distributions

tmax=10000;

muvec=collect(0:0.001:2);
xss=Array{Float64}(length(muvec));

pdist = Normal(0,0.0000000001);

for i=1:length(muvec)
  x = zeros(tmax);
  x[1] = 0.1;
  for t=1:tmax-1
    x[t+1] = x[t]+ muvec[i]*x[t] - x[t]^3 + rand(pdist);
  end
  burnin=0.9;
  xtrim = x[Int64(floor(burnin*tmax)):tmax];
  xss[i]=mean(xtrim);
end

R"""
plot($muvec,$xss,pch='.')
"""



x = zeros(tmax);
x[1] = 0.1;
for t=1:tmax-1
  x[t+1] = x[t]+ 0.9*x[t] - x[t]^3 + rand(pdist);
end
R"plot($x,type='l',ylim=c(0,2))"