function KevinEvolve_asym(
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
  perror,
  asym #Ranges from 0 to 0.5 (at which point rmax will be 100% part)
  )
  # sgn = [-1,1];
  rv1 = 0.0;
  rv2 = 0.0;
  if asym > 0.0
      # rv1 = abs(rand(Normal(0,asym)));
      # rv2 = -abs(rand(Normal(0,asym)));
      rv1 = asym;
      rv2 = -asym;
  end
  rmax1 = rmax*(1+rv1);
  rmax2 = rmax1*(1+rv2);
  
  # sb = rand(sgn);
  beta1 = beta*(1+rv1);
  beta2 = beta*(1+rv2);
  
  #Initialize variables
  n1 = zeros(Float64,tmax); n1[1]=2;
  n2 = zeros(Float64,tmax); n2[1]=2;
  
  x1 = zeros(Float64,tmax); x1[1]=theta1 + rand(Normal(0,0.0001));
  x2 = zeros(Float64,tmax); x2[1]=(theta1 + thetadiff)  + rand(Normal(0,0.0001));
  
  w1 = zeros(Float64,tmax-1);
  w2 = zeros(Float64,tmax-1);
  meanfit1 = zeros(Float64,tmax-1);
  meanfit2 = zeros(Float64,tmax-1);
  
  #Process Error Distribution
  if perror == 0
    pdist = [0,0];
  else
    pdist = Normal(0,perror);
  end
  
  rfunc = function(mu,theta,rmax)
    r = rmax*tau*(exp(-((theta-mu)^2)/(2*((sigmaE^2 + sigmaG^2)+tau^2))))/(sqrt(tau^2+(sigmaE^2 + sigmaG^2)))
  end
  
  
  for t=1:tmax-1
    
      
      w1[t] = (n1[t]-m*n1[t])/(n1[t]-m*n1[t]+m*n2[t]);
        
      w2[t] = (n2[t]-m*n2[t])/(n2[t]-m*n2[t]+m*n1[t]);
      
      # n1[t+1] = ((1-m)*n1[t] + m*n2[t])*exp(-z) +
      # ((rfunc(x1[t],theta1) + rand(pdist))*(n1[t]-m*n1[t]) +
      # (rfunc(x2[t],theta1) + rand(pdist))*m*n2[t]) * 
      # exp(-beta*((1-m)*n1[t]+m*n2[t]));
      # 
      # n2[t+1] = ((1-m)*n2[t] + m*n1[t])*exp(-z) +
      # ((rfunc(x2[t],(theta1 + thetadiff)) + rand(pdist))*(n2[t]-m*n2[t]) +
      # (rfunc(x1[t],(theta1 + thetadiff)) + rand(pdist))*m*n1[t]) * 
      # exp(-beta*((1-m)*n2[t]+m*n1[t]));
      
      # n1[t+1] = ((1-m)*n1[t] + m*n2[t])*exp(-z) +
      # ((rfunc(w1[t]*x1[t] + (1-w1[t])*x2[t],theta1,rmax1) + rand(pdist))*((1-m)*n1[t]+m*n2[t])) * 
      # exp(-beta1*((1-m)*n1[t]+m*n2[t]));
      # 
      # n2[t+1] = ((1-m)*n2[t] + m*n1[t])*exp(-z) +
      # ((rfunc(w2[t]*x2[t] + (1-w2[t])*x1[t],(theta1 + thetadiff),rmax2) + rand(pdist))*((1-m)*n2[t]+m*n1[t])) * 
      # exp(-beta2*((1-m)*n2[t]+m*n1[t]));
      
      n1[t+1] = ((rfunc(w1[t]*x1[t] + (1-w1[t])*x2[t],theta1,rmax1) + rand(pdist))*((1-m)*n1[t]+m*n2[t])) * 
      exp(-beta1*((1-m)*n1[t]+m*n2[t]));
      
      n2[t+1] = ((rfunc(w2[t]*x2[t] + (1-w2[t])*x1[t],(theta1 + thetadiff),rmax2) + rand(pdist))*((1-m)*n2[t]+m*n1[t])) * 
      exp(-beta2*((1-m)*n2[t]+m*n1[t]));
      

      
      # meanfit1[t] = (rmax*w1[t]*(theta1-x1[t])*tau*exp(-((theta1-x1[t])^2)/(2*((sigmaE^2 + sigmaG^2)+tau^2)))) / 
      # ((((sigmaE^2 + sigmaG^2) + tau^2)^(3/2))*(((exp(-((theta1-x2[t])^2)/(2*((sigmaE^2 + sigmaG^2)+tau^2)))*rmax*(1-w1[t])*tau)/
      # (sqrt((sigmaE^2 + sigmaG^2)+tau^2))) +  ((exp(-((theta1-x1[t])^2)/(2*((sigmaE^2 + sigmaG^2)+tau^2)))*rmax*w1[t]*tau)/(sqrt((sigmaE^2 + sigmaG^2)+tau^2)))));
      
      # meanfit2[t] = (rmax*w2[t]*((theta1 + thetadiff)-x2[t])*tau*exp(-(((theta1 + thetadiff)-x2[t])^2)/(2*((sigmaE^2 + sigmaG^2)+tau^2)))) / 
      # ((((sigmaE^2 + sigmaG^2) + tau^2)^(3/2))*(((exp(-(((theta1 + thetadiff)-x1[t])^2)/(2*((sigmaE^2 + sigmaG^2)+tau^2)))*rmax*(1-w2[t])*tau)/
      # (sqrt((sigmaE^2 + sigmaG^2)+tau^2))) +  ((exp(-(((theta1 + thetadiff)-x2[t])^2)/(2*((sigmaE^2 + sigmaG^2)+tau^2)))*rmax*w2[t]*tau)/(sqrt((sigmaE^2 + sigmaG^2)+tau^2)))));
      
      meanfit1[t] = (theta1 - w1[t]*x1[t] - (1 - w1[t])*x2[t])/((sigmaE^2 + sigmaG^2) + tau^2);
      
      meanfit2[t] = ((theta1 + thetadiff) - w2[t]*x2[t] - (1 - w2[t])*x1[t])/((sigmaE^2 + sigmaG^2) + tau^2);
      
      x1[t+1] = w1[t]*x1[t] + (1-w1[t])*x2[t] + (h)*(( sigmaG^2))*meanfit1[t];
      
      x2[t+1] = w2[t]*x2[t] + (1-w2[t])*x1[t] + (h)*(( sigmaG^2))*meanfit2[t];
      
    
  end
  
  return(
  n1,
  n2,
  x1,
  x2,
  w1,
  w2
  )
  
end
