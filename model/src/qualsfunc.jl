function qualsfunc(ssval1, ssval2, extinct_threshold, similarity_threshold)
  d1 = length(ssval1[:,1,1]);
  d2 = length(ssval1[1,:,1]);
  d3 = length(ssval1[1,1,:]);
  
  out = zeros(Int64,d1,d2,d3).-1;
  # extinct_threshold = 5.0;
  # similarity_threshold = 10.0;
  
  #both are extinct
  out[find(x->x==true,(ssval1.<extinct_threshold) .* (ssval2.<extinct_threshold))] = 0;
  
  #One is extinct and one is not
  #out[find(x->(x==0.0 || x==Inf),ssval1.<extinct_threshold ./ ssval2.<extinct_threshold)] = 1;
  
  #Trajectories are similar
  out[find(x->x==true,(abs(ssval1.-ssval2).<similarity_threshold) .* (ssval1 .> extinct_threshold) .* (ssval2 .> extinct_threshold))] = 1;
  
  #Trajectories are in alternative stable states
  out[find(x->x==true,(abs(ssval1.-ssval2) .> similarity_threshold) .* (ssval1 .> extinct_threshold) .* (ssval2 .> extinct_threshold))] = 2;
  
return out

end