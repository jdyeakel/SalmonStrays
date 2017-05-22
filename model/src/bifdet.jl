function bifdet(
  n1mean,
  n2mean,
  mvec,
  hvec
  )
  diffvec = abs(n1mean - n2mean);
  bifvalue = zeros(Float64,length(hvec),2);
  #over hvec
  for i = 1:length(hvec)
    vec = diffvec[i,:];
    aasvec = find(x->x>1,vec);
    slip = diff(diffvec[i,:]);
    if length(aasvec)==0
      bifvalue[i,:] = [0 0];
    else
      bifvalue[i,:] = [mvec[indmax(abs(slip))] hvec[i]];
    end
  end
  return(bifvalue)
end
