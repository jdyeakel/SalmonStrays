function timeSS(n1,n2,tstart)
  
  #Who went extinct? Calculate time to steady state for this one.
  if n1[tstart] == 0;
    n = n1;
  else
    n = n2;
  end
  
  #We know the future!
  tmax = length(n);
  finalmean = mean(n[tmax-100:tmax]);
  finalsd = std(n[tmax-100:tmax]);
  
  ss = false;
  vrange = 5;
  t = tstart+vrange;
  
  while ss == false
    #Steady state reached when n[t] is within the sd range of the final ss.
    if abs(diff([n[t],finalmean])[1]) <= finalsd
      ss = true;
    else
      ss = false;
      t = t+1;
    end
  end
  tss = t;
  
  relaxtime = tss - tstart;
  
  return(tss,relaxtime)
  
end