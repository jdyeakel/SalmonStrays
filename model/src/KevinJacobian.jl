function KevinJacobian(
  N1star,
  N2star,
  x1star,
  x2star,
  z,
  rmax,
  beta,
  theta1,
  thetadiff,
  tau,
  h,
  sigmaE,
  sigmaG,
  m
  )
  
  theta2 = theta1+thetadiff;
  
  #This function is only appropriate for sigmaE = 0;
  sigma = copy(sigmaG);
  
  #Define Jacobian Matrix
  Jac = Array{Float64}(4,4);

  Jac[1,1]= ((1 - m)*rmax*tau*e^(-(beta*((1 - m)*N1star + m*N2star)) - 
      (theta1 - x1star)^2/(2*(sigma^2 + tau^2))))/
   sqrt(sigma^2 + tau^2) - 
  (beta*(1 - m)*((rmax*tau*(N1star - m*N1star))/
      (e^((theta1 - x1star)^2/(2*(sigma^2 + tau^2)))*
       sqrt(sigma^2 + tau^2)) + (m*N2star*rmax*tau)/
      (e^((theta1 - x2star)^2/(2*(sigma^2 + tau^2)))*
       sqrt(sigma^2 + tau^2))))/e^(beta*((1 - m)*N1star + m*N2star)) + 
  (1 - m)/e^z

  Jac[1,2]= m/e^z + (e^((-beta)*((1 - m)*N1star + m*N2star) - (theta1 - x2star)^2/(2*(sigma^2 + tau^2)))*m*
    rmax*tau)/sqrt(sigma^2 + tau^2) - 
  (beta*m*(((N1star - m*N1star)*rmax*tau)/(e^((theta1 - x1star)^2/(2*(sigma^2 + tau^2)))*
       sqrt(sigma^2 + tau^2)) + (m*N2star*rmax*tau)/
      (e^((theta1 - x2star)^2/(2*(sigma^2 + tau^2)))*sqrt(sigma^2 + tau^2))))/
   e^(beta*((1 - m)*N1star + m*N2star))


  Jac[1,3]= (e^((-beta)*((1 - m)*N1star + m*N2star) - (theta1 - x1star)^2/(2*(sigma^2 + tau^2)))*
   (N1star - m*N1star)*rmax*tau*(theta1 - x1star))/(sigma^2 + tau^2)^(3/2)


  Jac[1,4]= (e^((-beta)*((1 - m)*N1star + m*N2star) - (theta1 - x2star)^2/(2*(sigma^2 + tau^2)))*m*N2star*
   rmax*tau*(theta1 - x2star))/(sigma^2 + tau^2)^(3/2)


  Jac[2,1]= m/e^z + (e^((-beta)*(m*N1star + (1 - m)*N2star) - (theta2 - x1star)^2/(2*(sigma^2 + tau^2)))*m*
    rmax*tau)/sqrt(sigma^2 + tau^2) - 
  (beta*m*((m*N1star*rmax*tau)/(e^((theta2 - x1star)^2/(2*(sigma^2 + tau^2)))*
       sqrt(sigma^2 + tau^2)) + ((N2star - m*N2star)*rmax*tau)/
      (e^((theta2 - x2star)^2/(2*(sigma^2 + tau^2)))*sqrt(sigma^2 + tau^2))))/
   e^(beta*(m*N1star + (1 - m)*N2star))


  Jac[2,2]= (1 - m)/e^z + (e^((-beta)*(m*N1star + (1 - m)*N2star) - (theta2 - x2star)^2/
       (2*(sigma^2 + tau^2)))*(1 - m)*rmax*tau)/sqrt(sigma^2 + tau^2) - 
  (beta*(1 - m)*((m*N1star*rmax*tau)/(e^((theta2 - x1star)^2/(2*(sigma^2 + tau^2)))*
       sqrt(sigma^2 + tau^2)) + ((N2star - m*N2star)*rmax*tau)/
      (e^((theta2 - x2star)^2/(2*(sigma^2 + tau^2)))*sqrt(sigma^2 + tau^2))))/
   e^(beta*(m*N1star + (1 - m)*N2star))


  Jac[2,3] = (e^((-beta)*(m*N1star + (1 - m)*N2star) - (theta2 - x1star)^2/(2*(sigma^2 + tau^2)))*m*N1star*
   rmax*tau*(theta2 - x1star))/(sigma^2 + tau^2)^(3/2)
        
  Jac[2,4] = (e^((-beta)*(m*N1star + (1 - m)*N2star) - (theta2 - x2star)^2/(2*(sigma^2 + tau^2)))*
   (N2star - m*N2star)*rmax*tau*(theta2 - x2star))/(sigma^2 + tau^2)^(3/2)
     
  Jac[3,1] = -((h*(N1star - m*N1star)*rmax*sigma^2*tau*(-(((1 - m)*(N1star - m*N1star)*rmax*tau)/
         (e^((theta1 - x1star)^2/(2*(sigma^2 + tau^2)))*((N1star - m*N1star + m*N2star)^2*
           sqrt(sigma^2 + tau^2)))) + ((1 - m)*rmax*tau)/
        (e^((theta1 - x1star)^2/(2*(sigma^2 + tau^2)))*((N1star - m*N1star + m*N2star)*
          sqrt(sigma^2 + tau^2))) + 
       ((((1 - m)*(N1star - m*N1star))/(N1star - m*N1star + m*N2star)^2 - 
          (1 - m)/(N1star - m*N1star + m*N2star))*rmax*tau)/
        (e^((theta1 - x2star)^2/(2*(sigma^2 + tau^2)))*sqrt(sigma^2 + tau^2)))*
      (theta1 - x1star))/e^((theta1 - x1star)^2/(2*(sigma^2 + tau^2)))/
    ((N1star - m*N1star + m*N2star)*(sigma^2 + tau^2)^(3/2)*
     (((N1star - m*N1star)*rmax*tau)/(e^((theta1 - x1star)^2/(2*(sigma^2 + tau^2)))*
         ((N1star - m*N1star + m*N2star)*sqrt(sigma^2 + tau^2))) + 
       ((1 - (N1star - m*N1star)/(N1star - m*N1star + m*N2star))*rmax*tau)/
        (e^((theta1 - x2star)^2/(2*(sigma^2 + tau^2)))*sqrt(sigma^2 + tau^2)))^2)) - 
  (h*(1 - m)*(N1star - m*N1star)*rmax*sigma^2*tau*(theta1 - x1star))/
    e^((theta1 - x1star)^2/(2*(sigma^2 + tau^2)))/((N1star - m*N1star + m*N2star)^2*
    (sigma^2 + tau^2)^(3/2)*(((N1star - m*N1star)*rmax*tau)/
      (e^((theta1 - x1star)^2/(2*(sigma^2 + tau^2)))*((N1star - m*N1star + m*N2star)*
        sqrt(sigma^2 + tau^2))) + ((1 - (N1star - m*N1star)/(N1star - m*N1star + m*N2star))*
       rmax*tau)/(e^((theta1 - x2star)^2/(2*(sigma^2 + tau^2)))*sqrt(sigma^2 + tau^2)))) + 
  (h*(1 - m)*rmax*sigma^2*tau*(theta1 - x1star))/e^((theta1 - x1star)^2/(2*(sigma^2 + tau^2)))/
   ((N1star - m*N1star + m*N2star)*(sigma^2 + tau^2)^(3/2)*
    (((N1star - m*N1star)*rmax*tau)/(e^((theta1 - x1star)^2/(2*(sigma^2 + tau^2)))*
       ((N1star - m*N1star + m*N2star)*sqrt(sigma^2 + tau^2))) + 
     ((1 - (N1star - m*N1star)/(N1star - m*N1star + m*N2star))*rmax*tau)/
      (e^((theta1 - x2star)^2/(2*(sigma^2 + tau^2)))*sqrt(sigma^2 + tau^2)))) - 
  ((1 - m)*(N1star - m*N1star)*x1star)/(N1star - m*N1star + m*N2star)^2 + 
  ((1 - m)*x1star)/(N1star - m*N1star + m*N2star) + 
  (((1 - m)*(N1star - m*N1star))/(N1star - m*N1star + m*N2star)^2 - 
    (1 - m)/(N1star - m*N1star + m*N2star))*x2star
          
  Jac[3,2] = -((h*(N1star - m*N1star)*rmax*sigma^2*tau*(-((m*(N1star - m*N1star)*rmax*tau)/
         (e^((theta1 - x1star)^2/(2*(sigma^2 + tau^2)))*((N1star - m*N1star + m*N2star)^2*
           sqrt(sigma^2 + tau^2)))) + (m*(N1star - m*N1star)*rmax*tau)/
        (e^((theta1 - x2star)^2/(2*(sigma^2 + tau^2)))*((N1star - m*N1star + m*N2star)^2*
          sqrt(sigma^2 + tau^2))))*(theta1 - x1star))/
     e^((theta1 - x1star)^2/(2*(sigma^2 + tau^2)))/((N1star - m*N1star + m*N2star)*
     (sigma^2 + tau^2)^(3/2)*(((N1star - m*N1star)*rmax*tau)/
        (e^((theta1 - x1star)^2/(2*(sigma^2 + tau^2)))*((N1star - m*N1star + m*N2star)*
          sqrt(sigma^2 + tau^2))) + ((1 - (N1star - m*N1star)/(N1star - m*N1star + m*N2star))*
         rmax*tau)/(e^((theta1 - x2star)^2/(2*(sigma^2 + tau^2)))*sqrt(sigma^2 + tau^2)))^2)) - 
  (h*m*(N1star - m*N1star)*rmax*sigma^2*tau*(theta1 - x1star))/
    e^((theta1 - x1star)^2/(2*(sigma^2 + tau^2)))/((N1star - m*N1star + m*N2star)^2*
    (sigma^2 + tau^2)^(3/2)*(((N1star - m*N1star)*rmax*tau)/
      (e^((theta1 - x1star)^2/(2*(sigma^2 + tau^2)))*((N1star - m*N1star + m*N2star)*
        sqrt(sigma^2 + tau^2))) + ((1 - (N1star - m*N1star)/(N1star - m*N1star + m*N2star))*
       rmax*tau)/(e^((theta1 - x2star)^2/(2*(sigma^2 + tau^2)))*sqrt(sigma^2 + tau^2)))) - 
  (m*(N1star - m*N1star)*x1star)/(N1star - m*N1star + m*N2star)^2 + 
  (m*(N1star - m*N1star)*x2star)/(N1star - m*N1star + m*N2star)^2
     
  Jac[3,3] = (N1star - m*N1star)/(N1star - m*N1star + m*N2star) - 
  (h*(N1star - m*N1star)*rmax*sigma^2*tau)/e^((theta1 - x1star)^2/(2*(sigma^2 + tau^2)))/
   ((N1star - m*N1star + m*N2star)*(sigma^2 + tau^2)^(3/2)*
    (((N1star - m*N1star)*rmax*tau)/(e^((theta1 - x1star)^2/(2*(sigma^2 + tau^2)))*
       ((N1star - m*N1star + m*N2star)*sqrt(sigma^2 + tau^2))) + 
     ((1 - (N1star - m*N1star)/(N1star - m*N1star + m*N2star))*rmax*tau)/
      (e^((theta1 - x2star)^2/(2*(sigma^2 + tau^2)))*sqrt(sigma^2 + tau^2)))) - 
  (h*(N1star - m*N1star)^2*rmax^2*sigma^2*tau^2*(theta1 - x1star)^2)/
    e^((theta1 - x1star)^2/(sigma^2 + tau^2))/((N1star - m*N1star + m*N2star)^2*
    (sigma^2 + tau^2)^3*(((N1star - m*N1star)*rmax*tau)/
       (e^((theta1 - x1star)^2/(2*(sigma^2 + tau^2)))*((N1star - m*N1star + m*N2star)*
         sqrt(sigma^2 + tau^2))) + ((1 - (N1star - m*N1star)/(N1star - m*N1star + m*N2star))*
        rmax*tau)/(e^((theta1 - x2star)^2/(2*(sigma^2 + tau^2)))*sqrt(sigma^2 + tau^2)))^2) + 
  (h*(N1star - m*N1star)*rmax*sigma^2*tau*(theta1 - x1star)^2)/
    e^((theta1 - x1star)^2/(2*(sigma^2 + tau^2)))/((N1star - m*N1star + m*N2star)*
    (sigma^2 + tau^2)^(5/2)*(((N1star - m*N1star)*rmax*tau)/
      (e^((theta1 - x1star)^2/(2*(sigma^2 + tau^2)))*((N1star - m*N1star + m*N2star)*
        sqrt(sigma^2 + tau^2))) + ((1 - (N1star - m*N1star)/(N1star - m*N1star + m*N2star))*
       rmax*tau)/(e^((theta1 - x2star)^2/(2*(sigma^2 + tau^2)))*sqrt(sigma^2 + tau^2))))
           
  Jac[3,4] = 1 - (N1star - m*N1star)/(N1star - m*N1star + m*N2star) - 
  (e^(-((theta1 - x1star)^2/(2*(sigma^2 + tau^2))) - (theta1 - x2star)^2/(2*(sigma^2 + tau^2)))*
    h*(N1star - m*N1star)*(1 - (N1star - m*N1star)/(N1star - m*N1star + m*N2star))*rmax^2*
    sigma^2*tau^2*(theta1 - x1star)*(theta1 - x2star))/((N1star - m*N1star + m*N2star)*
    (sigma^2 + tau^2)^3*(((N1star - m*N1star)*rmax*tau)/
       (e^((theta1 - x1star)^2/(2*(sigma^2 + tau^2)))*((N1star - m*N1star + m*N2star)*
         sqrt(sigma^2 + tau^2))) + ((1 - (N1star - m*N1star)/(N1star - m*N1star + m*N2star))*
        rmax*tau)/(e^((theta1 - x2star)^2/(2*(sigma^2 + tau^2)))*sqrt(sigma^2 + tau^2)))^2)
        
  Jac[4,1] = (m*(N2star - m*N2star)*x1star)/(m*N1star + N2star - m*N2star)^2 - 
  (h*(N2star - m*N2star)*rmax*sigma^2*tau*((m*(N2star - m*N2star)*rmax*tau)/
       (e^((theta2 - x1star)^2/(2*(sigma^2 + tau^2)))*((m*N1star + N2star - m*N2star)^2*
         sqrt(sigma^2 + tau^2))) - (m*(N2star - m*N2star)*rmax*tau)/
       (e^((theta2 - x2star)^2/(2*(sigma^2 + tau^2)))*((m*N1star + N2star - m*N2star)^2*
         sqrt(sigma^2 + tau^2))))*(theta2 - x2star))/
    e^((theta2 - x2star)^2/(2*(sigma^2 + tau^2)))/((m*N1star + N2star - m*N2star)*
    (sigma^2 + tau^2)^(3/2)*(((N2star - m*N2star)*rmax*tau)/
       (e^((theta2 - x2star)^2/(2*(sigma^2 + tau^2)))*((m*N1star + N2star - m*N2star)*
         sqrt(sigma^2 + tau^2))) + ((1 - (N2star - m*N2star)/(m*N1star + N2star - m*N2star))*
        rmax*tau)/(e^((theta2 - x1star)^2/(2*(sigma^2 + tau^2)))*sqrt(sigma^2 + tau^2)))^2) - 
  (h*m*(N2star - m*N2star)*rmax*sigma^2*tau*(theta2 - x2star))/
    e^((theta2 - x2star)^2/(2*(sigma^2 + tau^2)))/((m*N1star + N2star - m*N2star)^2*
    (sigma^2 + tau^2)^(3/2)*(((N2star - m*N2star)*rmax*tau)/
      (e^((theta2 - x2star)^2/(2*(sigma^2 + tau^2)))*((m*N1star + N2star - m*N2star)*
        sqrt(sigma^2 + tau^2))) + ((1 - (N2star - m*N2star)/(m*N1star + N2star - m*N2star))*
       rmax*tau)/(e^((theta2 - x1star)^2/(2*(sigma^2 + tau^2)))*sqrt(sigma^2 + tau^2)))) - 
  (m*(N2star - m*N2star)*x2star)/(m*N1star + N2star - m*N2star)^2
     
  Jac[4,2] = (((1 - m)*(N2star - m*N2star))/(m*N1star + N2star - m*N2star)^2 - 
    (1 - m)/(m*N1star + N2star - m*N2star))*x1star - 
  (h*(N2star - m*N2star)*rmax*sigma^2*tau*
     (-(((1 - m)*(N2star - m*N2star)*rmax*tau)/(e^((theta2 - x2star)^2/(2*(sigma^2 + tau^2)))*
         ((m*N1star + N2star - m*N2star)^2*sqrt(sigma^2 + tau^2)))) + 
      ((1 - m)*rmax*tau)/(e^((theta2 - x2star)^2/(2*(sigma^2 + tau^2)))*
        ((m*N1star + N2star - m*N2star)*sqrt(sigma^2 + tau^2))) + 
      ((((1 - m)*(N2star - m*N2star))/(m*N1star + N2star - m*N2star)^2 - 
         (1 - m)/(m*N1star + N2star - m*N2star))*rmax*tau)/
       (e^((theta2 - x1star)^2/(2*(sigma^2 + tau^2)))*sqrt(sigma^2 + tau^2)))*
     (theta2 - x2star))/e^((theta2 - x2star)^2/(2*(sigma^2 + tau^2)))/
   ((m*N1star + N2star - m*N2star)*(sigma^2 + tau^2)^(3/2)*
    (((N2star - m*N2star)*rmax*tau)/(e^((theta2 - x2star)^2/(2*(sigma^2 + tau^2)))*
        ((m*N1star + N2star - m*N2star)*sqrt(sigma^2 + tau^2))) + 
      ((1 - (N2star - m*N2star)/(m*N1star + N2star - m*N2star))*rmax*tau)/
       (e^((theta2 - x1star)^2/(2*(sigma^2 + tau^2)))*sqrt(sigma^2 + tau^2)))^2) - 
  (h*(1 - m)*(N2star - m*N2star)*rmax*sigma^2*tau*(theta2 - x2star))/
    e^((theta2 - x2star)^2/(2*(sigma^2 + tau^2)))/((m*N1star + N2star - m*N2star)^2*
    (sigma^2 + tau^2)^(3/2)*(((N2star - m*N2star)*rmax*tau)/
      (e^((theta2 - x2star)^2/(2*(sigma^2 + tau^2)))*((m*N1star + N2star - m*N2star)*
        sqrt(sigma^2 + tau^2))) + ((1 - (N2star - m*N2star)/(m*N1star + N2star - m*N2star))*
       rmax*tau)/(e^((theta2 - x1star)^2/(2*(sigma^2 + tau^2)))*sqrt(sigma^2 + tau^2)))) + 
  (h*(1 - m)*rmax*sigma^2*tau*(theta2 - x2star))/e^((theta2 - x2star)^2/(2*(sigma^2 + tau^2)))/
   ((m*N1star + N2star - m*N2star)*(sigma^2 + tau^2)^(3/2)*
    (((N2star - m*N2star)*rmax*tau)/(e^((theta2 - x2star)^2/(2*(sigma^2 + tau^2)))*
       ((m*N1star + N2star - m*N2star)*sqrt(sigma^2 + tau^2))) + 
     ((1 - (N2star - m*N2star)/(m*N1star + N2star - m*N2star))*rmax*tau)/
      (e^((theta2 - x1star)^2/(2*(sigma^2 + tau^2)))*sqrt(sigma^2 + tau^2)))) - 
  ((1 - m)*(N2star - m*N2star)*x2star)/(m*N1star + N2star - m*N2star)^2 + 
  ((1 - m)*x2star)/(m*N1star + N2star - m*N2star)
     
  Jac[4,3] = 1 - (N2star - m*N2star)/(m*N1star + N2star - m*N2star) - 
  (e^(-((theta2 - x1star)^2/(2*(sigma^2 + tau^2))) - (theta2 - x2star)^2/(2*(sigma^2 + tau^2)))*
    h*(N2star - m*N2star)*(1 - (N2star - m*N2star)/(m*N1star + N2star - m*N2star))*rmax^2*
    sigma^2*tau^2*(theta2 - x1star)*(theta2 - x2star))/((m*N1star + N2star - m*N2star)*
    (sigma^2 + tau^2)^3*(((N2star - m*N2star)*rmax*tau)/
       (e^((theta2 - x2star)^2/(2*(sigma^2 + tau^2)))*((m*N1star + N2star - m*N2star)*
         sqrt(sigma^2 + tau^2))) + ((1 - (N2star - m*N2star)/(m*N1star + N2star - m*N2star))*
        rmax*tau)/(e^((theta2 - x1star)^2/(2*(sigma^2 + tau^2)))*sqrt(sigma^2 + tau^2)))^2)
            
  Jac[4,4] = (N2star - m*N2star)/(m*N1star + N2star - m*N2star) - 
  (h*(N2star - m*N2star)*rmax*sigma^2*tau)/e^((theta2 - x2star)^2/(2*(sigma^2 + tau^2)))/
   ((m*N1star + N2star - m*N2star)*(sigma^2 + tau^2)^(3/2)*
    (((N2star - m*N2star)*rmax*tau)/(e^((theta2 - x2star)^2/(2*(sigma^2 + tau^2)))*
       ((m*N1star + N2star - m*N2star)*sqrt(sigma^2 + tau^2))) + 
     ((1 - (N2star - m*N2star)/(m*N1star + N2star - m*N2star))*rmax*tau)/
      (e^((theta2 - x1star)^2/(2*(sigma^2 + tau^2)))*sqrt(sigma^2 + tau^2)))) - 
  (h*(N2star - m*N2star)^2*rmax^2*sigma^2*tau^2*(theta2 - x2star)^2)/
    e^((theta2 - x2star)^2/(sigma^2 + tau^2))/((m*N1star + N2star - m*N2star)^2*
    (sigma^2 + tau^2)^3*(((N2star - m*N2star)*rmax*tau)/
       (e^((theta2 - x2star)^2/(2*(sigma^2 + tau^2)))*((m*N1star + N2star - m*N2star)*
         sqrt(sigma^2 + tau^2))) + ((1 - (N2star - m*N2star)/(m*N1star + N2star - m*N2star))*
        rmax*tau)/(e^((theta2 - x1star)^2/(2*(sigma^2 + tau^2)))*sqrt(sigma^2 + tau^2)))^2) + 
  (h*(N2star - m*N2star)*rmax*sigma^2*tau*(theta2 - x2star)^2)/
    e^((theta2 - x2star)^2/(2*(sigma^2 + tau^2)))/((m*N1star + N2star - m*N2star)*
    (sigma^2 + tau^2)^(5/2)*(((N2star - m*N2star)*rmax*tau)/
      (e^((theta2 - x2star)^2/(2*(sigma^2 + tau^2)))*((m*N1star + N2star - m*N2star)*
        sqrt(sigma^2 + tau^2))) + ((1 - (N2star - m*N2star)/(m*N1star + N2star - m*N2star))*
       rmax*tau)/(e^((theta2 - x1star)^2/(2*(sigma^2 + tau^2)))*sqrt(sigma^2 + tau^2))))
           
  return Jac
end
