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
  
  #This function is only appropriate for sigmae = 0;
  sigma = copy(sigmaG);
  
  #Define Jacobian Matrix
  Jac = Array{Float64}(4,4);

  Jac[1,1]= ((m - 1)*(-(rmax*
        tau*(m^2*
           N2star^2*(sigma^2*(1 - beta*m*N2star) - 
             beta*m*N2star*tau^2 + tau^2 - 
                      x2star*(theta1 + x1star) + theta1*x1star + 
             x2star^2) + beta*(m - 1)^3*N1star^3*
                   (sigma^2 + tau^2) - (m - 1)^2*
           N1star^2*(sigma^2 + tau^2)*(3*beta*m*N2star - 1) + 
                 (m - 1)*m*N1star*
           N2star*(sigma^2*(3*beta*m*N2star - 2) + 
             tau^2*(3*beta*m*N2star - 2) - 
                      (theta1 - x1star)*(x1star - x2star)))*
        e^(beta*(m - 1)*N1star - beta*m*N2star - 
                   ((m - 1)*N1star*(theta1 - x1star) + 
               m*N2star*(x2star - theta1))^2/
                     (2*(sigma^2 + 
                tau^2)*(-(m*N1star) + m*N2star + N1star)^2))) - 
         ((sigma^2 + tau^2)^(3/2)*(-(m*N1star) + m*N2star + 
           N1star)^2)/e^z))/
   ((sigma^2 + tau^2)^(3/2)*(-(m*N1star) + m*N2star + N1star)^2)

  Jac[1,2]= ((m*(N1star - m*N1star + m*N2star)^2*(sigma^2 + tau^2)^(3/2))/e^z + 
      e^(beta*(-1 + m)*N1star - beta*m*N2star - 
             ((-1 + m)*N1star*(theta1 - x1star) + 
           m*N2star*(-theta1 + x2star))^2/
               (2*(N1star - m*N1star + m*N2star)^2*(sigma^2 + 
            tau^2)))*m*rmax*tau*
        (beta*(-1 + m)^3*N1star^3*(sigma^2 + tau^2) - 
      m^2*N2star^2*(-1 + beta*m*N2star)*
             (sigma^2 + tau^2) - (-1 + m)^2*
       N1star^2*((-1 + 3*beta*m*N2star)*sigma^2 + 
                (-1 + 3*beta*m*N2star)*
          tau^2 + (theta1 - x1star)*(x1star - x2star)) + 
           (-1 + m)*m*N1star*
       N2star*((-2 + 3*beta*m*N2star)*
          sigma^2 + (-2 + 3*beta*m*N2star)*tau^2 + 
                (theta1 - x2star)*(x1star - x2star))))/((N1star - 
      m*N1star + m*N2star)^2*
      (sigma^2 + tau^2)^(3/2))


  Jac[1,3]= -((e^(beta*(-1 + m)*N1star - beta*m*N2star - 
              ((-1 + m)*N1star*(theta1 - x1star) + 
            m*N2star*(-theta1 + x2star))^2/
                (2*(N1star - m*N1star + m*N2star)^2*(sigma^2 + 
             tau^2)))*(-1 + m)*N1star*rmax*tau*
         ((-1 + m)*N1star*(theta1 - x1star) + 
       m*N2star*(-theta1 + x2star)))/
      (((-1 + m)*N1star - m*N2star)*(sigma^2 + tau^2)^(3/2)))


  Jac[1,4]= (e^(beta*(-1 + m)*N1star - beta*m*N2star - 
           ((-1 + m)*N1star*(theta1 - x1star) + 
          m*N2star*(-theta1 + x2star))^2/
             (2*(N1star - m*N1star + m*N2star)^2*(sigma^2 + tau^2)))*
   m*N2star*rmax*tau*
      ((-(-1 + m))*N1star*(theta1 - x1star) + 
     m*N2star*(theta1 - x2star)))/
   ((N1star - m*N1star + m*N2star)*(sigma^2 + tau^2)^(3/2))


  Jac[2,1]= (e^((-beta)*(m*(N1star - N2star) + N2star) - 
             (m*N1star*theta2 + N2star*theta2 - m*N2star*theta2 - 
           m*N1star*x1star + 
                    (-1 + m)*N2star*
            x2star)^2/(2*(m*(N1star - N2star) + N2star)^2*(sigma^2 + 
            tau^2)))*m*
        (m*(N1star - N2star) + N2star)^2*rmax*
    tau*(sigma^2 + tau^2) + 
      (m*(m*(N1star - N2star) + N2star)^2*(sigma^2 + tau^2)^(3/2))/
    e^z - 
      e^((-beta)*(m*(N1star - N2star) + N2star) - 
             (m*N1star*theta2 + N2star*theta2 - m*N2star*theta2 - 
           m*N1star*x1star + 
                    (-1 + m)*N2star*
            x2star)^2/(2*(m*(N1star - N2star) + N2star)^2*(sigma^2 + 
            tau^2)))*rmax*
        tau*(beta*
       m*(m*(N1star - N2star) + N2star)^3*(sigma^2 + tau^2) + 
           (-1 + m)*m*
       N2star*(x1star - x2star)*(m*N1star*theta2 + N2star*theta2 - 
         m*N2star*theta2 - 
                
         m*N1star*x1star + (-1 + m)*N2star*x2star)))/((m*(N1star - 
         N2star) + N2star)^2*
      (sigma^2 + tau^2)^(3/2))


  Jac[2,2]= ((-1 + m)*((-e^((-beta)*(m*(N1star - N2star) + N2star) - 
                   (m*N1star*theta2 + N2star*theta2 - 
               m*N2star*theta2 - m*N1star*x1star + 
                          (-1 + m)*N2star*
                x2star)^2/(2*(m*(N1star - N2star) + 
                 N2star)^2*(sigma^2 + tau^2))))*
           (m*(N1star - N2star) + N2star)^2*rmax*
      tau*(sigma^2 + tau^2) - 
         ((m*(N1star - N2star) + N2star)^2*(sigma^2 + tau^2)^(3/2))/
      e^z - 
         e^((-beta)*(m*(N1star - N2star) + N2star) - 
                (m*N1star*theta2 + N2star*theta2 - m*N2star*theta2 - 
             m*N1star*x1star + 
                       (-1 + m)*N2star*
              x2star)^2/(2*(m*(N1star - N2star) + 
               N2star)^2*(sigma^2 + tau^2)))*
           rmax*
      tau*((-beta)*(m*(N1star - N2star) + N2star)^3*(sigma^2 + 
           tau^2) + 
              
        m*N1star*(x1star - x2star)*(N2star*(-theta2 + x2star) + 
                   
           m*((-N1star)*theta2 + N2star*theta2 + N1star*x1star - 
              N2star*x2star)))))/
   ((m*(N1star - N2star) + N2star)^2*(sigma^2 + tau^2)^(3/2))


  Jac[2,3] = (e^((-beta)*(m*(N1star - N2star) + N2star) - 
           (m*N1star*theta2 + N2star*theta2 - m*N2star*theta2 - 
          m*N1star*x1star + 
                  (-1 + m)*N2star*
           x2star)^2/(2*(m*(N1star - N2star) + N2star)^2*(sigma^2 + 
           tau^2)))*m*
      N1star*rmax*
   tau*(m*N1star*theta2 + N2star*theta2 - m*N2star*theta2 - 
     m*N1star*x1star + 
         (-1 + m)*N2star*x2star))/((m*(N1star - N2star) + 
     N2star)*(sigma^2 + tau^2)^(3/2))
        
  Jac[2,4] = (e^((-beta)*(m*(N1star - N2star) + N2star) - 
           (m*N1star*theta2 + N2star*theta2 - m*N2star*theta2 - 
          m*N1star*x1star + 
                  (-1 + m)*N2star*
           x2star)^2/(2*(m*(N1star - N2star) + N2star)^2*(sigma^2 + 
           tau^2)))*
      (-1 + m)*N2star*rmax*
   tau*(m*(N1star*(-theta2 + x1star) + N2star*(theta2 - x2star)) + 
         N2star*(-theta2 + x2star)))/((m*(N1star - N2star) + 
     N2star)*(sigma^2 + tau^2)^(3/2))
     
  Jac[3,1] = -(((-1 + m)*m*
     N2star*(x1star - x2star)*(sigma^2 + tau^2 - h*sigma^2))/
      ((N1star - m*N1star + m*N2star)^2*(sigma^2 + tau^2)))
          
  Jac[3,2] = ((-1 + m)*m*N1star*(x1star - x2star)*(sigma^2 + tau^2 - h*sigma^2))/
   ((N1star - m*N1star + m*N2star)^2*(sigma^2 + tau^2))
     
  Jac[3,3] = ((-1 + m)*
   N1star*(sigma^2 + tau^2 - h*sigma^2))/(((-1 + m)*N1star - 
     m*N2star)*(sigma^2 + tau^2))
           
  Jac[3,4] = (m*N2star*(sigma^2 + tau^2 - h*sigma^2))/((N1star - m*N1star + 
     m*N2star)*(sigma^2 + tau^2))
        
  Jac[4,1] = -(((-1 + m)*m*
     N2star*(x1star - x2star)*(sigma^2 + tau^2 - h*sigma^2))/
      ((m*(N1star - N2star) + N2star)^2*(sigma^2 + tau^2)))
     
  Jac[4,2] = ((-1 + m)*m*N1star*(x1star - x2star)*(sigma^2 + tau^2 - h*sigma^2))/
   ((m*(N1star - N2star) + N2star)^2*(sigma^2 + tau^2))
     
  Jac[4,3] = (m*N1star*(sigma^2 + tau^2 - h*sigma^2))/((m*(N1star - N2star) + 
     N2star)*(sigma^2 + tau^2))
            
  Jac[4,4] = -(((-1 + m)*
     N2star*(sigma^2 + tau^2 - h*sigma^2))/((m*(N1star - N2star) + 
       N2star)*
         (sigma^2 + tau^2)))
           
  return Jac
end
