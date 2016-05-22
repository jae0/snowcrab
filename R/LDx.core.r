LDx.core = function ( obj, cf=c(1,2), p=0.5, method="default", covars=NULL ) {
    
   # print( "Expecting: y ~ aX + bT + cU + dXT + eXU ..." )
   # print( "Order of terms is important: main effects first, then first order effects")
   # print( "Where X is the variable of interest ... dose, size, date, etc)" )
   # print( "When there are interaction terms, the first interaction term must be the second main effect" )

    # modified from MASS::dose.p
    # cf = indices of vars of interest (intercept, xvar): 
    #   default :: no interaction terms
    
    eta <- family(obj)$linkfun(p)  #logit at associated Prob
    a <- coef(obj)
    out = NULL

    if (method == "default" ) {
      # identicle to MASS::dose.p
      x.p <- (eta - a[1])/a[2]
      pd <- -cbind(1, x.p)/a[cf[2]]
      SD <- sqrt(((pd %*% vcov(obj)[cf, cf]) * pd) %*% c(1, 1))
    }    

    if (method == "2main.0interaction" ) {
      # Solve[ 0  == eta - (a1 + a2*X + a3*T ), X] ;  x = main factor of interest
      T = covars
      x.p = (eta -  a[1] - a[3]* T)/  a[2] 
      vcv =  vcov(obj) 
      vcv[lower.tri(vcv)] = 0
      dy = vcv * 0
      Q = a[2] 
      dy1 = - 1/Q     # partial derivative wrt a[1] -- the intercept... etc
      dy2 = - ( eta - a[1] - a[3]*T )/Q^2 
      dy3 = - T / Q 
      diag(dy) = c( dy1^2, dy2^2, dy3^2 ) # partial derivatives associated with variances
      dy[1,2:3] = dy1 * c(dy2, dy3 )
      dy[2,3  ] = dy2 * c(dy3 )
      SD <- sqrt( sum(dy * vcv) )
      out <- data.frame(cbind( p=p, dose=x.p, SD=SD, Cov1=T ) )

    }


    if (method == "2main.1interaction" ) {
      # Solve[ 0  == eta - (a1 + a2*X + a3*T + a4*X*T), X] ;  x = main factor of interest
      T = covars
      x.p = (eta - a[1] - a[3]* T)/ ( a[2] + a[4]*T) 
      vcv =  vcov(obj) 
      vcv[lower.tri(vcv)] = 0
      dy = vcv * 0
      Q = a[2] + a[4] * T
      dy1 = - 1/Q     # partial derivative wrt a[1] -- the intercept... etc
      dy2 = - (-a[1] + eta -a[3]*T) / Q^2 
      dy3 = - T / Q 
      dy4 = - T * (-a[1] + eta -a[3]*T) / Q^2 
      diag(dy) = c( dy1^2, dy2^2, dy3^2, dy4^2 ) # partial derivatives associated with variances
      dy[1,2:4] = dy1 * c(dy2, dy3, dy4 )
      dy[2,3:4] = dy2 * c(dy3, dy4 )
      dy[3,4]   = dy3 * c(dy4 )
      SD <- sqrt( sum(dy * vcv) )
      out <- data.frame(cbind( p=p, dose=x.p, SD=SD, Cov1=T ) )

      
    }
  
    if (method=="3main.0interaction") {
      # Solve[eta == a1 + a2*X + a3*T + a4*U , X] ; x = main factor of interest
      
      T = covars[1]
      U = covars[2]
      x.p = ( eta -a[1]- a[3]*T - a[4]*U )/ a[2] 
      vcv =  vcov(obj) 
      vcv[lower.tri(vcv)] = 0
      dy = vcv * 0
      Q = a[2]
      dy1 = -1/Q
      dy2 = (-a[1] + eta - a[3]*T - a[4] *U)/Q^2
      dy3 = -T/Q
      dy4 = -U/Q
      
      diag(dy) = c( dy1^2, dy2^2, dy3^2, dy4^2 ) # partial derivatives associated with variances
      dy[1,2:4] = dy1 * c(dy2, dy3, dy4 )
      dy[2,3:4] = dy2 * c(     dy3, dy4 )
      dy[3,4  ] = dy3 * c(          dy4 )
      
      SD <- sqrt( sum(dy * vcv) )
      out <- data.frame(cbind( p=p, dose=x.p, SD=SD, Cov1=T, Cov2=U ) )

  
    }

    if (method=="3main.1interaction") {
      # Solve[eta == a1 + a2*X + a3*T + a4*U + a5*X*T, X] ; x = main factor of interest
      
      T = covars[1]
      U = covars[2]
      x.p = ( eta -a[1]- a[3]*T - a[4]*U )/ (a[2] + a[5]*T)
      vcv =  vcov(obj) 
      vcv[lower.tri(vcv)] = 0
      dy = vcv * 0
      Q = a[2] + a[5]*T 
      dy1 = -1/Q
      dy2 = -(-a[1] + eta - a[3]*T - a[4] *U)/Q^2
      dy3 = -T/Q
      dy4 = -U/Q
      dy5 = -T * (-a[1] + eta - a[3]*T - a[4]* U)/Q^2
      
      diag(dy) = c( dy1^2, dy2^2, dy3^2, dy4^2, dy5^2 ) # partial derivatives associated with variances
      dy[1,2:5] = dy1 * c(dy2, dy3, dy4, dy5 )
      dy[2,3:5] = dy2 * c(     dy3, dy4, dy5 )
      dy[3,4:5] = dy3 * c(          dy4, dy5 )
      dy[4,  5] = dy4 * c(               dy5 )
      
      SD <- sqrt( sum(dy * vcv) )
      out <- data.frame(cbind( p=p, dose=x.p, SD=SD, Cov1=T, Cov2=U ) )

   
    }
      
    if (method=="3main.2interaction.12.13") {
       # Solve[eta == a1 + a2*X + a3*T + a4*U + a5*X*T + a6*X*U, X] ; x = main factor of interest
     
  	  T = covars[1]
      U = covars[2]
      x.p = (eta - a[1] - a[3]*T - a[4]*U) / (a[2] + a[5]*T + a[6]*U)
      
      vcv =  vcov(obj) 
      vcv[lower.tri(vcv)] = 0
      dy = vcv * 0
      Q = a[2] + a[5]*T + a[6]* U
      
      dy1 = -1/Q
      dy2 = -(-a[1] + eta - a[3]*T - a[4] *U)/Q^2
      dy3 = -T/Q
      dy4 = -U/Q
      dy5 = -T * (-a[1] + eta - a[3]*T - a[4]* U)/Q^2
      dy6 = -U * (-a[1] + eta - a[3]*T - a[4]* U)/Q^2
      
      diag(dy) = c( dy1^2, dy2^2, dy3^2, dy4^2, dy5^2, dy6^2 ) # partial derivatives associated with variances
      dy[1,2:6] = dy1 * c(dy2, dy3, dy4, dy5, dy6 )
      dy[2,3:6] = dy2 * c(     dy3, dy4, dy5, dy6 )
      dy[3,4:6] = dy3 * c(          dy4, dy5, dy6 )
      dy[4,5:6] = dy4 * c(               dy5, dy6 )
      dy[5,  6] = dy5 * c(                    dy6 )
      
      SD <- sqrt( sum(dy * vcv) )
      out <- data.frame(cbind( p=p, dose=x.p, SD=SD, Cov1=T, Cov2=U ) )

    }
 
    if (method=="3main.2interaction.12.23") {
       # Solve[eta == a1 + a2*X + a3*T + a4*U + a5*X*T + a6*T*U, X] ; x = main factor of interest
     
  	  T = covars[1]
      U = covars[2]
      x.p = (eta - a[1] - a[3]*T - a[4]*U - a[6]*T * U) / (a[2] + a[5]*T )
      
      vcv =  vcov(obj) 
      vcv[lower.tri(vcv)] = 0
      dy = vcv * 0
      Q = a[2] + a[5]*T 
      
      dy1 = -1/Q
      dy2 = -(-a[1] + eta - a[3]*T - a[4] *U -a[6]*T*U)/Q^2
      dy3 = -T/Q
      dy4 = -U/Q
      dy5 = -T * (-a[1] + eta - a[3]*T - a[4]* U-a[6]*T*U)/Q^2
      dy6 = -U * (-a[1] + eta - a[3]*T - a[4]* U-a[6]*T*U)/Q^2
      
      diag(dy) = c( dy1^2, dy2^2, dy3^2, dy4^2, dy5^2, dy6^2 ) # partial derivatives associated with variances
      dy[1,2:6] = dy1 * c(dy2, dy3, dy4, dy5, dy6 )
      dy[2,3:6] = dy2 * c(     dy3, dy4, dy5, dy6 )
      dy[3,4:6] = dy3 * c(          dy4, dy5, dy6 )
      dy[4,5:6] = dy4 * c(               dy5, dy6 )
      dy[5,  6] = dy5 * c(                    dy6 )
      
      SD <- sqrt( sum(dy * vcv) )
      out <- data.frame(cbind( p=p, dose=x.p, SD=SD, Cov1=T, Cov2=U ) )

    }


    if (method=="3main.3interaction"  ) {
       # Solve[eta == a1 + a2*X + a3*T + a4*U + a5*X*T + a6*X*U + a7*T*U, X] ; x = main factor of interest
     
  	  T = covars[1]
      U = covars[2]
      x.p = (eta - a[1] - a[3]*T - a[4]*U -  a[7]*T*U) / (a[2] + a[5]*T + a[6]*U)
      
      vcv =  vcov(obj) 
      vcv[lower.tri(vcv)] = 0
      dy = vcv * 0
      Q = a[2] + a[5]*T + a[6]* U
      
      dy1 = -1/Q
      dy2 = -(-a[1] + eta - a[3]*T - a[4] *U - a[7]*T*U)/Q^2
      dy3 = -T/Q
      dy4 = -U/Q
      dy5 = -T * (-a[1] + eta - a[3]*T - a[4]* U - a[7]*T*U)/Q^2
      dy6 = -U * (-a[1] + eta - a[3]*T - a[4]* U - a[7]*T*U)/Q^2
      dy7 = - T * U  / Q
      
      diag(dy) = c( dy1^2, dy2^2, dy3^2, dy4^2, dy5^2, dy6^2, dy7^2 ) # partial derivatives associated with variances
      dy[1,2:7] = dy1 * c(dy2, dy3, dy4, dy5, dy6, dy7 )
      dy[2,3:7] = dy2 * c(     dy3, dy4, dy5, dy6, dy7  )
      dy[3,4:7] = dy3 * c(          dy4, dy5, dy6, dy7  )
      dy[4,5:7] = dy4 * c(               dy5, dy6, dy7  )
      dy[5,6:7] = dy5 * c(                    dy6, dy7  )
      dy[6,  7] = dy6 * c(                         dy7  )
      
      SD <- sqrt( sum(dy * vcv) )
      out <- data.frame(cbind( p=p, dose=x.p, SD=SD, Cov1=T, Cov2=U ) )

    }
    rownames(out) = NULL
    return(out)
}

