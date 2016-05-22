
  subset.biomass = function(Q, R, nodes, type="mature.only", fraction11 = 0.2) {
    out = NULL
    names(Q) = nodes
    names(R) = nodes

    Q[!is.finite(Q)] = 0
    R[!is.finite(R)] = 0

    if (type =="historic") {
      tot = fraction11 * (Q["imm.11"] + Q["imm.sm.11"] + Q["CC3to4.11"] + Q["CC5.11"]) +
                  (Q["imm.12"] + Q["imm.sm.12"] + Q["CC3to4.12"] + Q["CC5.12"]) +
                  (                               Q["CC3to4.13"] + Q["CC5.13"])

      err = sqrt( fraction11^2 * (R["imm.11"]^2 + R["imm.sm.11"]^2 + R["CC3to4.11"]^2 + R["CC5.11"]^2) +
                                  R["imm.12"]^2 + R["imm.sm.12"]^2 + R["CC3to4.12"]^2 + R["CC5.12"]^2 +
                                  R["CC3to4.13"]^2 + R["CC5.13"]^2
                )
    } else if (type =="mature.only") {
      tot = fraction11 * (  Q["CC3to4.11"] + Q["CC5.11"]) + 
                         (  Q["CC3to4.12"] + Q["CC5.12"]) +
                         (  Q["CC3to4.13"] + Q["CC5.13"])

      err = sqrt( fraction11^2 * (  R["CC3to4.11"]^2 + R["CC5.11"]^2) +
                                    R["CC3to4.12"]^2 + R["CC5.12"]^2 +
                                    R["CC3to4.13"]^2 + R["CC5.13"]^2
                )
    }
    
    names(tot) = NULL
    return( c(tot,err) )
  }


