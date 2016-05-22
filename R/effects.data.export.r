

effects.data.export = function( Y ) {
 
  Z = Y$x
  Z$predictedPr = inverse.logit( Y$fit )
  Z$se = Y$se # keep on logit scale as it is meaningless of the Pr. scale.
  Z$upper = inverse.logit( Y$upper )
  Z$lower = inverse.logit( Y$lower )

  return (Z)

}


