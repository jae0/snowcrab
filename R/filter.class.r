
  filter.class = function(x, type) {

    i=NULL

    if (type=="all") i = c(1:dim(x)[1])

    if (type=="mat")  i = which(x$mat==mature)
    if (type=="imm")  i = which(x$mat==immature)

    if (type=="male")  i = which(x$sex==male)
    if (type=="female")  i = which(x$sex==female)

    if (type=="m.mat")  i = which(x$sex==male & x$mat==mature)
    if (type=="m.imm")  i = which(x$sex==male & x$mat==immature)

    if (type=="m.com")  i = which(x$sex==male & x$cw>=95 & x$cw<200 )  # commerical size crab
    if (type=="m.ncom") i = which(x$sex==male & x$cw< 95 )

    if (type=="R0") { # "fishable" biomass (by sex, size, carapace condition and shell hardness) -- mature only
      i = which(x$sex==male & x$mat==mature & x$cw>=95 & x$cw<200 &  x$shell %in% c(3,4) )
      j = which(x$sex==male & x$mat==mature & x$cw>=95 & x$cw<200 &  x$shell==2 & x$durometer>=68 )
      k = which(x$sex==male & x$mat==mature & x$cw>=95 & x$cw<200 & !is.finite(x$shell) & x$durometer>=68)
      i = sort(unique(c(i,j,k)))
      # i = which(x$sex==male & x$mat==mature & x$cw>=95 & x$cw<200 & x$shell %in% c(3,4,5) )
    }
    
    if (type=="R1"){  # terminally moulted soft-shells (new recruits), mature only
      i = which(x$sex==male & x$mat==mature & x$cw>=95 & x$cw<200 &  x$shell==1 )
      j = which(x$sex==male & x$mat==mature & x$cw>=95 & x$cw<200 &  x$shell==2 & x$durometer<68 )
      k = which(x$sex==male & x$mat==mature & x$cw>=95 & x$cw<200 & !is.finite(x$shell) & x$durometer<68)
      i = sort(unique(c(i,j,k))) 
      # i = which(x$sex==male & x$mat==mature & x$cw>=95 & x$cw<200 & x$shell %in% c(1,2) )
    }
    
    if (type=="R2") i = which(x$sex==male & x$mat==immature & x$cw>=mb(8) & x$cw<95 )
    if (type=="R3") i = which(x$sex==male & x$mat==immature & x$cw>=mb(7) & x$cw<mb(8) )
    if (type=="R4") i = which(x$sex==male & x$mat==immature & x$cw<=mb(6) & x$cw<mb(7))
    if (type=="R5p") i = which(x$sex==male & x$mat==immature & x$cw<=mb(5) & x$cw<mb(6))

    if (type=="skip.moulter") i = which(x$sex==male & x$mat==immature & x$cw<=200 & (x$shell %in% c(3,4,5) ) )
    if (type=="m.dwarf") i = which(x$sex==male & x$mat==mature & x$cw<95 )

    if (type=="mi123") i = which(x$sex==male & x$mat==immature & x$cw<mb(1) )
    if (type=="mi4") i = which(x$sex==male & x$mat==immature & x$cw>=mb(1) & x$cw<mb(2))
    if (type=="mi5") i = which(x$sex==male & x$mat==immature & x$cw>=mb(2) & x$cw<mb(3)) # poor data beyond this point
    if (type=="mi6") i = which(x$sex==male & x$mat==immature & x$cw>=mb(3) & x$cw<mb(4)) # R7
    if (type=="mi7") i = which(x$sex==male & x$mat==immature & x$cw>=mb(4) & x$cw<mb(5)) # R6
    if (type=="mi8") i = which(x$sex==male & x$mat==immature & x$cw>=mb(5) & x$cw<mb(6)) # R5
    if (type=="mi9") i = which(x$sex==male & x$mat==immature & x$cw>=mb(6) & x$cw<mb(7)) # R4
    if (type=="mi10") i = which(x$sex==male & x$mat==immature & x$cw>=mb(7) & x$cw<mb(8)) # same as "R3"
    if (type=="mi11") i = which(x$sex==male & x$mat==immature & x$cw>=mb(8) & x$cw<mb(9)) # R0, R1 and R2
    if (type=="mi12") i = which(x$sex==male & x$mat==immature & x$cw>=mb(9) & x$cw<mb(10))# R0

    if (type=="mi9.skip.moulter") i = which(x$sex==male & x$mat==immature & x$cw>=mb(6) & x$cw<mb(7) & x$shell %in% c(3,4,5) ) # R4/5
    if (type=="mi10.skip.moulter") i = which(x$sex==male & x$mat==immature & x$cw>=mb(7) & x$cw<mb(8) & x$shell %in% c(3,4,5) ) # same as "R3/4"
    if (type=="mi11.skip.moulter") i = which(x$sex==male & x$mat==immature & x$cw>=mb(8) & x$cw<mb(9) & x$shell %in% c(3,4,5) ) # R0/1, R1/2 and R2/3
    if (type=="mi12.skip.moulter") i = which(x$sex==male & x$mat==immature & x$cw>=mb(9) & x$cw<mb(10) & x$shell %in% c(3,4,5) ) # R0/1

    if (type=="ma9") i = which(x$sex==male & x$mat==mature & x$cw>=mb(6) & x$cw<mb(7)) # dwarves
    if (type=="ma10") i = which(x$sex==male & x$mat==mature & x$cw>=mb(7) & x$cw<mb(8)) # dwarves
    if (type=="ma11") i = which(x$sex==male & x$mat==mature & x$cw>=mb(8) & x$cw<mb(9)) # mostly dwarves
    if (type=="ma12") i = which(x$sex==male & x$mat==mature & x$cw>=mb(9) & x$cw<mb(10))
    if (type=="ma13") i = which(x$sex==male & x$mat==mature & x$cw>=mb(10) )

    if (type=="ma9.CC1to2") i = which(x$sex==male & x$mat==mature & x$cw>=mb(6) & x$cw<mb(7) & (x$shell %in% c(1,2) ) )
    if (type=="ma10.CC1to2") i = which(x$sex==male & x$mat==mature & x$cw>=mb(7) & x$cw<mb(8) & (x$shell %in% c(1,2) ) )
    if (type=="ma11.CC1to2") i = which(x$sex==male & x$mat==mature & x$cw>=mb(8) & x$cw<mb(9) & (x$shell %in% c(1,2) ) )
    if (type=="ma12.CC1to2") i = which(x$sex==male & x$mat==mature & x$cw>=mb(9) & x$cw<mb(10) & (x$shell %in% c(1,2) ) )
    if (type=="ma13.CC1to2") i = which(x$sex==male & x$mat==mature & x$cw>=mb(10) & (x$shell %in% c(1,2) ) )

    if (type=="ma9.CC3to4") i = which(x$sex==male & x$mat==mature & x$cw>=mb(6) & x$cw<mb(7) & x$shell %in% c(3,4) )
    if (type=="ma10.CC3to4") i = which(x$sex==male & x$mat==mature & x$cw>=mb(7) & x$cw<mb(8) & x$shell %in% c(3,4) )
    if (type=="ma11.CC3to4") i = which(x$sex==male & x$mat==mature & x$cw>=mb(8) & x$cw<mb(9) & x$shell %in% c(3,4) )
    if (type=="ma12.CC3to4") i = which(x$sex==male & x$mat==mature & x$cw>=mb(9)  & x$cw<mb(10) & x$shell %in% c(3,4) )
    if (type=="ma13.CC3to4") i = which(x$sex==male & x$mat==mature & x$cw>=mb(10) & x$shell %in% c(3,4) )

    if (type=="ma9.CC5") i = which(x$sex==male & x$mat==mature & x$cw>=mb(6) & x$cw<mb(7) & x$shell %in% c(5) )
    if (type=="ma10.CC5") i = which(x$sex==male & x$mat==mature & x$cw>=mb(7) & x$cw<mb(8) & x$shell %in% c(5) )
    if (type=="ma11.CC5") i = which(x$sex==male & x$mat==mature & x$cw>=mb(8) & x$cw<mb(9) & x$shell %in% c(5) )
    if (type=="ma12.CC5") i = which(x$sex==male & x$mat==mature & x$cw>=mb(9) & x$cw<mb(10) & x$shell %in% c(5)  )
    if (type=="ma13.CC5") i = which(x$sex==male & x$mat==mature & x$cw>=mb(10) & x$shell %in% c(5)  )

    if (type=="fi1234") i = which(x$sex==female & x$mat==immature & x$cw<fb(1) )
    if (type=="fi5") i = which(x$sex==female & x$mat==immature & x$cw>=fb(1) & x$cw<fb(2))

    if (type=="fi6") i = which(x$sex==female & x$mat==immature & x$cw>=fb(2) & x$cw<fb(3))
    if (type=="fi7") i = which(x$sex==female & x$mat==immature & x$cw>=fb(3) & x$cw<fb(4))
    if (type=="fi8") i = which(x$sex==female & x$mat==immature & x$cw>=fb(4) & x$cw<fb(5))
    if (type=="fi9") i = which(x$sex==female & x$mat==immature & x$cw>=fb(5) & x$cw<fb(6))
    if (type=="fi10") i = which(x$sex==female & x$mat==immature & x$cw>=fb(6) )

    if (type=="fi6.adolescent") i = which(x$sex==female & x$mat==immature & x$cw>=fb(2) & x$cw<fb(3) & x$gonad<=2)
    if (type=="fi7.adolescent") i = which(x$sex==female & x$mat==immature & x$cw>=fb(3) & x$cw<fb(4) & x$gonad<=2)
    if (type=="fi8.adolescent") i = which(x$sex==female & x$mat==immature & x$cw>=fb(4) & x$cw<fb(5) & x$gonad<=2)
    if (type=="fi9.adolescent") i = which(x$sex==female & x$mat==immature & x$cw>=fb(5) & x$cw<fb(6) & x$gonad<=2)
    if (type=="fi10.adolescent") i = which(x$sex==female & x$mat==immature & x$cw>=fb(6) & x$gonad<=2 )

    if (type=="fi6.preprimiparous") i = which(x$sex==female & x$mat==immature & x$cw>=fb(2) & x$cw<fb(3) & x$gonad==3)
    if (type=="fi7.preprimiparous") i = which(x$sex==female & x$mat==immature & x$cw>=fb(3) & x$cw<fb(4) & x$gonad==3)
    if (type=="fi8.preprimiparous") i = which(x$sex==female & x$mat==immature & x$cw>=fb(4) & x$cw<fb(5) & x$gonad==3)
    if (type=="fi9.preprimiparous") i = which(x$sex==female & x$mat==immature & x$cw>=fb(5) & x$cw<fb(6) & x$gonad==3)
    if (type=="fi10.preprimiparous") i = which(x$sex==female & x$mat==immature & x$cw>=fb(6) & x$gonad==3 )

    if (type=="fa7") i = which(x$sex==female & x$mat==mature & x$cw>=fb(3) & x$cw<fb(4))
    if (type=="fa8") i = which(x$sex==female & x$mat==mature & x$cw>=fb(4) & x$cw<fb(5))
    if (type=="fa9") i = which(x$sex==female & x$mat==mature & x$cw>=fb(5) & x$cw<fb(6))
    if (type=="fa10") i = which(x$sex==female & x$mat==mature & x$cw>=fb(6) )

    if (type=="fa7.berried") i = which(x$sex==female & x$mat==mature & x$cw>=fb(3) & x$cw<fb(4) & x$eggPr %in% c(1:4) )
    if (type=="fa8.berried") i = which(x$sex==female & x$mat==mature & x$cw>=fb(4) & x$cw<fb(5) & x$eggPr %in% c(1:4) )
    if (type=="fa9.berried") i = which(x$sex==female & x$mat==mature & x$cw>=fb(5) & x$cw<fb(6) & x$eggPr %in% c(1:4) )
    if (type=="fa10.berried") i = which(x$sex==female & x$mat==mature & x$cw>=fb(6) & x$eggPr %in% c(1:4)  )

    if (type=="fa7.primiparous") i = which(x$sex==female & x$mat==mature & x$cw>=fb(3) & x$cw<fb(4) & x$shell<=2  )
    if (type=="fa8.primiparous") i = which(x$sex==female & x$mat==mature & x$cw>=fb(4) & x$cw<fb(5) & x$shell<=2   )
    if (type=="fa9.primiparous") i = which(x$sex==female & x$mat==mature & x$cw>=fb(5) & x$cw<fb(6) & x$shell<=2   )
    if (type=="fa10.primiparous") i = which(x$sex==female & x$mat==mature & x$cw>=fb(6) & x$shell<=2 )

    if (type=="fa7.multiparous") i = which(x$sex==female & x$mat==mature & x$cw>=fb(3) & x$cw<fb(4) & x$shell>=3  )
    if (type=="fa8.multiparous") i = which(x$sex==female & x$mat==mature & x$cw>=fb(4) & x$cw<fb(5) & x$shell>=3  )
    if (type=="fa9.multiparous") i = which(x$sex==female & x$mat==mature & x$cw>=fb(5) & x$cw<fb(6) & x$shell>=3  )
    if (type=="fa10.multiparous") i = which(x$sex==female & x$mat==mature & x$cw>=fb(6) & x$shell>=3   )

    if (type=="fa7.senile") i = which(x$sex==female & x$mat==mature & x$cw>=fb(3) & x$cw<fb(4) & x$shell >=4 & x$eggPr <=1 )
    if (type=="fa8.senile") i = which(x$sex==female & x$mat==mature & x$cw>=fb(4) & x$cw<fb(5) & x$shell >=4 & x$eggPr <=1 )
    if (type=="fa9.senile") i = which(x$sex==female & x$mat==mature & x$cw>=fb(5) & x$cw<fb(6) & x$shell >=4 & x$eggPr <=1 )
    if (type=="fa10.senile") i = which(x$sex==female & x$mat==mature & x$cw>=fb(6) & x$shell >=4 & x$eggPr <=1  )

    if (type=="pre.recruit") i = which(x$sex==male & x$mat==immature & x$cw>56 & x$cw<95)  # Koeller's request
     
    if (type=="f7" ) i = which(x$sex==female &x$cw>=fb(3) & x$cw<fb(4))
    if (type=="f8" ) i = which(x$sex==female &x$cw>=fb(4) & x$cw<fb(5))
    if (type=="f9" ) i = which(x$sex==female &x$cw>=fb(5) & x$cw<fb(6))
    if (type=="f10" ) i = which(x$sex==female &x$cw>=fb(6) )
    
    if (type=="m7" ) i = which(x$sex==male &x$cw>=mb(4) & x$cw<mb(5))
    if (type=="m8" ) i = which(x$sex==male &x$cw>=mb(5) & x$cw<mb(6))
    if (type=="m9" ) i = which(x$sex==male &x$cw>=mb(6) & x$cw<mb(7))
    if (type=="m10" ) i = which(x$sex==male &x$cw>=mb(7) & x$cw<mb(8))
    if (type=="m11" ) i = which(x$sex==male &x$cw>=mb(8) & x$cw<mb(9))
    if (type=="m12" ) i = which(x$sex==male &x$cw>=mb(9) & x$cw<mb(10))
    if (type=="m13" ) i = which(x$sex==male &x$cw>=mb(10) ) 

    if (type=="f.mat")   i = which(x$sex==female & x$mat==mature )
    if (type=="f.imm")   i = which(x$sex==female & x$mat==immature )

    if (type=="f.berried")  i = which(x$sex==female & x$eggPr %in% c(1:4) )
    if (type=="primiparous") i = which(x$sex==female & x$mat==mature & x$shell<=2 )
    if (type=="multiparous") i = which(x$sex==female & x$mat==mature & x$shell>=3 )
    if (type=="senile") i = which(x$sex==female & x$mat==mature & x$shell >=4 & x$eggPr <=1 )
    if (type=="adolescent") i = which(x$sex==female & x$mat==immature & x$gonad<=2)
    if (type=="preprimiparous") i = which(x$sex==female & x$mat==immature & x$gonad==3)

    if (type=="f.soft") {
      i = which(x$sex==female & x$shell==1)
      j = which(x$sex==female & !is.finite(x$shell))
      i = sort(unique(c(i,j)))
    }

    if (type=="f.hard") {
      i = which(x$sex==female & x$shell>1)
      j = which(x$sex==female & !is.finite(x$shell))
      i = sort(unique(c(i,j)))
    }

    if (type=="m.soft") {
      i = which(x$sex==male & x$shell==1)
      j = which(x$sex==male & !is.finite(x$shell))
      i = sort(unique(c(i,j)))
    }

    if (type=="m.hard") {
      i = which(x$sex==male & x$shell>1)
      j = which(x$sex==male & !is.finite(x$shell))
      i = sort(unique(c(i,j)))
    }
 
    if (type=="male.small") i = which(x$sex==male & x$cw<mb(8) )
    if (type=="male.large") i = which(x$sex==male & x$cw>=mb(8) )
    if (type=="female.small") i = which(x$sex==female & x$cw<fb(4) )
    if (type=="female.large") i = which(x$sex==female & x$cw>=fb(4) )

    if (type=="m.CC1to2") i = which(x$sex==male & x$cw>=95 & x$cw<=200 & (x$shell %in% c(1,2) ))
    if (type=="m.CC3to4") i = which(x$sex==male & x$cw>=95 & x$cw<=200 & (x$shell %in% c(3,4) ))
    if (type=="m.CC1") {  # legal comercial males, carapace condition 1
      j = which(x$sex==male & x$cw>=95 & x$cw<200 & x$shell ==1 )
      k = which(x$sex==male & x$cw>=95 & x$cw<200 & !is.finite(x$shell))
      i = sort(unique(c(j,k)))
    }
    if (type=="m.CC2") i = which(x$sex==male & x$cw>=95 & x$cw<=200 & x$shell==2 )
    if (type=="m.CC3") i = which(x$sex==male & x$cw>=95 & x$cw<=200 & x$shell==3 )
    if (type=="m.CC4") i = which(x$sex==male & x$cw>=95 & x$cw<=200 & x$shell==4 )
    if (type=="m.CC5") i = which(x$sex==male & x$cw>=95 & x$cw<=200 & x$shell==5 )

    if (type=="f.CC1to2") i = which(x$sex==female & (x$shell %in% c(1,2) ))
    if (type=="f.CC3to4") i = which(x$sex==female  & (x$shell %in% c(3,4) ))
    if (type=="f.CC3") i = which(x$sex==female & (x$shell %in% c(3) ))
    if (type=="f.CC4") i = which(x$sex==female & (x$shell %in% c(4) ))
    if (type=="f.CC5") i = which(x$sex==female & x$shell==5 )

    return(i)
  }



