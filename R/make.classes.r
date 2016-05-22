
  make.classes = function(sex) {

      if (sex==male) {
        yclass = as.character( c( 5:13 ) )
        cats = c("imm", "imm.sm", "CC1to2", "CC3to4", "CC5")
        varmap = list( imm = paste("mi", yclass, ".no", sep="" ),
                       imm.sm = paste("mi", yclass, ".skip.moulter.no", sep="" ),
                       CC1to2 = paste("ma", yclass, ".CC1to2.no", sep="" ),
                       CC3to4 = paste("ma", yclass, ".CC3to4.no", sep="" ),
                       CC5 = paste("ma", yclass, ".CC5.no", sep="" )
                    )
      }
      if (sex==female) {
        yclass = as.character( c( 5:10 ) )
        cats = c("imm", "adol", "preprimipar", "berried", "primipar", "multipar", "senile")
        varmap = list( imm = paste("fi", yclass, ".no", sep="" ),
                       adol = paste("fi", yclass, ".adolescent.no", sep="" ),
                       preprimipar = paste("fi", yclass, ".preprimiparous.no", sep="" ),
                       berried = paste("fa", yclass, ".berried.no", sep="" ),
                       primipar = paste("fa", yclass, ".primiparous.no", sep="" ),
                       multipar = paste("fa", yclass, ".multiparous.no", sep="" ),
                       senile = paste("fa", yclass, ".senile.no", sep="" )
                     )
      }

      return( list (yclass=yclass, cats=cats, varmap=varmap) )
  }


