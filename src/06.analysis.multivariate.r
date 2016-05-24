
 ...  incomplete

#   source(file.path(ordination.package, "ordination.package.r"))
	loadfunctions( "sorted.ordination") 
	loadfunctions( "snowcrab", functionname="initialise.local.environment.r") 
	


  set = snowcrab.db("set.with.cat")
  allvars = c(
      "totno.all", "totno.male", "totno.male.com", "totno.male.ncom", "totmass.male.imm", "totmass.male.mat",
      "totmass.all", "totmass.male", "totmass.male.com", "totmass.male.ncom", "totno.male.imm", "totno.male.mat",

      "totno.male.soft", "totno.male.hard",
      "totmass.male.soft", "totmass.male.hard",

      "totno.female.soft", "totno.female.hard",
      "totmass.female.soft", "totmass.female.hard",

      "R0.no", "R1.no", "R2.no", "R3.no", "R4.no", "R5p.no", "dwarf.no",
      "R0.mass", "R1.mass", "R2.mass", "R3.mass", "R4.mass", "R5p.mass", "dwarf.mass",

      "totmass.female", "totmass.female.mat", "totmass.female.imm",
      "totmass.female.berried", "totmass.female.primiparous", "totmass.female.multiparous",

      "totno.female", "totno.female.mat", "totno.female.imm",
      "totno.female.berried", "totno.female.primiparous", "totno.female.multiparous",

      "totmass.male.com.CC1to2", "totmass.male.com.CC3to4", "totmass.male.com.CC5",
      "totno.male.com.CC1to2", "totno.male.com.CC3to4", "totno.male.com.CC5",

      "cw.fem.mat.mean", "cw.fem.imm.mean",
      "cw.male.mat.mean", "cw.male.imm.mean",

      "cw.mean", "cw.comm.mean", "cw.notcomm.mean",

      "sexratio.imm", "sexratio.mat", "sexratio.all",

      "z", "t",

      "uniqueid"
    )

    x = set[, allvars]
    for (v in allvars)  x[,v] = recode.variable(x[,v], v)$x  # normalise the data where required

    rownames(x) = x$uniqueid
    x$uniqueid = NULL
    p = pca.modified(x)

 


