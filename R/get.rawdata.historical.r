# ----------------------------
# obtain raw biological data
# odb.historic = get.rawdata(outfile = file.path(p$current.assessment.year, "odb.historic.rdata") )
# save(odb.historic, file=file.path(p$annual.results, "odb.historic.rdata"), compress=T)

#  ss = get.rawdata (outfile="odb.historical.seasamples.datadump.rdata", type="seasamples")
#  ps = get.rawdata (outfile="odb.historical.portsamples.datadump.rdata", type="portsamples")



  get.rawdata.historical = function(type) {

    outfile=file.path( project.datadirectory("snowcrab"), "data", "observer", "odb.historical.datadump")
    directories = file.path( project.datadirectory("snowcrab"), "data", "observer", "archive", c(2001:2003))
    fl = list.files(path=directories, pattern="[*.txt]$", full.names=T, recursive=T)

    if (type=="seasamples") fl = fl[ grep( "sea sample", fl, ignore.case=T ) ]
    if (type=="portsamples") fl = fl[ grep( "port sample", fl, ignore.case=T ) ]

    temp.output = "tmp.dump"
    if (file.exists( temp.output ) )  file.remove(temp.output)

    for (i in fl) {
      print ( paste("reading",i) )
      out = NULL
      rawdata=NULL
      read0 = scan( file=i, what="", sep="\n", strip.white=F )  # there are control characters at the end of a file on occasion
      read0 = read0[which(nchar(read0)>140)]
      write.table( read0, file=temp.output, append=T, row.names=F, col.names=F, quote=F, sep=";" )
    }

    # data format is column-width structured
    w = c(9,1, 2, 1, 9, 3, 1, 1, 1, 3, 1,
          4, 2, 4, 4, 1,
          1, 2, 1, 6, 1, 14, 1, 2, 1, 6, 3, 6,
          3, 3, 2, 6, 1, 2, 1, 1, 1, 10, 5, 23, 8, 18, 1, 1, 1 )
    v = c("sdate", "d0", "trap", "dx1", "tripno", "crabno", "d4", "sex", "d0", "carapacewidth", "d5",
          "ndays", "trapno", "clawheightRight", "clawheightLeft", "shell_condition",
          "shell_condition_3M", "durometer", "dx2", "cfv", "dx3", "comment", "dx4", "ll", "dx5", "latitude", "dx6", "longitude",
          "d8", "depth", "clawlength", "weight", "d9", "zone", "d10", "data", "d11", "digits", "d12", "observer", "d13",
          "vessel", "d14", "eggmaturity", "d16" )
    keep = c( "sdate", "trap", "tripno", "crabno", "sex", "carapacewidth",
          "ndays", "trapno", "clawheightRight", "clawheightLeft", "shell_condition",
          "shell_condition_3M", "durometer", "cfv", "comment", "ll", "latitude", "longitude",
          "depth", "clawlength", "weight", "zone", "data", "digits", "observer",
          "vessel", "eggmaturity"
          )
    rawdata = read.fwf( file=temp.output, widths=w, as.is=T, colClasses="character", header=F, strip.white=F)
    colnames(rawdata) = v
    rawdata = rawdata[,keep]
    colnames(rawdata) = c( "sdate", "trap", "tripno", "crabno", "sex", "carapacewidth",
          "ndays", "trapno", "clawheightRight", "clawheightLeft", "shell_condition",
          "shell_condition_3M", "durometer", "cfv", "comment", "ll", "latitude", "longitude",
          "depth", "clawlength", "weight", "zone", "data", "digits", "observer",
          "vessel", "eggmaturity"
    )

    save( rawdata, file=file.path( paste( outfile, type, "rdata", sep=".") ) , compress=T)
    if (file.exists( temp.output ) )  file.remove(temp.output)

    return( "Done" )


    ## the following are not used and kept for reference

     #recode data
    recode.data =F
    if (recode.data) {
         # remove offending entries ...
    # most relate to the same entries for "FC157" in files: ECR02072.DE, ECR02073.DE, ECR02091.DE, ECR02117.DE
    baddata = grep("FC157",rawdata$sdate )   # 4 entries with no data
    baddata = c(baddata, which( rawdata$zone %in% c("2", "4", "5") ))
    baddata = c(baddata, grep("[438]", rawdata$subzone) )
    baddata = c(baddata, grep("[ ]", rawdata$datatype) )
    baddata = c(baddata, grep("[ ]", rawdata$trawlno) )
    baddata = c(baddata, grep("[ ]", rawdata$crabno) )

    baddata = c(baddata, grep("[ ]", rawdata$sex) )
    # one case date 17072001  zone24E   trawl009   crab0001 fileECR01305.txt

    baddata = unique(baddata)
    rawdata = rawdata[-baddata ,]



      rawdata$datatype[ which(rawdata$datatype=="3")] = "trawl"
      i.m = which(rawdata$sex=="1")
      i.f = which(rawdata$sex=="2")
      i.o = setdiff( 1:nrow(rawdata), c(i.m, i.f) )
      rawdata$sex[ i.m ] = male
      rawdata$sex[ i.f ] = female
      rawdata$sex[ i.o ] = sex.other

      i.mat = which(rawdata$maturity=="1")
      i.imm = which(rawdata$maturity=="2")
      i.other = setdiff( 1:nrow(rawdata), c(i.mat, i.imm) )

      rawdata$maturity[ i.mat ] = mature
      rawdata$maturity[ i.imm ] = immature
      rawdata$maturity[ i.other ] = mat.other

      rawdata$shell_condition[ which(rawdata$shell_condition %in% c("*", " ", "0")) ] = NA
      rawdata$shell_condition_3M = ifelse(rawdata$shell_condition_3M == "M" , T, F )
      rawdata$gonad_colour = gsub("[ ]", "", rawdata$gonad_colour)
      rawdata$gonad_colour[ which( rawdata$gonad_colour %in% c("**", " ", "0")) ] = NA  # (1=white, 2=beige, 3=orange)
      rawdata$egg_colour = gsub("[*]", "", rawdata$egg_colour)
      rawdata$egg_colour[ which( rawdata$egg_colour %in% c("**", " ", "0", "5")) ] = NA  # (1=light orange, 2=dark orange, 3=black, 4=coccoon?)
      rawdata$Pr_eggs_remaining[ which( rawdata$Pr_eggs_remaining %in% c("*", " ", "5")) ] = NA # (0=0, 1=1-49%, 2=50-74%, 3=75-99%, 4=100% )
      rawdata$tagno[ which(rawdata$tagno %in% c("00000000", "        ")) ] = NA # there are no tag data enterred into raw tables ... where are they?
      rawdata$legs[ which(rawdata$legs %in% c("          ")) ] = NA
      rawdata$position_type[ which(rawdata$position_type == "  ") ] = NA
      rawdata$loc0[ which(rawdata$loc0=="99999999") ] = NA  # explain ?
      rawdata$loc1[ which(rawdata$loc1=="99999999") ] = NA  # explain ?
      rawdata$depth[ which(rawdata$depth=="***") ] = NA  # in fathoms ?
      rawdata$soak_days[ which(rawdata$soak_days %in% c(" ", "*")) ] = NA  # no valid data in raw files ... unused?
      rawdata$durometer = gsub("[ ]", "", rawdata$durometer)
      rawdata$durometer[ which(rawdata$durometer %in% c("***", "")) ] = NA
      rawdata$trap_code = gsub("[* ]", "", rawdata$trap_code)
      rawdata$samplers =  gsub("[* ]", "",rawdata$samplers)
      rawdata$weight[ which(rawdata$weight=="*****") ] = NA  # no valid data ...
      rawdata$comments =  gsub("[*]", "", rawdata$comments)

      # additional variables
      rawdata$julian = chron( 
        dates.=paste( substr(rawdata$sdate,1,2), "-", substr(rawdata$sdate,3,4), "-", substr(rawdata$sdate,5,8), sep=""),
        format=c(dates="d-m-y"), out.format=c(dates="year-m-d") )

      numerics = c( "carapace_width", "female_abdomen", "chela_height", "loc0", "loc1",
                  "depth", "soak_days", "durometer", "weight")
      for(i in numerics) rawdata[,i] = as.numeric(rawdata[,i])
    }

    return (odb)

  }


