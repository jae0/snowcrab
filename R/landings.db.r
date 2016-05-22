  landings.db = function ( ) {

    # glue historical data with marfis data 
    
    threshold.year = 2004
    
    b = logbook.db( DS="logbook" ) # modern data:  mass in kg  -- must use all data and not postionally filtered data to get accurate totals
    b = b[b$yr >= threshold.year, ] 
    
    a = read.table(file.path(  project.datadirectory("snowcrab"), "data", "fisheries", "annual.landings.csv"), sep=",", as.is=T)
    names(a) = c("yr","cfa","nlicences","tac.tons","landings.tons","cpue.kg.trap","effort.100traps")
    numerics = c("yr","nlicences","tac.tons","landings.tons","cpue.kg.trap","effort.100traps")
    a = factor2number(a, numerics)
    a = a[a$yr < threshold.year,]
    
    a$year = a$yr
    a$lat = NA
    a$lon = NA
    a$depth = NA
    a$landings = a$landings.tons * 1000
    a$effort = a$effort.100traps * 100
    a$soak.time = NA
    a$cpue = a$cpue.kg.trap
    a$trap.type = NA
    a$cfv = NA
    a$status = NA
    a$licence = NA
    a$date.landed = NA
    a$date.fished = NA
    a$cfa0 = NA
    a$area_id = NA

    a$subarea = a$cfa
    a$cfa = NA
    a$cfa[which(a$subarea %in% c("cfa20", "cfa21", "cfa22", "cfanorth") )] = "cfanorth"
    a$cfa[which(a$subarea %in% c("cfa23", "cfa24", "cfasouth", "cfaslope") )] = "cfasouth"
    a$cfa[which(a$subarea %in% c("cfa4x") )] = "cfa4x"
    
    a$no.lic[ which(a$subarea=="cfaslope" & a$yr==2001) ] = 4
    a$no.lic[ which(a$subarea=="cfaslope" & a$yr==2002) ] = 4
    a$no.lic[ which(a$subarea=="cfaslope" & a$yr==2003) ] = 5

    d = a[,names(b)]

    c = rbind(d, b)

    return(c)
  }


