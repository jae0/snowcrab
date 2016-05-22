
  logbook.determine.region = function( logs ) {
      
      lic = logbook.db( DS="odbc.licence" )
      names(lic) = tolower( names(lic) )
     
      areas = logbook.db( DS="odbc.areas" )
      names(areas) = tolower( names(areas) )
      areas = areas[ , c("area_id", "area", "area_type_id", "desc_eng" ) ]  # reduce size
      areas$desc_eng = as.character(areas$desc_eng ) 

      snowcrab.area = c(27, 304, 305, 306, 307, 308, 615, 616, 619, 620, 790, 791, 792, 793, 794, 795, 796, 921, 922, 923, 1058, 1059  )
      
      areas = areas[ which( areas$area_id %in% snowcrab.area) , ] 
      lic = lic[ which( lic$area_id %in% snowcrab.area ) , ]  # reduce size to crab areas only
      
      lic = merge( lic, areas, by="area_id", all.x=T, all.y=F )
      lic$area = toupper( as.character( lic$area) )

# data dump from marfissci.areas table (2011)
# from 
# areas = areas[ grep ("crab", areas$desc_eng, ignore.case=T ) ,   ]  # reduce size
# AREA_ID   AREA AREA_TYPE_ID                      DESC_ENG                                DESC_FRE PARENT_AREA_ID ACT_FLAG      CUSER
#    304     20            9        CRAB FISHING AREA - 20  ZONE DE P?CHE DU CRABE DES NEIGES - 20             NA        Y CONVERSION
#    305     21            9        CRAB FISHING AREA - 21  ZONE DE P?CHE DU CRABE DES NEIGES - 21             NA        Y CONVERSION
#    306     22            9        CRAB FISHING AREA - 22  ZONE DE P?CHE DU CRABE DES NEIGES - 22             NA        Y CONVERSION
#    615    22I            9 CRAB FISHING AREA - 22I INNER ZONE DE P?CHE DU CRABE - 22I INT?RIEUR>            306        Y     MARFIS
#    619    22O            9 CRAB FISHING AREA - 22O OUTER ZONE DE P?CHE DU CRABE - 22O INT?RIEUR>            306        Y     MARFIS

#    307     23            9        CRAB FISHING AREA - 23  ZONE DE P?CHE DU CRABE DES NEIGES - 23             NA        Y CONVERSION
#    790    23B            9         CRAB FISHING AREA 23B ZONE DE P?CHE DU CRABE DES NEIGES - 23B             NA        Y   SCHLEITC
#     791    23C            9         CRAB FISHING AREA 23C ZONE DE P?CHE DU CRABE DES NEIGES - 23C             NA        Y   SCHLEITC
#    1058    23S            9       CRAB FISHING AREA - 23S            ZONE DE P?CHE DU CRABE - 23S             NA        Y     MARFIS
#     792    23D            9         CRAB FISHING AREA 23D ZONE DE P?CHE DU CRABE DES NEIGES - 23D             NA        Y   SCHLEITC

#     308     24            9        CRAB FISHING AREA - 24  ZONE DE P?CHE DU CRABE DES NEIGES - 24             NA        Y CONVERSION
#     616    24E            9  CRAB FISHING AREA - 24E EAST        ZONE DE P?CHE DU CRABE - 24E EST            308        Y     MARFIS
#     793    24B            9         CRAB FISHING AREA 24B ZONE DE P?CHE DU CRABE DES NEIGES - 24B             NA        Y   SCHLEITC
#     794    24C            9       CRAB FISHING AREA - 24C ZONE DE P?CHE DU CRABE DES NEIGES - 24C             NA        Y   SCHLEITC
#     795    24D            9         CRAB FISHING AREA 24D ZONE DE P?CHE DU CRABE DES NEIGES - 24D             NA        Y   SCHLEITC
#    1059    24S            9       CRAB FISHING AREA - 24S            ZONE DE P?CHE DU CRABE - 24S             NA        Y     MARFIS

#      27     4X            1            NAFO DIVISION - 4X             DIVISION DE L?OPANO - 4X             NA
#     620    24W            9  CRAB FISHING AREA - 24W WEST      ZONE DE P?CHE DU CRABE - 24W OUEST            308        Y     MARFIS
#     796    24H            9         CRAB FISHING AREA 24H              SNOW CRAB FISHING AREA 24H             NA        Y   SCHLEITC


#     922  CFA24            0 OBSERVER AREA-SNOW CRAB CFA24           OBSERVER AREA-SNOW CRAB CFA24             NA        Y     MARFIS
#     921  CFA23            0 OBSERVER AREA-SNOW CRAB cfa23           OBSERVER AREA-SNOW CRAB CFA23             NA        Y     MARFIS
#     923 SURVEY            0   OBSERVER AREA-SNOW CRAB SUR             OBSERVER AREA-SNOW CRAB SUR             NA        Y     MARFIS
    


      # licence information
      
      # determine 4x lic id's
      
      lic_CFA4X = unique( lic$licence_id[ which( lic$area %in% c( "4X", "24W" ) ) ])  # known to be in CFA 4X
      lic_CFA24 = unique( lic$licence_id[ which( lic$area %in% c( "24A", "24B", "24C", "24D", "24E", "24S","24H","CFA24" ) ) ]) # known to be in CFA 24

      north = which( lic$area %in% c( "20", "21" ,"22", "22I", "22O" )   )  # use upper case
      cfa23 = which( lic$area %in% c( "23", "23A", "23B", "23C", "23D", "23S","CFA23") )
      
      # "area=24" contains both CFA24 and CFA4X .. by default consider all to be part of CFA24
      #  and then recode after the fact those that belong to 4X---Feb 2015 add in the area 308 ie '24' to cfa24 only this was the big difference between ben's and my data
      #cfa24 = which( lic$area %in% c( "24A", "24B", "24C", "24D", "24E", "24S" ,'24') | ( (lic$area=="24") & ( lic$licence_id %in% lic_CFA24 ) ) ) 
      #cfa4x = which( lic$area %in% c( "4X", "24W" ) | ( (lic$area=="24") & ( lic$licence_id %in% lic_CFA4X ) 
    cfa24 = which( lic$area %in% c( "24A", "24B", "24C", "24D", "24E", "24S" ,'24') ) 
     cfa4x = which( lic$area %in% c( "4X", "24W" )) 
      lic$cfa0 = NA
      lic$cfa0 [north] = "cfanorth"
      lic$cfa0 [cfa23] = "cfa23"
      lic$cfa0 [cfa24] = "cfa24"
      lic$cfa0 [cfa4x] = "cfa4x"

      lic0 = lic[ which( !is.na(lic$cfa0)) , ]
      lic0 = lic0[ - which(duplicated( lic0$licence_id) ) ,]  # finalised licencing information
      lic0$licence = as.character(lic0$licence_id) 
      lic0 = lic0[, c("licence", "area_id", "cfa0", "area", "desc_eng") ]
      #il = which(lic0$licence =='100198')
      #lic0[il,c('cfa0','area','desc_eng')] <- c('cfa4x','4X','NAFO DIVISION - 4X')
      QQ = merge(logs, lic0, by="licence", all.x=T, all.y=F, sort=F)
     

      debug = F
      if (debug) {
        
        
        plot( QQ$lon, QQ$lat )
        ss = which( QQ$area_id ==796 )  # area "24H" is in area 24
        s4x = which( QQ$cfa0 =="cfa4x" )  # area "24" spans 4X and area 24
        s23 = which( QQ$cfa0 =="cfa23" )  # area "24" spans 4X and area 24
        s24 = which( QQ$cfa0 =="cfa24" )  # area "24" spans 4X and area 24
        ssn = which( QQ$cfa0 =="cfanorth" )  # area "24" spans 4X and area 24
        points( QQ$lon[s4x], QQ$lat[s4x], pch=16, col="blue" )
        points( QQ$lon[s23], QQ$lat[s23], pch=16, col="green" )
        points( QQ$lon[s24], QQ$lat[s24], pch=16, col="red" )
        points( QQ$lon[ssn], QQ$lat[ssn], pch=16, col="yellow" )

        unique( QQ$licence[ s4x] )  
        length(unique( QQ$licence[ ssn] ))
        length(unique( QQ$licence[ s23] ))
        length(unique( QQ$licence[ s24] ))
        length(unique( QQ$licence[ s4x] ))

      }

 
      i.missing = which( is.na( QQ$cfa0 ) )
      if (length( i.missing) > 1) {
      # try to determine via geographics:
        FF = QQ[ i.missing , c("lon", "lat")]
      # cfa 4X has a fishing season that spans two years recode "yr" to accomodate this
        icfa4x = filter.region.polygon(FF, "cfa4x")
        icfanorth = filter.region.polygon(FF, "cfanorth")
        icfa23 = filter.region.polygon(FF, "cfa23")
        icfa24 = filter.region.polygon(FF, "cfa24")
        
        G = rep( NA, length(i.missing) )

        G[icfa4x] = "cfa4x"
        G[icfanorth] = "cfanorth"
        G[icfa23] = "cfa23"
        G[icfa24] = "cfa24"

        QQ$cfa0 [ i.missing ] = G
    }

    # overwrite cfa with internally consistent areas
    QQ$cfa = QQ$cfa0
    QQ$cfa[ which( QQ$cfa0 %in% c("cfa23", "cfa24"))] = "cfasouth"

    # remove vars no more needed (for debugging and verification)
    QQ$area_id = NULL
    QQ$area = NULL
    QQ$area_type_id
    QQ$desc_eng = NULL

    return ( QQ )
  }

 

