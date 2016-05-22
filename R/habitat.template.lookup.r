

	habitat.template.lookup = function( V ) {

		out = V
    
		if ( grepl( "R0", V )   | grepl( "male.com", V ) )  out = "R0.mass"

    if ( grepl( "male.large", V ) |
         grepl( "ma11", V ) | grepl( "ma12", V ) | grepl( "ma13", V )     | 
         grepl( "mi11", V ) | grepl( "mi12", V ) | grepl( "mi13", V )     | 
         grepl( "m11", V )  | grepl( "m12", V )  | grepl( "m13", V )       | 
         grepl( "R1", V )   | grepl( "R2", V ) )  out = "male.large.mass"
			
    if ( grepl( "male.small", V )    |
         grepl( "R3", V )   | grepl( "R4", V )   | grepl( "R5", V )       | 
         grepl( "ma9", V )  | grepl( "ma10", V ) | 
         grepl( "mi6", V )  | grepl( "mi7", V )  | grepl( "mi7", V )      | 
         grepl( "mi8", V )  | grepl( "mi9", V )  | grepl( "mi10", V )     | 
         grepl( "m6", V )   | grepl( "m7", V )   | grepl( "m8", V )       | 
         grepl( "m9", V )   | grepl( "m10", V )   )  out = "male.large.mass"

    if ( grepl( "female.large", V )    |
         grepl( ".female.prim", V )   | grepl( ".female.multi", V )   | 
         grepl( ".female.berried", V )   | grepl( ".female.mat", V )       | 
         grepl( "fa8", V )  | grepl( "fa9", V )  | grepl( "fa10", V )     |
         grepl( "fi8", V )  | grepl( "fi9", V )  | grepl( "fi10", V )      | 
         grepl( "f8", V )   | grepl( "f9", V )   | grepl( "f10", V )       | 
         grepl( "f9", V )   | grepl( "f10", V )   )  out = "female.large.mass"

    if ( grepl( "female.small", V )    |
         grepl( "R3", V )   | grepl( ".female.small", V )   | grepl( ".female.imm", V )       | 
         grepl( "fa9", V )  | grepl( "fa10", V ) | 
         grepl( "fi6", V )  | grepl( "fi7", V )  | grepl( "fi7", V )      | 
         grepl( "fi8", V )  | grepl( "fi9", V )  | grepl( "fi10", V )     | 
         grepl( "f6", V )   | grepl( "f7", V )   | grepl( "f8", V )       | 
         grepl( "f9", V )   | grepl( "f10", V )   )  out = "female.small.mass"

		return (out )
	}



