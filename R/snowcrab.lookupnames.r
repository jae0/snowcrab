
  snowcrab.lookupnames = function() {
    
    lookuptable = c(
    "longitude",    "lon", 
    "latitude",     "lat", 
    "minilog_file", "minilog", 
    "comarea_id",   "cfa", 
    "vr_number", "cfv", 
    "vessel_name", "vessel", 
    "captain_name", "captain", 
    "licence_id", "licence", 
    "date_landed", "date.landed", 
    "date_fished", "date.fished", 
    "num_of_traps", "effort", 

    "set_no", "set", 
    "nafarea_id", "nafo", 
    "depth", "Zx",                     # m
    "water_temperature", "Tx", 
    "distance_gps", "dist", 
    "area_swept", "sa",                # m2
    "fish_no", "crabno", 
    "sexcd_id", "sex", 
    "fish_length", "cw", 
    "measured_wgt", "mass", 
    "female_abdomen", "abdomen", 
    "chela_height", "chela", 
    "maturity_cd", "mat", 
    "shellcond_cd", "shell", 
    "gonade_cd", "gonad", 
    "eggcolor_cd", "eggcol", 
    "eggpercent_cd", "eggPr", 
    "durometre", "durometer", 
    "missing_legs", "legs", 
    "speccd_id", "spec", 
    "est_catch", "totmass", 
    "soak_time", "soak.time", 
    "start_settime", "time0", 
    "start_setdate", "date0", 
    "end_settime", "time1", 
    "end_setdate", "date1", 
    "gearcd_id", "gear", 

    "netmind_file", "netmind", 
    "haulccd_id", "towquality", 
    "est_num_caught", "totno", 
    "est_discard_wt", "totmass.discarded",       # kg
    "tow_date", "date", 
    "start_long", "lon", 
    "start_lat", "lat", 
    "end_lat", "lat1",
    "end_long", "lon1",
    "start_time", "stime", 
    "board_date", "sdate",
	"setcd_id","set_type") 

    lookuptable = matrix( lookuptable, byrow=T, ncol=2 )

    return(lookuptable)
 }


