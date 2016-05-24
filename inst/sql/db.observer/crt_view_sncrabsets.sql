REM
REM     File: 		crt_view_sncrabsets.sql
REM 	Creation Date: 	Aug 02, 2004
REM	Author:		Kohila Thana
REM	Purpose:	This script creates the view sncrabsets 
REM			The view includes the crab details  
REM	
REM			   
REM	Parameters: 	None
REM
REM	Path:		
REM	
REM	Called by: 	$ORACLE_HOME/sqlplus		
REM	
REM


CREATE OR REPLACE view sncrabsets AS
SELECT trip.trip_id, trip.trip, tt.trip_type, vess.vessel_name, trip.license_no,
trip.board_date, trip.landing_date, prof.settime start_time, obs.observer,
 gr.gearcd_id, st.haulccd_id, 
st.set_no,   st.comarea_id, st.nafarea_id,  
prof.depth, prof.water_temperature,
prof.latitude start_lat, prof.longitude start_long, 
(SELECT latitude 
 FROM issetprofile prof
 WHERE prof.fishset_id = st.fishset_id
 AND prof.pntcd_id = 3) end_lat, 
(SELECT longitude 
 FROM issetprofile prof
 WHERE prof.fishset_id = st.fishset_id
 AND prof.pntcd_id = 3) end_long, sn.distance_gps, sn.area_swept, 
 st.est_catch, st.comments
FROM istrips trip, sntows sn, istriptypecodes tt, isvessels vess, isobservercodes obs,
    isgears gr, isfishsets st, issetprofile prof 
WHERE trip.tripcd_id = 7061
AND trip.trip_id = gr.trip_Id
AND trip.tripcd_id = tt.tripcd_id
AND trip.vess_id = vess.vess_id(+)
AND trip.obscd_id = obs.obscd_id
AND (trip.trip_id = st.trip_Id
    AND gr.gear_id = st.gear_id)
AND (st.fishset_id = prof.fishset_id
     AND prof.pntcd_id = 2)
AND (trip.trip = sn.trip
AND sn.set_number = st.set_no)
order by board_date, set_no;



