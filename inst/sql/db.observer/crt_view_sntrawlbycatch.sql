REM
REM     File: 		crt_view_sntrawlbycatch.sql
REM 	Creation Date: 	Aug 02, 2004
REM	Author:		Kohila Thana
REM	Purpose:	This script creates the view sntrawlbycatch 
REM			The view includes the by catch info 
REM	
REM			   
REM	Parameters: 	None
REM
REM	Path:		
REM	
REM	Called by: 	$ORACLE_HOME/sqlplus		
REM	
REM

CREATE OR REPLACE view sntrawlbycatch AS
SELECT trip.trip_id, trip.trip, 
trip.board_date, 
st.set_no,  st.est_catch, ca.speccd_id, ca.est_num_caught, ca.est_discard_wt
FROM istrips trip, isgears gr, isfishsets st, iscatches ca
WHERE trip.tripcd_id = 7061
AND trip.trip_id = gr.trip_Id
AND (trip.trip_id = st.trip_Id
    AND gr.gear_id = st.gear_id)
AND (st.fishset_id = ca.fishset_id
     AND ca.speccd_id !=2526)
order by board_date, set_no, speccd_id;



