REM 	Creation Date: 	Aug 02, 2004
REM	Author:		Kohila Thana
REM	Purpose:	observer bycatches: data requested, modified from: crt_view_sncrabdetails.sql
REM			   
REM	Parameters: 	snowcrab fishery == 2526

SELECT trip.trip_id, trip.trip, trip.board_date, st.set_no, 
       ca.speccd_id, ca.est_num_caught,  ca.est_kept_wt, ca.est_discard_wt, 
       fish.fish_no , fish.fish_length, fish.fish_weight 
FROM istrips trip, isgears gr, isfishsets st, iscatches ca, isfish fish
WHERE trip.tripcd_id = 2509
  AND trip.trip_id = gr.trip_Id
  AND (trip.trip_id = st.trip_Id AND gr.gear_id = st.gear_id)
  AND st.fishset_id = ca.fishset_id(+)
  AND ca.speccd_id(+) = 204
  AND ca.catch_id = fish.catch_id(+)
order by board_date, set_no, fish_no;


