REM
REM     File: 		crt_view_sncrabdetails.sql
REM 	Creation Date: 	Aug 02, 2004
REM	Author:		Kohila Thana
REM	Purpose:	This script creates the view sncrabdetails 
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

CREATE OR REPLACE view sncrabdetails AS
SELECT trip.trip_id, trip.trip,
trip.board_date,
st.set_no,  
ca.est_num_caught,  ca.est_discard_wt, 
fish.fish_no , fish.sexcd_id, fish.fish_length,
fish.fish_weight measured_wgt,
DECODE(fish.fish_weight, NULL, ROUND(0.0002665 * POWER(fish_length, 3.098), 3)) calc_wgt,
SUM(DECODE(mrph.mrphcd_id, 54, mrph.quant_value, NULL)) female_abdomen,
SUM(DECODE(mrph.mrphcd_id, 37, mrph.quant_value, NULL)) chela_height,
SUM(DECODE(mrph.mrphcd_id, 56, mrph.mrphvcd_id, NULL)) maturity_cd,
SUM(DECODE(mrph.mrphcd_id, 52, mrph.mrphvcd_id, NULL)) shellcond_cd,
SUM(DECODE(mrph.mrphcd_id, 55, mrph.mrphvcd_id, NULL)) gonade_cd,
SUM(DECODE(mrph.mrphcd_id, 53, mrph.mrphvcd_id, NULL)) eggcolor_cd,
SUM(DECODE(mrph.mrphcd_id, 50, mrph.mrphvcd_id, NULL)) eggpercent_cd,
SUM(DECODE(mrph.mrphcd_id, 51, mrph.quant_value, NULL)) durometre,
NVL(TO_CHAR(SUM(DECODE(mrph.mrphcd_id, 61, mrph.mrphvcd_id, NULL))), 
DECODE(fish.fish_no, NULL, NULL,'*'))||
NVL(TO_CHAR(SUM(DECODE(mrph.mrphcd_id, 62, mrph.mrphvcd_id, NULL))), 
DECODE(fish.fish_no, NULL, NULL,'*'))||
NVL(TO_CHAR(SUM(DECODE(mrph.mrphcd_id, 63, mrph.mrphvcd_id, NULL))), 
DECODE(fish.fish_no, NULL, NULL,'*'))||
NVL(TO_CHAR(SUM(DECODE(mrph.mrphcd_id, 64, mrph.mrphvcd_id, NULL))), 
DECODE(fish.fish_no, NULL, NULL,'*'))||
NVL(TO_CHAR(SUM(DECODE(mrph.mrphcd_id, 65, mrph.mrphvcd_id, NULL))), 
DECODE(fish.fish_no, NULL, NULL,'*'))||
NVL(TO_CHAR(SUM(DECODE(mrph.mrphcd_id, 71, mrph.mrphvcd_id, NULL))), 
DECODE(fish.fish_no, NULL, NULL,'*'))||
NVL(TO_CHAR(SUM(DECODE(mrph.mrphcd_id, 72, mrph.mrphvcd_id, NULL))), 
DECODE(fish.fish_no, NULL, NULL,'*'))||
NVL(TO_CHAR(SUM(DECODE(mrph.mrphcd_id, 73, mrph.mrphvcd_id, NULL))), 
DECODE(fish.fish_no, NULL, NULL,'*'))||
NVL(TO_CHAR(SUM(DECODE(mrph.mrphcd_id, 74, mrph.mrphvcd_id, NULL))), 
DECODE(fish.fish_no, NULL, NULL,'*'))||
NVL(TO_CHAR(SUM(DECODE(mrph.mrphcd_id, 75, mrph.mrphvcd_id, NULL))), 
DECODE(fish.fish_no, NULL, NULL,'*')) missing_legs
FROM istrips trip, isgears gr, isfishsets st, 
      iscatches ca, isfish fish,
      isfishmorphs mrph
WHERE trip.tripcd_id = 7061
AND trip.trip_id = gr.trip_Id
AND (trip.trip_id = st.trip_Id
    AND gr.gear_id = st.gear_id)
AND st.fishset_id = ca.fishset_id(+)
AND ca.speccd_id(+) = 2526
AND ca.catch_id = fish.catch_id(+)
AND fish.fish_id = mrph.fish_id(+)
GROUP BY trip.trip_id, trip.trip, 
trip.board_date, 
st.set_no , st.specscd_id, 
ca.est_num_caught, ca.est_kept_wt, ca.est_discard_wt, ca.est_reduction_wt,
ca.est_combined_wt, 
fish.fish_no , fish.sexcd_id, fish.fish_length, fish.fish_weight
order by board_date, set_no, fish_no;
