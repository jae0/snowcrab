CREATE or REPLACE VIEW MARFIS_CRAB2 AS 
(	SELECT 
	-- any mon docs with a crab slip record, with any associated effort.  slip weight is prorated based on estimated log weights if exists else overall count of effort recs
		'MD-' || slips.mon_doc_id doc_id, --mon doc id 
		slips.vr_number,  --CFV 
		(select vessel_name from vessels where vr_number = slips.vr_number) vessel_name, --VESSEL NAME 
		to_number(to_char(slips.landing_date_time,'yyyy')) year,
		slips.captain_name, -- CAPTAIN NAME 
		slips.LFA, --lobster fishing area 
		nvl(logs.CFA,slips.CFA) CFA,
		logs.licence_id, --LICENCE 
		(select licence_id_old from licences_sf_old where licence_id = logs.licence_id) licence_id_old,
		slips.community_code, --PORT 
		slips.landing_date_time  date_landed, --DATE LANDED 
		slips.ssf_species_code species_code, --SPECIES 
		round( slips.slip_weight_lbs, 0) slip_weight_lbs,  -- LIVE WEIGHT 
		decode( /* IF  estimated log weights is 0*/
			(select nvl(sum(ls.weight),0) from log_spc_std_info ls, log_efrt_std_info le 
											where ls.ssf_species_code = slips.ssf_species_code and slips.mon_doc_id = le.mon_doc_id 
											and le.log_efrt_std_info_id = ls.log_efrt_std_info_id) ,0,
			/*THEN  there are no estimated log weights - have to use a count*/
				slips.slip_weight_lbs / (select decode(count(le.log_efrt_std_info_id),0,1, null,1,count(le.log_efrt_std_info_id))
											from log_spc_std_info ls, log_efrt_std_info le 
											where  slips.mon_doc_id = le.mon_doc_id 
												and le.log_efrt_std_info_id = ls.log_efrt_std_info_id(+) 
												and nvl(ls.ssf_species_code, slips.ssf_species_code) = slips.ssf_species_code) ,
			 /* ELSE pro-rate based on log estimated weights*/
			(logs.est_weight_log_lbs/
				(select sum(decode(ls.unit_of_measure_id,10,ls.weight*2.2046,ls.weight) )
				 from log_spc_std_info ls, log_efrt_std_info le 
				 where ls.ssf_species_code = slips.ssf_species_code and slips.mon_doc_id = le.mon_doc_id and le.log_efrt_std_info_id = ls.log_efrt_std_info_id)
			) * 	 slips.slip_weight_lbs
		) pro_rated_slip_wt_lbs,
		logs.log_efrt_std_info_id, --shouldn't need to include this since the UNION ALL should not eliminate log records that look the same, but....
		logs.fv_fished_datetime  date_fished, -- DATE FISHED 
		logs.est_weight_log_lbs,--ESTIMATE-QTY-LBS 
		logs.ent_latitude latitude, --LATITUDE 
		logs.ent_longitude longitude, --LONGITUDE 
		to_number(logs.depth_fm) depth_fm,--DEPTH 
		to_number(logs.num_of_traps) num_of_traps,--NBR-TRAPS 
		to_number(logs.soak_days) soak_days, --NBR-HOURS -soak days 
		logs.catch_usage_code, --KEPT/DISCARD 
		logs.target_spc, -- license target species 
		'MD-' || slips.mon_doc_defn_id doc_type--  mon doc type 
	FROM 
		(-- slip info -all crab slips 
				select 
						mon_doc_id, 
						vr_number, 
						captain_name, 
						LFA, CFA,
						community_code, 
						landing_date_time, 
						ssf_species_code, 
						sum(live_wt_kgs*2.2046 ) slip_weight_lbs, 
						-- use sum for case of multiple slip_spc_info recs for same species 
						mon_doc_defn_id 
				from( 
						select 
								md.mon_doc_id, md.vr_number, 
								(select data_value from mon_doc_entrd_dets where mon_doc_id = md.mon_doc_id and column_defn_id = 17) captain_name, 
								(select data_value from mon_doc_entrd_dets where mon_doc_id = md.mon_doc_id and column_defn_id = 147) LFA, 
								(select area from areas where area_id = md.fv_fishing_area_id and area_type_id = 9 /* crab fishing area type */) CFA,
								so.community_code, 
								so.landing_date_time, 
								ss.slip_spc_std_info_id, 
								(nvl(ss.weight,0) / nvl(   (select factor from unit_of_measures where unit_of_measure_id = ss.unit_of_measure_id)   ,1))  * 
											(select form_factor from ssf_factors where ss.ssf_species_code = ssf_species_code 
											and nvl(ss.ssf_species_size_code,9) = ssf_species_size_code 
											and nvl(ss.ssf_landed_form_code,1) = ssf_landed_form_code 
								) live_wt_kgs,--account for uom and ssf_factor for round weight in kgs 
								ss.ssf_species_code, 
								md.mon_doc_defn_id 
						from 
								mon_docs md, slip_offld_std_info so, slip_buyr_std_info sb, slip_spc_std_info ss 
						where 
								 md.mon_doc_id = so.mon_doc_id 
								and so.slip_offld_std_info_id = sb.slip_offld_std_info_id 
								and sb.slip_buyr_std_info_id = ss.slip_buyr_std_info_id 
								and (ss.ssf_species_code = 704 /*rock*/or ss.ssf_species_code = 703 /*jonah*/or 
								ss.ssf_species_code = 708 /*stone*/or ss.ssf_species_code = 701 /*toad*/or ss.ssf_species_code = 705/*snow*/) 
				) 
				group by 
						mon_doc_id, 
						vr_number, 
						captain_name, 
						LFA, CFA,
						community_code, 
						landing_date_time, 
						ssf_species_code, 
						mon_doc_defn_id 
		)SLIPS, 
		(-- log info -all logs and effort associated with a mon_doc_id that had a crab slip.  If no catch for effort, the species is set to the licence targeted species and weight = 0 
				select 
						mds.mon_doc_id, mds.licence_id, 
						le.log_efrt_std_info_id, 
						le.fv_fished_datetime, -- DATE FISHED 
						nvl((select area from areas where area_id = le.fv_fishing_area_id and area_type_id = 9 /* crab fishing area type */),
									(select area from areas where area_id = le.det_fishing_area_id and area_type_id = 9) ) CFA,
						nvl (ls.ssf_species_code, mds.target_spc) species_code,  -- species caught 
						decode(ls.unit_of_measure_id,10,ls.weight*2.2046,ls.weight) est_weight_log_lbs,--ESTIMATE-QTY-LBS 
						le.ent_latitude, --LATITUDE 
						le.ent_longitude, --LONGITUDE 
						ls.catch_usage_code, --KEPT/DISCARD 
						mds.target_spc, 
						(select data_value from log_efrt_entrd_dets where log_efrt_std_info_id = le.log_efrt_std_info_id and column_defn_id = 158) depth_fm,--DEPTH 
						(select data_value from log_efrt_entrd_dets where log_efrt_std_info_id = le.log_efrt_std_info_id and column_defn_id in(702,112)) num_of_traps,--NBR-TRAPS 
						(select data_value from log_efrt_entrd_dets where log_efrt_std_info_id = le.log_efrt_std_info_id and column_defn_id in (703,111)) soak_days --NBR-HOURS -SOAK DAYS
				from 
						(	select distinct 
								mdLic.mon_doc_id, mdLic.licence_id, 
								(select species_code from licences where licence_id = mdLic.licence_id) target_spc
							from slip_spc_std_info ss, mon_doc_lics mdLic 
							where 
								(ss.ssf_species_code = 704 /*rock*/or ss.ssf_species_code = 703 /*jonah*/or 
											ss.ssf_species_code = 708 /*stone*/or ss.ssf_species_code = 701 /*toad*/or ss.ssf_species_code = 705/*snow*/) 
								and ss.mon_doc_lic_id = mdLic.mon_doc_lic_id 
						)mds, 
						log_efrt_std_info le, log_spc_std_info ls 
				where 
						mds.mon_doc_id = le.mon_doc_id(+) 
						and le.log_efrt_std_info_id = ls.log_efrt_std_info_id(+) 
		)LOGS
	WHERE 
		slips.mon_doc_id = logs.mon_doc_id 
		and slips.ssf_species_code = logs.species_code
UNION ALL
	SELECT 
		-- get the rest of the effort associated with mon_don_defn_id's = 10 and 19 (crab mon. doc. and offshore lobster mon. doc.) that did not have a crab slip record
		'MD-' || md.mon_doc_id doc_id, 
		md.vr_number,
		(select vessel_name from vessels where vr_number = md.vr_number) vessel_name,
		to_number(to_char(le.fv_fished_datetime,'yyyy')) year,
		(select data_value from mon_doc_entrd_dets where mon_doc_id = md.mon_doc_id and column_defn_id = 17) captain_name, 
		nvl((select data_value from mon_doc_entrd_dets where mon_doc_id = md.mon_doc_id and column_defn_id = 147),
				(select area from areas where area_id = md.fv_fishing_area_id and area_type_id = 5)
			) LFA, 
		nvl(
				(select area from areas where area_id = md.fv_fishing_area_id and area_type_id = 9 /* crab fishing area type */) ,
				(nvl(	(select area from areas where area_id = le.fv_fishing_area_id and area_type_id = 9 ),	(select area from areas where area_id = le.det_fishing_area_id and area_type_id = 9))	) 
			)CFA,
		to_char(null) /*licence_id*/, to_char(null) /*licence_id_old*/, to_number(null) /*community_code*/, to_date(null) /*date_landed*/, 
		to_number(null) /*species_code*/, to_number(null) /*slip_weight_lbs*/, to_number(null)/*pro_rated_slip_wt_lbs*/,
		le.log_efrt_std_info_id,
		le.fv_fished_datetime  date_fished,
		to_number(null) /*est_weight*/,
		le.ent_latitude latitude,
		le.ent_longitude longitude,
		(select to_number(data_value) from log_efrt_entrd_dets where log_efrt_std_info_id = le.log_efrt_std_info_id and column_defn_id = 158) depth_fm,--DEPTH 
		(select to_number(data_value) from log_efrt_entrd_dets where log_efrt_std_info_id = le.log_efrt_std_info_id and column_defn_id in(702,112)) num_of_traps,--NBR-TRAPS 
		(select to_number(data_value) from log_efrt_entrd_dets where log_efrt_std_info_id = le.log_efrt_std_info_id and column_defn_id in (703,111)) soak_days, --NBR-HOURS -SOAK DAYS
		to_number(null) /*catch_usage_code*/, to_number(null) /*target_spc*/,
		'MD-' || md.mon_doc_defn_id doc_type
	FROM mon_docs md, log_efrt_std_info le
	WHERE  (md.mon_doc_defn_id = 10 or md.mon_doc_defn_id = 19)
		and md.mon_doc_id = le.mon_doc_id
		and md.mon_doc_id not in	(	select distinct mdLic.mon_doc_id
										from slip_spc_std_info ss, mon_doc_lics mdLic 
										where (ss.ssf_species_code = 704 /*rock*/or ss.ssf_species_code = 703 /*jonah*/or 
											ss.ssf_species_code = 708 /*stone*/or ss.ssf_species_code = 701 /*toad*/or ss.ssf_species_code = 705/*snow*/) 
										and ss.mon_doc_lic_id = mdLic.mon_doc_lic_id 	)		
UNION ALL
	SELECT 
	-- all sum docs with a crab slip rec with any associated effort, slip weight prorated by number of effort records
			'SD-' || slips.sum_doc_id doc_id, -- sum doc id 
			slips.vr_number,  --CFV 
			slips.vessel_name, --VESSEL NAME 
			to_number(to_char(slips.date_landed,'yyyy')) year,
			to_char(null) captain_name, --CAPTAIN 
			(select area from areas where area_id = le.fv_ent_fishing_area_id and area_type_id = 5)  LFA,  -- area_type 5 is lobster fishing area 
			(select area from areas where area_id = le.fv_ent_fishing_area_id and area_type_id = 9 /* crab fishing area type */) CFA,
			slips.licence_id, --LICENCE 
			(select licence_id_old from licences_sf_old where licence_id = slips.licence_id) licence_id_old,
			slips.fv_community_code, --PORT 
			slips.date_landed, --DATE LANDED 
			slips.fv_ssf_species_code, --SPECIES 
			round(slips.live_wt * 2.2046, 0) slip_weight_lbs,-- LIVE WEIGHT 
		/*  there are no estimated log weights -  use a count*/
			((slips.live_wt * 2.2046)  / (select decode(count(le.sd_log_eff_std_info_id),0,1, null,1,count(le.sd_log_eff_std_info_id))
															from sd_log_eff_std_info le 
															where  slips.sum_doc_id = le.sum_doc_id) 
			) pro_rated_slip_wt_lbs,			
			le.sd_log_eff_std_info_id, 
			le.fv_date_fished, -- DATE FISHED 
			to_number(null) est_weight_log_lbs,--ESTIMATE-QTY-LBS 
			nvl(le.fv_ent_latitude,le.fv_det_latitude) latitude, --LATITUDE 
			nvl(le.fv_ent_longitude, le.fv_det_longitude) longitude, --LONGITUDE 
			to_number(null)  depth_fm,--DEPTH 
			decode(le.fv_gear_code, 62, le.fv_num_of_gear_units, null) num_of_traps, --gear_code 62 is TRAP 
			to_number(null)  soak_days, --NBR-HOURS -soak days 
			to_number(null) catch_usage_code, --KEPT/DISCARD 
			(select species_code from licences where licence_id = slips.licence_id) target_spc, -- license target species 
			'SD-' || slips.sum_doc_defn_id doc_type-- sum doc type 
	FROM 
		(-- SD slip info 
			select 
					sd.sum_doc_id, 
					sd.vr_number,  --CFV 
					sd.vessel_name, --VESSEL NAME 
					sd.submitr_fname || sd.submitr_sname captain_name, --CAPTAIN 
					sd.licence_id, --LICENCE 
					lnd.fv_community_code, --PORT 
					lnd.fv_landed_date date_landed, --DATE LANDED 
					spc.fv_ssf_species_code, --SPECIES 
					sum((nvl(spc.fv_weight,0) / nvl(   (select factor from unit_of_measures where unit_of_measure_id = spc.fv_unit_of_measure_id)   ,1))  * 
														(select form_factor from ssf_factors where spc.fv_ssf_species_code = ssf_species_code 
														and nvl(spc.fv_ssf_species_size_code,9) = ssf_species_size_code 
														and nvl(spc.fv_ssf_landed_form_code,1) = ssf_landed_form_code) 
					)live_wt, 
					sd.sum_doc_defn_id -- sum doc type 
			from 
					sum_docs sd, sd_slp_lnd_std_info lnd, sd_slp_byr_std_info byr, sd_slp_spc_std_info spc 
			where 
					(spc.fv_ssf_species_code = 704 /*rock*/or spc.fv_ssf_species_code = 703 /*jonah*/or 
								spc.fv_ssf_species_code = 708 /*stone*/or spc.fv_ssf_species_code = 701 /*toad*/or spc.fv_ssf_species_code = 705/*snow*/) 
					and sd.sum_doc_id = lnd.sum_doc_id 
					and lnd.sd_slp_lnd_std_info_id = byr.sd_slp_lnd_std_info_id 
					and byr.sd_slp_byr_std_info_id = spc.sd_slp_byr_std_info_id 
			group by 
					sd.sum_doc_id, 
					sd.vr_number,  --CFV 
					sd.vessel_name, --VESSEL NAME 
					sd.submitr_fname , sd.submitr_sname, --CAPTAIN 
					sd.licence_id, --LICENCE 
					lnd.fv_community_code, --PORT 
					lnd.fv_landed_date , --DATE LANDED 
					spc.fv_ssf_species_code, --SPECIES 
					sd.sum_doc_defn_id 
		)SLIPS, 
		sd_log_eff_std_info le 
	WHERE 
	slips.sum_doc_id = le.sum_doc_id(+) 
);

