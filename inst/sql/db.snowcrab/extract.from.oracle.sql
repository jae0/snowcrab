
set pagesize 10000
set arraysize 5000
set linesize 500
set colsep ;;
set pause off
set tab off
set term off
set flush off
spool datadump.dat
select * from scs_main;

scs_main
scs_area
scs_count







# merging in oracle






set pagesize 50000
set arraysize 5000
set linesize 500
set colsep ,
set pause off
set tab off
set term off
set flush off
spool /users/choij/q.out

select  c.year, c.series, c.sdate, c.strat, c.slong, c.slat, c.totwgt, c.sampwgt, c.dist,
        d.mission, d.setno,
        d.spec species,
        d.flen length,
        d.fwt mass,
	d.age, 
	d.fmat mat,
	d.fsex sex
from    groundfish.gsdet  d,
        ( select i.year, i.series, i.sdate, i.strat, i.type, i.slong, i.slat,
                 cat.mission, cat.setno, cat.spec,
                 cat.totwgt, cat.sampwgt, i.dist
          from
                groundfish.gscat cat,
                ( select mission.year, mission.fk_series_id series, inf.setno,
                         inf.mission, inf.sdate, inf.strat, inf.dist, inf.type,
                         inf.slong, inf.slat
                  from  groundfish.gsmission_list mission,
                        groundfish.gsinf inf
                  where inf.mission(+) = mission.pk_mission
                ) i
          where
                cat.mission(+) = i.mission
          and   cat.setno(+) = i.setno
        ) c
where   d.mission(+)=c.mission
and     d.setno(+)=c.setno
and     d.spec(+)=c.spec
and     type=1
and     year = 1970
and     d.spec = 400;

