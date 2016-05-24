

# convert dbf files to csv with ";" delimiters
dbf --separator "|" --csv main.tmp main.dbf

# seems to be a copy of main
# dbf --separator "|" --csv main.year.tmp main.year.dbf

# calculated values
# dbf --separator "|" --csv mean.carapace.width.tmp mean.carapace.width.dbf
# dbf --separator "|" --csv mean.weight.tmp mean.weight.dbf

dbf --separator "|" --csv minilog.temp.alldata.tmp minilog.temp.alldata.dbf
dbf --separator "|" --csv minilog.temp.tmp minilog.temp.dbf
dbf --separator "|" --csv netmid.surf.tmp netmid.surf.dbf 
dbf --separator "|" --csv species.tmp species.dbf

# remove the first line (the header) and quotations and dos carriage returns
sed -e '1d' -e 's/\"//g' -e 's/\r//g' main.tmp > main.csv

# sed -e '1d' -e 's/\"//g' -e 's/\r//g' main.year.tmp > main.year.csv
# sed -e '1d' -e 's/\"//g' -e 's/\r//g' mean.carapace.width.tmp > mean.carapace.width.csv
# sed -e '1d' -e 's/\"//g' -e 's/\r//g' mean.weight.tmp > mean.weight.csv

sed -e '1d' -e 's/\"//g' -e 's/\r//g' minilog.temp.alldata.tmp > minilog.temp.alldata.csv
sed -e '1d' -e 's/\"//g' -e 's/\r//g' minilog.temp.tmp > minilog.temp.csv
sed -e '1d' -e 's/\"//g' -e 's/\r//g' netmid.surf.tmp > netmid.surf.csv
sed -e '1d' -e 's/\"//g' -e 's/\r//g' species.tmp > species.csv


# --- create snowcrab database

mysql -u user -p
# pass

create database snowcrab;
use snowcrab;


# 1 -----------
drop table main;
CREATE TABLE main(
datecode date,
zone varchar(2),
sampling_t int,
trawl_number int,
crab_number int,
sex varchar(2),
cw numeric(7,3),
abdomen_wi varchar(5),
ch numeric(7,3),
maturity varchar(1),
sc varchar(1),
scm varchar(2),
gonad varchar(2),
femage_egg varchar(1),
percent_egg varchar(1),
stage int,
missinge_legs varchar(12),
durometer int,
code varchar(8),
samplers varchar(21),
weight varchar(15),
comments varchar(29),
maturity_c int,
zone_year varchar(11),
zone_date_year varchar(20),
weight_cal numeric(10,4)
);

load data infile "/home/jae/data/snowcrab/datadump/main.csv" into table main fields terminated by '|'; 
show warnings;



# 2 -----
drop table netmind;
CREATE TABLE netmind(
scanmarfile varchar(16),
setdate date,
starttime varchar(8),
endtime varchar(8),
distance numeric(10,4),
pre_valida numeric(10,4),
doorspread numeric(8,4),
dummy1 varchar(1),
latstart numeric(10,4),
lonstart numeric(10,4),
latend numeric(10,4),
lonend numeric(10,4),
commentsnetmind varchar(73),
dummy2 varchar(1),
sequence int,
recoreds int,
currentrec int,
year int,
midlat numeric(10,4),
midlon numeric(10,4),
datetext varchar(1),
zone varchar(3),
totalarea numeric(12,4),
avergae_di numeric(12,4),
average_wi numeric(12,4),
count_netm int,
average_ar numeric(12,4),
problem_wi1 varchar(3),
problem_wi2 varchar(11),
global_var numeric(12,4),
elasped_ti varchar(9),
uncertaint numeric(10,4),
did_not_pa varchar(6),
validated numeric(10,4),
trimmed_p1 varchar(3),
trimmed_p2 varchar(11)
);

load data infile "/home/jae/data/snowcrab/datadump/netmid.surf.csv" into table netmind fields terminated by '|' ;

show warnings;


# 3 ----
drop table scanmar;
CREATE TABLE scanmar(

year integer,
zone_date_tow varchar(20),
zone varchar(6),
setdate integer,
tow_number integer,
scanmar_file varchar(18),
tow_quality varchar(4),
sequence integer,
verify varchar(2),
crab_presence varchar(3),
total_crab_counts integer,
female_crab_counts integer,
male_crab_counts integer,
records integer,
curent_record integer,
bad_tow_indice integer,
count_by_zone integer,
depth numeric(5, 1),
bad_tow_status varchar(1),
temperature_logbook numeric(12, 6),
bad_records integer,
comment varchar(199),
dummy1 varchar(1), 
logbook_time varchar(8),
lat_start_logbook numeric(12, 6),
lon_start_log_book numeric(12, 6),
lat_end_logbook numeric(12, 6),
lon_end_logbook numeric(12, 6),
lat_mid_logbook numeric(12, 6),
lon_mid_logbook numeric(12, 6),
loranc_x_start numeric(12, 6),
loranc_x_end numeric(12, 6),
loranc_y_start numeric(12, 6),
loranc_y_end numeric(12, 6),
distance_logbook numeric(12, 6),
surface_logbook numeric(12, 6),
total_surface_logbook numeric(12, 6),
avg_width_logbook numeric(12, 6),
projected_logbook integer,
dummy2 varchar(1),
gps_start_time varchar(10),
gps_end_time varchar(10),
gps_elapsed_time varchar(10),
lat_start_gps numeric(11, 8),
lon_start_gps numeric(11, 8),
lat_end_gps numeric(11, 8),
lon_end_gps numeric(11, 8),
mid_lat_gps numeric(11, 8),
mid_lon_gps numeric(11, 8),
distance_gps numeric(11, 8),
scanmar_width numeric(11, 8),
scanmar_total_surface numeric(11, 8),
scanmar_area numeric(11, 8),
projected_gps integer,
dummy3 varchar(1),
lat_start_final numeric(16, 13),
lon_start_final numeric(16, 13),
lat_end_final numeric(16, 13),
lon_end_final numeric(16, 13),
lat_mid_boat_final numeric(16, 13),
lon_mid_boat_final numeric(16, 13),
dummy4 varchar(1),
alligator_count integer,
american_plaice_count integer,
cod_count integer,
dogfish_count integer,
hake_count integer,
herring_count integer,
hyas_count integer,
jonah_crab_count integer,
lumpfish_count integer,
northern_stone_crab_count integer,
redfish_count integer,
rock_cod_count integer,
scallop_count integer,
sculpin_count integer,
skate_count integer,
squid_count integer,
rock_crab_count integer,
thorny_skate_count integer,
whelks_count integer,
winter_flounder_count integer,
witch_flounder_count integer,
dummy5 varchar(1),
algae_pres integer,
alligator_pres integer,
american_plaice_pres integer,
barnacles_pres integer,
basket_star_pres integer,
brittle_star_pres integer,
capelin_pres integer,
clam_pres integer,
cod_pres integer,
dogfish_pres integer,
eel_pouts_count integer,
eel_pres integer,
eel_pouts_pres integer,
gaspereau_pres integer,
gastropod_pres integer,
haddock_pres integer,
hake_pres integer,
halibut_pres integer,
hermit_crab_pres integer,
herring_pres integer,
hyas_pres integer,
jelly_fish_pres integer,
jonah_crab_pres integer,
laminaire_pres integer,
limaces_pres integer,
lumpfish_pres integer,
mackerel_pres integer,
monkfish_pres integer,
mussles_pres integer,
northern_stone_crab_pres integer,
octapus_pres integer,
pollock_pres integer,
poux_de_mer_pres integer,
redfish_pres integer,
reeds_pres integer,
rock_crab_pres integer,
rock_cod_pres integer,
roundnose_grenadir_pres integer,
sand_dollar_pres integer,
scallop_pres integer,
sculpin_pres integer,
sea_anemone_pres integer,
sea_cucumber_pres integer,
sea_mouse_pres integer,
sea_potatoe_pres integer,
sea_raven_pres integer,
sea_snake_pres integer,
sea_star_pres integer,
sea_urchin_pres integer,
sea_worm_pres integer,
shrimp_pres integer,
skate_pres integer,
smelt_pres integer,
sponge_pres integer,
squid_pres integer,
starfish_pres integer,
trout_pres integer,
thorny_skate_pres integer,
tommy_cod_pres integer,
tuna_pres integer,
turbot_pres integer,
whelks_pres integer,
winter_flounder_pres integer,
witch_flounder_pres integer,
wolf_fish_pres integer,
yellowtail_pres integer,
picture_pres varchar(5),
species_id varchar(3),
date_tow varchar(13),
male_count_crab_database integer,
female_count_crab_database integer,
test2 varchar(1),
minilog_file varchar(20),
start_time varchar(8),
end_time varchar(8),
mean_surface numeric(15,10),
count_surface numeric(16, 9),
giant_scallop_pres integer,
giant_scallop_count integer,
clams_count integer,
sea_urchin_count integer,
sand_dollars_count integer,
starfish_count integer,
sea_star_count integer,
poux_de_mer_count integer,
brittle_star_count integer,
shrimps_count integer,
anemones_count integer,
basket_star_count integer,
sea_cucumbers_count integer,
barnacles_count integer,
sponges_count integer,
sea_potato_count integer,
jellyfish_count integer,
snailfish_count integer,
gastropod_count integer,
mussels_count integer,
laminaire_count integer,
reeds_count integer,
tommy_cod_count integer,
halibut_count integer,
greenland_halibut_count integer,
yellowtail_count integer,
hermit_crab_count integer,
capelin_count integer,
wolf_fish_count integer,
roundnose_grenadier_count integer,
sea_raven_count integer,
monkfish_count integer,
smelt_count integer,
haddock_count integer,
pollock_count integer,
mackeral_count integer,
tuna_count integer,
gaspareau_count integer,
trout_count integer,
eel_count integer,
sea_snake_count integer,
octopus_count integer,
sea_worm_count integer,
sea_mouse_count integer,
coque_count integer,
coque_pres integer,
region varchar(3),
subzone varchar(8),
sector varchar(10),
region_date_tow varchar(17),
region_year varchar(8),
subzone_year varchar(13),
sector_year varchar(15),
project_category varchar(14),
    
);

load data infile "/home/jae/data/snowcrab/datadump/scan.mar.tab" into table scanmar fields terminated by '\t' ;

show warnings;



# 4 -----------
drop table species;
CREATE TABLE species(
habitat varchar(254),
distribution varchar(254),
commonnames varchar(196),
description varchar(254),
speciesname varchar(34),
current_re int,
records int,
sciname1 varchar(49),
sciname2 varchar(33),
idnumber int,
habitat_la varchar(7),
distribution2 varchar(12),
description2 varchar(11),
speciesid int
);

load data infile "/home/jae/data/snowcrab/datadump/species.csv" into table species fields terminated by '|';

show warnings;



# 5 ---  minilog.temp.alldata.dbf
drop table minilogall;
CREATE TABLE minilogall(
DAT date,
TIME varchar(11),
TEMPERATUR numeric(6,2),
FILENAME varchar(14),
REFERENCE varchar(17),
START_TIME varchar(7),
END_TIME_T varchar(8),
REC_NUMBER int,
SERIAL_NUM int,
SEQUENCER int,
TOTAL_TIME varchar(8),
COMMENT varchar(30),
YEAR int,
ZONE varchar(4),
DEPTH numeric(12,4),
POS_FILENA varchar(3)
);

load data infile "/home/jae/data/snowcrab/datadump/minilog.temp.alldata.csv" into table minilogall fields terminated by '|';

show warnings;



# 6 --- minilog.temp.dbf
drop table minilog;
CREATE TABLE minilog(
DATE date,
TOTAL_TIME varchar(8),
NUMBER_OF int,
REFERENCE varchar(17),
TEMP_CHART varchar(3),
SOURCE varchar(1),
MIN_TEMPER numeric(5,2)
);

load data infile "/home/jae/data/snowcrab/datadump/minilog.temp.csv" into table minilog fields terminated by '|';

show warnings;




