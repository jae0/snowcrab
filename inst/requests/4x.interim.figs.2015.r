
p = bio.snowcrab::load.environment()

a = snowcrab.timeseries.db( DS="biologicals.direct", regions='cfa4x', vn=c('R0.mass','t'), trim=0 )
figure.timeseries.R0( outdir=file.path(p$annual.results, "timeseries", "survey"),infile = a,specific.area='cfa4x' )
figure.timeseries.t( outdir=file.path(p$annual.results, "timeseries", "survey"),infile = a,specific.area='cfa4x' )
