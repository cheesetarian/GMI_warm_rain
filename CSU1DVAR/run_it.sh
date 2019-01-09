#!/bin/sh

dat="150601"
ppfile="../pp/output/geos5.GMIg.$dat.pp" 
#ppfile2="../pp/output/geos5.GMIg.$dat.pp" #if using separate cal for either GPROF/NT2
out="../sample_output/integrated.vX.$dat.nc"
gout="../sample_output/gprof_out.vX.$dat.BIN"
csu1out="../sample_output/1dvar.vX.$dat.nc" #code currently outputs this only

goe $ppfile $ppfile $out $gout $csu1out

#gzip -f $gout
