# Script to run fk
# For green function, it only depending on focal depth and epicentral distance
# So, here we will calculate once for repeated depth and distance.

######################### i/o paths ###########################
# velocity model
# vp in third value(off k) and no need to flatten (off f)
model=sd
# path for Glib and model
model_dir=./glib/$model
if [ ! -d $model_dir ]; then
  mkdir -p $model_dir
fi
cd $model_dir
# path for input file of sources and receivers
input_file=../../input/fk.in
######################### i/o paths ###########################

########################## params #############################
# focal depth(km)
dep_array=(16 17 18 19)
# number of points(2^n) and sample interval(s)
# T = npts * dt, fmax=1/2dt, df=1/T
npts=512; dt=0.1
# smooth and low pass filter
# npts_inter=npts*smth; fc=(1-taper)*fmax
smth=1; taper=0.3
# integration interval
# dk=dk*pi/max(dist, delta_h)
dk=0.2
# high pass for narrow frequency calculation range
# from [0, df, fmax] to [f1, df, fmax] (Hz)
f1=0.01; f2=0.02
# limit of integration
# kmin=wpmin; kmax=sqrt(kmax^2+wpmax)
pmin=0; pmax=1; kmax=15
# type of sources(0-explosion, 1-single force, 2-DC)
# for explosion, need to synthesize DC first
S=2
# list of epicentral distances
dist_array=$(awk '{print $3}' $input_file) # exract depth array from input file
dist_array=($(printf "%s\n" "${dist_array[@]}" | awk '!seen[$0]++')) # remove the repeated elements
########################## params ###########################

# run fk.pl
for dep in ${dep_array[@]}; do
    fk.pl -M$model/$dep -N$npts/$dt/$smth/$dk/$taper -H$f1/$f2 -P$pmin/$pmax/$kmax -S$S ${dist_array[@]}
done

