# Script to run gcap3D

######################### i/o paths ###########################
# velocity model
model=cus
# path for Glib and model
glib=./glib
data_dir=./event_data_rtz
######################### i/o paths ###########################

########################## params #############################
# library for green function
G=$glib
# sample interval for data(s) (-H)
dt=0.1
H=$dt
# corner frequencies for bandpass filter of pnl and surface wave (-C)
f1_pnl=0.05; f2_pnl=0.2; f1_sw=0.02; f2_sw=0.2
C=$f1_pnl/$f2_pnl/$f1_sw/$f2_sw
# weight for Pnl wave and distance scaling powers for Pnl and surface wave (-D)
w1=2; p1=1; p2=0.5
D=$w1/$p1/$p2
# max time window for Pnl and surface wave (-T)
len_pnl=35; len_sw=70
T=$len_pnl/$len_sw
# max shift for pnl and surface wave, and tie between SH shift and SV shift(-S)
pnl_shift=5; sw_shift=10; tie=0
S=$pnl_shift/$sw_shift/$tie
# search range for strike, dip rake (-R)
strike1=0; strike2=360; dip1=0; dip2=90; rake1=-90; rake2=90
R=$strike1/$strike2/$dip1/$dip2/$rake1/$rake2
# search interval for strike, dip rake (-I)
ddegree=5
I=$ddegree
# initial value and search step for ISO and CLVD components (-J)
iso=0; diso=0; clvd=0; dclvd=0
J=$iso/$diso/$clvd/$dclvd
# data type for seismic record (-W)
# 1-velocity; 2-displacement
data_type=1
W=$data_type
# output other local minimums whose error no more than uncertainty range (-X)
# misfit-min<n*sigma
nsigma=10
X=$nsigma
# generate waveform-fit plot with plotting scale.
# yscale: amplitude in inch for the first trace of the page
yscale=100000 
P=$yscale
# list of depths for inverse
dep_list=(5 10 15 20)
########################## params #############################

for source in `ls $data_dir`; do
    seis_dir=$data_dir/$source
    out_dir=./results/$source
    if [ ! -d $out_dir ]; then
    mkdir -p $out_dir
    fi
    for dep in "${dep_list[@]}"; do
        # run gcap3D
        cmd="-M${model}_$dep -G$G -H$H -C$C -D$D -T$T -S$S -R$R -I$I -J$J -W$W -X$X -P$P $seis_dir "
        cap3D.pl $cmd
        # Convert ps to pdf and move to outdir
        mv $seis_dir/${model}_${dep}* $out_dir
        ps2pdf $out_dir/${model}_${dep}.ps $out_dir/${model}_${dep}.pdf
    done
    # # output the results for different depth
    # grep -h Event $out_dir/${model}_*.out >  $out_dir/depth.out
    # perl depth.pl $out_dir/depth.out $seis_dir >  $out_dir/depth.ps
    # ps2pdf $out_dir/depth.ps $out_dir/depth.pdf
    # # remove 
    rm -r .gmt* $seis_dir/.gmt* $out_dir/*ps
done

