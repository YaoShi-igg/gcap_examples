# Plot the earthquake and stations by gmt.
# Plot observation system using Geographic Coordinate System.

######################################### i/o paths ######################################### 
source_file="./input/source.in"
receiver_file="./input/receiver.in" 
out_file="./figures/Observation_System_Geo"
# info for plot
figure_type="png,pdf"
title="Observation-System"
########################################## i/o paths ######################################### 

########################################### Region ############################################ 
# X,Y,Z range; 
lon1=114.84; lon2=117.84
lat1=35.66; lat2=38.66
dep1=0; dep2=30
dlon=`echo "$lon2 - $lon1" | bc`; absdlon=${dlon#-}
dlat=`echo "$lat2 - $lat1" | bc`; absdlat=${dlat#-}
ddep=`echo "$dep2 - $dep1" | bc`; absddep=${ddep#-}
# scale lon based on lat
absdlon=`echo "$absdlon * c(($lat1 + $lat2)/2*3.1415926/180)" | bc -l`

# length unit: inch
lon_range=5
dep_zoom=2.5  # zoom the dep for plot
lat_range=$(printf "%.3f" `echo "scale=4;$absdlat*$lon_range/$absdlon"|bc`)
dep_range=$(printf "%.3f" `echo "scale=4;$absddep*$lon_range*$dep_zoom/($absdlon*111.19)"|bc`)

lon_range=`echo $lon_range|awk '{print $1"i"}'`  
lat_range=`echo $lat_range|awk '{print $1"i"}'`  
dep_range=`echo $dep_range|awk '{print $1"i"}'`  

# echo $lon_range $lat_range $dep_range
########################################### Region ############################################ 
gmt makecpt -Cgray -T-8000/2000 > topo.cpt
topo=./topo.cpt

gmt begin $out_file $figure_type
gmt subplot begin 2x2 -Fs$lon_range,$dep_range/$lat_range,$dep_range -A -M0.2c/0.1c -T$title

# set the figure and fonts
gmt set PS_PAGE_ORIENTATION Portrait PS_MEDIA 9ix10.5i
gmt set FONT_TAG 20p FONT_HEADING 25p MAP_HEADING_OFFSET 0p FONT_LABEL 20p FONT_ANNOT 16p

    gmt subplot set 0,0
    region="$lon1/$lon2/$lat1/$lat2"
    projection="X$lon_range/$lat_range"
    # gmt basemap -R$region -J$projection -Bya0.5+l"Latitude" -BWsrt
    gmt grdimage @earth_relief_15s -J$projection -R$region -Bya0.5+l"Latitude (deg.)" -BWsen -I+a-45+nt1+m0 -C$topo
    # plot sources and receivers
    awk '{print $2,$3}' $source_file      | gmt plot -Sa0.5c -W0.5p,black -Gred
    awk '{print $2,$3,$1}' $source_file   | gmt text -F+f8p,1,black+jTC -D0c/-0.3c
    awk '{print $2,$3}' $receiver_file    | gmt plot -St0.5c -W0.5p,black -Gblue
    awk '{print $2,$3,$1}' $receiver_file | gmt text -F+f6p,1,black+jTC -D0c/-0.2c

    gmt subplot set 0,1
    region="$dep1/$dep2/$lat1/$lat2"
    projection="X$dep_range/$lat_range"
    gmt basemap -R$region -J$projection -Bxa10+l"Depth" -BwSen
    # plot sources
    awk '{print $4,$3}' $source_file   | gmt plot -Sa0.5c -W0.5p,black -Gred
    awk '{print $4,$3,$1}' $source_file   | gmt text -F+f8p,1,black+jTC -D0c/-0.3c

    gmt subplot set 1,0
    region="$lon1/$lon2/$dep1/$dep2"
    projection="X$lon_range/-$dep_range"
    gmt basemap -R$region -J$projection -Bxa0.5+l"Longitude" -Bya10+l"Depth" -BWSrt
    # plot sources
    awk '{print $2,$4}' $source_file   | gmt plot -Sa0.5c -W0.5p,black -Gred
    awk '{print $2,$4,$1}' $source_file   | gmt text -F+f8p,1,black+jTC -D0c/-0.3c

gmt subplot end
gmt end
rm -r *cpt
