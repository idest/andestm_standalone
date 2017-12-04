#!/bin/bash

#  Qs.sh
#  
#
#  Created by Andrés Tassara Oddo on 16-04-14.
#

for dir in `ls ../output_termal`

do

cd ../output_termal
cd $dir

mkdir Post

cd Post

T=../estructura_termal_lon_lat_Qs_sigma_area_Z_Tz

#### crea grilla Qs

awk '{if ($3!=NaN) print $1, $2, $3*-1000}' $T > Qs.xyz

R=-R-80/-60/-45/-10

xyz2grd -GQs.grd Qs.xyz -I0.2 $R

#### plot Qs

B=-B5g5WSen

P=Qs.ps

makecpt -Chot -T0/120/5 > Q.cpt

grdimage -JM8 $R Qs.grd $B -CQ.cpt -K > $P

pscoast -JM8 $R -W1 -N1 -O -K >> $P

awk '{if ($5==3) print $1, $2, $3}' ../../../Post/ObsQs/QsObs.txt | psxy -JM8 $R -St0.3 -CQ.cpt -W0.7 -K -O >> $P

awk '{if ($5==4) print $1, $2, $3}' ../../../Post/ObsQs/QsObs.txt | psxy -JM8 $R -Ss0.3 -CQ.cpt -W0.7 -K -O >> $P

awk '{if ($5==1) print $1, $2, $3}' ../../../Post/ObsQs/QsObs.txt | psxy -JM8 $R -Sc0.3 -CQ.cpt -W0.7 -K -O >> $P

awk '{if ($5==2) print $1, $2, $3}' ../../../Post/ObsQs/QsObs.txt | psxy -JM8 $R -Sh0.3 -CQ.cpt -W0.7 -K -O >> $P


psscale -CQ.cpt -D0.1/6/8/0.5 -B20 -K -O >> $P




### plot diferencia Obs-Pred

grdtrack ../../../Post/ObsQs/QsObs.txt -GQs.grd > tmp

awk '{if ($5==1) print $1, $2, $3, $6, "c0.3"}' tmp > tmp1
awk '{if ($5==2) print $1, $2, $3, $6, "h0.3"}' tmp >> tmp1
awk '{if ($5==3) print $1, $2, $3, $6, "t0.3"}' tmp >> tmp1
awk '{if ($5==4) print $1, $2, $3, $6, "s0.3"}' tmp >> tmp1



makecpt -Cgray -T-6000/6000/200 -Z > TB.cpt

grdimage -JM8 $R ../../../Post/GEBCO/gebco_08_-80_-45_60_-10.nc -B5g5wSen -CTB.cpt -K -O -X8.5 >> $P

#psbasemap -JM8 $R $B -X8.5 -K -O >> $P

awk '{if (($3-$4)>=80) print $1, $2, $5}' tmp1 | psxy -JM8 $R -K -O -S -W1 -Gblack >> $P

awk '{if ((($3-$4)>=20) && (($3-$4)<=40)) print $1, $2, $5}' tmp1 | psxy -JM8 $R -K -O -S -W1 -Glightblue >> $P

awk '{if ((($3-$4)<=-20) && (($3-$4)>=-40)) print $1, $2, $5}' tmp1 | psxy -JM8 $R -K -O -S -W1 -Glightred >> $P

awk '{if ((($3-$4)>=40) && (($3-$4)<=80)) print $1, $2, $5}' tmp1 | psxy -JM8 $R -K -O -S -W1 -Gblue >> $P

awk '{if (($3-$4)<=-40) print $1, $2, $5}' tmp1 |psxy -JM8 $R -K -O -S -W1 -Gred >> $P

awk '{if ((($3-$4)>=-20) && (($3-$4)<=20)) print $1, $2, $5}' tmp1 |psxy -JM8 $R -K -O -S -W1 -Gwhite >> $P

pscoast -JM8 $R -W1 -N1 -O -K >> $P


##### plot obs vs pred

psbasemap -JX8 -R0/100/0/100 -B10g10WSen -K -O -x9.5 >> $P

psxy -JX8 -R0/100/0/100 -W4 -K -O << END >> $P
0 0
100 100
END

psxy -JX8 -R0/100/0/100 -Glightred -K -O -L << END >> $P
20 0
100 0
100 80
END

psxy -JX8 -R0/100/0/100 -Glightblue -K -O -L << END >> $P
0 20
80 100
0 100
END

psxy -JX8 -R0/100/0/100 -L -Gblue -K -O << END >> $P
0 100
0 40
60 100
END


psxy -JX8 -R0/100/0/100 -Gred -K -O -L << END >> $P
40 0
100 0
100 60
END

awk '{print $4, $3, $5}' tmp1 | psxy -JX8 -R0/100/0/100 -W1 -S -K -O >> $P


psbasemap -JX8 -R0/100/0/100 -B10g10WSen -K -O >> $P

gs $P


### vuelve a directorio output_termal

cd ../../

done
