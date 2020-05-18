#!/bin/tcsh

set ifile = "/media/cwharris/Transcend/CC/Xcorr/fta.csv"
set lines1 = "lines1.gmt"
set lines2 = "lines2.gmt"
set lines3 = "lines3.gmt"
set lines4 = "lines4.gmt"
set lines5 = "lines5.gmt"
set lines6 = "lines6.gmt"
rm -f $lines1 $lines2 $lines3 $lines4 $lines5 $lines6
sort -k8 -g -t',' $ifile | awk -F, '{if($1=="S"&&$8!="nan")print (">-Z"$8"\n"$4,$5"\n"$6,$7)}' >> $lines1
sort -k9 -g -t',' $ifile | awk -F, '{if($1=="S"&&$9!="nan")print (">-Z"$9"\n"$4,$5"\n"$6,$7)}' >> $lines2
sort -k10 -g -t',' $ifile | awk -F, '{if($1=="S"&&$10!="nan")print (">-Z"$10"\n"$4,$5"\n"$6,$7)}' >> $lines3
sort -k11 -g -t',' $ifile | awk -F, '{if($1=="S"&&$11!="nan")print (">-Z"$11"\n"$4,$5"\n"$6,$7)}' >> $lines4
sort -k13 -g -t',' $ifile | awk -F, '{if($1=="S"&&$13!="nan")print (">-Z"$13"\n"$4,$5"\n"$6,$7)}' >> $lines5
sort -k15 -g -t',' $ifile | awk -F, '{if($1=="S"&&$15!="nan")print (">-Z"$15"\n"$4,$5"\n"$6,$7)}' >> $lines6
set ps='gv_lines.ps'
set land='sienna'
set sea='lightsteelblue'
set minlon=-12
set maxlon=33
set minlat=31
set maxlat=64
set misc = "-R -J -Cvel.cpt -K -O -W0.1p"
gmt makecpt -T1.5/4.5 -D -Cseis > vel.cpt
echo 6.2 -0.5 "Group Velocity km/s" | gmt pstext -R0/10/0/10 -JX8i -K -N > $ps
echo 4 8.75 "8s" | gmt pstext -R -J -O -K -N >> $ps
echo 8.45 8.75 "15s" | gmt pstext -R -J -O -K -N >> $ps
echo 12.95 8.75 "25s" | gmt pstext -R -J -O -K -N >> $ps
echo 4 4.05 "45s" | gmt pstext -R -J -O -K -N >> $ps
echo 8.45 4.05 "85s" | gmt pstext -R -J -O -K -N >> $ps
echo 12.95 4.05 "125s" | gmt pstext -R -J -O -K -N >> $ps
gmt pscoast -R$minlon/$maxlon/$minlat/$maxlat -A100 -Bx10 -By10 -BNseW -G$land -S$sea -JM7.5 -X0c -Y9c -K -O>> $ps
gmt psxy $lines1 $misc >> $ps
gmt pscoast -R -A100 -Bx10 -By10 -Bnsew -G$land -S$sea -J -X9c -Y0c -K -O >> $ps
gmt psxy $lines2 $misc >> $ps
gmt pscoast -R -A100 -Bx10 -By10 -Bnsew -G$land -S$sea -J -X9c -Y0c -K -O >> $ps
gmt psxy $lines3 $misc >> $ps
gmt pscoast -R -A100 -Bx10 -By10 -Bnsew -G$land -S$sea -J -X-18c -Y-9.5c -K -O >> $ps
gmt psxy $lines4 $misc >> $ps
gmt pscoast -R -A100 -Bx10 -By10 -Bnsew -G$land -S$sea -J -X9c -Y0c -K -O >> $ps
gmt psxy $lines5 $misc >> $ps
gmt pscoast -R -A100 -Bx10 -By10 -Bnsew -G$land -S$sea -J -X9c -Y0c -K -O >> $ps
gmt psxy $lines6 $misc >> $ps
gmt psscale -D-5.15/-0.85/4.6/0.3h -B1g1 -Cvel.cpt -O -K >> $ps
ps2pdf $ps
rm gv_lines.ps vel.cpt $lines1 $lines2 $lines3 $lines4 $lines5 $lines6
xdg-open gv_lines.pdf
