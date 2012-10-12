#!/bin/bash

echo $1

mencoder mf://./$1*.png  -fps 5 -ofps 5 -o video.flv -of lavf -oac mp3lame -lameopts abr:br=64 -srate 22050 -ovc lavc -lavcopts vcodec=flv:keyint=50:vbitrate=700:mbd=2:mv0:trell:v4mv:cbp:last_pred=3 -vf yadif -sws 9
