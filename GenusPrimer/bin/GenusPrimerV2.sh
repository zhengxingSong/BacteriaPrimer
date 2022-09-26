#!/usr/bin/env bash

data=$(dirname $PWD)
:> $data/DownLoadc.list
:> $data/PrimerDc.list
for i in `ls $data/result2`;
do
	echo $i
	echo -ne "bash $data/bin/Pipeline_StepI.sh $data/result2/$i\n" >> $data/DownLoadc.list
	echo -ne "bash $data/bin/Pipeline_StepII_GenusV2.sh $data/result2/$i\n" >> $data/PrimerDc.list
done 
#ParaFly -c $data/Downloadc.list -CPU 10
