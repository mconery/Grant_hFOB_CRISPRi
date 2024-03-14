#!/bin/bash
#for each SNPs, this script selects the first 3 (in order of pick order) non overlapping guides#

module load BEDTools

for j in `cat snp88.list`;do
grep ${j} guides_filt_short.bed | sed "s/_/	/g"|sort -k4 -k5 -n > list_${j}.txt
while read -a a$idx; do let idx++;done <list_${j}.txt


readarray -t myArray < list_${j}.txt
N=$(wc -l list_${j}.txt|cut -d" " -f1)
k=$(($N-1))

#initialize
echo ${myArray[0]}|awk '{print $4"_"$5}' > ${j}_selected.txt
echo ${myArray[0]}|sed "s/ /\t/g" > ${j}_selected_1.bed
echo ${myArray[0]}|sed "s/ /\t/g" > ${j}_selected_2.bed

#for i in 0 1 2; do

for i in `seq 0 ${k}`;do

#create first bed file
echo ${myArray[${i}]} | sed "s/ /\t/g" > ${j}_$(($i+1)).bed

#create second bed file, except for last step
if (($i<$k));then
echo ${myArray[${i}+1]} | sed "s/ /\t/g" > ${j}_$(($i+2)).bed
fi


#create intersection file of new guide with previous guides
intersectBed -a ${j}_$(($i+2)).bed -b ${j}_selected_1.bed > int_$(($i+2))-1.txt
intersectBed -a ${j}_$(($i+2)).bed -b ${j}_selected_2.bed > int_$(($i+2))-2.txt

#if intersection files are empty, add guide 
	if ! [ -s "int_$(($i+2))-1.txt" ] && ! [ -s "int_$(($i+2))-2.txt" ];then 
	echo `awk '{print $4"_"$5}' ${j}_$(($i+2)).bed` >> ${j}_selected.txt
	cp ${j}_$(($i+2)).bed ${j}_selected_2.bed
	fi

#count how many guides are in selected file
a=`wc -l ${j}_selected.txt|cut -d" " -f1`
#echo $a

#if there are 3 guides in selected file, exit
if (($a>=3)); then break;fi

done
done

