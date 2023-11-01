#!/bin/sh
method=$1 # Block construction methods: GAB, GAM or SPI
Haploview=~/software/Haploview.jar # Path to your Haploview.jar #

if [ $method = "GAB" ]
then
 block='GABRIELblocks'
elif [ $method = "GAM" ]
then
 block='4GAMblocks'
else
 block='SPINEblocks'
fi

### Subset dataset by plink ###
plink_1.9 --recode HV \
--bfile ../Data/geno.DHGC.plink \
--keep ../Data/Sample.list \
--extract ../Data/Marker.list \
--out ../Data/geno_forHaploView &

wait

### Run HaploView ###
for chr in {1..10}
do
nohup java -jar $Haploview \
 -nogui \
 -skipcheck \
 -memory 500000 \
 -pedfile ../Data/geno_forHaploView.chr-${chr}.ped \
 -info ../Data/geno_forHaploView.chr-${chr}.info \
 -out ../Data/HVBlock.chr-${chr} \
 -blockoutput ${method} > ../Data/HVBlock-${chr}.log 2>&1 &

wait
done

### Modify format ###
for chr in {1..10}
do

grep BLOCK ../Data/HVBlock.chr-${chr}.${block}  | awk -v chr="$chr" 'BEGIN{FS=" "}{for(j=4;j<=NF;j++){print chr,NR,$j}}' >> ../Data/HVBlock.${block}.list

done

