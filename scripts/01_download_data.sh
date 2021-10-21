#!/bin/bash

# alorenzetti 20211014

# description ####
# this script will download
# raw data from SRA
# using the default runInfo
# file selected from SRA search

# requires sra-tools from bioconda
# version should be > 2.10.0
# to use fasterq-dump

# requires pigz from conda
# requires rename from bioconda
 
# setting up number of threads
threads=16

# initializing conda 
eval "$(conda shell.bash hook)"

# getting started ####
conda activate sratools_env

# parsing data in a user friendly way
# source was the full metadata
# file from SRA Run Selector

# PRJNA609253
sed 's/\".*\",//' data/PRJNA609253.csv | \
tail -n +2 | \
cut -d',' -f1,13 | \
sort -t, -Vk2,2 > data/PRJNA609253_usr_friendly.csv

# PRJNA630692
sed 's/\".*\",//' data/PRJNA630692.csv | \
tail -n +2 | \
cut -d',' -f1,20,18,24,27,28 | \
sed 's/CO2: //;s/Coffea canephora/CC/;s/Coffea arabica/CA/;s/biological replicate /BR/' | \
sed 's/\(^.*C[AC],\).* = \(.*$\)/\1\2/' | \
sort -t, -Vk6,6 -Vk3,3 > data/PRJNA630692_usr_friendly.csv

# downloading experiments
# PRJNA609253 and PRJNA630692
for i in data/PRJNA609253.csv data/PRJNA630692.csv ; do

    prjna=${i##*/}
    prjna=${prjna%.*}

    if [[ ! -d raw/$prjna ]] ; then
        mkdir -p raw/$prjna

        # getting accessions
        acc=`cut -d',' -f1,1 data/${prjna}.csv | tail -n +2`

        # downloading
        for j in $acc ; do
            if [[ ! -f raw/$prjna/${j}_1.fastq.gz ]] ; then
            fasterq-dump --split-3 \
                         --threads $threads \
                         -O raw/$prjna \
                         $j
            fi
        done
    fi

done

# for the next step we are going
# to need an efficient compressing tool
conda activate pigz_env

# giving alias to libraries
# and specially unifying the 
# libraries of experiment PRJNA630692
# PRJNA609253
while IFS="," read srr sample ; do
    if [[ -f raw/PRJNA609253/${srr}_1.fastq && ! -f raw/PRJNA609253/${srr}_1.fastq.gz ]] ; then pigz -p $threads -c raw/PRJNA609253/${srr}_1.fastq > raw/PRJNA609253/${srr}_1.fastq.gz ; fi
    if [[ -f raw/PRJNA609253/${srr}_2.fastq && ! -f raw/PRJNA609253/${srr}_2.fastq.gz ]] ; then pigz -p $threads -c raw/PRJNA609253/${srr}_2.fastq > raw/PRJNA609253/${srr}_2.fastq.gz ; fi
    if [[ -f raw/PRJNA609253/${srr}_1.fastq.gz && ! -f raw/PRJNA609253/${srr}-${sample}_R1.fastq.gz ]] ; then mv raw/PRJNA609253/${srr}_1.fastq.gz raw/PRJNA609253/${srr}-${sample}_R1.fastq.gz ; fi
    if [[ -f raw/PRJNA609253/${srr}_2.fastq.gz && ! -f raw/PRJNA609253/${srr}-${sample}_R2.fastq.gz ]] ; then mv raw/PRJNA609253/${srr}_2.fastq.gz raw/PRJNA609253/${srr}-${sample}_R2.fastq.gz ; fi
done < data/PRJNA609253_usr_friendly.csv

# PRJNA630692
# getting list of samples
slist=`cut -d',' -f6 data/PRJNA630692_usr_friendly.csv | sort -V | uniq`

for i in $slist ; do
    filelist=`awk -v s=$i -v FS=',' '{if($6 == s){print "raw/PRJNA630692/"$1".fastq"}}' data/PRJNA630692_usr_friendly.csv | xargs`
    atts=`awk -v s=$i -v FS=',' '{if($6 == s){print}}' data/PRJNA630692_usr_friendly.csv | head -n +1`
    srr=`echo $atts | cut -d',' -f1`
    ppm=`echo $atts | cut -d',' -f2`
    code=`echo $atts | cut -d',' -f3 | sed 's/_R1.*$//'`
    spp=`echo $atts | cut -d',' -f4`
    br=`echo $atts | cut -d',' -f5`
    sample=`echo $atts | cut -d',' -f6`

    if [[ ! -f raw/PRJNA630692/${srr}-${spp}_${ppm}_${i}_R1.fastq.gz ]] ; then
        cat $filelist | pigz -p $threads -c > raw/PRJNA630692/${srr}-${spp}_${ppm}_${i}_R1.fastq.gz

        # removing uncompressed files if the compression was successful
        if [[ $? == 0 ]] ; then
                rm $filelist
                echo "$i done!"
        else
                echo "Error in $i compression."
        fi
    else
        echo "raw/PRJNA630692/${srr}-${spp}_${ppm}_${i}_R1.fastq.gz already exists."
    fi
done

# adjusting names of samples according to supp. material table s1 from
# https://www.mdpi.com/1422-0067/22/6/3125/s1
# removing SRR from names
rename 's/SRR.*-//' raw/PRJNA630692/*.fastq.gz

if [[ ! -f data/dictPRJNA630692.txt ]] ; then 
    touch data/dictPRJNA630692.txt

    for i in "1,CA_380ppm,CA_380ppm_25C" \
             "2,CA_380ppm,CA_380ppm_42C" \
             "3,CA_700ppm,CA_700ppm_25C" \
             "4,CA_700ppm,CA_700ppm_42C" \
             "5,CC_380ppm,CC_380ppm_25C" \
             "6,CC_380ppm,CC_380ppm_42C" \
             "7,CC_700ppm,CC_700ppm_25C" \
             "8,CC_700ppm,CC_700ppm_42C" \
             "9,CA_380ppm,CA_380ppm_37C" \
             "10,CA_700ppm,CA_700ppm_37C" \
             "11,CC_380ppm,CC_380ppm_37C" \
             "12,CC_700ppm,CC_700ppm_37C" ; do

        echo $i >> data/dictPRJNA630692.txt
    done
fi

while IFS="," read sample oldname newname ; do
    for i in {A..C} ; do
        completesample=`echo ${sample}${i}`
        mv raw/PRJNA630692/${oldname}_${completesample}_R1.fastq.gz raw/PRJNA630692/${newname}_${completesample}_R1.fastq.gz 
    done
done < data/dictPRJNA630692.txt