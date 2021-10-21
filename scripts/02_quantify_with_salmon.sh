#!/bin/bash

# alorenzetti 20211019

# description ####
# this script will
# use salmon to generate bootstrapped
# quantification of trasncriptomes

# requires salmon from bioconda
# conda install -c conda-forge boost-cpp=1.74.0
# conda install -c bioconda salmon=1.5.2

# setting up number of threads
threads=16

# initializing conda 
eval "$(conda shell.bash hook)"

# getting started ####
conda activate salmon_env

# getting transcriptome and annotation files
# from ncbi refseq
if [[ ! -f data/Carabica.fa ]] ; then
    wget -O data/Carabica.fa.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/713/225/GCF_003713225.1_Cara_1.0/GCF_003713225.1_Cara_1.0_rna.fna.gz
    gunzip data/Carabica.fa.gz
fi

if [[ ! -f data/Carabica.gtf ]] ; then
    wget -O data/Carabica.gtf.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/713/225/GCF_003713225.1_Cara_1.0/GCF_003713225.1_Cara_1.0_genomic.gtf.gz
    gunzip data/Carabica.gtf.gz
fi

if [[ ! -f data/Carabica.gff ]] ; then
    wget -O data/Carabica.gff.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/713/225/GCF_003713225.1_Cara_1.0/GCF_003713225.1_Cara_1.0_genomic.gff.gz
    gunzip data/Carabica.gff.gz
fi

if [[ ! -d data/Carabica ]] ; then
    # making ref genome index
    salmon index -p $threads -t data/Carabica.fa -i data/Carabica > data/salmon_index_out.log 2> data/salmon_index_err.log
fi

# salmon quant for
# PRJNA609253
if [[ ! -d salmon_quant/PRJNA609253 ]] ; then
    mkdir -p salmon_quant/PRJNA609253

    # running salmon quant for each sample
    for i in raw/PRJNA609253/*.fastq.gz ; do
        prefix=${i%%_R[12]*}
        prefix=${prefix##*/}

            if [[ ! -d salmon_quant/PRJNA609253/${prefix} ]] ; then
                salmon quant -i data/Carabica \
                            -l A \
                            -1 raw/PRJNA609253/${prefix}_R1.fastq.gz \
                            -2 raw/PRJNA609253/${prefix}_R2.fastq.gz \
                            -p $threads \
                            -o salmon_quant/PRJNA609253/${prefix} \
                            --softclip \
                            --numBootstraps 20 > salmon_quant/PRJNA609253/${prefix}.log 2>&1
            fi
    done

    # getting the names of all output dirs
    outdirs=`ls raw/PRJNA609253/*.fastq.gz | sed 's/_R[12].*//;s/raw/salmon_quant/' | sort | uniq | xargs`

    # merging results from samples
    salmon quantmerge --quants $outdirs -o salmon_quant/PRJNA609253/PRJNA609253_merged_tpm.sf -c tpm
    salmon quantmerge --quants $outdirs -o salmon_quant/PRJNA609253/PRJNA609253_merged_numreads.sf -c numreads
fi

# PRJNA630692
if [[ ! -d salmon_quant/PRJNA630692 ]] ; then
    mkdir -p salmon_quant/PRJNA630692

    # running salmon quant for each sample
    for i in raw/PRJNA630692/*.fastq.gz ; do
        prefix=${i%%_R[12]*}
        prefix=${prefix##*/}

            if [[ ! -d salmon_quant/PRJNA630692/${prefix} ]] ; then
                salmon quant -i data/Carabica \
                            -l A \
                            -r raw/PRJNA630692/${prefix}_R1.fastq.gz \
                            -p $threads \
                            -o salmon_quant/PRJNA630692/${prefix} \
                            --softclip \
                            --numBootstraps 20 > salmon_quant/PRJNA630692/${prefix}.log 2>&1
            fi
    done

    # getting the names of all output dirs
    outdirs=`ls raw/PRJNA630692/*.fastq.gz | sed 's/_R[12].*//;s/raw/salmon_quant/' | sort | uniq | xargs`

    # merging results from samples
    salmon quantmerge --quants $outdirs -o salmon_quant/PRJNA630692/PRJNA630692_merged_tpm.sf -c tpm
    salmon quantmerge --quants $outdirs -o salmon_quant/PRJNA630692/PRJNA630692_merged_numreads.sf -c numreads
fi
