# COSMO

To cite: [![DOI](https://zenodo.org/badge/182076507.svg)](https://zenodo.org/badge/latestdoi/182076507)

Condition-specific Mapping of Operons (COSMO) is a tool for detecting operons in microorganisms

![Github Actions Cosmo Testing](https://github.com/sanbi-sa/COSMO/actions/workflows/cosmo.yml/badge.svg)

An **Operon** consists of a group of structural genes that codes for enzymes involved in a metabolic pathway.
The genes of an operon are located contiguously on a stretch DNA and are under the control of one promoter.
A single mRNA unit is transcribed and translated into separate proteins.

*COSMO* takes as input, a reference genome aligned BAM file and a GTF file with gene coordinates. 
COSMO predicts operons by calculating average coverages of the genes/CDSs and their intergenic regions (IGRs). 
The user provides four cut-offs which determine whether genes/CDSs form part of the same operon. 
If they fail any one of these cut-offs, then they are not part of the operon.

#### Requirements
1. Python >= 3.7
2. Pysam >=0.15.0

#### Installation (Run from a python script)
```
python setup.py install  
```

#### Usage (Run from command line)
```
user_input.py [-h] [-D GDEPTH] [-d IDEPTH] [-F GFACTOR] [-f IFACTOR]
                     [-p PREFIX]
                     ref length bam gtf

Detect possible genome operons using RNA expression coverages

positional arguments:
  ref                   Name of reference sequence in BAM file
  length                Length of reference genome according to GTF file
  bam                   Bam input file
  gtf                   Gtf input file

optional arguments:
  -h, --help            show this help message and exit
  -D GDEPTH, --gdepth GDEPTH
                        Average number of reads per base required to consider
                        a gene expressed
  -d IDEPTH, --idepth IDEPTH
                        Average number of reads per base required to consider
                        an IGR expressed
  -F GFACTOR, --gfactor GFACTOR
                        Maximum fold difference (FD) allowed between two genes/CDSs of an operon
  -f IFACTOR, --ifactor IFACTOR
                        Maximum FD allowed between an intergenic region (IGR) and its flanking genes/CDSs
  -o OUTPUT, --output OUTPUT
                        CSV output filename (output file is written to "output" folder)
                        

Example code for running a file:

python user_input.py -D 1 -d 2 -F 5.0 -f 10.0 -o Strain_5_control.csv NC_000962 4411532 /home/Sally/Input_bam_files/Strain_5_control.bam /home/Sally/Input_bam_files/NC_000962_tuberculosis.gtf

 Where the:

1) minimum CDS coverage =  1
2) min IGR coverage = 2
3) maximum FD between adjacent CDSs/genes = 5.0
4) max FD between IGR and flanking CDSs = 10.0
5) output file name given by user = Strain_5_control.csv 
6) Genome name in GTF file = NC_000962 
7) Size of genome in GTF file = 4411532
8) Input bam file = Strain_5_control.bam 
9) GTF file name = NC_000962_tuberculosis.gtf

```    
