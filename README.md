# CRISPR sgRNA finding #

CRISPR sgRNA finding identify all available sgRNAs in given genomic regions.

### Prerequisite ###

* Documents: 2bit file
```
For example, hg38.2bit can be downloaded from http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/
```

* Python package
```
twobitreader: only works for python2
To use python3, the 'long' in __init__.py of twobitreader needs to be replaced with 'int'.
```

### Input ###

* genomic region in bed file
```
chr1	67093579	67093604	C1orf141
chr1	67096251	67096321	C1orf141
```


### Output ###

* Non_repeat_19_bp_sgRNA_file
```
chr1	67093597	67093616	TTTAATCTACAGGAACAAA	C1orf141	-
chr1	67096253	67096272	AGAAGGTGATAGAAATAAA	C1orf141	-
```
* Non_repeat_40_bp_sgRNA_file
```
chr1	67093587	67093627	ttattcttttcTTTAATCTACAGGAACAAAaggacacaaa	C1orf141	-
chr1	67096243	67096283	gttggtctcttAGAAGGTGATAGAAATAAAaggtaattta	C1orf141	-
```
* Repeat_19_bp_sgRNA_file
```
chr1	8358585	8358604	GAGCGGGAGATCCGAGAGC	RERE	-
chr1	8358600	8358619	GAGCGGGAGATCCGAGAGC	RERE	-
```
* Repeat_40_bp_sgRNA_file
```
chr1	8358575	8358615	gggagatccgaGAGCGGGAGATCCGAGAGCgggagctgcg	RERE	-
chr1	8358590	8358630	gggagctccggGAGCGGGAGATCCGAGAGCgggagatccg	RERE	-
```
### Arguments ###
```
usage: CRISPR_sgRNA_finding.py [-h] -f INPUT_FILE -b BIT_PATH -v BIT_VERSION
                               -l SPACER_LENGTH -n {SpCas9,SaCas9,Cpf1}
                               [-t TARGET_PATH] [-p PROCESSING] [-m PAM_MOTIF]
optional arguments:
  -h, --help            show this help message and exit

Required arguments:

  -f INPUT_FILE, --input_file INPUT_FILE
                        The genomic regions targeted by sgRNAs. The format is chromosome, start, end, gene_symbol.
  -b BIT_PATH, --bit_path BIT_PATH
                        Directory where 2bit file is stored
  -v BIT_VERSION, --bit_version BIT_VERSION
                        Bit version used to extract the sequence. Please choose from:
                         1) hg38
                         2) hg19
                         3) mm10
                         4) mm9
  -l SPACER_LENGTH, --spacer_length SPACER_LENGTH
                        Spacer length. Ex: Commonly used spacer length for SpCas9 is 20.
  -n {SpCas9,SaCas9,Cpf1}, --PAM_name {SpCas9,SaCas9,Cpf1}
                        The name for PAM motif.

Optional arguments:

  -t TARGET_PATH, --target_path TARGET_PATH
                        Directory in which all the output stored, default is current path/input_file_sgRNA_folder
  -p PROCESSING, --processing PROCESSING
                        Number of multiprocessing! Default=100
  -m PAM_MOTIF, --PAM_motif PAM_MOTIF
                        PAM motif. Default motifs include:
                         1) SpCas9: NGG
                         2) SaCas9: NNGRR
                         3) Cpf1: TTN
                         You could modify it if you have updated motif.
```
### Example ###
```
sgRNA_finding.py -f demo_hg38_refseq_CDS.txt -v hg38 -l 19 -b /path/of/2bit/file/ -n SpCas9
```