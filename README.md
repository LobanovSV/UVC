# UVC
Universal Variant Caller

This software performs in-depth read alignment and variant calling. It calls any type of variant including SNP, INDEL, STR, QTR, and any their combination.

It was intially written to analyse bam files and than re-written to analyse sam files with specific region because this way allowed me to process the data in the quickest way.

Please be sure that your folder structure satisfy the requirements written below in order to run the script `AlignRegion.m`.

#### 1. Software
Download the software from https://github.com/LobanovSV/UVC.git.

#### 2. Reference genome
If you use `hg19` and `hg38` genomes, you can download the required files from
https://drive.google.com/drive/folders/11Iu7mCbnSSuJc8jnCswQRLQBn01WqDyW?usp=sharing.
If you wish to use your own genome, look into structure of the `hg19` and `hg38` genomes and create your own in the same way:
- The folder with software must contain the sub-folder `RefGenome`.
- The sub-folder `RefGenome` must contain the sub-sub-folder with reference genome name (for example, `hg19` or `hg38`).
- The sub-sub-folder must contain mat files `chr.mat` (`chr` = 1, 2, etc.) each having only one variable `FASTA` with nucleotide sequence of the respective chromosome.
- The nucleotide sequence must contain only capital letters `A`, `C`, `G`, `T`, `N`.

#### 3. The script `AlignRegion.m`
Specify the variable `DirWork` with path to the folder you wish to analyse.
This folder must contain the `SAMs` sub-folder with

(a) `Pos.mat` file containing two variables: `chr` and `posW`. For example, `chr = 5` and `posW = [145838397, 145838950]` used for TCERG1 QTR calling. `chr` is chromosome and `posW` contains coordinates of the specific region you wish to be analysed.

(b) `ID.mat` file containing the variable `ID` with names of individuals (or samples). For example, `ID = {'HD322-45'; 'HD524-75'; 'SZ62'}`.

(c) sam files with reads overlaping the region `chr:posW(1)-posW(2)`. Sam files must be named with integer numbers, for example, `1.sam`, `2.sam`, `3.sam`, etc. and correspond to the names of individuals in the variable `ID`. For example, `HD322-45 <---> 1.sam`, `HD524-75 <---> 2.sam`, `SZ62 <---> 3.sam`, etc.
For example, you can use `samtools` to create such files:
```
samtools view BAMs/HD322-45.bam 5:145838397-145838950 > SAMs/1.sam
```
**NOTE:** sam files should contain specific region. If sam files contain all reads or reads for the whole chromosome, the script might not work.


Please cite the following publication if you use this code:
Lobanov et. al., bioRxiv 2021.07.16.452643.


This software is written by S.V. Lobanov and is openly available at
https://github.com/LobanovSV/UVC




<a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-nd/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/4.0/">Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License</a>.
