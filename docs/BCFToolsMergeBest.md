# BCFToolsMergeBest

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Scan a VCF file generated by 'bcftools merge --force-samples', identify duplicate samples, keep the best


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar bcftoolsmergebest  [options] Files

Usage: bcftoolsmergebest [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --out
      Output file. Optional . Default: stdout
    --rename-samples
      TSV file for manually assigning samples . Syntax: old-name(tab)new-name. 
      Both names must be present in the VCF file
    --version
      print version and exit
    -c
      comma separated list of criteria used to compare to genotype: GQ:best 
      GQ, DP: highest depth, AD : best AD ratio according to diploy genotype, 
      PL: highest value, GQ: better called than no-call, FT: PASS variant are 
      better, RND: random
      Default: GQ,DP,AD,PL,FT,RND

```


## Keywords

 * merge
 * vcf
 * bcftools



## See also in Biostars

 * [https://www.biostars.org/p/9594639](https://www.biostars.org/p/9594639)



## Creation Date

20240604

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bcftoolsmergebest/BCFToolsMergeBest.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bcftoolsmergebest/BCFToolsMergeBest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bcftoolsmergebest** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Motivation

can a VCF file generated by 'bcftools merge --force-samples', identify duplicate samples, keep the best"

## Example

```
$ bcftools merge --force-samples src/test/resources/rotavirus_rf.vcf.gz src/test/resources/rotavirus_rf.freebayes.vcf.gz  |\
	grep "#CHROM" -A 3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5	2:S5	2:S2	2:S4	2:S3	2:S1
RF01	243	.	A	C	0	.	DPB=28;EPPR=3.37221;GTI=0;MQMR=60;NS=5;NUMALT=1;ODDS=7.64661;PAIREDR=1;PQR=0;PRO=0;QR=408;RO=24;RPPR=3.37221;SRF=24;SRP=55.1256;SRR=0;DP=28;AB=0.5;ABP=3.0103;AF=0.1;AO=4;CIGAR=1X;DPRA=1.66667;EPP=11.6962;LEN=1;MEANALT=1;MQM=60;PAIRED=1;PAO=0;PQA=0;QA=68;RPL=4;RPP=11.6962;RPR=0;RUN=1;SAF=4;SAP=11.6962;SAR=0;TYPE=snp;AN=10;AC=1	GT:AD:AO:DP:PL:QA:QR:RO	./.:.:.:.:.:.:.:.	./.:.:.:.:.:.:.:.	./.:.:.:.:.:.:.:.	./.:.:.:.:.:.:.:.	./.:.:.:.:.:.:.:.	0/0:4,0:0:4:0,12,63:0:68:4	0/0:7,1:1:8:0,7,92:17:119:7	0/1:2,2:2:4:20,0,20:34:34:2	0/0:7,1:1:8:0,7,92:17:119:7	0/0:4,0:0:4:0,12,63:0:68:4
RF01	280	.	A	C	0	.	DPB=24;EPPR=4.58955;GTI=0;MQMR=60;NS=5;NUMALT=1;ODDS=6.90616;PAIREDR=1;PQR=0;PRO=0;QR=374;RO=22;RPPR=4.58955;SRF=22;SRP=50.7827;SRR=0;DP=24;AB=0.4;ABP=3.44459;AF=0.1;AO=2;CIGAR=1X;DPRA=1.05263;EPP=3.0103;LEN=1;MEANALT=1;MQM=60;PAIRED=1;PAO=0;PQA=0;QA=34;RPL=1;RPP=3.0103;RPR=1;RUN=1;SAF=2;SAP=7.35324;SAR=0;TYPE=snp;AN=10;AC=1	GT:AD:AO:DP:PL:QA:QR:RO	./.:.:.:.:.:.:.:.	./.:.:.:.:.:.:.:.	./.:.:.:.:.:.:.:.	./.:.:.:.:.:.:.:.	./.:.:.:.:.:.:.:.	0/1:3,2:2:5:17,0,33:34:51:3	0/0:7,0:0:7:0,21,109:0:119:7	0/0:2,0:0:2:0,6,32:0:34:2	0/0:7,0:0:7:0,21,109:0:119:7	0/0:3,0:0:3:0,9,48:0:51:3
RF01	351	.	T	A	0	.	DPB=25;EPPR=3.94093;GTI=1;MQMR=60;NS=5;NUMALT=1;ODDS=8.27411;PAIREDR=1;PQR=0;PRO=0;QR=357;RO=21;RPPR=3.94093;SRF=21;SRP=48.6112;SRR=0;DP=25;AB=0.4;ABP=3.87889;AF=0.2;AO=4;CIGAR=1X;DPRA=1;EPP=11.6962;LEN=1;MEANALT=1;MQM=60;PAIRED=1;PAO=0;PQA=0;QA=68;RPL=0;RPP=11.6962;RPR=4;RUN=1;SAF=4;SAP=11.6962;SAR=0;TYPE=snp;AN=10;AC=2	GT:AD:AO:DP:PL:QA:QR:RO	./.:.:.:.:.:.:.:.	./.:.:.:.:.:.:.:.	./.:.:.:.:.:.:.:.	./.:.:.:.:.:.:.:.	./.:.:.:.:.:.:.:.	0/0:7,0:0:7:0,21,109:0:119:7	0/1:3,2:2:5:17,0,33:34:51:3	0/0:4,0:0:4:0,12,63:0:68:4	0/1:3,2:2:5:17,0,33:34:51:3	0/0:4,0:0:4:0,12,63:0:68:4


$ bcftools merge --force-samples src/test/resources/rotavirus_rf.vcf.gz src/test/resources/rotavirus_rf.freebayes.vcf.gz  |\
	java -jar dist/jvarkit.jar bcftoolsmergebest |\
	grep "#CHROM" -A 3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S3	S4	S5	S1	S2
RF01	243	.	A	C	0	.	AB=0.5;ABP=3.0103;AC=1;AF=0.1;AN=10;AO=4;CIGAR=1X;DP=28;DPB=28;DPRA=1.66667;EPP=11.6962;EPPR=3.37221;GTI=0;LEN=1;MEANALT=1;MQM=60;MQMR=60;NS=5;NUMALT=1;ODDS=7.64661;PAIRED=1;PAIREDR=1;PAO=0;PQA=0;PQR=0;PRO=0;QA=68;QR=408;RO=24;RPL=4;RPP=11.6962;RPPR=3.37221;RPR=0;RUN=1;SAF=4;SAP=11.6962;SAR=0;SRF=24;SRP=55.1256;SRR=0;TYPE=snp	GT:AD:AO:DP:PL:QA:QR:RO	0/0:7,1:1:8:0,7,92:17:119:7	0/1:2,2:2:4:20,0,20:34:34:2	0/0:4,0:0:4:0,12,63:0:68:4	./.	./.
RF01	280	.	A	C	0	.	AB=0.4;ABP=3.44459;AC=1;AF=0.1;AN=10;AO=2;CIGAR=1X;DP=24;DPB=24;DPRA=1.05263;EPP=3.0103;EPPR=4.58955;GTI=0;LEN=1;MEANALT=1;MQM=60;MQMR=60;NS=5;NUMALT=1;ODDS=6.90616;PAIRED=1;PAIREDR=1;PAO=0;PQA=0;PQR=0;PRO=0;QA=34;QR=374;RO=22;RPL=1;RPP=3.0103;RPPR=4.58955;RPR=1;RUN=1;SAF=2;SAP=7.35324;SAR=0;SRF=22;SRP=50.7827;SRR=0;TYPE=snp	GT:AD:AO:DP:PL:QA:QR:RO	0/0:7,0:0:7:0,21,109:0:119:7	0/0:2,0:0:2:0,6,32:0:34:2	0/1:3,2:2:5:17,0,33:34:51:3	./.	./.
RF01	351	.	T	A	0	.	AB=0.4;ABP=3.87889;AC=2;AF=0.2;AN=10;AO=4;CIGAR=1X;DP=25;DPB=25;DPRA=1;EPP=11.6962;EPPR=3.94093;GTI=1;LEN=1;MEANALT=1;MQM=60;MQMR=60;NS=5;NUMALT=1;ODDS=8.27411;PAIRED=1;PAIREDR=1;PAO=0;PQA=0;PQR=0;PRO=0;QA=68;QR=357;RO=21;RPL=0;RPP=11.6962;RPPR=3.94093;RPR=4;RUN=1;SAF=4;SAP=11.6962;SAR=0;SRF=21;SRP=48.6112;SRR=0;TYPE=snp	GT:AD:AO:DP:PL:QA:QR:RO	0/1:3,2:2:5:17,0,33:34:51:3	0/0:4,0:0:4:0,12,63:0:68:4	0/0:7,0:0:7:0,21,109:0:119:7	./.	./.
```

