# VcfLiftOver

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Lift-over a VCF file


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar vcfliftover  [options] Files

Usage: vcfliftover [options] Files
  Options:
    --adaptivematch
      Use adapative liftover minmatch using the ratio between the min allele 
      size and the longest allele size
      Default: false
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
  * -f, --chain
      LiftOver chain file. Can be a local chain file, a URL 'https://hgdownload.soe.ucsc.edu/goldenpath/hg19/liftOver/hg19ToCriGri1.over.chain.gz', 
      or a chain identifier like 'hg19ToHg38'.
    -check, --check
      Check variant allele sequence is the same on REF
      Default: false
    -x, --failed
      (file.vcf) write variants failing the liftOver here. Optional.
    -failtag, --failtag
      failed INFO tag
      Default: LIFTOVER_FAILED
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --ignore-chain-validation
      Ignore LiftOver chain validation
      Default: false
    --indel, --indels
      do not LiftOver indels
      Default: false
    --info
      remove attribute from INFO on the fly
      Default: []
    -m, --minmatch
      lift over min-match.
      Default: 0.95
    -o, --out
      Output file. Optional . Default: stdout
  * -D, -R, -r, --reference
      Destination Reference. Indexed fasta Reference file. This file must be 
      indexed with samtools faidx and with picard/gatk 
      CreateSequenceDictionary or samtools dict
    -T, --tag
      INFO tag
      Default: LIFTOVER
    --version
      print version and exit

```


## Keywords

 * vcf
 * liftover



## Creation Date

20240114

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/liftover/VcfLiftOver.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/liftover/VcfLiftOver.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/liftover/VcfLiftOverTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/liftover/VcfLiftOverTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfliftover** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Warning

you'd better have a look at gatk/picard LiftOverVcf


## Example

Lift over a hg19 file to hg18

Before liftover:

```bash
$  curl -s "https://raw.github.com/arq5x/gemini/master/test/test5.vep.snpeff.vcf" | grep -v "#

chr1	145273345	.	T	C	289.85	.	.	GT:AD:DP:GQ:PL	0/0:226,22:250:99:0,158,4259	0/1:224,24:250:5.77:6,0,5314	0/1:219,28:249:57.30:57,0,5027	0/1:215,34:250:99:269,0,3796
chr1	156011444	.	T	C	2523.46	.	.	GT:AD:DP:GQ:PL	0/1:24,15:40:99:214,0,443	0/1:32,36:68:99:702,0,794	1/1:1,59:61:99:1656,132,0	0/0:34,1:35:69.10:0,69,717
chr5	64982321	.	T	C	61.12	.	.	GT:AD:DP:GQ:PL	1/1:0,2:2:6:58,6,0	1/1:0,1:1:3.01:37,3,0	./.	./.
chr10	1142208	.	T	C	3404.30	.	.	GT:AD:DP:GQ:PL	1/1:1,37:39:87.16:940,87,0	1/1:0,29:29:78.20:899,78,0	1/1:0,24:24:66.14:729,66,0	1/1:0,30:30:75.18:836,75,0
chr10	126678092	.	G	A	89.08	.	.	GT:AD:DP:GQ:PL	0/0:64,3:67:99:0,165,1505	0/0:11,1:12:7.31:0,7,240	0/0:52,10:62:54.97:0,55,1263	0/1:35,9:44:99:125,0,693
chr10	135210791	.	T	C	65.41	.	.	GT:AD:DP:GQ:PL	0/0:4,0:4:9:0,9,84	1/1:0,3:3:6.02:74,6,0	1/1:0,1:1:3.01:37,3,0	0/0:3,0:3:9.02:0,9,100
chr13	48873835	.	G	A	58.95	.	.	GT:AD:DP:GQ:PL	./.	./.	1/1:0,2:2:6.01:62,6,0	1/1:0,1:1:3.01:31,3,0
chr20	36779424	.	G	A	128.76	.	.	GT:AD:DP:GQ:PL	0/0:49,1:52:63.68:0,64,969	0/0:17,0:17:30.05:0,30,320	0/0:93,0:94:99:0,216,2384	0/1:24,9:33:99:165,0,505
chrX	17819377	.	T	C	7515.25	.	.	GT:AD:DP:GQ:PL	1/1:0,125:126:99:2343,237,0	1/1:0,26:26:78.14:837,78,0	1/1:0,90:92:99:2640,244,0	1/1:0,74:75:99:1695,171,0
```

Running the liftover:

```bash
$  curl -s "https://raw.github.com/arq5x/gemini/master/test/test5.vep.snpeff.vcf" |\
  java -jar dist/jvarkit.jar vcfliftover \
      -f hg19ToHg18 \
      -X failing.vcf |\
grep -v "#"

chr1	143984702	.	T	C	289.85	.	LIFTOVER=chr1|145273345	GT:AD:DP:GQ:PL0/0:226,22:250:99:0,158,4259	0/1:224,24:250:5.77:6,0,5314	0/1:219,28:249:57.30:57,0,5027	0/1:215,34:250:99:269,0,3796
chr1	154278068	.	T	C	2523.46	.	LIFTOVER=chr1|156011444	GT:AD:DP:GQ:PL0/1:24,15:40:99:214,0,443	0/1:32,36:68:99:702,0,794	1/1:1,59:61:99:1656,132,0	0/0:34,1:35:69.10:0,69,717
chr5	65018077	.	T	C	61.12	.	LIFTOVER=chr5|64982321	GT:AD:DP:GQ:PL1/1:0,2:2:6:58,6,0	1/1:0,1:1:3.01:37,3,0	./.	./.
chr10	1132208	.	T	C	3404.30	.	LIFTOVER=chr10|1142208	GT:AD:DP:GQ:PL	1/1:1,37:39:87.16:940,87,0	1/1:0,29:29:78.20:899,78,0	1/1:0,24:24:66.14:729,66,0	1/1:0,30:30:75.18:836,75,0
chr10	126668082	.	G	A	89.08	.	LIFTOVER=chr10|126678092	GT:AD:DP:GQ:PL	0/0:64,3:67:99:0,165,1505	0/0:11,1:12:7.31:0,7,240	0/0:52,10:62:54.97:0,55,1263	0/1:35,9:44:99:125,0,693
chr10	135060781	.	T	C	65.41	.	LIFTOVER=chr10|135210791	GT:AD:DP:GQ:PL	0/0:4,0:4:9:0,9,84	1/1:0,3:3:6.02:74,6,0	1/1:0,1:1:3.01:37,3,0	0/0:3,0:3:9.02:0,9,100
chr13	47771836	.	G	A	58.95	.	LIFTOVER=chr13|48873835	GT:AD:DP:GQ:PL./.	./.	1/1:0,2:2:6.01:62,6,0	1/1:0,1:1:3.01:31,3,0
chr20	36212838	.	G	A	128.76	.	LIFTOVER=chr20|36779424	GT:AD:DP:GQ:PL0/0:49,1:52:63.68:0,64,969	0/0:17,0:17:30.05:0,30,320	0/0:93,0:94:99:0,216,2384	0/1:24,9:33:99:165,0,505
chrX	17729298	.	T	C	7515.25	.	LIFTOVER=chrX|17819377	GT:AD:DP:GQ:PL1/1:0,125:126:99:2343,237,0	1/1:0,26:26:78.14:837,78,0	1/1:0,90:92:99:2640,244,0	1/1:0,74:75:99:1695,171,0
```





