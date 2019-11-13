# VcfRebase

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Restriction sites overlaping variations in a vcf


## Usage

```
Usage: vcfrebase [options] Files
  Options:
    -A, --attribute
      VCF INFO attribute
      Default: ENZ
    -E, -enzyme, --enzyme
      restrict to that enzyme name. Default: use all enzymes
      Default: []
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --out
      Output file. Optional . Default: stdout
    -R, -reference, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    --version
      print version and exit
    -w, -weight, --weight
      min enzyme weight 6 = 6 cutter like GAATTC, 2 = 2 cutter like ATNNNNNNAT
      Default: 5.0

```


## Keywords

 * vcf
 * rebase
 * restriction
 * enzyme


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfrebase
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20131115

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfrebase/VcfRebase.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfrebase/VcfRebase.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfrebase/VcfRebaseTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfrebase/VcfRebaseTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfrebase** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

 
 ## Example
 

 ```
$ java -jar dist/vcfrebase.jar -w 6 -R ~/data/human_g1k_v37.fasta src/test/resources/test_vcf01.vcf | bcftools annotate -x '^INFO/ENZ' | bcftools view --drop-genotypes | grep ENZ

##INFO=<ID=ENZ,Number=.,Type=String,Description="Enzyme overlapping: Format: (Name|Site|Sequence|pos|strand)">
##bcftools_annotateCommand=annotate -x ^INFO/ENZ; Date=Wed Nov 13 10:38:39 2019
1	852063	.	G	A	387	PASS	ENZ=PflMI|CCANNNN^NTGG|CCAGGCCCTGG|852064|+
1	866893	.	T	C	431	PASS	ENZ=SacI|GAGCT^C|GAGCtC|866889|+
1	875770	.	A	G	338	PASS	ENZ=ClaI|AT^CGAT|ATCGaT|875766|+
1	909238	.	G	C	229	PASS	ENZ=PmaCI|CAC^GTG|CACgTG|909235|+
1	913889	.	G	A	372	PASS	ENZ=BsaXI|(9/12)ACNNNNNCTCC(10/7)|GGAGGCCCCgT|913880|-
1	918384	.	G	T	489	PASS	ENZ=DraIII|CACNNN^GTG|CACgCCGTG|918381|+
1	933790	.	G	A	436	PASS	ENZ=BsaXI|(9/12)ACNNNNNCTCC(10/7)|GGAGGAGGGgT|933781|-
1	940005	.	A	G	188	PASS	ENZ=GsuI|CTGGAG(16/14)|CTGGAG|940006|+,BaeI|(10/15)ACNNNNGTAYC(12/7)|GGTaCTGGAGT|940002|-
1	940096	.	C	T	487	PASS	ENZ=BcgI|(10/12)CGANNNNNNTGC(12/10)|cGAGGTGGGTGC|940096|+
1	950113	.	GAAGT	G	1427	PASS	ENZ=Eco57I|CTGAAG(16/14)|CTgaag|950111|+
1	950243	.	A	C	182	PASS	ENZ=BclI|T^GATCA|TGaTCA|950241|+
1	951283	.	C	T	395	PASS	ENZ=NarI|GG^CGCC|GGcGCC|951281|+
1	951564	.	A	G	105	PASS	ENZ=BstXI|CCANNNNN^NTGG|CCaAGTAGTTGG|951562|+
1	952003	.	G	A	177	PASS	ENZ=Bpu10I|CCTNAGC(-5/-2)|CCTCAGC|952004|+,BbvCI|CCTCAGC(-5/-2)|CCTCAGC|952004|+
1	952428	.	G	A	456	PASS	ENZ=EciI|GGCGGA(11/9)|TCCgCC|952425|-
1	953952	.	G	A	490	PASS	ENZ=BsrDI|GCAATG(2/0)|CATTgC|953948|-
1	959155	.	G	A	370	PASS	ENZ=BarI|(7/12)GAAGNNNNNNTAC(12/7)|gAAGCCGCTCTAC|959155|+
1	959231	.	G	A	350	PASS	ENZ=BsaXI|(9/12)ACNNNNNCTCC(10/7)|GGAGGGTCCgT|959222|-
1	960409	.	G	C	357	PASS	ENZ=BseYI|CCCAGC(-5/-1)|CCCAgC|960405|+
1	962210	.	A	G	300	PASS	ENZ=NcoI|C^CATGG|CCaTGG|962208|+
1	964389	.	C	T	32	LowGQXHetSNP;LowGQXHomSNP	ENZ=BseYI|CCCAGC(-5/-1)|cCCAGC|964389|+
1	967658	.	C	T	515	PASS	ENZ=StuI|AGG^CCT|AGGCcT|967654|+
1	970215	.	G	C	379	PASS	ENZ=DrdI|GACNNNN^NNGTC|GACCCCTCGGTC|970216|+
1	972180	.	G	A	403	PASS	ENZ=AgeI|A^CCGGT|ACCgGT|972177|+
1	1004957	.	G	A	316	PASS	ENZ=BsgI|GTGCAG(16/14)|gTGCAG|1004957|+
1	1004980	.	G	A	292	PASS	ENZ=BsePI|G^CGCGC|gCGCGC|1004980|+
1	1011087	.	CG	C	1052	PASS	ENZ=Eam1105I|GACNNN^NNGTC|GACTCTCAGTc|1011077|+
1	1017170	.	C	G	507	PASS	ENZ=AloI|(7/12)GAACNNNNNNTCC(12/7)|GAACAGAGcATCC|1017162|+
```
 
