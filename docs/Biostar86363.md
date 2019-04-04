# Biostar86363

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Set genotype of specific sample/genotype comb to unknown in multisample vcf file. See http://www.biostars.org/p/86363/


## Usage

```
Usage: biostar86363 [options] Files
  Options:
    --disable-vc-attribute-recalc
      When genotypes are removed/changed, Dd not recalculate variant 
      attributes like DP, AF, AC, AN...
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    --vc-attribute-recalc-ignore-filtered
      When recalculating variant attributes like DP AF, AC, AN, ignore 
      FILTERed **Genotypes**
      Default: false
    --vc-attribute-recalc-ignore-missing
      Ignore missing VCF headers (DP, AF, AC, AN). Default behavior: adding 
      VCF header if they're missing
      Default: false
    --version
      print version and exit
  * -G
      genotypes to reset. Format :CHROM(tab)POS(tab)ref(tab)SAMPLE. REQUIRED.

```


## Keywords

 * sample
 * genotype
 * vcf



## See also in Biostars

 * [https://www.biostars.org/p/86363](https://www.biostars.org/p/86363)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew biostar86363
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar86363.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar86363.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **biostar86363** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example
```bash
$ cat reset.txt
20	14370	NA00001
20	1234567	NA00003
20	1110696	NA00002

$ curl "https://raw.github.com/jamescasbon/PyVCF/master/vcf/test/example-4.1.vcf" |\
  java -jar dist/biostar86363.jar -G reset.txt 

##fileformat=VCFv4.1
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GR,Number=1,Type=Integer,Description="(1) = Genotype was reset by Biostar86363:Set genotype of specific sample/genotype comb to unknown in multisample vcf fi
le. See http://www.biostars.org/p/86363/">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
##fileDate=20090805
##phasing=partial
##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
##source=myImputationProgramV3.1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002	NA00003
20	14370	rs6054257	G	A	29	PASS	AF=0.5;DB;DP=14;H2;NS=3	GT:DP:GQ:GR:HQ	.|.:1:48:1:51,51	1|0:8:48:0:51,51	1/1:5:43:0
20	17330	.	T	A	3	q10	AF=0.017;DP=11;NS=3	GT:GQ:DP:HQ	0|0:49:3:58,50	0|1:3:5:65,3	0/0:41:3
20	1110696	rs6040355	A	G,T	67	PASS	AA=T;AF=0.333,0.667;DB;DP=10;NS=2	GT:DP:GQ:GR:HQ	1|2:6:21:0:23,27	.|.:0:2:1:18,2	2/2:4:35:0
20	1230237	.	T	.	47	PASS	AA=T;DP=13;NS=3	GT:GQ:DP:HQ	0|0:54:7:56,60	0|0:48:4:51,51	0/0:61:2
20	1234567	microsat1	GTC	G,GTCT	50	PASS	AA=G;DP=9;NS=3	GT:DP:GQ:GR	0/1:4:35:0	0/2:2:17:0	./.:3:40:1
```

