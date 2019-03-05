# Biostar130456

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Split individual VCF files from multisamples VCF file


## Usage

```
Usage: biostar130456 [options] Files
  Options:
    --disable-vc-attribute-recalc
      When genotypes are removed/changed, Dd not recalculate variant 
      attributes like DP, AF, AC, AN...
      Default: false
    -f, --filtered
      remove filtered Genotype
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -z, --homref
      remove homzygote REF/REF
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
  * -p, --pattern
      output file pattern. Must contain the word __SAMPLE__
    -x, --uncalled
      remove uncalled genotypes
      Default: false
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

```


## Keywords

 * vcf
 * samples
 * sample



## See also in Biostars

 * [https://www.biostars.org/p/130456](https://www.biostars.org/p/130456)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew biostar130456
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar130456.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar130456.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/biostar/Biostar130456Test.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/biostar/Biostar130456Test.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **biostar130456** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
bash
$   curl -sL "https://raw.githubusercontent.com/arq5x/bedtools2/bc2f97d565c36a82c1a0b12f570fed4398001e5f/test/map/test.vcf" |\
    java -jar dist/biostar130456.jar -x -z -p "sample.__SAMPLE__.vcf.gz" 
sample.NA00003.vcf.gz
sample.NA00001.vcf.gz
sample.NA00002.vcf.gz

$ gunzip -c sample.NA00003.vcf.gz
(...)
##source=myImputationProgramV3.1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00003
chr1	10	rs6054257	G	A	29	PASS	AF=0.5;DB;DP=14;H2;NS=3	GT:DP:GQ:HQ	1/1:5:43
chr1	20	rs6040355	A	G,T	67	PASS	AA=T;AF=0.333,0.667;DB;DP=10;NS=2	GT:DP:GQ	2/2:4:35
chr1	130	microsat1	GTC	G,GTCT	50	PASS	AA=G;DP=9;NS=3	GT:DP:GQ	1/1:3:40
chr2	130	microsat1	GTC	G,GTCT	50	PASS	AA=G;DP=9;NS=3	GT:DP:GQ	1/1:3:40
```

## See also

 * GATK SelectVariants with option -sn , or bcftools view --samples-file

