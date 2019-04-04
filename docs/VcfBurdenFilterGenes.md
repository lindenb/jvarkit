# VcfBurdenFilterGenes

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Filter VEP/SnpEff Output from a list of genes.


## Usage

```
Usage: vcfburdenfiltergenes [options] Files
  Options:
    -a, --add
      [20180627] Gene Names: Add this gene, multiple separated by 
      comma,spaces,semicolon 
      Default: <empty string>
    -filter, --filter
      If empty: remove the variants from the VCF. If not empty, add a token in 
      the column FILTER.
      Default: <empty string>
    -g, --genes
      Gene/transcript file: one name per line
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * gene
 * vcf
 * vep
 * snpeff



## See also in Biostars

 * [https://www.biostars.org/p/353011](https://www.biostars.org/p/353011)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfburdenfiltergenes
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burden/VcfBurdenFilterGenes.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burden/VcfBurdenFilterGenes.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfburdenfiltergenes** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


### Example

```
echo "IL2" > genes.txt
echo "NOCTH2" >>  genes.txt
gunzip -c input.vcf.gz |\
	java -jar dit/vcfburdenfiltergenes.jar -g genes.txt
```

### Example

```
$ wget -O - -q "https://github.com/immune-health/antigen.garnish/raw/f0453336a4859c83d27640c286d3960c1672f164/inst/extdata/testdata/antigen.garnish_hg19anno_example.vcf"  |\
	grep -w 216371854  | cut -f 8 | tr ";" "\n"   | grep ^ANN= | cut -c5- | tr "," "\n"
T|missense_variant|MODERATE|USH2A|USH2A|transcript|NM_206933.2|protein_coding|18/72|c.3884G>A|p.Arg1295Gln|4271/18883|3884/15609|1295/5202||
T|missense_variant|MODERATE|USH2A|USH2A|transcript|NM_007123.5|protein_coding|18/21|c.3884G>A|p.Arg1295Gln|4271/6316|3884/4641|1295/1546||
```

```
$ wget -O - -q "https://github.com/immune-health/antigen.garnish/raw/f0453336a4859c83d27640c286d3960c1672f164/inst/extdata/testdata/antigen.garnish_hg19anno_example.vcf"  |\
	sed 's/PASS\t[A-Z0-9]*;/PASS\t/' |\ ## the vcf above is malformed, quick hack
	java -jar dist/vcfburdenfiltergenes.jar -a "NM_206933.2" |\
	grep -w 216371854  | cut -f 8 | tr ";" "\n"   | grep ^ANN= | cut -c5- | tr "," "\n"
T|missense_variant|MODERATE|USH2A|USH2A|transcript|NM_206933.2|protein_coding|18/72|c.3884G>A|p.Arg1295Gln|4271/18883|3884/15609|1295/5202||
```



## History

  * 20181205 : snpeff scan transcriptID
  * 20180617 : for SNpEFF, now looks into GeneName OR GeneId (was only GeneName)

