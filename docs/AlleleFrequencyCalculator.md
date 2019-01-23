# AlleleFrequencyCalculator

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Allele Frequency Calculator


## DEPRECATED

Use bioalcidae

## Usage

```
Usage: allelefreqcalc [options] Files
  Options:
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

 * vcf
 * af


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew allelefreqcalc
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/AlleleFrequencyCalculator.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/AlleleFrequencyCalculator.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/AlleleFrequencyCalculatorTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/AlleleFrequencyCalculatorTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **allelefreqcalc** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)





```

 curl -s  "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/integrated_call_sets/ALL.chr22.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz" | gunzip -c | \
java -jar dist/allelefreqcalc.jar | head

CHR	POS	ID	REF	ALT	TOTAL_CNT	ALT_CNT	FRQ
22	16050408	rs149201999	T	C	2184	134	0.06135531
22	16050612	rs146752890	C	G	2184	184	0.08424909
22	16050678	rs139377059	C	T	2184	113	0.051739927
22	16050984	rs188945759	C	G	2184	5	0.0022893774
22	16051107	rs6518357	C	A	2184	127	0.058150183
22	16051249	rs62224609	T	C	2184	157	0.07188645
22	16051347	rs62224610	G	C	2184	650	0.29761904
22	16051453	rs143503259	A	C	2184	160	0.07326008
22	16051477	rs192339082	C	A	2184	2	9.157509E-4

```




