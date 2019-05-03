# ConvertBedChromosomes

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Convert the names of the chromosomes in a Bed file


## Usage

```
Usage: bedrenamechr [options] Files
  Options:
    -c, --column
      1-based chromosome column(s), multiple separated by commas
      Default: 1
    -convert, --convert
      What should I do when  a converstion is not found
      Default: RAISE_EXCEPTION
      Possible Values: [RAISE_EXCEPTION, SKIP, RETURN_ORIGINAL]
    -d, --delim
      field delimiter.
      Default: 	
    -s, --header
      Ignore lines starting with this java regular expression
      Default: (#|browser|track)
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
  * -f, --mapping, -m, -R
      load a custom name mapping.Chromosome mapping file. If the file looks 
      like a NGS file (fasta, vcf, bam...) the mapping is extracted from a 
      dictionary; Otherwise, it is interpreted as a mapping file ( See 
      https://github.com/dpryan79/ChromosomeMappings )
    -o, --out
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * bed
 * chromosome
 * contig
 * convert


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew bedrenamechr
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/ConvertBedChromosomes.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/ConvertBedChromosomes.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bedrenamechr** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



## Example

```bash
$   curl -s http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz |\
    gunzip -c | \
    java -jar dist/bedrenamechr.jar -f src/main/resources/chromnames/hg19_to_g1kv37.tsv -c 2 |\
   tail


uc011nca.2	Y	+	59213948	59276439	59230880	59274621	7	59213948,59222126,59230781,59233166,59252482,59272370,59274552,	59214117,59222281,59230919,59233257,59252550,59272463,59276439,	P51809	uc011nca.2
uc004fxl.3	Y	+	59213948	59276439	59222135	59274621	7	59213948,59222126,59230781,59233166,59252482,59272370,59274552,	59214117,59222216,59230919,59233257,59252550,59272463,59276439,	P51809-3	uc004fxl.3
uc004fxk.3	Y	+	59213948	59276439	59222135	59274809	7	59213948,59222126,59228291,59230781,59233166,59272370,59274552,	59214117,59222281,59228349,59230919,59233257,59272463,59276439,	P51809-2	uc004fxk.3
uc011ncb.2	Y	+	59213948	59276439	59222274	59274621	7	59213948,59222126,59230781,59233166,59252482,59272370,59274552,	59214117,59222277,59230919,59233257,59252550,59272463,59276439,	B4DE96	uc011ncb.2
uc010nxr.2	Y	+	59330251	59340490	59335611	59340461	9	59330251,59334000,59335576,59336119,59336347,59337090,59337948,59338753,59340193,	59330458,59334179,59335690,59336231,59336526,59337236,59338150,59338859,59340490,	B4E011	uc010nxr.2
uc004fxm.1	Y	+	59330251	59343488	59330414	59342523	9	59330251,59334078,59335552,59336119,59336354,59337119,59337948,59338753,59342486,	59330458,59334179,59335690,59336231,59336526,59337236,59338150,59338859,59343488,	B9ZVT0	uc004fxm.1
uc004fxn.1	Y	+	59330251	59343488	59330430	59343080	9	59330251,59335576,59336119,59336347,59337090,59337948,59338753,59340193,59342486,	59330458,59335690,59336231,59336526,59337236,59338150,59338859,59340278,59343488,	Q01113	uc004fxn.1
uc004fxo.1	Y	+	59352972	59356131	59353631	59356130	7	59352972,59354351,59354669,59354993,59355369,59355682,59355972,	59353819,59354463,59354816,59355130,59355505,59355884,59356131,	I3L0A4	uc004fxo.1
uc022cpg.1	Y	+	59354984	59358336	59355427	59358045	7	59354984,59355369,59355682,59355972,59356790,59357702,59357911,	59355130,59355505,59355884,59356131,59356943,59357771,59358336,	Q9NQA3	uc022cpg.1
uc011ncc.1	Y	-	59358328	59360854	59358328	59358328	3	59358328,59360006,59360500,	59359508,59360115,59360854,	uc011ncc.1

```

