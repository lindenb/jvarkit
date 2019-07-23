# VcfAlleleDepthFraction

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

filter VCF for strange FORMAT:AD fraction


## Usage

```
Usage: vcfadfraction [options] Files
  Options:
    -filter, --filter
      Variant FILTER
      Default: AD_RATIO
    -f, --filtered
      ignore FILTER-ed **GENOTYPES**
      Default: false
    -gtf, --gtf
      Genotype FILTER
      Default: AD_RATIO
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -het, --het
      AD ratio for **HET** genotypes. HET genotype should have x <= 
      AD[1]/(AD[0]+AD[1])<= (1-x)
      Default: 0.2
    -hom, --hom
      AD ratio for **HOM_REF** or **HOM_VAR** genotypes. HOM_REF genotype 
      should have x <= AD[1]/(AD[0]+AD[1]). HOM_VAR genotype should have  
      AD[1]/(AD[0]+AD[1]) >= (1-x).
      Default: 0.05
    -maxFilteredGenotypes, --maxFilteredGenotypes
      Set Variant FILTER if number of BAD genotype is greater than 'x'. 
      Negative is ignore.
      Default: -1
    -maxFractionFilteredGenotypes, --maxFractionFilteredGenotypes
      Set Variant FILTER if percent of BAD genotype is greater than 'x'. 
      Negative is ignore.
      Default: -1
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * vcf
 * allele-balance
 * depth


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfadfraction
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20190723

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfAlleleDepthFraction.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfAlleleDepthFraction.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfadfraction** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example


