# MantaMerger

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Merge Vcf from Manta VCF.


## Usage

```
Usage: mantamerger [options] Files
  Options:
    --bnd-distance
      Two BND variants are the same if their bounds are distant by less than 
      xxx bases. A distance specified as a positive integer.Commas are 
      removed. The following suffixes are interpreted : b,bp,k,kb,m,mb
      Default: 100
    --check-bnd-mate
      When comparing two BND, check that their mate (using the ALT allele) are 
      the same too
      Default: false
    -c, --contig
      limit to this contig
    -x, --exclude
      Exclude regions in this bed file
    --force-svtype
      When comparing two SV variants, their INFO/SVTYPE should be the same. 
      Default is to just use coordinates to compare non-BND variants.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --no-bnd
      discar BND
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    --sv-alleles-bases
      When comparing two non-BND SV variants, use their ALT alleles to adjust 
      the interval. It solves the problem of  
      'chr2:10556028:AATTATATAATTAAATTAATTATATAATT:A'  vs 
      'chr2:10556028:A:AATTATATAATTAAATTAATTATATAATT'. See 
      https://twitter.com/yokofakun/status/1169182491606999046 
      Default: false
    --sv-overlap-fraction
      Two CNV/DEL/.. variants are the same if they share 'x' fraction of their 
      size. 
      Default: 0.75
    --sv-small-overlap
      Two non-BND variants are the same if they overlap and both have a 
      length<= 'x'. A distance specified as a positive integer.Commas are 
      removed. The following suffixes are interpreted : b,bp,k,kb,m,mb
      Default: 10
    --version
      print version and exit

```


## Keywords

 * sv
 * burden
 * manta
 * vcf


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew mantamerger
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20190916

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/manta/MantaMerger.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/manta/MantaMerger.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/structvar/manta/MantaMergerTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/structvar/manta/MantaMergerTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **mantamerger** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

 
 # Input
 
 input is a list of indexed vcf files or one file with the '.list' suffix containing the path to the vcfs
 
 
 # Example
 
 ```
 $ find src -name "manta*z" > jeter.list
 $ java -jar dist/mantamerger.jar jeter.list 2> /dev/null
 
 (...)
 
 
 
 ```
