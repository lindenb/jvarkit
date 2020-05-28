# ScanStructuralVariants

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Scan structural variants for case/controls data


## Usage

```
Usage: scansv [options] Files
  Options:
    --all
      Print all original variants from each file instead of printing just one.
      Default: false
    --bed
      A source of intervals. The following suffixes are recognized: vcf, 
      vcf.gz bed, bed.gz, gtf, gff, gff.gz, gtf.gz.Otherwise it could be an 
      empty string (no interval) or a list of plain interval separated by '[ 
      \t\n;,]' 
    --bnd-distance
      Two BND variants are the same if their bounds are distant by less than 
      xxx bases. A distance specified as a positive integer.Commas are 
      removed. The following suffixes are interpreted : b,bp,k,kb,m,mb
      Default: 100
    --check-bnd-mate
      When comparing two BND, check that their mate (using the ALT allele) are 
      the same too
      Default: false
    -c, --controls
      Controls indexed VCF files. a file endings with the suffix '.list' is 
      interpretted as a list of path.
      Default: []
    --force-svtype
      When comparing two SV variants, their INFO/SVTYPE should be the same. 
      Default is to just use coordinates to compare non-BND variants.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -L, --large
      Large number of controls: By default, all VCF readers for controls are 
      opened and are kept opened. It's fast but requires a lot of resources. 
      This option open+close the controls if needed but it makes things 
      slower. It's the number of VCF that should be keept open, So '0' = 
      ignore/all re-open+close (slow)
      Default: 0
    --maf
      Max frequency of variants found in controls. 0:no control should carry 
      the variant
      Default: 0.0
    -o, --out
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

 * cnv
 * indel
 * sv
 * pedigree


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew scansv
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20190815

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/ScanStructuralVariants.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/ScanStructuralVariants.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/structvar/ScanStructuralVariantsTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/structvar/ScanStructuralVariantsTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **scansv** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
find CONTROLS/ -name "*.vcf.gz" > controls.list

java -Xmx3g -Djava.io.tmpdir=. -jar scansv.jar --controls controls.list -d2 25 --fraction 0.6 cases1.vcf cases2.vcf > out.vcf

```


