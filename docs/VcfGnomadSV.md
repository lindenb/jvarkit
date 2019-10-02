# VcfGnomadSV

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Peek annotations from gnomad structural variants


## Usage

```
Usage: vcfgnomadsv [options] Files
  Options:
    --any-overlap-filter
      If not empty, set this FILTER if any variant in gnomad is found 
      overlaping the variant BUT we didn't find a correct match
      Default: <empty string>
    --bnd-distance
      Two BND variants are the same if their bounds are distant by less than 
      xxx bases. A distance specified as a positive integer.Commas are 
      removed. The following suffixes are interpreted : b,bp,k,kb,m,mb
      Default: 100
    --check-bnd-mate
      When comparing two BND, check that their mate (using the ALT allele) are 
      the same too
      Default: false
    --discordant_svtype
      If not empty, set this FILTER if SVTYPE are discordants
      Default: <empty string>
    --force-svtype
      When comparing two SV variants, their INFO/SVTYPE should be the same. 
      Default is to just use coordinates to compare non-BND variants.
      Default: false
  * -g, --gnomad
      Gnomad-SV VCF file. see 
      https://gnomad.broadinstitute.org/downloads#structural-variants 
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --in-gnomad-filter
      If not empty, set this FILTER is variant was found in gnomad
      Default: <empty string>
    -o, --output
      Output file. Optional . Default: stdout
    -p, --prefix
      INFO field prefix
      Default: GNOMAD_
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

 * vcf
 * annotation
 * gnomad
 * sv


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfgnomadsv
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20190814

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/gnomad/VcfGnomadSV.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/gnomad/VcfGnomadSV.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/gnomad/VcfGnomadSVTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/gnomad/VcfGnomadSVTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfgnomadsv** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


# Example:

```
java -jar dist/vcfgnomadsv.jar \
	-g src/test/resources/gnomad_v2_sv.sites.vcf.gz \
	./src/test/resources/manta.B00GWGD.vcf.gz
```

