# VcfToSvg

write a vcf to svg , with gene context


## Usage

```
Usage: vcf2svg [options] Files
  Options:
    --alphaFILTER
      Variant having FILTER!=PASS opacity (0== hide FILTERED variants)
      Default: 1.0
    --alphaINDEL
      Variant INDEL opacity (0== hide INDEL variants)
      Default: 1.0
    --exon, --exons
      Only keep variants in exons
      Default: false
    -gw, --genotypeWidth
      Genotype square width
      Default: 10
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
  * -k, --knownGenes
      UCSC knownGene URI. Beware chromosome names are formatted the same as 
      your REFERENCE. A typical KnownGene file is 
      http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz 
      .If you only have a gff file, you can try to generate a knownGene file 
      with [http://lindenb.github.io/jvarkit/Gff2KnownGene.html](http://lindenb.github.io/jvarkit/Gff2KnownGene.html)
    -m, --manifest
      Manifest file containing the names of the files.
    --nonCoding
      Ignore Non-coding genes
      Default: false
    -o, --out
      Output SVG file. If defined, MUST Contains the word '__SEGMENT__'
    --stopAfterFirst
      when writing multiple SVG docs, stop after the first one. It avoids 
      writing multiple concatenated SVG documents when writing to stdout
      Default: false
    -trim2ctx, --trimToVariant
      Don't use gene interval for graphics but min/max of variants
      Default: false
    --version
      print version and exit

```


## Keywords

 * vcf
 * svg
 * xlm
 * visualization


## Compilation

### Requirements / Dependencies

* java compiler SDK 1.8 http://www.oracle.com/technetwork/java/index.html (**NOT the old java 1.7 or 1.6**) . Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )
* GNU Make >= 3.81
* curl/wget
* git
* xsltproc http://xmlsoft.org/XSLT/xsltproc2.html (tested with "libxml 20706, libxslt 10126 and libexslt 815")


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ make vcf2svg
```

The *.jar libraries are not included in the main jar file, so you shouldn't move them (https://github.com/lindenb/jvarkit/issues/15#issuecomment-140099011 ).
The required libraries will be downloaded and installed in the `dist` directory.

### edit 'local.mk' (optional)

The a file **local.mk** can be created edited to override/add some definitions.

For example it can be used to set the HTTP proxy:

```
http.proxy.host=your.host.com
http.proxy.port=124567
```
## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfToSvg.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfToSvg.java
)
## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcf2svg** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
# download and gzip refGene from UCSC
$ curl -s "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/refGene.txt.gz" |\
   gunzip -c | LC_ALL=C sort -t '	' -k3,3 -k5,5n | bgzip > refGene.txt.gz

# index the refGene file
$ tabix  -0 -s  3 -b 5 -e 6 -f refGene.txt.gz

# run vcf2svg 
$ java -jar dist/vcf2svg.jar \
   --stopAfterFirst \
   -k refGene.txt.gz input.vcf > out.svg
```

## Screenshot

https://twitter.com/yokofakun/status/851875435948462080

![screenshot](https://pbs.twimg.com/media/C9J4LeoXkAEqvIN.jpg)





