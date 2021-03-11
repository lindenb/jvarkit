# VcfStrechToSvg

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

another VCF to SVG


## Usage

```
Usage: vcfstrech2svg [options] Files
  Options:
    --af
      How to extract the AlleleFrequencies from a variant. Multiple separated 
      with comma or semicolon. e.g: 
      "AC/AN;exome_CEU_*;genome_NFE_AF;another_AC/another/AN". Input is a set 
      of AC/AN field pairs or/and AF field separated by semicolon. The special 
      value 'FORMAT/GT' will re-compute the Frequencies from the 
      genotypes.'x/y' means AC/AN fields. '*' will be replaced with AC and AN, 
      hence, 'exome_CEU_*' will be interpreted as exome_CEU_AC/exome_CEU_AN. 
      Other field will be interpreted as an AF field.
      Default: FORMAT/GT
    --bam-list, --bam, --bams
      Add BAM/CRAM for plotting the depth. A file with the suffix '.list' is a 
      list of path to the bams/crams.
      Default: []
    --color-tag
      specify the optional INFO/tag defining a named svg color. If defined for 
      a variant, a vertical line with the color will be painted.
    --dp
      minimum sum(FORMAT/AD[0]+FORMAT/AD[1]). Doesn't work with FORMAT/PL.
      Default: 1
    --extend
      Extend each area with 'x' bp. A distance specified as a positive 
      integer.Commas are removed. The following suffixes are interpreted : 
      b,bp,k,kb,m,mb 
      Default: 1000
    --format
      wich format to use to calculate the allele depth ratio.
      Default: AD
      Possible Values: [PL, AD]
    --gq
      minimum FORMAT/GQ
      Default: 1
    --gtf
      Plot gene structure using this GTF file.
    --gzip
      Generate gzipped compressed svg files.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -hr, --hom-ref
      Hide HOM_REF genotypes (0/0)
      Default: false
    -u, --url, --hyperlink
      creates a hyperlink an area is clicked. creates a hyperlink when 'click' 
      in an area. The URL must contains __CHROM__, __START__ and __END__ that 
      will be replaced by their values. Predefined values are 
      'hg19','hg38','igv'. IGV : 
      "http://localhost:60151/goto?locus=__CHROM__%3A__START__-__END__" , 
      UCSC: "http://genome.ucsc.edu/cgi-bin/hgTracks?org=Human&db=hg19&position=__CHROM__%3A__START__-__END__"
      Default: <empty string>
    --keep-filtered
      keep FILTERed variants
      Default: false
    --manifest
      Output BED manifest
    --max-af
      Discard variant with an internal AF > 'x' A decimal number between 0.0 
      and 1.0. If the value ends with '%' it is interpretted as a percentage 
      eg. '1%' => '0.01'. A slash '/' is interpretted as a ratio. e.g: '1/100' 
      => '0.01'.
      Default: 1.0
    --min-af
      Discard variant with an internal AF < 'x' A decimal number between 0.0 
      and 1.0. If the value ends with '%' it is interpretted as a percentage 
      eg. '1%' => '0.01'. A slash '/' is interpretted as a ratio. e.g: '1/100' 
      => '0.01'.
      Default: 0.0
    --no-tooltip
      remove contextual tooltips (reduce the size of the svg)
      Default: false
  * -o, --output
      An existing directory or a filename ending with the '.zip' or '.tar' or 
      '.tar.gz' suffix.
    --pack-distance
      pack variant in the same area if they're close to 'x' bp. A distance 
      specified as a positive integer.Commas are removed. The following 
      suffixes are interpreted : b,bp,k,kb,m,mb
      Default: 10000
    --param
      Other parameters. Undocumented
      Syntax: --paramkey=value
      Default: {gt.r1=1, gt.r2=7, sample.height=50}
    --reference, -R
      For reading/writing CRAM files. Indexed fasta Reference file. This file 
      must be indexed with samtools faidx and with picard 
      CreateSequenceDictionary 
  * -r, --region, --bed
      BED File
    --samples
      Limit to those samples. (comma or space separated). If the list starst 
      with '^' then the samples are excluded.
    --version
      print version and exit
    -w, --width
      image width.
      Default: 1000

```


## Keywords

 * vcf
 * deletion
 * cnv
 * svg


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfstrech2svg
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20210304

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/VcfStrechToSvg.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/VcfStrechToSvg.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfstrech2svg** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Motivation

another VCF to SVG

plot  allele ratio for each sample for each diallelic variant in a given region

colors define the genotype type : HOM_REF, HET, HOM_VAR

the size the dots represent the inverse internal allele frequency.

## Example

```
java -jar dist/vcfstrech2svg.jar --bed intervals.bed  -o TMP indexed.vcf.gz
```

## Screenshot

![twitter](https://pbs.twimg.com/media/EvtRyyXWEAEqpz7?format=jpg&name=small "Screenshot")

[https://twitter.com/yokofakun/status/1367778079813341185](https://twitter.com/yokofakun/status/1367778079813341185)


