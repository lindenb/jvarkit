# VCFBedSetFilter

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Set FILTER for VCF if intersects with BED.


## Usage

```
Usage: vcfbedsetfilter [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    -e, --exclude, --blacklist
      Tribble or Tabix bed file containing the regions to be FILTERED. Must be 
      indexed with tribble or tabix, or use '--fast' to load in memory.
    -x, --extend
      Extend the variant coordinates per 'x' bases. A distance specified as a 
      positive integer.Commas are removed. The following suffixes are 
      interpreted : b,bp,k,kb,m,mb
      Default: 0
    -f, --filter
      FILTER name. If `--filter` is empty, FILTERED variant will be discarded.
      Default: VCFBED
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --fast, --memory
      Load the bed in memory: faster than tribble/tabix but memory consumming)
      Default: false
    --min-bed-fraction
      Min BED fraction overlap after extension. Only consider BED records if 
      variant overlap >= 'x' percent of bed length. A decimal number between 
      0.0 and 1.0. If the value ends with '%' it is interpretted as a 
      percentage eg. '1%' => '0.01'. A slash '/' is interpretted as a ratio. 
      e.g: '1/100' => '0.01'.
    --min-vc-fraction
      Min Variant fraction overlap after extension. Only consider BED records 
      if bed overlap >= 'x' percent of vc length. A decimal number between 0.0 
      and 1.0. If the value ends with '%' it is interpretted as a percentage 
      eg. '1%' => '0.01'. A slash '/' is interpretted as a ratio. e.g: '1/100' 
      => '0.01'.
    -o, --out
      Output file. Optional . Default: stdout
    --version
      print version and exit
    -i, --include, --whitelist
      Tribble or Tabix bed file containing the regions to be accepted. Regions 
      NOT overlapping those regions will be FILTERED. Must be indexed with 
      tribble or tabix, or use '--fast' to load in memory.

```


## Keywords

 * vcf
 * bed
 * filter



## See also in Biostars

 * [https://www.biostars.org/p/9465226](https://www.biostars.org/p/9465226)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfbedsetfilter
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20150415

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfbed/VCFBedSetFilter.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfbed/VCFBedSetFilter.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfbed/VCFBedSetFilterTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfbed/VCFBedSetFilterTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfbedsetfilter** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


### Examples

```
$java -jar dist/vcfbedsetfilter.jar -f MYFILTER - -B in.bed in.vcf 
```

## history:

2191104: changed the logic which was wrongly defined in the documentation. :-/

