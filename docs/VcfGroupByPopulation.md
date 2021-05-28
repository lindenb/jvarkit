# VcfGroupByPopulation

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Group VCF data by population, creates a VCF  where each 'SAMPLE' is a population


## Usage

```
Usage: vcfgroupbypop [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
  * -p, --mapping
      mapping file: each line is (SAMPLE)\t(POP)\n
    -M, --max-fisher
      max inclusive value of fisher test. A decimal number between 0.0 and 
      1.0. If the value ends with '%' it is interpretted as a percentage eg. 
      '1%' => '0.01'. A slash '/' is interpretted as a ratio. e.g: '1/100' => 
      '0.01'. 
      Default: 1.0
    -m, --min-fisher
      min inclusive value of fisher test. A decimal number between 0.0 and 
      1.0. If the value ends with '%' it is interpretted as a percentage eg. 
      '1%' => '0.01'. A slash '/' is interpretted as a ratio. e.g: '1/100' => 
      '0.01'. 
      Default: 0.0
    -o, --out
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * vcf
 * pedigree
 * population


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfgroupbypop
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20190319

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfGroupByPopulation.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfGroupByPopulation.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfgroupbypop** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

