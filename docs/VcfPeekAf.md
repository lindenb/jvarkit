# VcfPeekAf

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Peek the AF from another VCF


## Usage

```
Usage: vcfpeekaf [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    -b, --buffer-size
      When we're looking for variant in a lare VCF file, load the variants in 
      an interval of 'N' bases instead of doing a random access for each 
      variant. A distance specified as a positive integer.Commas are removed. 
      The following suffixes are interpreted : b,bp,k,kb,m,mb
      Default: 10000
  * -F, --database, --tabix, --resource
      An indexed VCF file. Source of the annotations
    -f, --filter
      soft FILTER the variant of this data if AF is not found or it greater > 
      max-af or lower than min-af. If empty, just DISCARD the variant
      Default: <empty string>
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -l, --list
      List available AF peekers and exit.
    --min-af
      AF min treshold. Variant is accepted is computed AF >= treshold.
      Default: 0.0
    --no-alt
      Do not look at the alternate alleles concordance
      Default: false
    -o, --out
      Output file. Optional . Default: stdout
    --peek-id
      Peek database variant ID if it is missing in the processed VCF.
      Default: false
    -P, --peek-info
      Name of INFO tag in the vcf database to extract the AF value for 
      exractor .'Custom'
  * -p, --peeker
      AF Peeker name. Use option --list to get a list of peekers.
    -T, --tag
      INFO tag to put found frequency. empty: no extra tag.
      Default: <empty string>
  * -t, --treshold, --max-af
      AF max treshold. Variant is accepted is computed AF <= treshold.
      Default: 1.0
    --version
      print version and exit

```


## Keywords

 * vcf
 * annotation
 * af


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfpeekaf
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20200624

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfvcf/VcfPeekAf.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfvcf/VcfPeekAf.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfvcf/VcfPeekAfTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfvcf/VcfPeekAfTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfpeekaf** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

ALLELE_FREQUENCY_KEY
