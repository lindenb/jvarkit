# VcfAfInfoFilter

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Filter VCF annotated with external (AF or AC/AN) frequency information like vcfgnomad


## Usage

```
Usage: vcfafinfofilter [options] Files
  Options:
    --disable-vc-attribute-recalc
      When genotypes are removed/changed, Dd not recalculate variant 
      attributes like DP, AF, AC, AN...
      Default: false
    -F, --fields
      [20180905]How to extract the AlleleFrequencies from a variant. Multiple 
      separated with comma or semicolon. e.g: 
      "AC/AN;exome_CEU_*;genome_NFE_AF;another_AC/another/AN". Input is a set 
      of AC/AN field pairs or/and AF field separated by semicolon. 'x/y' means 
      AC/AN fields. '*' will be replaced with AC and AN, hence, 'exome_CEU_*' 
      will be interpreted as exome_CEU_AC/exome_CEU_AN. Other field will be 
      interpreted as an AF field.
      Default: <empty string>
    --filter, -f
      set this filter if all ALT fails the treshold. If empty :remove the 
      variant 
      Default: <empty string>
    --gtfilter, -gtf
      set this *GENOTYPE* filter if all ALT for a Genotype fail the treshold. 
      If empty :set genotype to NO_CALL
      Default: <empty string>
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -nfe, --nfe
      Add INFO fields for the 'NFE' population created by vcfgnomad: gnomad_exome_AC_NFE,gnomad_exome_AF_NFE,gnomad_exome_AN_NFE,gnomad_genome_AC_NFE,gnomad_genome_AF_NFE,gnomad_genome_AN_NF
      Default: false
    -i, --no-valid
      Ignore INFO Field Validation. (e.g INFO field not declarated in VCF 
      header) 
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    --treshold, -t
      Treshold for allele Frequency. ALT alleles above this AF value will be 
      subject to filtration.
      Default: 0.001
    --vc-attribute-recalc-ignore-filtered
      When recalculating variant attributes like DP AF, AC, AN, ignore 
      FILTERed **Genotypes**
      Default: false
    --vc-attribute-recalc-ignore-missing
      Ignore missing VCF headers (DP, AF, AC, AN). Default behavior: adding 
      VCF header if they're missing
      Default: false
    --version
      print version and exit

```


## Keywords

 * vcf
 * annotation
 * gnomad


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfafinfofilter
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfAfInfoFilter.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfAfInfoFilter.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/VcfAfInfoFilterTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/VcfAfInfoFilterTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfafinfofilter** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Motivation

I'm often asked to filter out variant that are too frequent in gnomad, but I must keep the data if any ALT allele is NOT in gnomad.

This tool filters VCF containing external allele frequency information (AF or AC/AN). Used as a  complement of VcfGnomadPext.

## Example

```
$ java -jar dist/vcfafinfofilter.jar -nfe input.vcf
$ java -jar dist/vcfafinfofilter.jar -af 'gnomad_exome_AF_NFE,gnomad_genome_AF_NFE'   input.vcf
$ java -jar dist/vcfafinfofilter.jar -acn 'gnomad_genome_AC_NFE,gnomad_genome_AN_NFE'   input.vcf
$ java -jar dist/vcfafinfofilter.jar -acn 'gnomad_genome_*_NFE'   input.vcf

```



