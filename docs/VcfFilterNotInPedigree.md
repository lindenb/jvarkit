# VcfFilterNotInPedigree

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Adds a FILTER 'NotInPedigree' if the only not(homref) genotypes are not in a pedigree


## Usage

```
Usage: vcffilternotinpedigree [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    -f, --filter
      FILTER name. Will be set for variant where the only genotypes non-homref 
      are NOT in the pedigree
      Default: NoGenotypeInPedigree
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -gtf, --ignore-filtered-gt
      [20180406] Do not consider a *genotype* if it is FILTERED.
      Default: false
    -o, --out
      Output file. Optional . Default: stdout
  * -p, --pedigree
      A pedigree file. tab delimited. Columns: family,id,father,mother, 
      sex:(0:unknown;1:male;2:female), phenotype 
      (-9|?|.:unknown;1|affected|case:affected;0|unaffected|control:unaffected) 
    -r, --remove
      remove the variant instead of setting the FILTER (hard filtering)
      Default: false
    -sf, --sfilter
      FILTER name for option --singleton
      Default: SingletonAlt
    -s, --singleton
      Variant is flagged/FILTERed as SingletonAlt if the ALT if found in less 
      or equal times 'singleton-times' in the genotypes. -1:ignore
      Default: 1
    --version
      print version and exit

```


## Keywords

 * burden
 * vcf
 * pedigree


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcffilternotinpedigree
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burden/VcfFilterNotInPedigree.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burden/VcfFilterNotInPedigree.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/burden/VcfFilterNotInPedigreeTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/burden/VcfFilterNotInPedigreeTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcffilternotinpedigree** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example



