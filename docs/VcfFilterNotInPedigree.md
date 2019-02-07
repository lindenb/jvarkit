# VcfFilterNotInPedigree

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Adds a FILTER 'NotInPedigree' if the only not(homref) genotypes are not in a pedigree


## Usage

```
Usage: vcffilternotinpedigree [options] Files
  Options:
    -f, --filter
      FILTER name. Will be set for variant where the only genotypes non-homref 
      are NOT in the pedigree
      Default: NoGenotypeInPedigree
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -gtf, --ignore-filtered-gt
      [20180406] Do not consider a *genotype* if it is FILTERED.
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    -p, --pedigree
      [20180406]A pedigree is a text file delimited with tabs. No header. 
      Columns are (1) Family (2) Individual-ID (3) Father Id or '0' (4) Mother 
      Id or '0' (5) Sex : 1 male/2 female / 0 unknown (6) Status : 0 
      unaffected, 1 affected,-9 unknown  Default is to try to read the 
      pedigree in the VCF header
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



