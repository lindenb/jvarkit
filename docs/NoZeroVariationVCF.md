# NoZeroVariationVCF

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

cat a whole VCF, or, if there is no variant, creates a fake one


## Usage

```
Usage: nozerovariationvcf [options] Files
  Options:
    -f, --filter
      FILTER name
      Default: FAKESNP
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    -R, -r, --reference
      The parameter is the path to an Indexed fasta Reference file. This fasta 
      file must be indexed with samtools faidx and with picard 
      CreateSequenceDictionary. The parameter can also be a 'key' (matching 
      the regular expression `[A-Za-z][A-Za-z0-9_\\-]*`) in a catalog file. A 
      'catalog' file is a java property file ( 
      https://docs.oracle.com/javase/tutorial/essential/environment/properties.html 
      ) where the values are the path to the fasta file.  Catalogs are 
      searched in that order : `${PWD}/fasta-ref.properties`, 
      `${HOME}/.fasta-ref.properties`, `/etc/jvarkit/fasta-ref.properties`.  
      If the key or the path are not defined by the user, they will be 
      searched in that order 1) the java property 
      -Djvarkit.fasta.reference=pathTofastaOrCatalogKey . 2) the linux 
      environement variable $FASTA_REFERENCE=pathTofastaOrCatalogKey 3) The 
      catalogs. 
      Default: <<Default Fasta Reference Supplier>>
    --version
      print version and exit

```


## Keywords

 * vcf


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew nozerovariationvcf
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/NoZeroVariationVCF.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/NoZeroVariationVCF.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/NoZeroVariationVCFTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/NoZeroVariationVCFTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **nozerovariationvcf** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```bash
# we use grep to create an empty VCF
$ gunzip -c file.vcf.gz | \
 grep  "#" |\
 java -jar dist/nozerovariationvcf.jar -r human_g1k_v37.fasta

##fileformat=VCFv4.1
##FILTER=<ID=FAKESNP,Description="Fake SNP created because vcf input was empty. See https://github.com/lindenb/jvarkit">
(...)
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1
GL000207.1	1	.	C	A	1	FAKESNP	.	GT:DP:GQ	0/1:1:1
```

