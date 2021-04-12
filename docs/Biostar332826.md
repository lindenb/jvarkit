# Biostar332826

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Fast Extraction of Variants from a list of IDs


## Usage

```
Usage: biostar332826 [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    -d, --delete
      When found , remove the ID from the list of identifiers unless it's a 
      '.'. Should be faster but don't use it if two variants have the same ID.
      Default: false
    -f, --filter
      if not blank soft filter the variants that are NOT in the list. If 
      '--inverse' is specified then soft-filter the variants IN the list.
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -r, -i, --ids
      A list of identifiers, one per line
    --inverse
      Inverse: don't print the variants containing the IDS.
      Default: false
    -o, --out
      Output file. Optional . Default: stdout
    --version
      print version and exit
    -R, -I
      A semicolon/comma/space separated list of identifiers
      Default: <empty string>

```


## Keywords

 * vcf
 * rs
 * id



## See also in Biostars

 * [https://www.biostars.org/p/332826](https://www.biostars.org/p/332826)
 * [https://www.biostars.org/p/433062](https://www.biostars.org/p/433062)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew biostar332826
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20180817

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar332826.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar332826.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/biostar/Biostar332826Test.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/biostar/Biostar332826Test.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **biostar332826** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

 
 ## Example
 
 ```
 $ wget -O - -q "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" |\
 	gunzip -c |\
 	java -jar dist/biostar332826 --ids ids.txt > out.vcf 
 ```
 
