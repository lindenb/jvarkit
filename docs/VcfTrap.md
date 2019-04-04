# VcfTrap

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

annotate vcf with trap database http://trap-score.org/


## Usage

```
Usage: vcftrap [options] Files
  Options:
    -A, --attribute
      VCF INFO attribute Format:(ALT|GENE|SCORE)
      Default: TRAP
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --ignore-filtered
      Ignore FILTERed variants (faster)
      Default: false
  * -m, --manifest
      Manifest file. A tab delimited file with two columns : 
      chromosome(tab)path-to-file-indexed-with-trapindex. 
    -o, --out
      Output file. Optional . Default: stdout
    --outputbcf
      Output bcf (for streams)
      Default: false
    --vcfcreateindex
      VCF, create tribble or tabix Index when writing a VCF/BCF to a file.
      Default: false
    --vcfmd5
      VCF, create MD5 checksum when writing a VCF/BCF to a file.
      Default: false
    --version
      print version and exit

```


## Keywords

 * vcf
 * trap
 * annotation


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcftrap
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/trap/VcfTrap.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/trap/VcfTrap.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcftrap** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example:

```
java -jar dist/trapindexer.jar  -o chr22.dat chr22.TraPv2.txt.gz
echo -e "22\tchr22.dat" > out.manifest
java -jar dist/vcftrap.jar -m out.manifest input.vcf
```


## See also

* TrapIndexer

