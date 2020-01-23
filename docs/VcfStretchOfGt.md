# VcfStretchOfGt

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Try to finds deletion by searching strech of HOM_REF/HOM_VAR/NO_CALL Genotypes.


## Usage

```
Usage: vcfstrechofgt [options] Files
  Options:
    -a, --affected
      Same as option --pedigree but provide the name of the samples using a 
      comma/space/semicolon separated string
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -nc, --no-call
      Do not accept NO_CALL genotypes.
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    -p, --pedigree
      If defined, the tool will use the affected sample and find strech where 
      all affected could be a DEL.A pedigree file. tab delimited. Columns: 
      family,id,father,mother, sex:(0:unknown;1:male;2:female), phenotype 
      (-9|?|.:unknown;1|affected|case:affected;0|unaffected|control:unaffected) 
    --version
      print version and exit

```


## Keywords

 * vcf
 * deletion
 * cnv


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfstrechofgt
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20190103

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/VcfStretchOfGt.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/VcfStretchOfGt.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/structvar/VcfStretchOfGtTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/structvar/VcfStretchOfGtTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfstrechofgt** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
$ java -jar dist/vcfstrechofgt.jar -p src/test/resources/test_vcf01.ped src/test/resources/test_vcf01.vcf

#chrom	start0	end0	length	count.affected.variants	average.affected.depth	count.other.variants
1	870316	870317	1	1	5.0	0
1	919500	919501	1	1	0.0	0
1	963703	963704	1	1	0.0	0
1	1004201	1004202	1	1	0.0	0
```

## See also

```
bcftools roh
```

