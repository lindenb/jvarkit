# VcfInjectPedigree

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Injects a pedigree (.ped) file in the VCF header


## Usage

```
Usage: vcfinjectpedigree [options] Files
  Options:
    -clean, --clean
      Remove all previous data about pedigree in the VCF header before adding 
      the new one.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -imih, --ignoreMissingInHeader
      Ignore errors if a sample is declared in the pedigree but is missing in 
      the VCF header
      Default: false
    -imip, --ignoreMissingInPedigree
      Ignore errors if a sample is declared in the VCF header but is missing 
      in the pedigree
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    -p, --pedigree
      A pedigree is a text file delimited with tabs. No header. Columns are 
      (1) Family (2) Individual-ID (3) Father Id or '0' (4) Mother Id or '0' 
      (5) Sex : 1 male/2 female / 0 unknown (6) Status : 0 unaffected, 1 
      affected,-9 unknown
    -valid, --valid
      Ignore pedigree validation
      Default: false
    --version
      print version and exit

```


## Keywords

 * vcf
 * pedigree
 * burden


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfinjectpedigree
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burden/VcfInjectPedigree.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burden/VcfInjectPedigree.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfinjectpedigree** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



This tools reads a pedigree file and inject it in the VCF header  



```
$ java -jar dist/vcfinjectpedigree.jar \
	-imih -imip -p input.ped \
	input.vcf.gz > out.vcf

$ grep Sample out.vcf
(...)
##Sample=<Family=F1,ID=INDI1,Father=0,Mother=0,Sex=1,Status=1>
##Sample=<Family=F2,ID=INDI2,Father=0,Mother=0,Sex=2,Status=1>
##Sample=<Family=F3,ID=INDI3,Father=INDI1,Mother=INDI2,Sex=1,Status=1>
(...)

```






This tools reads a pedigree file and inject it in the VCF header  


```
$ java -jar dist/vcfinjectpedigree.jar \
	-imih -imip -p input.ped \
	input.vcf.gz > out.vcf

$ grep Sample out.vcf
(...)
##Sample=<Family=F1,ID=INDI1,Father=0,Mother=0,Sex=1,Status=1>
##Sample=<Family=F2,ID=INDI2,Father=0,Mother=0,Sex=2,Status=1>
##Sample=<Family=F3,ID=INDI3,Father=INDI1,Mother=INDI2,Sex=1,Status=1>
(...)

```

