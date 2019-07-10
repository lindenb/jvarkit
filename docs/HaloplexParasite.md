# HaloplexParasite

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

for @SolenaLS : remove artctifacts from haloplex that gives indels in GATK hapcaller 


## Usage

```
Usage: haloplexparasite [options] Files
  Options:
    -B, --bams
      A list of path to indexed BAM files
    -m, --clipsize
      Min. Soft Clipping size
      Default: 10
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    -t, --treshold
      treshold
      Default: 5.0E-4
    --version
      print version and exit

```


## Keywords

 * vcf
 * haloplex


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew haloplexparasite
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/haloplex/HaloplexParasite.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/haloplex/HaloplexParasite.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **haloplexparasite** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


### Examples

```

echo "input.bam" > all.list
gunzip -c input.vcf.gz |
  java -jar dist/haloplexparasite.jar -B all.list
rm all.list


```



### Examples

```

echo "input.bam" > all.list
gunzip -c input.vcf.gz |
  java -jar dist/haloplexparasite.jar -B all.list
rm all.list



```




