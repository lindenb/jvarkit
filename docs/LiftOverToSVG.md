# LiftOverToSVG

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Convert LiftOver chain files to animated SVG


## Usage

```
Usage: liftover2svg [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --version
      print version and exit
    -b
      add this build name
      Default: []
    -t
       seconds per step
      Default: 20.0
    -w
      width
      Default: 1000

```


## Keywords

 * svg
 * liftover
 * ucsc
 * xml


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew liftover2svg
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/liftover/LiftOverToSVG.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/liftover/LiftOverToSVG.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **liftover2svg** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Screenshot

https://twitter.com/yokofakun/status/525014376387215360

![svg ](https://pbs.twimg.com/media/B0k5ytFCIAAiBPJ.png)

## Examples

### Examples 1

hg16 to hg38 http://cardioserve.nantes.inserm.fr/~lindenb/liftover2svg/hg16ToHg38.svg

BIG FILE : hg38 to panTro4 http://cardioserve.nantes.inserm.fr/~lindenb/liftover2svg/hg38ToPanTro4.svg



### Example 2

create a **SVG** file for **hg16** to **hg38**.

```makefile
.PHONY:all 
all: test.svg


test.svg: hg16ToHg17.over.chain hg17ToHg18.over.chain hg18ToHg19.over.chain hg19ToHg38.over.chain dist/liftover2svg.jar
	java -jar dist/liftover2svg.jar -b hg16 -b hg17 -b hg18 -b hg19 -b hg38 \
			 hg16ToHg17.over.chain hg17ToHg18.over.chain hg18ToHg19.over.chain hg19ToHg38.over.chain > $@

hg16ToHg17.over.chain : 
	curl -s "http://hgdownload-test.cse.ucsc.edu/goldenPath/hg16/liftOver/$@.gz" | gunzip -c |\
	grep chain > $@

hg17ToHg18.over.chain : 
	curl -s "http://hgdownload-test.cse.ucsc.edu/goldenPath/hg17/liftOver/$@.gz" | gunzip -c |\
	grep chain > $@

hg18ToHg19.over.chain : 
	curl -s "http://hgdownload-test.cse.ucsc.edu/goldenPath/hg18/liftOver/$@.gz" | gunzip -c |\
	grep chain > $@
	
hg38ToPanTro4.over.chain :
	curl -s "http://hgdownload-test.cse.ucsc.edu/goldenPath/hg38/liftOver/$@.gz" | gunzip -c |\
	grep chain > $@

hg19ToHg38.over.chain hg19ToPanTro4.over.chain :
	curl -s "http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/liftOver/$@.gz" | gunzip -c |\
	grep chain > $@


dist/liftover2svg.jar: ./src/main/java/com/github/lindenb/jvarkit/tools/liftover/LiftOverToSVG.java
	$(MAKE)  liftover2svg

```

