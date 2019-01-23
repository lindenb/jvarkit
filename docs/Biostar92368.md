# Biostar92368

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Binary interactions depth.


## Usage

```
Usage: biostar92368 [options] Files
  Options:
  * -D, --bdbhome
      berkeleydb home
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -M, --maxdepth
      Max depth
      Default: 3
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * protein
 * interaction
 * interactome



## See also in Biostars

 * [https://www.biostars.org/p/92368](https://www.biostars.org/p/92368)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew biostar92368
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar92368.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar92368.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **biostar92368** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


##Example

```bash
$ cat input.txt
A1	A2
A2	A3
A3	A4
A4	A5
A1	A6


$ mkdir -p tmp
$  java -jar dist/biostar92368.jar -D tmp input.txt

A1	A2	0
A1	A3	1
A1	A4	2
A1	A6	0
A2	A3	0
A2	A4	1
A2	A5	2
A2	A6	1
A3	A4	0
A3	A5	1
A3	A6	2
A4	A5	0

```
