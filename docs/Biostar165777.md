# Biostar165777

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Split a XML file


## Usage

```
Usage: biostar165777 [options] Files
  Options:
    -N, --count
      Number of files to be created
      Default: 100
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
  * -o, --output
      Output file. Must contains __SPLIT__
    -T, --tag
      XML tag to be split.e.g 'Hit' in blast
    --version
      print version and exit

```


## Keywords

 * xml


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew biostar165777
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar165777.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar165777.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **biostar165777** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```bash
$  java -jar dist-1.139/biostar165777.jar -o out__SPLIT__.xml -T Hit -N 5 ~/blastn.xml

$ ls -la ~/blastn.xml out*.xml
-rw-rw-r-- 1 lindenb lindenb 422606 nov.  14 12:47 /home/lindenb/blastn.xml
-rw-rw-r-- 1 lindenb lindenb  86319 nov.  14 16:17 out001.xml
-rw-rw-r-- 1 lindenb lindenb  83570 nov.  14 16:17 out002.xml
-rw-rw-r-- 1 lindenb lindenb  85096 nov.  14 16:17 out003.xml
-rw-rw-r-- 1 lindenb lindenb  88297 nov.  14 16:17 out004.xml
-rw-rw-r-- 1 lindenb lindenb  87123 nov.  14 16:17 out005.xml


$ grep -cF "<Hit>" ~/blastn.xml out*.xml
/home/lindenb/blastn.xml:100
out001.xml:20
out002.xml:20
out003.xml:20
out004.xml:20
out005.xml:20

```
