# KnimeToText

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

converts a Knime Workflow to a html representation.


## Usage

```
Usage: knime2txt [options] Files
  Options:
    -g, --dot
      graphiz dot file.
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * knime
 * workflow
 * convert


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew knime2txt
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/KnimeToText.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/KnimeToText.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **knime2txt** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

 
 ## Example
 
```
 $ java -jar dist/knime2txt.jar -g out.dot  knime_3.3.2/workspace/TEST/  | xmllint --format - > out.html
 [INFO][KnimeToText]read XML: knime_3.3.2/workspace/TEST/workflow.knime
[INFO][KnimeToText]read XML: knime_3.3.2/workspace/TEST/Interactive Table (#3)/settings.xml
[INFO][KnimeToText]read XML: knime_3.3.2/workspace/TEST/File Reader (#4)/settings.xml
[INFO][KnimeToText]read XML: knime_3.3.2/workspace/TEST/Java Snippet Row Filter (#5)/settings.xml
[INFO][KnimeToText]read XML: knime_3.3.2/workspace/TEST/Extract Column Header (#6)/settings.xml
[INFO][KnimeToText]read XML: knime_3.3.2/workspace/TEST/Unpivoting (#7)/settings.xml
[INFO][KnimeToText]read XML: knime_3.3.2/workspace/TEST/Java Snippet _simple_ (#8)/settings.xml
[INFO][KnimeToText]read XML: knime_3.3.2/workspace/TEST/Column Filter (#10)/settings.xml
[INFO][KnimeToText]read XML: knime_3.3.2/workspace/TEST/CSV Writer (#9)/settings.xml
```
 
```
$ xmllint out.html  | head
<?xml version="1.0"?>
<div>
  <div>
    <p>Interactive Table (#3)/settings.xml</p>
    <ul>
      <li><b>factory</b>  value:"<i>org.knime.base.node.viz.table.TableNodeFactory</i>"</li>
      <li><b>customDescription</b>  value:"<i/>"</li>
    </ul>
  </div>
  <div>
(...)
```
 
 
 
 
