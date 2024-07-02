# MergeBlastXml

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

merge XML blast results (same Iteration/Iteration_query-def in multiple xml files


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar mergeblastxml  [options] Files

Usage: mergeblastxml [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --maxRecordsInRam
      Max Records in RAM
      Default: 50000
    -o, --out
      Output SVG file or stdout
    --tmpDir
      Tmp Directory
      Default: /tmp
    --version
      print version and exit

```


## Keywords

 * blast
 * xml



## See also in Biostars

 * [https://www.biostars.org/p/246958](https://www.biostars.org/p/246958)


## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/blast/MergeBlastXml.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/blast/MergeBlastXml.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **mergeblastxml** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
$ java -jar dist/jvarkit.jar mergeblastxml input1.blastn.xml  input2.blastn.xml  input2.blastn.xml > out.xml
``` 



