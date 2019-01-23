# GexfTransformer

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Gexf file manipulation


## Usage

```
Usage: gexftr [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -c, --max-cluster
      limit to this number of clusters. Negative: no limit
      Default: 1
    -M, --max-nodes
      max number of nodes per cluster inclusive. Negative: no limit
      Default: -1
    -m, --min-nodes
      min number of nodes per cluster inclusive. Negative: no limit
      Default: -1
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * gexf
 * graph
 * network


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew gexftr
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/gephi/GexfTransformer.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/gephi/GexfTransformer.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **gexftr** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Cluster 

a **Cluster** is a set of nodes that are all connected with some edges.

## Example:

```
$ java -jar dist/pubmeddump.jar "NSP1 Rotavirus" | \
  java -jar dist/pubmedauthorgraph.jar -D BDB |\
  java -jar dist/gexftr.jar | xmllint --format - 
  
(...)
      <!--Skip edge  from HOWARD~C_R to pmid:9614866-->
      <!--Skip edge  from BRIDGER~J_C to pmid:9614866-->
      <!--Skip edge  from LÓPEZ~S to pmid:9645203-->
      <!--Skip edge  from GONZÁLEZ~R_A to pmid:9645203-->
      <!--Skip edge  from ARIAS~C_F to pmid:9645203-->
      <!--Skip edge  from TORRES-VEGA~M_A to pmid:9645203-->
    </edges>
  </graph>
  <!--Number of nodes 692-->
</gexf>


```


