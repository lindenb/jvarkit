# WorldMapGenome

Genome/Map of the world. Input is a BED file: chrom/start/end/country.


## Usage

```
Usage: worldmapgenome [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    --version
      print version and exit
    -R
       (ref) fasta reference file indexed with samtools. Required.
    -T
       just list countries and exit.
    -o
       (file.jpg) ouput file. Required
    -u
       (uri) URL world SVG map.
      Default: http://upload.wikimedia.org/wikipedia/commons/7/76/World_V2.0.svg
    -w
       (int: size) square size
      Default: 1000

```


## Keywords

 * gis


## Compilation

### Requirements / Dependencies

* java compiler SDK 1.8 http://www.oracle.com/technetwork/java/index.html (**NOT the old java 1.7 or 1.6**) . Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )
* GNU Make >= 3.81
* curl/wget
* git
* xsltproc http://xmlsoft.org/XSLT/xsltproc2.html (tested with "libxml 20706, libxslt 10126 and libexslt 815")


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ make worldmapgenome
```

The *.jar libraries are not included in the main jar file, so you shouldn't move them (https://github.com/lindenb/jvarkit/issues/15#issuecomment-140099011 ).
The required libraries will be downloaded and installed in the `dist` directory.

### edit 'local.mk' (optional)

The a file **local.mk** can be created edited to override/add some definitions.

For example it can be used to set the HTTP proxy:

```
http.proxy.host=your.host.com
http.proxy.port=124567
```
## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/circular/WorldMapGenome.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/circular/WorldMapGenome.java
)
## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **worldmapgenome** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
head map.bed

1	141535	143008	taiwan
1	564463	570304	china
1	564463	570304	china
1	564463	564813	canada
1	564510	569779	estonia
1	564510	569779	estonia
1	564510	569170	germany
1	564633	569139	italy
1	564660	565200	iran
1	564733	569130	brazil
1	564776	569431	mexico
```


```bash
$  cat map.bed |\
     java -jar dist/worldmapgenome.jar \
     -u World_V2.0.svg \
     -w 2000 -o ~/ouput.jpg \
     -R human_g1k_v37.fasta

```
![worldmapgenome](https://pbs.twimg.com/media/BfGE0X4CMAAfRAR.jpg)

## See also

* https://twitter.com/yokofakun/status/428269474869817344
* https://twitter.com/GenomeBrowser/status/426398651103997953




