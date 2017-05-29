# Biostar172515

Convert BAI to XML


## Usage

```
Usage: biostar172515 [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * bai
 * bam
 * xml



## See also in Biostars

 * [https://www.biostars.org/p/172515](https://www.biostars.org/p/172515)


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
$ make biostar172515
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar172515.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar172515.java
)
## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **biostar172515** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)




### Example




```
 
 find DIR -name "*.bam" | xargs java -jar dist/biostar172515.jar  | xmllint --format -

<?xml version="1.0" encoding="UTF-8"?>
<bai-list>
<bam bam="DIR/exampleBAM.bam" has-index="true" n_ref="1">
    <reference ref-id="0" ref-name="chr1" ref-length="100000" n_aligned="33" n_bin="12" n_no_coor="0">
      <bin first-locus="1" last-locus="536870912" level="0" first-offset="0" n_chunk="0"/>
      <bin first-locus="1" last-locus="67108864" level="1" first-offset="0" n_chunk="0"/>
      <bin first-locus="1" last-locus="8388608" level="2" first-offset="0" n_chunk="0"/>
      <bin first-locus="1" last-locus="1048576" level="3" first-offset="0" n_chunk="0"/>
      <bin first-locus="1" last-locus="131072" level="4" first-offset="0" n_chunk="0"/>
      <bin first-locus="1" last-locus="16384" level="5" first-offset="828" n_chunk="1">
        <chunk chunk_beg="828" chunk_end="1963"/>
      </bin>
      <bin first-locus="16385" last-locus="32768" level="5" first-offset="1963" n_chunk="1">
        <chunk chunk_beg="1963" chunk_end="3323"/>
      </bin>
      <bin first-locus="32769" last-locus="49152" level="5" first-offset="3323" n_chunk="1">
        <chunk chunk_beg="3323" chunk_end="4687"/>
      </bin>
      <bin first-locus="49153" last-locus="65536" level="5" first-offset="4687" n_chunk="1">
        <chunk chunk_beg="4687" chunk_end="6501"/>
      </bin>
      <bin first-locus="65537" last-locus="81920" level="5" first-offset="0" n_chunk="0"/>
      <bin first-locus="81921" last-locus="98304" level="5" first-offset="6501" n_chunk="1">
        <chunk chunk_beg="6501" chunk_end="238223360"/>
      </bin>
      <bin first-locus="98305" last-locus="114688" level="5" first-offset="0" n_chunk="0"/>
    </reference>
  </bam>
  <bam bam="/DIR/exampleBAM2.bam" has-index="true" n_ref="1">
    <reference ref-id="0" ref-name="chr1" ref-length="100000" n_aligned="33" n_bin="12" n_no_coor="0">
      <bin first-locus="1" last-locus="536870912" level="0" first-offset="0" n_chunk="0"/>
      <bin first-locus="1" last-locus="67108864" level="1" first-offset="0" n_chunk="0"/>
      <bin first-locus="1" last-locus="8388608" level="2" first-offset="0" n_chunk="0"/>
      <bin first-locus="1" last-locus="1048576" level="3" first-offset="0" n_chunk="0"/>
      <bin first-locus="1" last-locus="131072" level="4" first-offset="0" n_chunk="0"/>
      <bin first-locus="1" last-locus="16384" level="5" first-offset="828" n_chunk="1">
        <chunk chunk_beg="828" chunk_end="1963"/>
      </bin>
      <bin first-locus="16385" last-locus="32768" level="5" first-offset="1963" n_chunk="1">
        <chunk chunk_beg="1963" chunk_end="3323"/>
      </bin>
      <bin first-locus="32769" last-locus="49152" level="5" first-offset="3323" n_chunk="1">
        <chunk chunk_beg="3323" chunk_end="4687"/>
      </bin>
      <bin first-locus="49153" last-locus="65536" level="5" first-offset="4687" n_chunk="1">
        <chunk chunk_beg="4687" chunk_end="6501"/>
      </bin>
      <bin first-locus="65537" last-locus="81920" level="5" first-offset="0" n_chunk="0"/>
      <bin first-locus="81921" last-locus="98304" level="5" first-offset="6501" n_chunk="1">
        <chunk chunk_beg="6501" chunk_end="238223360"/>
      </bin>
      <bin first-locus="98305" last-locus="114688" level="5" first-offset="0" n_chunk="0"/>
    </reference>
  </bam>
</bai-list>
```


