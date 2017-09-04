# VcfBurdenSplitter

Split VCF Using a Burden Splitter (by gene, transcript, etc..)


## Usage

```
Usage: vcfburdensplitter [options] Files
  Options:
    -all_enst, --all_enst
      enable grouping by ALL_ENST: gene starting with ENST
      Default: false
    -all_filtered, --all_filtered
      If defined, the group where ALL the variants are FILTERED will be saved 
      here. 
    -all_genes, --all_genes
      enable grouping by all transcript for a gene using transcript name (e.g 
      PRKCB1) 
      Default: false
    -all_nm, --all_nm
      enable grouping by ALL_NM : gene not empty and transcript starting with 
      NM_ 
      Default: false
    -all_refseq, --all_refseq
      enable grouping by ALL_REFSEQ: gene not empty and transcript NOT 
      starting with ENST
      Default: false
    -all_transcripts, --all_transcripts
      enable grouping by all transcript for a gene using gene name (e.g 
      Nm_12345) 
      Default: false
    -gh, --galaxyhtml
      When used with galaxy, the files will be expanded in that path.
      Default: <empty string>
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -if, --ignorefilter
      accept variants having a FILTER column. Default is ignore variants with 
      a FILTER column
      Default: false
    -ls, --listsplitters
      List available splitters and exit
    --maxRecordsInRam
      When writing  files that need to be sorted, this will specify the number 
      of records stored in RAM before spilling to disk. Increasing this number 
      reduces the number of file  handles needed to sort a file, and increases 
      the amount of RAM needed
      Default: 50000
    -o, --output
      Output file. Optional . Default: stdout
    -sp, --splitterName
      Splitter Name
      Default: vepso
    --tmpDir
      tmp working directory. Default: java.io.tmpDir
      Default: []
    -vepEnsg, --vepEnsg
      enable VEP 'ENSG'
      Default: false
    -vepEnsp, --vepEnsp
      enable VEP 'ENSP'
      Default: false
    -vepEnst, --vepEnst
      enable VEP 'FEATURE' starting with 'ENST'
      Default: false
    -vepFeature, --vepFeature
      enable VEP 'FEATURE' (transcript)
      Default: false
    -vepHgnc, --vepHgnc
      enable VEP 'HGNC'
      Default: false
    -vepRefSeq, --vepRefSeq
      enable VEP 'SYMBOL'= XM_ or NM_
      Default: false
    -vepSymbol, --vepSymbol
      enable VEP 'SYMBOL'
      Default: false
    --version
      print version and exit

```

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
$ make vcfburdensplitter
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burden/VcfBurdenSplitter.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burden/VcfBurdenSplitter.java
)
## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfburdensplitter** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


### Description

This tools reads a VCF  (it should be sorted on chrom/POS and annotated with Ensembl Variation Predictor) and split data into genomic area of interest (gene, transcripts...).
For each area, a small VCF is produced and a Fished test is computed.
The final output is a set of concatenated VCF files. You could insert in a database using VcfDerby01



