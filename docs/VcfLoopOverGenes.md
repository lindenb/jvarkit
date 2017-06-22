# VcfLoopOverGenes

generate a BED file of the Genes in a VCF, loop over those genes


## Usage

```
Usage: vcfloopovergenes [options] Files
  Options:
    -g, --gene, -gene, --genes
      Loop over the gene file. If not defined VCF will be scanned for SnpEff 
      annotations and output will be a BED file with the gene names and 
      provenance.If defined, I will create a VCF for each Gene.
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    --maxRecordsInRam
      When writing  files that need to be sorted, this will specify the number 
      of records stored in RAM before spilling to disk. Increasing this number 
      reduces the number of file  handles needed to sort a file, and increases 
      the amount of RAM needed
      Default: 50000
    -o, --output
      Output file. Optional . Default: stdout
    -p, -prefix, --prefix
      File prefix when saving the individual VCF files.
      Default: <empty string>
    --tmpDir
      tmp working directory. Default: java.io.tmpDir
      Default: []
    --version
      print version and exit

```


## Keywords

 * vcf
 * gene
 * burden


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
$ make vcfloopovergenes
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burden/VcfLoopOverGenes.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burden/VcfLoopOverGenes.java
)
## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfloopovergenes** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

Generate the bed file from a VCF annotated with SnpEff

```
$ java -jar dist/vcfloopovergenes.jar -p KARAKA input.vcf.gz > genes.bed 
$ head jeter.bed
13	124462807	124462808	KARAKA.000000002	2V30:A_440-A_477:ENST00000232607	ANN_FeatureID	1
13	124411689	124420735	KARAKA.000000004	AC080008.1	ANN_GeneName	30
13	124475803	124490595	KARAKA.000000006	ENSG00000082781	ANN_GeneID	284
13	124444306	124468961	KARAKA.000000008	ENSG00000114491	ANN_GeneID	1545
(...)
```

Generate the VCFs:


```
 $ java -jar dist/vcfloopovergenes.jar -p KARAKA -g genes.bed -o tmp input.vcf.gz
 

 $ head tmp/KARAKA.manifest.txt 
KARAKA.3_KARAKA.000000001.vcf
KARAKA.3_KARAKA.000000002.vcf
KARAKA.3_KARAKA.000000003.vcf
KARAKA.3_KARAKA.000000004.vcf
(..)
 ```

 ```
$ ls tmp/*.vcf | head
tmp/KARAKA.3_KARAKA.000000001.vcf
tmp/KARAKA.3_KARAKA.000000002.vcf
tmp/KARAKA.3_KARAKA.000000003.vcf
tmp/KARAKA.3_KARAKA.000000004.vcf
(...)
```





