# VcfRenameSamples

Rename the Samples in a VCF


## Usage

```
Usage: vcfrenamesamples [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -o, --output
      Output file. Optional . Default: stdout
    --outputbcf
      Output bcf (for streams)
      Default: false
    --vcfcreateindex
      VCF, create tribble or tabix Index when writing a VCF/BCF to a file.
      Default: false
    --vcfmd5
      VCF, create MD5 checksum when writing a VCF/BCF to a file.
      Default: false
    --version
      print version and exit
    -E
      error like src sample missing in VCF
      Default: false
  * -f
      Tab delimited file containing old-name\tnew-name

```


## Keywords

 * vcf
 * sample


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
$ make vcfrenamesamples
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfRenameSamples.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfRenameSamples.java
)
## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfrenamesamples** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


```bash
$ curl -s "https://raw.github.com/arq5x/gemini/master/test/test1.snpeff.vcf" |\
grep -A 2 CHROM  | cut -f 1-5,10-

#CHROM	POS	ID	REF	ALT	1094PC0005	1094PC0009	1094PC0012	1094PC0013
chr1	30860	.	G	C	0/0:7,0:7:15.04:0,15,177	0/0:2,0:2:3.01:0,3,39	0/0:6,0:6:12.02:0,12,143	0/0:4,0:4:9.03:0,9,119
chr1	69270	.	A	G	./.	./.	1/1:0,3:3:9.03:106,9,0	1/1:0,6:6:18.05:203,18,0

curl -s "https://raw.github.com/arq5x/gemini/master/test/test1.snpeff.vcf" |\
java -jar dist/vcfrenamesamples.jar -f <(echo -e "1094PC0005\tALPHA\n1094PC0012\tBETA\nSAMPLE\tGAMMA\n1094PC0009\tEPSILON") -E  |\
grep -A 2 "#CHROM" | cut -f 1-5,10-

#CHROM	POS	ID	REF	ALT	ALPHA	EPSILON	BETA	1094PC0013
chr1	30860	.	G	C	0/0:7,0:7:15:0,15,177	0/0:2,0:2:3:0,3,39	0/0:6,0:6:12:0,12,143	0/0:4,0:4:9:0,9,119
chr1	69270	.	A	G	./.	./.	1/1:0,3:3:9:106,9,0	1/1:0,6:6:18:203,18,0
```


