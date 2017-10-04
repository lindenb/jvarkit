# VcfConcat

Concatenante sorted VCF with same sample, does NOT merge genotypes


## Usage

```
Usage: vcfconcat [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -o, --out
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * vcf


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
$ make vcfconcat
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfconcat/VcfConcat.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfconcat/VcfConcat.java)


<details>
<summary>Git History</summary>

```
Tue Jun 6 18:06:17 2017 +0200 ; postponed vcf ; https://github.com/lindenb/jvarkit/commit/bcd52318caf3cd76ce8662485ffaacaabde97caf
Sun Jun 4 21:53:22 2017 +0200 ; writing bcf ; https://github.com/lindenb/jvarkit/commit/784fdac37cd7e6eca04e35d0a3ddad8637826b4a
Mon May 29 16:53:42 2017 +0200 ; moved to docs ; https://github.com/lindenb/jvarkit/commit/6c0535d7add884e75b424af89a4f00aff6fae75f
Mon May 22 17:20:59 2017 +0200 ; moving to jcommaner ; https://github.com/lindenb/jvarkit/commit/60cbfa764f7f5bacfdb78e48caf8f9b66e53a6a0
Sun Apr 16 12:35:50 2017 +0200 ; knime logger ; https://github.com/lindenb/jvarkit/commit/d16e6079bfa7cf6a338365eed90542f9fa551995
Tue Mar 17 16:59:10 2015 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/b74a74277f72f240dab3360a49fdb8357f7bfbbd
Mon Mar 16 21:28:39 2015 +0100 ; fix this tomorrow ; https://github.com/lindenb/jvarkit/commit/6597489a5c1f6c04280f9e98f93264e732564b44
Thu Jan 15 17:30:17 2015 +0100 ; multi to one ; https://github.com/lindenb/jvarkit/commit/4a5399a30c7780fadee3a54109ba22940140d01f
Mon Jun 23 12:34:44 2014 +0200 ; find-a-variation + using abstractcodec instead of vcfcodec ; https://github.com/lindenb/jvarkit/commit/da621ba8326d56da8f6907c845c539e4ea785284
Mon May 12 10:28:28 2014 +0200 ; first sed on files ; https://github.com/lindenb/jvarkit/commit/79ae202e237f53b7edb94f4326fee79b2f71b8e8
Sun Feb 2 18:55:03 2014 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/abd24b56ec986dada1e5162be5bbd0dac0c2d57c
Mon Dec 30 19:33:23 2013 +0100 ; vcf concat ; https://github.com/lindenb/jvarkit/commit/b5ebf67dd2926d8a6afadb4d1e36a4959508057f
```

</details>

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfconcat** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

### From stdin

```bash
$ find ./ -name "*.vcf" | grep Sample1 | java -jar dist/vcfconcat.jar > out.vcf
```

### From files

```bash
$ java -jar dist/vcfconcat.jar Sample1.samtools.vcf Sample1.gatk.vcf > out.vcf
```



