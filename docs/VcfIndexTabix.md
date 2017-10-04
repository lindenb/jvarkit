# VcfIndexTabix

Index and sort a VCF on the fly with Tabix


## Usage

```
Usage: vcfindextabix [options] Files
  Options:
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
    -s, --sort
      sort VCF prior to saving
      Default: false
    --tmpDir
      tmp working directory. Default: java.io.tmpDir
      Default: []
    --version
      print version and exit

```


## Keywords

 * vcf
 * tabix


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
$ make vcfindextabix
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfIndexTabix.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfIndexTabix.java)


<details>
<summary>Git History</summary>

```
Mon May 29 16:53:42 2017 +0200 ; moved to docs ; https://github.com/lindenb/jvarkit/commit/6c0535d7add884e75b424af89a4f00aff6fae75f
Thu May 11 16:20:27 2017 +0200 ; move to jcommander ; https://github.com/lindenb/jvarkit/commit/15b6fabdbdd7ce0d1e20ca51e1c1a9db8574a59e
Tue May 9 10:40:20 2017 +0200 ; moving to jcommander ; https://github.com/lindenb/jvarkit/commit/88cfdecb60c1f193ae8b3176ad86181c4a15256b
Thu May 4 14:27:36 2017 +0200 ; fix make ; https://github.com/lindenb/jvarkit/commit/044933b2091b68f3241c60c85de5dd605ae02756
Mon Jan 25 18:43:22 2016 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/f49222d5a112d79fc4148f2a2a56e46a7ee5f517
Tue Sep 15 10:51:01 2015 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/ac8707d91503a75fabbc6bd2a2035092177249f8
Wed May 20 12:10:16 2015 +0200 ; vcf index tabix ; https://github.com/lindenb/jvarkit/commit/b73558a884ab1429f2c20721ff963d650faa12f6
Tue May 19 17:40:33 2015 +0200 ; vcf index tabix ; https://github.com/lindenb/jvarkit/commit/9c0d532b734ce8f5d9e2a2b67c4b2e3a764ca1d6
```

</details>

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfindextabix** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```bash
$ java -jar dist/vcfshuffle.jar in.vcf |\
  java -jar dist/vcfindextabix.jar -s -o out.vcf.gz

```

