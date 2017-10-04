# VCFMerge2

Merge VCF Files


## DEPRECATED

use GATK combineVariants.

## Usage

```
Usage: vcfmerge [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -homref, --homref
      Use HomRef 0/0 for unknown variant
      Default: false
    --maxRecordsInRam
      When writing  files that need to be sorted, this will specify the number 
      of records stored in RAM before spilling to disk. Increasing this number 
      reduces the number of file  handles needed to sort a file, and increases 
      the amount of RAM needed
      Default: 50000
    -m, --nomerge
      Do NOT merge VariantContext lines, but create multiple lines
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    -region, --region
      Merge in that region: An interval as the following syntax : 
      "chrom:start-end" or "chrom:middle+extend"  or 
      "chrom:start-end+extend".A program might use a Reference sequence to fix 
      the chromosome name (e.g: 1->chr1)
      Default: <empty string>
    -s, --sorted
      files are known to be ROD sorted
      Default: false
    --tmpDir
      tmp working directory. Default: java.io.tmpDir
      Default: []
    --version
      print version and exit

```


## Keywords

 * vcf
 * sort


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
$ make vcfmerge
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfmerge/VCFMerge2.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfmerge/VCFMerge2.java)


<details>
<summary>Git History</summary>

```
Mon Sep 18 17:29:05 2017 +0200 ; adding test for ConvertVcfChromosomes ; https://github.com/lindenb/jvarkit/commit/39c750097a13c007c850449a4586de5b51962242
Thu Jun 22 13:16:05 2017 +0200 ; vcfloopovergenes ; https://github.com/lindenb/jvarkit/commit/aa4a6f29c853efddcee5678f9441d9994a2deee6
Thu Jun 22 09:36:31 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/918c9929cfbfa5fd2b915776bc95422297615c4b
Wed Jun 21 17:31:49 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/556a9b4ed5d047a215e160c0a480ea241cea83d9
Wed Jun 21 15:27:13 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/034f57d0e8d0399c12b290385d89e498e6138e1d
Tue Jun 20 15:07:17 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/8e712ff2b8d4b73d71cd2035cfac57381d3e9d4b
Tue Jun 6 18:06:17 2017 +0200 ; postponed vcf ; https://github.com/lindenb/jvarkit/commit/bcd52318caf3cd76ce8662485ffaacaabde97caf
Sun Jun 4 21:53:22 2017 +0200 ; writing bcf ; https://github.com/lindenb/jvarkit/commit/784fdac37cd7e6eca04e35d0a3ddad8637826b4a
Mon May 29 16:53:42 2017 +0200 ; moved to docs ; https://github.com/lindenb/jvarkit/commit/6c0535d7add884e75b424af89a4f00aff6fae75f
Wed May 17 14:09:36 2017 +0200 ; fix typo bioalcidae ; https://github.com/lindenb/jvarkit/commit/9db2344e7ce840df02c5a7b4e2a91d6f1a5f2e8d
Fri Apr 21 18:16:07 2017 +0200 ; scan sv ; https://github.com/lindenb/jvarkit/commit/49b99018811ea6a624e3df556627ebdbf3f16eab
Thu Jun 2 09:49:17 2016 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/2ae46b7df29c6f1b66ce5104ea03bf6390db120d
Thu Feb 11 16:18:06 2016 +0100 ; vcfmerge ; https://github.com/lindenb/jvarkit/commit/54cb61fc38df634279aadcabc3dc7bb0c0765314
Fri Jan 22 23:49:23 2016 +0100 ; vcfiterator is now an interface ; https://github.com/lindenb/jvarkit/commit/9f9b9314c4b31b21044c5911a7e79e1b3fb0af7a
Mon Jun 8 17:24:41 2015 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/f9a941d604f378ff40a32666c8381cb2450c7cfa
Fri Jun 5 12:42:21 2015 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/cc909f9f4ceea181bb65e4203e3fdbde176c6f2f
Mon Jun 1 15:27:11 2015 +0200 ; change getChrom() to getContig() ; https://github.com/lindenb/jvarkit/commit/5abd60afcdc2d5160164ae6e18087abf66d8fcfe
Fri Feb 20 17:50:08 2015 +0100 ; continue integration in knime ; https://github.com/lindenb/jvarkit/commit/5693e07342a30f21c807a4e3a655e3446019458f
Fri Nov 7 17:11:35 2014 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/17b77998c42387988386f5a696bba464d130cf86
Wed Nov 5 10:45:45 2014 +0100 ; vcf merge fixed ; https://github.com/lindenb/jvarkit/commit/7ddcffc73f823f9e377ffd2a3644cbf50cf26581
Mon Oct 13 18:29:16 2014 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/c83f20cde867920870918ee6eb5e5406f554e2bb
Thu Sep 11 09:36:01 2014 +0200 ; problem with java dataInputSTream: writeUTF requires line.length < SHORt_MAX ; https://github.com/lindenb/jvarkit/commit/19eac4ee36909a730903546b50461de3c19a5c1f
Mon Jun 23 12:34:44 2014 +0200 ; find-a-variation + using abstractcodec instead of vcfcodec ; https://github.com/lindenb/jvarkit/commit/da621ba8326d56da8f6907c845c539e4ea785284
Mon May 12 15:27:08 2014 +0200 ; moving to htsjdk ; https://github.com/lindenb/jvarkit/commit/fd30a81154a16835b5bab3d8e1ef90c9fee6bdcb
Mon May 12 14:06:30 2014 +0200 ; continue moving to htsjdk ; https://github.com/lindenb/jvarkit/commit/011f098b6402da9e204026ee33f3f89d5e0e0355
Mon May 12 10:28:28 2014 +0200 ; first sed on files ; https://github.com/lindenb/jvarkit/commit/79ae202e237f53b7edb94f4326fee79b2f71b8e8
Sun Feb 2 18:55:03 2014 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/abd24b56ec986dada1e5162be5bbd0dac0c2d57c
Sat Jan 4 16:27:34 2014 +0100 ; vcf merge, peek vcf, samgrep cmd line, sw factory ; https://github.com/lindenb/jvarkit/commit/bd9a33f08c6e32efe9b54e1909527dc57dc55060
```

</details>

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfmerge** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


 
## Example




```bash
$  find ./ -name "*.vcf.gz" | xargs java -jar dist/vcfmerge.jar   > out.vcf
```


