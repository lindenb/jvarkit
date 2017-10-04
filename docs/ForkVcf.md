# ForkVcf

Fork a VCF.


## Usage

```
Usage: forkvcf [options] Files
  Options:
    -n, --count
      number of vcf files to generate
      Default: 2
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -m, --manifest
      optional save produced vcf filenames in this file.
    -maxRecordsInRam, --maxRecordsInRam
      Max records in RAM
      Default: 50000
  * -o, --output
      Output file Must contains __GROUPID__
    -c, --splitbychunk
      When this option is used, the variant are first saved in a temporary 
      file, the number of variant is dividided by 'count' and the output files 
      are lineray produced. The default is to dispatch the variants as they 
      are coming in the stream.
      Default: false
    -T, --tmpDir
      mp directory
      Default: /tmp
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
$ make forkvcf
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/ForkVcf.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/ForkVcf.java)


<details>
<summary>Git History</summary>

```
Sun May 7 13:21:47 2017 +0200 ; rm xml ; https://github.com/lindenb/jvarkit/commit/f37088a9651fa301c024ff5566534162bed8753d
Tue Apr 25 15:40:45 2017 +0200 ; cont jcommander ; https://github.com/lindenb/jvarkit/commit/16aeb209fda502b60dd75689b85d1304f469775b
Wed Jun 8 12:51:03 2016 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/3a139dad3aa0c899b4a84c9a0d2908d47ecccd58
```

</details>

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **forkvcf** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)





### Output

Output filename (option -o) MUST contain the word __GROUPID__.



### Example



```
$ 

```






```


```

cat input.vcf | java -jar dist/forkvcf.jar -n 3 -o "_tmp.__GROUPID__.vcf"
[main] INFO jvarkit - opening VCF file "_tmp.00001.vcf" for writing
[main] INFO jvarkit - opening VCF file "_tmp.00002.vcf" for writing
[main] INFO jvarkit - opening VCF file "_tmp.00003.vcf" for writing

$ wc _tmp.0000*
   226   6819 143947 _tmp.00001.vcf
   226   6819 140792 _tmp.00002.vcf
   225   6161 125219 _tmp.00003.vcf
   
   
   


### See also


 *  https://github.com/lindenb/jvarkit/wiki/SplitVcf






