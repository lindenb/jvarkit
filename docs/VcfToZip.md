# VcfToZip

Reads a stream of concatenated VCFs and insert them into a Zip file


## Usage

```
Usage: vcf2zip [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -o, --output
      Output zip file.
    -p, --prefix
      Prefix all zip entries with this prefix
      Default: VCF
    -t, --title
      Try to find ##(TITLE)=abcdefghijk in the VCF header and use it as the 
      name of the inserted VCF file
      Default: <empty string>
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
$ make vcf2zip
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfconcat/VcfToZip.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfconcat/VcfToZip.java
)
## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcf2zip** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)





### Motivation

This tool was used to create a zip from the output of VCFburdensplitter which is a stream of VCFs.



### Example



```

$ cat ~/input.vcf ~/input.vcf ~/input.vcf | java -jar dist/vcf2zip.jar -o jeter.zip
[main] INFO jvarkit - Command Line args : -o jeter.zip
[main] INFO jvarkit - Executing as lindenb@kaamelot-master01 on Linux 2.6.32-431.17.1.el6.x86_64 amd64; Java HotSpot(TM) 64-Bit Server VM 1.8.0_60-b27
[main] INFO jvarkit - reading concatenated vcf from stdin
[main] INFO jvarkit - VCF/vcf2zip.00001.vcf
[main] INFO jvarkit - Count: 499 Elapsed: 10 seconds(0.05%) Remains: 6 hours(99.95%) Last: 1:1431105
[main] INFO jvarkit - done: N=870
[main] INFO jvarkit - VCF/vcf2zip.00002.vcf
[main] INFO jvarkit - Count: 530 Elapsed: 10 seconds(0.05%) Remains: 5 hours(99.95%) Last: 1:1510577
[main] INFO jvarkit - done: N=870
[main] INFO jvarkit - VCF/vcf2zip.00003.vcf
[main] INFO jvarkit - Count: 530 Elapsed: 10 seconds(0.05%) Remains: 5 hours(99.95%) Last: 1:1510577
[main] INFO jvarkit - done: N=870
[main] INFO jvarkit - done. Number of VCFs:3
[main] INFO jvarkit - End JOB  [Mon May 02 12:30:24 CEST 2016] VcfToZip done. Elapsed time: 0.85 minutes.
$ unzip -t jeter.zip 
Archive:  jeter.zip
    testing: VCF/vcf2zip.00001.vcf    OK
    testing: VCF/vcf2zip.00002.vcf    OK
    testing: VCF/vcf2zip.00003.vcf    OK
No errors detected in compressed data of jeter.zip.

```






