# SplitVcf

split a vcf...


## Usage

```
Usage: splitvcf [options] Files
  Options:
    -g, --groupfile
      Chromosome group file. Intervals are 1-based. If undefined, splitvcf 
      will use the sequence dictionary to output one vcf per contig.
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -m, --multi
      if set, allow one variant to be mapped on multiple chromosome group (the 
      record is duplicated)
      Default: false
  * -o, --out
      Output file (or stdout). Name must contain '__GROUPID__'
    -u, --unmapped
      unmapped interval name
      Default: OTHER
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
$ make splitvcf
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/SplitVcf.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/SplitVcf.java)


<details>
<summary>Git History</summary>

```
Fri May 19 17:10:13 2017 +0200 ; cont doc ; https://github.com/lindenb/jvarkit/commit/d2aea1eaa554d0498b197fb8fac01893b10ceb83
Tue May 9 20:36:16 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/517cc3660251857061fa955cce5c8e07362c5bee
Fri Mar 31 17:08:11 2017 +0200 ; moving to jcommander ; https://github.com/lindenb/jvarkit/commit/f78937d19c4b038e69a32fbcfa2aeab8fd8417c6
Wed Jun 8 12:51:03 2016 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/3a139dad3aa0c899b4a84c9a0d2908d47ecccd58
Fri Jun 3 19:44:07 2016 +0200 ; splitvcf ; https://github.com/lindenb/jvarkit/commit/802aa5e2bc96303bf0b64e3188f1a277764db453
```

</details>

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **splitvcf** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
$ cat groups.txt
G1	10:112583204-112583210
G2	11
```


```
$ java -jar dist/splitvcf.jar  -o tmp__GROUPID__.vcf.gz -g groups.txt in.vcf
$ ls tmp*
tmpG1.vcf.gz
tmpG2.vcf.gz
tmpOTHER.vcf.gz
```

## See also

* https://github.com/lindenb/jvarkit/wiki/SplitBam3


