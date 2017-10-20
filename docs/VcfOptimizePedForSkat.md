# VcfOptimizePedForSkat

Optimize ped file for SKAT


## Usage

```
Usage: vcfoptimizeped4skat [options] Files
  Options:
    --bootstrap
      bootstrap samples. Multiple list of sample separated with space, comma 
      or semicolons
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    --max-iter
      max number of iterations. -1 == infinite
      Default: -1
    --max-results
      max number of results.
      Default: 10
    -o, --output
      Output file. Optional . Default: stdout
    -ped, --pedigree
      A pedigree is a text file delimited with tabs. No header. Columns are 
      (1) Family (2) Individual-ID (3) Father Id or '0' (4) Mother Id or '0' 
      (5) Sex : 1 male/2 female / 0 unknown (6) Status : 0 unaffected, 1 
      affected,-9 unknown  If not defined, I will try to extract the pedigree 
      from  the VCFheader.
    -n, --remove
      max number of samples to remove
      Default: 1
    -seed, --seed
      random seed; -1=currentTimeMillis
      Default: 0
    --skat-accept-filtered
      accept variants FILTER-ed
      Default: false
    --skat-adjusted
      SKAT adjusted
      Default: false
    --skat-num-retry
      compute n-times the p-value
      Default: 1
    --skat-optimized
      SKAT optimized (SKATO)/ davies method.
      Default: false
    --skat-random-seed
      Rstats value for `set.seed`. -1 == use random
      Default: -1
    --version
      print version and exit

```


## Keywords

 * vcf
 * pedigree
 * skat
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
$ make vcfoptimizeped4skat
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/skat/VcfOptimizePedForSkat.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/skat/VcfOptimizePedForSkat.java)


<details>
<summary>Git History</summary>

```
Fri Oct 20 16:36:05 2017 +0200 ; skat continue ; https://github.com/lindenb/jvarkit/commit/54e62cdc08a38d1685b3842d300ec30740f2788a
Thu Oct 19 15:53:48 2017 +0200 ; skat continue ; https://github.com/lindenb/jvarkit/commit/5c71e1cbcacfd5b034a49580655db7066d83c50e
Wed Oct 18 19:20:26 2017 +0200 ; skat optimize vcf ; https://github.com/lindenb/jvarkit/commit/75da9b6ddd1f2daaf04365ecec4f712ee79851a6
```

</details>

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfoptimizeped4skat** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



