# Biostar105754

bigwig : peak distance from specific genomic region


## Usage

```
Usage: biostar105754 [options] Files
  Options:
    -B, --bigwig
      Big Wig file
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * wig
 * bigwig



## See also in Biostars

 * [https://www.biostars.org/p/105754](https://www.biostars.org/p/105754)


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
$ make biostar105754
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar105754.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar105754.java)


<details>
<summary>Git History</summary>

```
Mon Aug 7 09:53:19 2017 +0200 ; fixed unicode problems after https://github.com/lindenb/jvarkit/issues/82 ; https://github.com/lindenb/jvarkit/commit/68254c69b027a9ce81d8b211447f1c0bf02dc626
Wed May 24 17:27:28 2017 +0200 ; lowres bam2raster & fix doc ; https://github.com/lindenb/jvarkit/commit/6edcfd661827927b541e7267195c762e916482a0
Sun May 21 17:11:09 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/aa4f02194fe00a1a842949e448661e227f16fe9f
Fri Apr 21 18:16:07 2017 +0200 ; scan sv ; https://github.com/lindenb/jvarkit/commit/49b99018811ea6a624e3df556627ebdbf3f16eab
Thu Apr 20 17:17:22 2017 +0200 ; continue transition jcommander ; https://github.com/lindenb/jvarkit/commit/fcf5def101925bea9ddd001d8260cf65aa52d6a0
Thu Apr 6 18:34:56 2017 +0200 ; moving to jcommander ; https://github.com/lindenb/jvarkit/commit/883b4ba4b693661663694256f16b137e371147fa
Wed Mar 30 17:45:00 2016 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/ef77dce3c82c470017916555304df6a470fbdad4
Thu Nov 26 17:41:15 2015 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/899c60335614350d463be66ec21e994b34dc55be
Fri Oct 2 16:58:59 2015 +0200 ; bioalci ; https://github.com/lindenb/jvarkit/commit/5c80e5e544c70bcfcd50fbf8cbb4020513894e98
Tue Jul 8 12:50:47 2014 +0200 ; checkErr in out ; https://github.com/lindenb/jvarkit/commit/235a18b083ea15c4ad94060de512f1edc74cec42
Tue Jul 8 12:27:19 2014 +0200 ; biostar105754 ; https://github.com/lindenb/jvarkit/commit/9e7774bca435cbf7bf0b5e1a9b203608044b40eb
```

</details>

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **biostar105754** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
$  echo -e "1\t1000\t20000\n3\t100\t200\nUn\t10\t11"  |\
  java -jar dist/biostar105754.jar -B path/to//All_hg19_RS_noprefix.b


#no data found for	Un	10	11
1	1000	1001	0.0	1	1000	20000
3	100	101	0.0	3	100	200
```


