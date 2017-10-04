# Biostar78285

Extract regions of genome that have 0 coverage See http://www.biostars.org/p/78285/


## Usage

```
Usage: biostar78285 [options] Files
  Options:
    -f, --filter
      A filter expression. Reads matching the expression will be filtered-out. 
      Empty String means 'filter out nothing/Accept all'. See https://github.com/lindenb/jvarkit/blob/master/src/main/resources/javacc/com/github/lindenb/jvarkit/util/bio/samfilter/SamFilterParser.jj 
      for a complete syntax.
      Default: mapqlt(1) || MapQUnavailable() || Duplicate() || FailsVendorQuality() || NotPrimaryAlignment() || SupplementaryAlignment()
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

 * sam
 * bam
 * depth
 * coverage



## See also in Biostars

 * [https://www.biostars.org/p/78285](https://www.biostars.org/p/78285)


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
$ make biostar78285
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar78285.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar78285.java)


<details>
<summary>Git History</summary>

```
Wed May 24 17:27:28 2017 +0200 ; lowres bam2raster & fix doc ; https://github.com/lindenb/jvarkit/commit/6edcfd661827927b541e7267195c762e916482a0
Thu May 11 16:20:27 2017 +0200 ; move to jcommander ; https://github.com/lindenb/jvarkit/commit/15b6fabdbdd7ce0d1e20ca51e1c1a9db8574a59e
Wed Apr 19 10:40:28 2017 +0200 ; rm-xml ; https://github.com/lindenb/jvarkit/commit/971b090382a1b0b96e250030a5c8e7be500593b7
Mon Dec 28 20:23:04 2015 +0100 ; sam2axt ; https://github.com/lindenb/jvarkit/commit/a2edef74730256e93d244e440a79e7362d647795
Mon Mar 9 14:47:06 2015 +0100 ; moving vcf2sql to mysql ; https://github.com/lindenb/jvarkit/commit/f2813fc2fbf434da37526f038b60181564881c8e
Mon Mar 9 10:52:57 2015 +0100 ; rewrote biostar78285  (regions with 0 coverage) with htsjdk #tweet ; https://github.com/lindenb/jvarkit/commit/3b1521878efdbf6b5966b461438e4344633966a3
Fri May 23 15:00:53 2014 +0200 ; cont moving to htsjdk ; https://github.com/lindenb/jvarkit/commit/81f98e337322928b07dfcb7a4045ba2464b7afa7
Mon May 12 14:06:30 2014 +0200 ; continue moving to htsjdk ; https://github.com/lindenb/jvarkit/commit/011f098b6402da9e204026ee33f3f89d5e0e0355
Mon May 12 10:28:28 2014 +0200 ; first sed on files ; https://github.com/lindenb/jvarkit/commit/79ae202e237f53b7edb94f4326fee79b2f71b8e8
Tue Nov 26 12:29:03 2013 +0100 ; unclipped start -> align start ; https://github.com/lindenb/jvarkit/commit/3944b21281c2b4afc1ef682f0abe020b26940e37
Tue Aug 6 18:54:07 2013 +0200 ; biostar + cigar ; https://github.com/lindenb/jvarkit/commit/218d1fa11e545c30b1b0a93198a7f5ec701c3c88
Tue Aug 6 17:31:26 2013 +0200 ; samlocusiterator for Biostar78285 ; https://github.com/lindenb/jvarkit/commit/25fad045dc0c4a118aa3b59049fa6b1c2c46880c
Tue Aug 6 15:04:00 2013 +0200 ; ops ; https://github.com/lindenb/jvarkit/commit/96ada2e69fb2a5c3b51cabcf849768610d614d91
Tue Aug 6 14:51:11 2013 +0200 ; biostar78285 ; https://github.com/lindenb/jvarkit/commit/43f1fe3d2f6ee4c1ec159034ca552f2839074611
```

</details>

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **biostar78285** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```bash
 $ java -jar dist/biostar78285.jar  sorted.bam 
 	

seq1	1569	1575
seq2	1567	1584
```


