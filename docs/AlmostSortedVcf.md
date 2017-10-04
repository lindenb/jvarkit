# AlmostSortedVcf

Sort an 'almost' sorted VCF. Most variants should be sorted but a few  consecutive lines might have been switched by a caller.


## Usage

```
Usage: almostsortedvcf [options] Files
  Options:
    -N, --buffersize
      max records in ram
      Default: 1000
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
$ make almostsortedvcf
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/sortvcfonref/AlmostSortedVcf.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/sortvcfonref/AlmostSortedVcf.java)


<details>
<summary>Git History</summary>

```
Mon May 29 16:53:42 2017 +0200 ; moved to docs ; https://github.com/lindenb/jvarkit/commit/6c0535d7add884e75b424af89a4f00aff6fae75f
Mon May 22 17:20:59 2017 +0200 ; moving to jcommaner ; https://github.com/lindenb/jvarkit/commit/60cbfa764f7f5bacfdb78e48caf8f9b66e53a6a0
Wed Apr 26 17:26:23 2017 +0200 ; cont jcommander ; https://github.com/lindenb/jvarkit/commit/ab6c7b760cd5376e08da24426cede7f84a6b3ae2
Fri Jan 27 16:16:46 2017 +0100 ; gatkwalker with range display ; https://github.com/lindenb/jvarkit/commit/7e6cde60fa6f90883a47df594769fbb4e6f1161b
Mon Jun 30 13:37:16 2014 +0200 ; almost sorted ; https://github.com/lindenb/jvarkit/commit/2375fb817a80f9ed9b7688421ebdb722a214a052
```

</details>

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **almostsortedvcf** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

## Example

```bash
$  java -jar dist/almostsortedvcf.jar < in.vcf > out.vcf
```

