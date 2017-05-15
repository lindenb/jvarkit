# VariantsInWindow


## Usage

```
Usage: vcfwindowvariants [options] Files
  Options:
    -filter, --filter
      if --treshold is != -1 and the number of matches is greater than 
      threshold, set this FILTER
      Default: TOO_MANY_CLOSE_VARIANTS
    -h, --help
      print help and exits
    -o, --output
      Output file. Optional . Default: stdout
    -treshold, --treshold
      Number of variants to set the FILTER
      Default: -1
    --version
      print version and exits
    -best
      Only print the window with the hightest number of matches
      Default: false
    -noemptywin
      Don't print Windows in INFO having zero match.
      Default: false
    -select
      Optional Jexl expression to use when selecting the adjacent variants
      Default: []
    -windowName
      INFO Attribute name that will be added
      Default: WINDOW
    -shift, -windowShift
      Window shift (in bp.)
      Default: 50
    -windowSize
      Window Size (in bp.)
      Default: 150

```


## Description

Annotate Variants using a sliding window.


## Keywords

 * vcf
 * annotation


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
$ make vcfwindowvariants
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

https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VariantsInWindow.java

## Contribute

- Issue Tracker: http://github.com/lindenb/jvarkit/issues
- Source Code: http://github.com/lindenb/jvarkit

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfwindowvariants** ? https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md

The current reference is:

http://dx.doi.org/10.6084/m9.figshare.1425030

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> http://dx.doi.org/10.6084/m9.figshare.1425030


