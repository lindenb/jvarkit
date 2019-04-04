# NgsFilesScanner

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Build a persistent database of NGS file. Dump as XML.


## Usage

```
Usage: ngsfilesscanner [options] Files
  Options:
  * -B, --bdb-home
      berkeleydb home directory
    -D, --dump
      dump as XML to stdout and exit
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * ngs
 * bam
 * sam
 * vcf
 * xml


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew ngsfilesscanner
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/ngsfiles/NgsFilesScanner.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/ngsfiles/NgsFilesScanner.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **ngsfilesscanner** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

### Example 01 : scanning all files:

we scan all files under /common/data/projects.
all missing files/dir previously inserted will be deleted from the berkeley database

```bash
$ java -jar dist/ngsfilesscanner.jar -B /path/to/bdbdir /commun/data/projects/
```
### Example 02 : dumping the database

we scan all files under /common/data/projects.

```bash
$ java -jar dist/ngsfilesscanner.jar -B /path/to/bdbdir -D 
```
```xml
<?xml version="1.0" encoding="UTF-8"?>
<ngs-files>
  <fastq-dir directory="/commun/data/projects/Project_P1/Samplex1">
    <samples>
      <sample size="2834419615">x1</sample>
    </samples>
  </fastq-dir>
  (...)
 <vcf timestamp="1398643093000" file="/commun/data/projects/path/Samples/S2/S2.varscan.annotations.vcf.gz" filename="S2.varscan.annotations.vcf.gz" modified="Mon Apr 28 01:58:13 CEST 2014" size="21053412">
    <samples>
      <sample>S2</sample>
    </samples>
  </vcf>
</ngs-files>
```
## setting-up a CRON job

content of `/etc/cron.daily/ngsfilesscanner.cron`
```bash
#!/bin/bash
/usr/bin/java -cp /path/to/picard-tools-1.100/picard-1.100.jar:/path/to/picard-tools-1.100/sam-1.100.jar:/path/to/picard-tools-1.100/tribble-1.100.jar:/path/to/picard-tools-1.100/variant-1.100.jar:/path/to/picard-tools-1.100/commons-jexl-2.1.1.jar:/path/to/picard-tools-1.100/commons-logging-1.1.1.jar:/path/to/je-5.0.34/lib/je-5.0.34.jar:/path/to/jvarkit-git/ngsfilesscanner.jar com.github.lindenb.jvarkit.tools.ngsfiles.NgsFilesScanner -B /var/www/cgi-bin/ngsfiles /commun/data/
```
set permissions
```
$ sudo chown root:root /etc/cron.daily/ngsfilesscanner.cron
$ sudo chmod u+x ngsfilesscanner.cron
```
## creating a CGI dumping the results:

content of `/var/www/cgi-bin/ngsfiles.cgi`
```bash
#!/bin/bash

echo "Content-Type: text/xml"
echo ""
/usr/bin/java -cp /path/to/picard-tools-1.100/picard-1.100.jar:/path/to/picard-tools-1.100/sam-1.100.jar:/path/to/picard-tools-1.100/tribble-1.100.jar:/path/to/picard-tools-1.100/variant-1.100.jar:/path/to/picard-tools-1.100/commons-jexl-2.1.1.jar:/path/to/picard-tools-1.100/commons-logging-1.1.1.jar:/path/to/je-5.0.34/lib/je-5.0.34.jar:/path/to/jvarkit-git/ngsfilesscanner.jar com.github.lindenb.jvarkit.tools.ngsfiles.NgsFilesScanner -B /var/www/cgi-bin/ngsfiles  -D
```

set permissions
```
$ sudo chown root:root /var/www/cgi-bin/ngsfiles.cgi
$ sudo chmod u+x /var/www/cgi-bin/ngsfiles.cgi
```

check URL: `curl -s http://localhost/cgi-bin/ngsfiles.cgi | xmllint --format -`
```xml
<?xml version="1.0" encoding="UTF-8"?>
<ngs-files>
  <fastq-dir directory="/commun/data/projects/Project_P1/Samplex1">
    <samples>
      <sample size="2834419615">x1</sample>
    </samples>
  </fastq-dir>
  (...)
 <vcf timestamp="1398643093000" file="/commun/data/projects/path/Samples/S2/S2.varscan.annotations.vcf.gz" filename="S2.varscan.annotations.vcf.gz" modified="Mon Apr 28 01:58:13 CEST 2014" size="21053412">
    <samples>
      <sample>S2</sample>
    </samples>
  </vcf>
</ngs-files>
```

