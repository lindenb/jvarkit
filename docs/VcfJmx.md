# VcfJmx

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Monitor/interrupt/break a VCF stream with java JMX http://www.oracle.com/technetwork/articles/java/javamanagement-140525.html


## Usage

```
Usage: vcfjmx [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --out
      Output vcf , ot stdin
    --version
      print version and exit
    -p
      Stream identifier

```


## Keywords

 * java
 * jmx
 * vcf


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfjmx
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/jmx/VcfJmx.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/jmx/VcfJmx.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfjmx** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

## Example

```bash
$   java -jar dist/vcfjmx.jar -p MyWorkflow1 input.vcf > /dev/null
```

while the stream is running, open a new jconsole https://docs.oracle.com/javase/7/docs/technotes/guides/management/jconsole.html . here you can get the number of records, the elapsed time. Two operation are available:

* doBreak: interrupt current streaming , exit with success (0)
* doAbort: interrupt current streaming , exit with failure (-1)

```
$ java -jar dist/vcfjmx.jar -p 1000G  ALL.wgs.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf.gz | \
    java -jar dist/vcfjmx.jar  -p 1000G-2  |\
    java -jar dist/vcfjmx.jar  -p 1000G-3  |\
    java -jar dist/vcfjmx.jar  -p 1000G-4 > /dev/null

[INFO/VcfJmx] 2015-07-10 14:11:08 "Starting JOB at Fri Jul 10 14:11:08 CEST 2015 com.github.lindenb.jvarkit.tools.jmx.VcfJmx version=4f797a9fbf2c3ceac9cec3c431c719ad794953c2  built=2015-07-10:13-07-05"
[INFO/VcfJmx] [INFO/VcfJmx] 2015-07-10 14:11:08 "Command Line args : -p 1000G-2"
2015-07-10 14:11:08 "Starting JOB at Fri Jul 10 14:11:08 CEST 2015 com.github.lindenb.jvarkit.tools.jmx.VcfJmx version=4f797a9fbf2c3ceac9cec3c431c719ad794953c2  built=2015-07-10:13-07-05"
[INFO/VcfJmx] 2015-07-10 14:11:08 "Command Line args : -p 1000G-3"
[INFO/VcfJmx] 2015-07-10 14:11:08 "Starting JOB at Fri Jul 10 14:11:08 CEST 2015 com.github.lindenb.jvarkit.tools.jmx.VcfJmx version=4f797a9fbf2c3ceac9cec3c431c719ad794953c2  built=2015-07-10:13-07-05"
[INFO/VcfJmx] 2015-07-10 14:11:08 "Command Line args : -p 1000G ALL.wgs.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf.gz"
[INFO/VcfJmx] 2015-07-10 14:11:08 "Executing as lindenb@kaamelot-master01 on Linux 2.6.32-431.17.1.el6.x86_64 amd64; Java HotSpot(TM) 64-Bit Server VM 1.7.0_60-b19"
[INFO/VcfJmx] 2015-07-10 14:11:08 "Executing as lindenb@kaamelot-master01 on Linux 2.6.32-431.17.1.el6.x86_64 amd64; Java HotSpot(TM) 64-Bit Server VM 1.7.0_60-b19"
[INFO/VcfJmx] 2015-07-10 14:11:08 "Executing as lindenb@kaamelot-master01 on Linux 2.6.32-431.17.1.el6.x86_64 amd64; Java HotSpot(TM) 64-Bit Server VM 1.7.0_60-b19"
[INFO/VcfJmx] 2015-07-10 14:11:08 "reading from stdin"
[INFO/VcfJmx] 2015-07-10 14:11:08 "reading from stdin"
[INFO/VcfJmx] 2015-07-10 14:11:08 "reading from ALL.wgs.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf.gz"

[SEVERE/VcfJmx] 2015-07-10 14:11:46 "#### Process "1000G-3" received message BREAK"

[INFO/VcfJmx] 2015-07-10 14:11:46 "Number of Variants:1130774"
[INFO/VcfJmx] 2015-07-10 14:11:46 "End JOB status=0 [Fri Jul 10 14:11:46 CEST 2015] com.github.lindenb.jvarkit.tools.jmx.VcfJmx done. Elapsed time: 0.64 minutes."
[INFO/VcfJmx] 2015-07-10 14:11:46 "Number of Variants:1120774"
[INFO/VcfJmx] 2015-07-10 14:11:46 "Number of Variants:1110975"
[INFO/VcfJmx] 2015-07-10 14:11:46 "End JOB status=0 [Fri Jul 10 14:11:46 CEST 2015] com.github.lindenb.jvarkit.tools.jmx.VcfJmx done. Elapsed time: 0.64 minutes."
[INFO/VcfJmx] 2015-07-10 14:11:46 "End JOB status=0 [Fri Jul 10 14:11:46 CEST 2015] com.github.lindenb.jvarkit.tools.jmx.VcfJmx done. Elapsed time: 0.64 minutes."
```
