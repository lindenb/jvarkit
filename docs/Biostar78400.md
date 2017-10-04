# Biostar78400

add the read group info to the sam file on a per lane basis


## Usage

```
Usage: biostar78400 [options] Files
  Options:
    --bamcompression
      Compression Level.
      Default: 5
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -o, --output
      Output file. Optional . Default: stdout
    -p, --regex
      Regular expression that can be used to parse read names in the incoming 
      SAM file. Flowcell: (group 1)and the lane (group 2). Another pattern 
      could be 
      '[a-zA-Z0-9\-]+:[0-9]+:([a-zA-Z0-9]+):([0-9]):[0-9]+:[0-9]+:[0-9]+.*.' 
      (Highseq) 
      Default: ([a-zA-Z0-9]+):([0-9]):[0-9]+:[0-9]+:[0-9]+.*
    --samoutputformat
      Sam output format.
      Default: TypeImpl{name='SAM', fileExtension='sam', indexExtension='null'}
    --version
      print version and exit
  * -x, --xmlFile
      XML description of the groups.

```


## Keywords

 * sam
 * bam
 * xml
 * read-group


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
$ make biostar78400
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar78400.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar78400.java)


<details>
<summary>Git History</summary>

```
Wed May 24 17:27:28 2017 +0200 ; lowres bam2raster & fix doc ; https://github.com/lindenb/jvarkit/commit/6edcfd661827927b541e7267195c762e916482a0
Sun May 21 20:02:10 2017 +0200 ; instanceMain -> instanceMainWithExit ; https://github.com/lindenb/jvarkit/commit/4fa41d198fe7e063c92bdedc333cbcdd2b8240aa
Thu May 11 16:20:27 2017 +0200 ; move to jcommander ; https://github.com/lindenb/jvarkit/commit/15b6fabdbdd7ce0d1e20ca51e1c1a9db8574a59e
Sat Apr 29 18:45:47 2017 +0200 ; partition ; https://github.com/lindenb/jvarkit/commit/7d72633d50ee333fcad0eca8aaa8eec1a475cc4d
Wed Apr 19 10:40:28 2017 +0200 ; rm-xml ; https://github.com/lindenb/jvarkit/commit/971b090382a1b0b96e250030a5c8e7be500593b7
Tue Jul 26 08:53:18 2016 +0200 ; fix bug https://github.com/lindenb/jvarkit/issues/59 ; https://github.com/lindenb/jvarkit/commit/8b78fc53d6b40eb0b32264b05fd9d465467eaa94
Tue Jul 19 08:41:16 2016 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/cb9a8f435beb0ff3db6ae0e45459d5a9eac5d4c3
Wed Jan 6 17:39:57 2016 +0100 ; VcfMultiToOneInfo ; https://github.com/lindenb/jvarkit/commit/4bea71a6d15bdb12288d2c34ffc43892004476c9
Mon Dec 14 17:18:02 2015 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/9b271459821d8061aa07e98bc7f30232597f47c9
Fri Jun 5 12:42:21 2015 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/cc909f9f4ceea181bb65e4203e3fdbde176c6f2f
Fri May 23 15:00:53 2014 +0200 ; cont moving to htsjdk ; https://github.com/lindenb/jvarkit/commit/81f98e337322928b07dfcb7a4045ba2464b7afa7
Mon May 12 14:06:30 2014 +0200 ; continue moving to htsjdk ; https://github.com/lindenb/jvarkit/commit/011f098b6402da9e204026ee33f3f89d5e0e0355
Mon May 12 10:28:28 2014 +0200 ; first sed on files ; https://github.com/lindenb/jvarkit/commit/79ae202e237f53b7edb94f4326fee79b2f71b8e8
Wed Aug 7 19:39:25 2013 +0200 ; biostar78400 ; https://github.com/lindenb/jvarkit/commit/6bac83632d0646999f4fca2dba75fa83c91add99
```

</details>

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **biostar78400** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)





### Read names

Reads' name should start with the following signature:




### XML

the XML should look like this:


```

<read-groups>
<flowcell name="HS2000-1259_127">
 <lane index="1">
   <group ID="X1">
     <library>L1</library>
     <platform>P1</platform>
     <sample>S1</sample>
     <platformunit>PU1</platformunit>
     <center>C1</center>
     <description>blabla</description>
   </group>
 </lane>
</flowcell>
<flowcell name="HS2000-1259_128">
 <lane index="2">
   <group ID="x2">
     <library>L2</library>
     <platform>P2</platform>
     <sample>S2</sample>
     <platformunit>PU1</platformunit>
     <center>C1</center>
     <description>blabla</description>
   </group>
 </lane>
</flowcell>
</read-groups>

```




### Example



```

$ cat input.sam 
@SQ SN:ref  LN:45
@SQ SN:ref2 LN:40
HS2000-1259_127:1:1210:15640:52255  163 ref 7   30  8M4I4M1D3M  =   37  39  
TTAGATAAAGAGGATACTG *   XX:B:S,12561,2,20,112
HS2000-1259_128:2:1210:15640:52255  0   ref 9   30  1S2I6M1P1I1P1I4M2I  *   0   
0   AAAAGATAAGGGATAAA   *

$java -jar dist/biostar78400.jar \
    -x groups.xml \
    input.sam \
   

@HD VN:1.4  SO:unsorted
@SQ SN:ref  LN:45
@SQ SN:ref2 LN:40
@RG ID:X1   PL:P1   PU:P1   LB:L1   DS:blabla   SM:S1   CN:C1
@RG ID:x2   PL:P2   PU:P2   LB:L2   DS:blabla   SM:S2   CN:C1
@PG ID:Biostar78400 PN:Biostar78400 PP:Biostar78400 VN:1.0  (...)
HS2000-1259_127:1:1210:15640:52255  163 ref 7   30  8M4I4M1D3M  =   37  39  TTAGATAAAGAGGATACTG *   RG:Z:X1 XX:B:S,12561,2,20,112
HS2000-1259_128:2:1210:15640:52255  0   ref 9   30  1S2I6M1P1I1P1I4M2I  *   0   0AAAAGATAAGGGATAAA  *   RG:Z:x2

```





