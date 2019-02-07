# Biostar78400

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

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
      What kind of help. One of [usage,markdown,xml].
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
      Default: SAM
      Possible Values: [BAM, SAM, CRAM]
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



## See also in Biostars

 * [https://www.biostars.org/p/78400](https://www.biostars.org/p/78400)
 * [https://www.biostars.org/p/302798](https://www.biostars.org/p/302798)
 * [https://www.biostars.org/p/202358](https://www.biostars.org/p/202358)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew biostar78400
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar78400.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar78400.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/biostar/Biostar78400Test.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/biostar/Biostar78400Test.java)


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

