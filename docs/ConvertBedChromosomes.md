# ConvertBedChromosomes

Convert the names of the chromosomes in a Bed file


## Usage

```
Usage: bedrenamechr [options] Files
  Options:
    -c, --column
      1-based chromosome column
      Default: 1
    -convert, --convert
      What should I do when  a converstion is not found
      Default: RAISE_EXCEPTION
      Possible Values: [RAISE_EXCEPTION, SKIP, RETURN_ORIGINAL]
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
  * -f, --mapping, -m
      load a custom name mapping. Format (chrom-source\tchrom-dest\n)+
    -o, --out
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * bed
 * chromosome
 * contig
 * convert


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
$ make bedrenamechr
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/ConvertBedChromosomes.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/ConvertBedChromosomes.java)


<details>
<summary>Git History</summary>

```
Mon May 15 10:41:51 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/c13a658b2ed3bc5dd6ade57190e1dab05bf70612
Fri May 12 19:41:30 2017 +0200 ; fix make, empty doc ; https://github.com/lindenb/jvarkit/commit/52fcf6d46a779fd7153ebc032fae643d2e266e7e
Wed Apr 5 13:49:50 2017 +0200 ; cont, fix bug in findallcovatpos ; https://github.com/lindenb/jvarkit/commit/7db18c7fe90fd5bf64d3ff3a4505607a1974ce6b
Thu Jun 2 09:49:17 2016 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/2ae46b7df29c6f1b66ce5104ea03bf6390db120d
Mon May 12 10:28:28 2014 +0200 ; first sed on files ; https://github.com/lindenb/jvarkit/commit/79ae202e237f53b7edb94f4326fee79b2f71b8e8
Sun Feb 2 18:55:03 2014 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/abd24b56ec986dada1e5162be5bbd0dac0c2d57c
Tue Jan 28 13:07:40 2014 +0100 ; bed rename chr ; https://github.com/lindenb/jvarkit/commit/3d1fbea5935084195d0b854089efcf571e42e0c6
```

</details>

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bedrenamechr** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



## Example

```bash
$   curl -s http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz |\
    gunzip -c | \
    java -jar dist/bedrenamechr.jar -f src/main/resources/chromnames/hg19_to_g1kv37.tsv -c 2 |\
   tail


uc011nca.2	Y	+	59213948	59276439	59230880	59274621	7	59213948,59222126,59230781,59233166,59252482,59272370,59274552,	59214117,59222281,59230919,59233257,59252550,59272463,59276439,	P51809	uc011nca.2
uc004fxl.3	Y	+	59213948	59276439	59222135	59274621	7	59213948,59222126,59230781,59233166,59252482,59272370,59274552,	59214117,59222216,59230919,59233257,59252550,59272463,59276439,	P51809-3	uc004fxl.3
uc004fxk.3	Y	+	59213948	59276439	59222135	59274809	7	59213948,59222126,59228291,59230781,59233166,59272370,59274552,	59214117,59222281,59228349,59230919,59233257,59272463,59276439,	P51809-2	uc004fxk.3
uc011ncb.2	Y	+	59213948	59276439	59222274	59274621	7	59213948,59222126,59230781,59233166,59252482,59272370,59274552,	59214117,59222277,59230919,59233257,59252550,59272463,59276439,	B4DE96	uc011ncb.2
uc010nxr.2	Y	+	59330251	59340490	59335611	59340461	9	59330251,59334000,59335576,59336119,59336347,59337090,59337948,59338753,59340193,	59330458,59334179,59335690,59336231,59336526,59337236,59338150,59338859,59340490,	B4E011	uc010nxr.2
uc004fxm.1	Y	+	59330251	59343488	59330414	59342523	9	59330251,59334078,59335552,59336119,59336354,59337119,59337948,59338753,59342486,	59330458,59334179,59335690,59336231,59336526,59337236,59338150,59338859,59343488,	B9ZVT0	uc004fxm.1
uc004fxn.1	Y	+	59330251	59343488	59330430	59343080	9	59330251,59335576,59336119,59336347,59337090,59337948,59338753,59340193,59342486,	59330458,59335690,59336231,59336526,59337236,59338150,59338859,59340278,59343488,	Q01113	uc004fxn.1
uc004fxo.1	Y	+	59352972	59356131	59353631	59356130	7	59352972,59354351,59354669,59354993,59355369,59355682,59355972,	59353819,59354463,59354816,59355130,59355505,59355884,59356131,	I3L0A4	uc004fxo.1
uc022cpg.1	Y	+	59354984	59358336	59355427	59358045	7	59354984,59355369,59355682,59355972,59356790,59357702,59357911,	59355130,59355505,59355884,59356131,59356943,59357771,59358336,	Q9NQA3	uc022cpg.1
uc011ncc.1	Y	-	59358328	59360854	59358328	59358328	3	59358328,59360006,59360500,	59359508,59360115,59360854,	uc011ncc.1

```


