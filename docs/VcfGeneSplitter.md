# VcfGeneSplitter

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Split VCF+VEP by gene/transcript.


## Usage

```
Usage: vcfgenesplitter [options] Files
  Options:
    -e, -E, --extractors
      Gene Extractors Name. Space/semicolon/Comma separated
      Default: ANN/GeneId VEP/GeneId
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --ignore-filtered
      Ignore FILTERED variant
      Default: false
    -l, --list
      list all available extractors
    -m, --manifest
      Manifest Bed file output containing chrom/start/end of each gene
    -M, --max-variant
      Maximum number of variants required to write a vcf. don't write if 
      num(variant) > 'x' . '<=0' is ignore
      Default: -1
    --maxRecordsInRam
      When writing  files that need to be sorted, this will specify the number 
      of records stored in RAM before spilling to disk. Increasing this number 
      reduces the number of file  handles needed to sort a file, and increases 
      the amount of RAM needed
      Default: 50000
    -n, --min-variant
      Minimum number of variants required to write a vcf. don't write if 
      num(variant) < 'x'
      Default: 1
  * -o, --output
      An existing directory or a filename ending with the '.zip' suffix.
    --tmpDir
      tmp working directory. Default: java.io.tmpDir
      Default: []
    --version
      print version and exit

```


## Keywords

 * genes
 * vcf


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfgenesplitter
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfGeneSplitter.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfGeneSplitter.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/VcfGeneSplitterTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/VcfGeneSplitterTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfgenesplitter** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


### Example

```
java -jar dist/vcfgenesplitter.jar -o jeter.zip src/test/resources/rotavirus_rf.ann.vcf.gz -m jeter.mf

$ unzip -l jeter.zip 
Archive:  jeter.zip
  Length      Date    Time    Name
---------  ---------- -----   ----
     1565  2019-05-27 11:26   2c/8fb9d2539e3f30d1d9b06f9ec54c4c/Gene_18_3284.vcf.gz
     2278  2019-05-27 11:26   4e/4897c51fe2dd067a8b75c19f111477/Gene_1621_1636.vcf.gz
     2278  2019-05-27 11:26   74/ca4273c3d5803c5865891c808234da/UniProtKB_Swiss-Prot:P12472.vcf.gz
     2264  2019-05-27 11:26   23/6b59cfe4fdd33a5f4feeb55521dd34/Gene_50_2557.vcf.gz
     2169  2019-05-27 11:26   b3/4bda8d8502e64e442fce077e45ded6/Gene_9_2339.vcf.gz
     2106  2019-05-27 11:26   b7/83f96c410c7cd75bc732d44a1522a7/Gene_32_1507.vcf.gz
     2023  2019-05-27 11:26   6f/8472e9f192c92bf46e4893b2367b7e/Gene_23_1216.vcf.gz
     1862  2019-05-27 11:26   3c/513d82eaea18447dd5f621f92b40e6/Gene_0_1073.vcf.gz
     1655  2019-05-27 11:26   84/977eac8cdef861cbd3109209675d21/Gene_0_1058.vcf.gz
     1754  2019-05-27 11:26   db/aee9cc8f5c9c3d39c7af4cec63b7a5/Gene_0_1061.vcf.gz
     1746  2019-05-27 11:26   b0/133c483f0ea676f8d29ab1f2daee5d/Gene_41_568.vcf.gz
     1664  2019-05-27 11:26   59/0fd5c1e8d6d60a986a0021fe357514/Gene_20_616.vcf.gz
     1663  2019-05-27 11:26   83/bc905cf311428ab80ce59aaf503838/Gene_78_374.vcf.gz
---------                     -------
    25027                     13 files

$ cat jeter.mf
#chrom	start	end	key	path	Count_Variants
RF01	969	970	ANN/GeneId	Gene_18_3284	2c/8fb9d2539e3f30d1d9b06f9ec54c4c/Gene_18_3284.vcf.gz	1
RF02	250	1965	ANN/GeneId	Gene_1621_1636	4e/4897c51fe2dd067a8b75c19f111477/Gene_1621_1636.vcf.gz	5
RF02	250	1965	ANN/GeneId	UniProtKB/Swiss-Prot:P12472	74/ca4273c3d5803c5865891c808234da/UniProtKB_Swiss-Prot:P12472.vcf.gz	5
RF03	1220	2573	ANN/GeneId	Gene_50_2557	23/6b59cfe4fdd33a5f4feeb55521dd34/Gene_50_2557.vcf.gz	8
RF04	886	1920	ANN/GeneId	Gene_9_2339	b3/4bda8d8502e64e442fce077e45ded6/Gene_9_2339.vcf.gz	7
RF05	40	1339	ANN/GeneId	Gene_32_1507	b7/83f96c410c7cd75bc732d44a1522a7/Gene_32_1507.vcf.gz	6
RF06	516	1132	ANN/GeneId	Gene_23_1216	6f/8472e9f192c92bf46e4893b2367b7e/Gene_23_1216.vcf.gz	5
RF07	97	952	ANN/GeneId	Gene_0_1073	3c/513d82eaea18447dd5f621f92b40e6/Gene_0_1073.vcf.gz	4
RF08	925	992	ANN/GeneId	Gene_0_1058	84/977eac8cdef861cbd3109209675d21/Gene_0_1058.vcf.gz	2
RF09	293	414	ANN/GeneId	Gene_0_1061	db/aee9cc8f5c9c3d39c7af4cec63b7a5/Gene_0_1061.vcf.gz	3
RF10	45	175	ANN/GeneId	Gene_41_568	b0/133c483f0ea676f8d29ab1f2daee5d/Gene_41_568.vcf.gz	3
RF11	73	79	ANN/GeneId	Gene_20_616	59/0fd5c1e8d6d60a986a0021fe357514/Gene_20_616.vcf.gz	1
RF11	73	79	ANN/GeneId	Gene_78_374	83/bc905cf311428ab80ce59aaf503838/Gene_78_374.vcf.gz	1


```

