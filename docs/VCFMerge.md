# VCFMerge

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Merge VCF Files


## DEPRECATED

use GATK combineVariants.

## Usage

```
Usage: vcfmerge [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -homref, --homref
      Use HomRef 0/0 for unknown variant
      Default: false
    --maxRecordsInRam
      When writing  files that need to be sorted, this will specify the number 
      of records stored in RAM before spilling to disk. Increasing this number 
      reduces the number of file  handles needed to sort a file, and increases 
      the amount of RAM needed
      Default: 50000
    -m, --nomerge
      Do NOT merge VariantContext lines, but create multiple lines
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    -region, --region
      Merge in that region: An interval as the following syntax : 
      "chrom:start-end" or "chrom:middle+extend"  or "chrom:start-end+extend" 
      or "chrom:start-end+extend-percent%".A program might use a Reference 
      sequence to fix the chromosome name (e.g: 1->chr1)
      Default: <empty string>
    -s, --sorted
      files are known to be ROD sorted
      Default: false
    --tmpDir
      tmp working directory. Default: java.io.tmpDir
      Default: []
    --version
      print version and exit

```


## Keywords

 * vcf
 * sort


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfmerge
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfmerge/VCFMerge.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfmerge/VCFMerge.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfmerge/VCFMergeTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfmerge/VCFMergeTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfmerge** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


 
## Example


```bash
$  find ./ -name "*.vcf.gz" | xargs java -jar dist/vcfmerge.jar   > out.vcf
```

