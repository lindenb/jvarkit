# SetFileTools

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Utilities for the setfile format


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar setfiletools  [options] Files

Usage: setfiletools [options] Files
  Options:
    --bed
      Restrict input to this bed file.
    --disable-uniq
      disable unique-name checker.
      Default: false
    --extend
      Extends each interval. Extending interval. The following syntaxes are 
      supported: 1000; 1kb; 1,000; 30%(shrink); 150% (extend); 0.5 (shrink); 
      1.5 (extend)
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --min-variants-per-setfile
      when using vcf, only keep the setfile is there are at least 'x' 
      overlapping variants.
      Default: 1
    -o, --out
      Output file. Optional . Default: stdout. For action=cluster, output is: 
      An existing directory or a filename ending with the '.zip' or '.tar' or 
      '.tar.gz' suffix.
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard/gatk CreateSequenceDictionary or samtools dict
    -U, --remove-unused-interval
      Remove
      Default: false
    --stringency
      Validation Stringency
      Default: LENIENT
      Possible Values: [STRICT, LENIENT, SILENT]
    -t, --trim-chr
      Remove chr prefix in chromosome names on output.
      Default: false
    --vcf, --vcfs
      Restrict input to thoses vcf file(s). A file with the '.list' suffix is 
      interpreted as a list of paths to the vcfs.
      Default: []
    --version
      print version and exit

```


## Keywords

 * setfile
 * bed



## Creation Date

20210125

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/setfile/SetFileTools.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/setfile/SetFileTools.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **setfiletools** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



## action = combine
```
$ echo -e "w RF01:150-200\nx RF01:1-100,RF02:1-100\ny RF01:90-200,RF03:1-100\nz RF03:50-150,RF04:100-200" | java -jar dist/setfiletools.jar -R src/test/resources/rotavirus_rf.fa combine
w_x	RF01:1-100,RF01:150-200,RF02:1-100
w_y	RF01:90-200,RF03:1-100
w_z	RF01:150-200,RF03:50-150,RF04:100-200
x_y	RF01:1-200,RF02:1-100,RF03:1-100
x_z	RF01:1-100,RF02:1-100,RF03:50-150,RF04:100-200
y_z	RF01:90-200,RF03:1-150,RF04:100-200
```

## action = bedbed

creates a setfile from two bed files. First bed is "peaks.bed" second bed is "genes.bed".
The output is a set file . Each record in the output setFile is a 'gene' where all items are the overlapping peaks. 
```
gunzip -c in.gtf.gz|\
awk -F '\t' '($3=="gene") {B=int($4)-1;X=100;printf("%s\t%d\t%s\n",$1,(B<X?0:B-X),int($5)+X);}' |\
sort -t $'\t' -k1,1 -k2,2n  |\
java -jar dist/setfiletools.jar -R ref.dict  in.bed.gz -
```

## action = intersectbed

print  whole setrecords overlapping bed file, there is to trimming

```
java -jar dist/setfiletools.jar -R ref.dict  in.bed.gz in.setfile
```



