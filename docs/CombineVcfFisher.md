# CombineVcfFisher

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Combine multiple VCF to perform a 'vertical' fisher test.


## Usage

```
Usage: combinevcffisher [options] Files
  Options:
    --buffer
      buffer all genes/vcf in memory. Avoid to read each VCF twice or more but 
      consummes mmemory.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --line
      switch to 'line' mode. Default is 'list' mode. See doc
      Default: false
    --max-fisher
      max fisher value to be printed.
      Default: 1.0
    -o, --output
      Output file. Optional . Default: stdout
  * -p, --pedigree
      A pedigree file. tab delimited. Columns: family,id,father,mother, 
      sex:(0:unknown;1:male;2:female), phenotype 
      (-9|?|.:unknown;1|affected|case:affected;0|unaffected|control:unaffected) 
    --recursion
      pool size in 'list mode'. 1, 2 or 3.
      Default: 2
    --version
      print version and exit

```


## Keywords

 * vcf
 * burden
 * fisher


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew combinevcffisher
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20200910

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burden/CombineVcfFisher.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burden/CombineVcfFisher.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **combinevcffisher** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Input 

input is a tab delimited text file with the list of paths to the indexed VCF files.

### List mode

All vcf are stored in a list and we run the all combinations of 'recursion' VCFs.

### Line mode

all the VCFs on the line are combined and tested.

## Example

```
$ find dir -type f -name "*.vcf.gz" > paths.txt
$ java -jar ${JVARKIT_HOME}/dist/combinevcffisher.jar --buffer -p the.ped paths.txt

#genes        n-variants      cases-R cases-A controls-R      controls-A      fisher
G1:G2 13      238     4       90      10      0.0010580863010766208
G1:G3 6       241     1       94      6       0.0029829802538355616
G1:G4 6       241     1       94      6       0.0029829802538355616
G1:G5 6       241     1       94      6       0.0029829802538355616
G2:G3 6       241     1       94      6       0.0029829802538355616
G2:G4 6       241     1       94      6       0.0029829802538355616
G2:G5 6       241     1       94      6       0.0029829802538355616
G3:G4 6       241     1       94      6       0.0029829802538355616
```

