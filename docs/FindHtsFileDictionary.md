# FindHtsFileDictionary

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Scan a set of HTS files (VCF, BAM, CRAM), return a tab delimited file (path-of-file,path-to-fasta)


## Usage

```
Usage: findhtsfiledict [options] Files
  Options:
  * -D, -R, --repository, --dictionaries
      A repository is a set of record separated by a white space. Each record 
      is a line of 'key : value' or 'key=value'. Two keys are required: name 
      and fasta
    -e, --errors
      What shall we do on error
      Default: raise
      Possible Values: [ignore, raise, empty]
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

 * sam
 * bam
 * cram
 * vcf
 * dict


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew findhtsfiledict
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20190912

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/FindHtsFileDictionary.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/FindHtsFileDictionary.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **findhtsfiledict** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example



```
$ cat references.txt

name : rf
fasta : /home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.fa

name : toy
fasta : /home/lindenb/src/jvarkit-git/src/test/resources/toy.fa


$ find . -name "*.bam" -o -name "*.vcf" -o -name "*.vcf.gz" -o -name "*.cram" |\
	java -jar dist/findhtsfiledict.jar  -R references.txt

[SEVERE][FindHtsFileDictionary]No Reference found for ./e90235dca2ff2d151ef796ad0be98915/test01.vcf
java.io.IOException: No Reference found for ./e90235dca2ff2d151ef796ad0be98915/test01.vcf
	at com.github.lindenb.jvarkit.tools.misc.FindHtsFileDictionary.scan(FindHtsFileDictionary.java:193)
	at com.github.lindenb.jvarkit.tools.misc.FindHtsFileDictionary.doWork(FindHtsFileDictionary.java:213)
	at com.github.lindenb.jvarkit.util.jcommander.Launcher.instanceMain(Launcher.java:756)
	at com.github.lindenb.jvarkit.util.jcommander.Launcher.instanceMainWithExit(Launcher.java:919)
	at com.github.lindenb.jvarkit.tools.misc.FindHtsFileDictionary.main(FindHtsFileDictionary.java:241)
#hts-file	dict.name	fasta
/home/lindenb/src/jvarkit-git/./b7/83f96c410c7cd75bc732d44a1522a7/Gene_32_1507.vcf.gz	rf	/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.fa
[INFO][Launcher]findhtsfiledict Exited with failure (-1)


$ find ${PWD}/ -name "*.bam" -o -name "*.vcf" -o -name "*.vcf.gz" -o -name "*.cram" |\
	java -jar dist/findhtsfiledict.jar  -R references.txt -e ignore

#hts-file	dict.name	fasta
/home/lindenb/src/jvarkit-git/b7/83f96c410c7cd75bc732d44a1522a7/Gene_32_1507.vcf.gz	rf	/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.fa
/home/lindenb/src/jvarkit-git/db/aee9cc8f5c9c3d39c7af4cec63b7a5/Gene_0_1061.vcf.gz	rf	/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.fa
/home/lindenb/src/jvarkit-git/jeter.cram	toy	/home/lindenb/src/jvarkit-git/src/test/resources/toy.fa
/home/lindenb/src/jvarkit-git/src/test/resources/S5.bam	rf	/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.fa
/home/lindenb/src/jvarkit-git/src/test/resources/toy.vcf.gz	toy	/home/lindenb/src/jvarkit-git/src/test/resources/toy.fa
/home/lindenb/src/jvarkit-git/src/test/resources/S1.vcf.gz	rf	/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.fa
/home/lindenb/src/jvarkit-git/src/test/resources/S2.bam	rf	/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.fa
/home/lindenb/src/jvarkit-git/src/test/resources/S5.vcf.gz	rf	/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.fa
/home/lindenb/src/jvarkit-git/src/test/resources/S4.bam	rf	/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.fa
/home/lindenb/src/jvarkit-git/src/test/resources/S3.bam	rf	/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.fa
/home/lindenb/src/jvarkit-git/src/test/resources/toy.cram	toy	/home/lindenb/src/jvarkit-git/src/test/resources/toy.fa
/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.unifiedgenotyper.vcf.gz	rf	/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.fa
/home/lindenb/src/jvarkit-git/src/test/resources/S2.vcf.gz	rf	/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.fa
/home/lindenb/src/jvarkit-git/src/test/resources/S4.vcf.gz	rf	/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.fa
/home/lindenb/src/jvarkit-git/src/test/resources/S1.bam	rf	/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.fa
/home/lindenb/src/jvarkit-git/src/test/resources/toy.bam	toy	/home/lindenb/src/jvarkit-git/src/test/resources/toy.fa
/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.ann.vcf.gz	rf	/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.fa
/home/lindenb/src/jvarkit-git/src/test/resources/S3.vcf.gz	rf	/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.fa
/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.vcf.gz	rf	/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.fa
/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.freebayes.vcf.gz	rf	/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.fa
```



