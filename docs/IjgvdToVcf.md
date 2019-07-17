# IjgvdToVcf

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Convert zips of Integrative Japanese Genome Variation to VCF file.


## Usage

```
Usage: ijgv2vcf [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -F, --no-filtered
      ignore 'filtered' entries
      Default: false
    -M, --no-multiallelic
      ignore 'multiallelic' entries
      Default: false
    -o, --out
      Output file. Optional . Default: stdout
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    --version
      print version and exit

```


## Keywords

 * vcf
 * jgvd
 * japan
 * tommo


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew ijgv2vcf
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20190717

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/ijgvd/IjgvdToVcf.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/ijgvd/IjgvdToVcf.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **ijgv2vcf** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


Integrative Japanese Genome Variation (iJGVD https://ijgvd.megabank.tohoku.ac.jp/ ) provides data of genomic variations obtained by whole-genome sequencing of Japanese individuals, who participate in the genome cohort study by ToMMo, IMM and other cohort projects in Japan.

> Rare variant discovery by deep whole-genome sequencing of 1,070 Japanese individuals, Nagasaki M, Yasuda J, Katsuoka F, Nariai N, Kojima K, Kawai Y, Yamaguchi-Kabat
a Y, Yokozawa J, Danjoh I, Saito S, Sato Y, Mimori T, Tsuda K, Saito R, Pan X, Nishikawa S, Ito S, Kuroki Y, Tanabe O, Fuse N, Kuriyama S, Kiyomoto H, Hozawa A, Minegi
shi N, Douglas Engel J, Kinoshita K, Kure S, Yaegashi N, ToMMo Japanese Reference Panel Project and Yamamoto M, Nat Commun, 21; 6:8018 (2015) 


##Example


```
java -jar dist/ijgv2vcf.jar -R ~/src/jvarkit-git/src/test/resources/human_b37.dict ~/Downloads/*.zip > out.vcf
```

