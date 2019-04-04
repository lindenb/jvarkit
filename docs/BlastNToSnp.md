# BlastNToSnp

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

print indel/mismatch in a blastn stream


## Usage

```
Usage: blastn2snp [options] Files
  Options:
    -n, --gapsize
      min gap
      Default: 3
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

 * blast
 * snp



## See also in Biostars

 * [https://www.biostars.org/p/89151](https://www.biostars.org/p/89151)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew blastn2snp
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/blast/BlastNToSnp.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/blast/BlastNToSnp.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/blast/BlastNToSnpTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/blast/BlastNToSnpTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **blastn2snp** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


### Example

```
$  java -jar dist/blastn2snp.jar < blast.xml | column -t
```

```
#query              hit                                                                                              hit-index  hsp-index  query-POS  hit-POS    STRAND  REF(hit)  ALT(query)  blast.align_length  blast.hit.var  blast.query.var  blast.mid.var
No definition line  Homo sapiens chromosome 6, alternate assembly CHM1_1.1                                           1          9          21         74567818   -       T         A           18                  A              T                .
No definition line  Homo sapiens chromosome 6, alternate assembly HuRef                                              2          9          21         71600901   -       T         A           18                  A              T                .
No definition line  Homo sapiens chromosome 6, GRCh37.p13 Primary Assembly                                           3          9          21         74401398   -       T         A           18                  A              T                .
No definition line  Homo sapiens chromosome 5, alternate assembly CHM1_1.1                                           4          1          7          107821121  -       A         G           28                  T              C                .
No definition line  Homo sapiens chromosome 5, alternate assembly CHM1_1.1                                           4          9          16         14262358   +       G         C           18                  G              C                .
No definition line  Homo sapiens chromosome 5, alternate assembly CHM1_1.1                                           4          13         8          132662461  -       T         C           18                  A              G                .
No definition line  Homo sapiens chromosome 5, alternate assembly CHM1_1.1                                           4          20         14         170329095  -       G         C           18                  C              G                .
No definition line  Homo sapiens chromosome 5, alternate assembly HuRef                                              5          1          7          103561224  -       A         G           28                  T              C                .
No definition line  Homo sapiens chromosome 5, alternate assembly HuRef                                              5          9          16         14234054   +       G         C           18                  G              C                .
No definition line  Homo sapiens chromosome 5, alternate assembly HuRef                                              5          13         8          128416747  -       T         C           18                  A              G                .
No definition line  Homo sapiens chromosome 5, alternate assembly HuRef                                              5          19         14         165993804  -       G         C           18                  C              G                .
No definition line  Homo sapiens chromosome 5, GRCh37.p13 Primary Assembly                                           6          1          7          108388040  -       A         G           28                  T              C                .
No definition line  Homo sapiens chromosome 5, GRCh37.p13 Primary Assembly                                           6          9          16         14262514   +       G         C           18                  G              C                .
No definition line  Homo sapiens chromosome 5, GRCh37.p13 Primary Assembly                                           6          13         8          133231874  -       T         C           18                  A              G                .
No definition line  Homo sapiens chromosome 5, GRCh37.p13 Primary Assembly                                           6          19         14         170896537  -       G         C           18                  C              G                .
No definition line  Homo sapiens chromosome 19, alternate assembly CHM1_1.1                                          8          13         18         54878835   -       A         C           18                  T              G                .
No definition line  Homo sapiens chromosome 19, alternate assembly CHM1_1.1                                          8          14         18         54999463   +       T         G           18                  T              G                .
No definition line  Homo sapiens chromosome 19, alternate assembly CHM1_1.1                                          8          15         21         55131801   -       G         A           18                  C              T                .
No definition line  Homo sapiens chromosome 19, alternate assembly HuRef                                             9          14         18         51207279   -       A         C           18                  T              G                .
No definition line  Homo sapiens chromosome 19, alternate assembly HuRef                                             9          15         18         51329261   +       T         G           18                  T              G                .
No definition line  Homo sapiens chromosome 19, alternate assembly HuRef                                             9          16         21         51461318   -       G         A           18                  C              T                .

```

