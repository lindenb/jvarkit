# BlastNToSnp

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
      What kind of help
      Possible Values: [usage, markdown, xml]
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * blast
 * snp


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
$ make blastn2snp
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/blast/BlastNToSnp.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/blast/BlastNToSnp.java)


<details>
<summary>Git History</summary>

```
Fri Aug 11 15:52:02 2017 +0200 ; git history, formatting ; https://github.com/lindenb/jvarkit/commit/cf2eb57ad251cd15ae1332db9dcd062cae607d38
Wed May 24 17:27:28 2017 +0200 ; lowres bam2raster & fix doc ; https://github.com/lindenb/jvarkit/commit/6edcfd661827927b541e7267195c762e916482a0
Tue May 16 12:40:09 2017 +0200 ; doc ; https://github.com/lindenb/jvarkit/commit/ce1caf182662dc4690ec9c90e8fdd567fafa7a1e
Tue May 9 12:56:11 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/9bb79d41ffeb58983b93209b7b66484fd35da515
Tue May 17 12:25:33 2016 +0200 ; bam format was ignored ; https://github.com/lindenb/jvarkit/commit/947b48244f25bc7bedafd3ab833daec8ed4034cb
Thu Dec 10 12:31:38 2015 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/2242c107247d288754dfb47b373e3ebbd9da07f6
Mon Apr 13 10:32:12 2015 +0200 ; added option -n for https://github.com/lindenb/jvarkit/issues/26#issuecomment-92193999 ; https://github.com/lindenb/jvarkit/commit/19072426981ab0f13e61755012a4a35da591cc95
Mon May 12 10:28:28 2014 +0200 ; first sed on files ; https://github.com/lindenb/jvarkit/commit/79ae202e237f53b7edb94f4326fee79b2f71b8e8
Sun Feb 2 18:55:03 2014 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/abd24b56ec986dada1e5162be5bbd0dac0c2d57c
Tue Dec 17 11:19:39 2013 +0100 ; blastn2var ; https://github.com/lindenb/jvarkit/commit/2ed17dd4072565d55e028efb3b8d55b9e1ed66eb
Thu Nov 28 14:54:21 2013 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/6bd741fe898f5d735e5ada6b59222f8818c08baf
Wed Nov 27 20:00:16 2013 +0100 ; abstract bam filter ; https://github.com/lindenb/jvarkit/commit/6da95f7c2f27ea15634c8f3504cdc71495020248
Tue Nov 26 18:16:23 2013 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/bdb815b62e7d646360779bc136be36ebcf57a89b
```

</details>

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





### See also


 *  http://www.biostars.org/p/89151






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




