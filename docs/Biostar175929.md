# Biostar175929

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Construct a combination set of fasta sequences from a vcf


## Usage

```
Usage: biostar175929 [options] Files
  Options:
    -b, --bracket
      Surround variant with '[' and ']'
      Default: false
    -x, --extend
      extend FASTA sequence by 'n' bases
      Default: 100
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
  * -R, --reference
      indexed Fasta reference
    --version
      print version and exit

```


## Keywords

 * fasta
 * vcf



## See also in Biostars

 * [https://www.biostars.org/p/175929](https://www.biostars.org/p/175929)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew biostar175929
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar175929.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar175929.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **biostar175929** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



### Example

```

$ java -jar dist-2.0.1/biostar175929.jar -x 2 -R ~/src/gatk-ui/testdata/ref.fa -b ~/src/gatk-ui/testdata/S1.vcf.gz  | more

>rotavirus:127|rotavirus:130-130(T)|rotavirus:232-232(T)|rotavirus:267-267(C)|rotavirus:286-286(G)|rotavirus:536-536(A)|rotavirus:693-693(T)|rotavirus:833-833(G)|rotavirus:916-916
(A)|rotavirus:961-961(T)
at[T]caatatgattacaatgaagtatttaccagagttaaaagtaaatttgattatgtga
tggatgactctggtgttaaaaacaatcttttgggtaaagctataac[T]attgatcaggc
gttaaatggaaagtttagctcag[C]tattagaaatagaaattg[G]atgactgattcta
aaacggttgctaaattagatgaagacgtgaataaacttagaatgactttatcttctaaag
ggatcgaccaaaagatgagagtacttaatgcttgttttagtgtaaaaagaataccaggaa
aatcatcatcaataattaaatgcactagacttatgaaggataaaatagaacgtggagaag
ttgaggttgatgattcatatgttgatgagaaaatggaaattgatactattgattgg[A]a
atctcgttatgatcagttagaaaaaagatttgaatcactaaaacagagggttaatgagaa
atacaatacttgggtacaaaaagcgaagaaagtaaatgaaaatatgtactctcttcagaa
tgttatctcacaacagcaaaaccaaatagcagatc[T]tcaacaatattgtagtaaattg
gaagctgatttgcaaggtaaatttagttcattagtgtcatcagttgagtggtatctaagg
tctatggaattaccagatgatgtaaagaatgacattgaacagcagttaaattcaatt[G]
atttaattaatcccattaatgctatagatgatatcgaatcgctgattagaaatttaattc
aagattatgacagaacattttt[A]atgttaaaaggactgttgaagcaatgcaactatga
atatgcata[T]tg
>rotavirus:127|rotavirus:130-130(T)|rotavirus:232-232(T)|rotavirus:267-267(C)|rotavirus:286-286(G)|rotavirus:536-536(A)|rotavirus:693-693(T)|rotavirus:833-833(G)|rotavirus:916-916
(A)|rotavirus:961-961(A)
at[T]caatatgattacaatgaagtatttaccagagttaaaagtaaatttgattatgtga
tggatgactctggtgttaaaaacaatcttttgggtaaagctataac[T]attgatcaggc
gttaaatggaaagtttagctcag[C]tattagaaatagaaattg[G]atgactgattcta
aaacggttgctaaattagatgaagacgtgaataaacttagaatgactttatcttctaaag
ggatcgaccaaaagatgagagtacttaatgcttgttttagtgtaaaaagaataccaggaa
aatcatcatcaataattaaatgcactagacttatgaaggataaaatagaacgtggagaag
ttgaggttgatgattcatatgttgatgagaaaatggaaattgatactattgattgg[A]a
atctcgttatgatcagttagaaaaaagatttgaatcactaaaacagagggttaatgagaa
atacaatacttgggtacaaaaagcgaagaaagtaaatgaaaatatgtactctcttcagaa
tgttatctcacaacagcaaaaccaaatagcagatc[T]tcaacaatattgtagtaaattg
gaagctgatttgcaaggtaaatttagttcattagtgtcatcagttgagtggtatctaagg
tctatggaattaccagatgatgtaaagaatgacattgaacagcagttaaattcaatt[G]
atttaattaatcccattaatgctatagatgatatcgaatcgctgattagaaatttaattc
aagattatgacagaacattttt[A]atgttaaaaggactgttgaagcaatgcaactatga
atatgcata[A]tg
>rotavirus:127|rotavirus:130-130(T)|rotavirus:232-232(T)|rotavirus:267-267(C)|rotavirus:286-286(G)|rotavirus:536-536(A)|rotavirus:693-693(T)|rotavirus:833-833(G)|rotavirus:916-916
(T)|rotavirus:961-961(T)
at[T]caatatgattacaatgaagtatttaccagagttaaaagtaaatttgattatgtga
tggatgactctggtgttaaaaacaatcttttgggtaaagctataac[T]attgatcaggc
gttaaatggaaagtttagctcag[C]tattagaaatagaaattg[G]atgactgattcta
aaacggttgctaaattagatgaagacgtgaataaacttagaatgactttatcttctaaag
ggatcgaccaaaagatgagagtacttaatgcttgttttagtgtaaaaagaataccaggaa
aatcatcatcaataattaaatgcactagacttatgaaggataaaatagaacgtggagaag
ttgaggttgatgattcatatgttgatgagaaaatggaaattgatactattgattgg[A]a
atctcgttatgatcagttagaaaaaagatttgaatcactaaaacagagggttaatgagaa
atacaatacttgggtacaaaaagcgaagaaagtaaatgaaaatatgtactctcttcagaa
tgttatctcacaacagcaaaaccaaatagcagatc[T]tcaacaatattgtagtaaattg
gaagctgatttgcaaggtaaatttagttcattagtgtcatcagttgagtggtatctaagg
tctatggaattaccagatgatgtaaagaatgacattgaacagcagttaaattcaatt[G]
atttaattaatcccattaatgctatagatgatatcgaatcgctgattagaaatttaattc
aagattatgacagaacattttt[T]atgttaaaaggactgttgaagcaatgcaactatga
atatgcata[T]tg
>rotavirus:127|rotavirus:130-130(T)|rotavirus:232-232(T)|rotavirus:267-267(C)|rotavirus:286-286(G)|rotavirus:536-536(A)|rotavirus:693-693(T)|rotavirus:833-833(G)|rotavirus:916-916
(T)|rotavirus:961-961(A)
at[T]caatatgattacaatgaagtatttaccagagttaaaagtaaatttgattatgtga
tggatgactctggtgttaaaaacaatcttttgggtaaagctataac[T]attgatcaggc
gttaaatggaaagtttagctcag[C]tattagaaatagaaattg[G]atgactgattcta
aaacggttgctaaattagatgaagacgtgaataaacttagaatgactttatcttctaaag
ggatcgaccaaaagatgagagtacttaatgcttgttttagtgtaaaaagaataccaggaa
aatcatcatcaataattaaatgcactagacttatgaaggataaaatagaacgtggagaag
ttgaggttgatgattcatatgttgatgagaaaatggaaattgatactattgattgg[A]a
atctcgttatgatcagttagaaaaaagatttgaatcactaaaacagagggttaatgagaa
atacaatacttgggtacaaaaagcgaagaaagtaaatgaaaatatgtactctcttcagaa
tgttatctcacaacagcaaaaccaaatagcagatc[T]tcaacaatattgtagtaaattg
gaagctgatttgcaaggtaaatttagttcattagtgtcatcagttgagtggtatctaagg
tctatggaattaccagatgatgtaaagaatgacattgaacagcagttaaattcaatt[G]
atttaattaatcccattaatgctatagatgatatcgaatcgctgattagaaatttaattc
aagattatgacagaacattttt[T]atgttaaaaggactgttgaagcaatgcaactatga
atatgcata[A]tg

```



