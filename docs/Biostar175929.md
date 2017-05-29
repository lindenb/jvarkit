# Biostar175929

Construct a combination set of fasta sequences from a vcf see also https://www.biostars.org/p/175929/


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
      What kind of help
      Possible Values: [usage, markdown, xml]
    -o, --output
      Output file. Optional . Default: stdout
  * -R, --reference
      indexed Fasta reference
    --version
      print version and exit

```

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
$ make biostar175929
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar175929.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar175929.java
)
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




