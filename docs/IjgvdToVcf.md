# IjgvdToVcf

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Convert zips of Integrative Japanese Genome Variation to VCF file.


## DEPRECATED

Deprecated since data are now available as VCF

## Usage

```
Usage: ijgv2vcf [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
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


## DEPRECATED

Deprecated since data are now available as VCF

## DESCRIPTION

Integrative Japanese Genome Variation (iJGVD https://ijgvd.megabank.tohoku.ac.jp/ ) provides data of genomic variations obtained by whole-genome sequencing of Japanese individuals, who participate in the genome cohort study by ToMMo, IMM and other cohort projects in Japan.

> Rare variant discovery by deep whole-genome sequencing of 1,070 Japanese individuals, Nagasaki M, Yasuda J, Katsuoka F, Nariai N, Kojima K, Kawai Y, Yamaguchi-Kabat
a Y, Yokozawa J, Danjoh I, Saito S, Sato Y, Mimori T, Tsuda K, Saito R, Pan X, Nishikawa S, Ito S, Kuroki Y, Tanabe O, Fuse N, Kuriyama S, Kiyomoto H, Hozawa A, Minegi
shi N, Douglas Engel J, Kinoshita K, Kure S, Yaegashi N, ToMMo Japanese Reference Panel Project and Yamamoto M, Nat Commun, 21; 6:8018 (2015) 


##Example


```
$ ls ~/Downloads/chr* ~/Downloads/multiallelic_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip 
~/Downloads/chr10_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr11_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr1_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr12_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr13_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr14_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr15_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr16_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr17_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr18_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr19_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr20_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr21_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr2_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr22_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr3_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr4_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr5_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr6_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr7_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr8_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr9_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/multiallelic_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip


$ java -jar dist/ijgv2vcf.jar -R ~/src/jvarkit-git/src/test/resources/human_b37.dict ~/Downloads/*.zip > out.vcf

$ java -jar dist/ijgv2vcf.jar -F -M -R ~/src/jvarkit-git/src/test/resources/human_b37.dict ~/Downloads/*.zip | bcftools view -O z -o ~/Downloads/3.5KJPN_tommo_2019071.vcf.gz
$ ls -lah ~/Downloads/3.5KJPN_tommo_2019071.vcf.gz
-rw-r--r-- 1 lindenb lindenb 514M juil. 17 15:04 ~/Downloads/3.5KJPN_tommo_2019071.vcf.gz
```

