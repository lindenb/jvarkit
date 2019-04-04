# Biostar234230

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Sliding Window : discriminate partial and fully contained fragments (from a bam file)


## Usage

```
Usage: biostar234230 [options] Files
  Options:
    -filter, --filter
      A filter expression. Reads matching the expression will be filtered-out. 
      Empty String means 'filter out nothing/Accept all'. See https://github.com/lindenb/jvarkit/blob/master/src/main/resources/javacc/com/github/lindenb/jvarkit/util/bio/samfilter/SamFilterParser.jj 
      for a complete syntax.
      Default: mapqlt(1) || MapQUnavailable() || Duplicate() || FailsVendorQuality() || NotPrimaryAlignment() || SupplementaryAlignment()
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit
    -s, --winshift
      Shift each window by 's' bases
      Default: 50
    -w, --winsize
      Window size
      Default: 100

```


## See also in Biostars

 * [https://www.biostars.org/p/234230](https://www.biostars.org/p/234230)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew biostar234230
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar234230.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar234230.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **biostar234230** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



Example:


```
$ curl -s "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20110915_CEUtrio_b37_decoy_alignment/CEUTrio.HiSeq.WGS.b37_decoy.NA12892.clean.dedup.recal.bam" |  java -jar dist/biostar234230.jar 
#contig	start	end	pairs_in_window	pairs_over_window	pairs_partial_overlap
1	10000	10100	0	2	240
1	10050	10150	4	615	274
1	10100	10200	0	800	276
1	10150	10250	0	216	649
1	10200	10300	0	2982	809
1	10250	10350	0	2918	207
1	10300	10400	0	1923	2851
1	10350	10450	0	227	4498
1	10400	10500	0	31	1971
(...)

```





