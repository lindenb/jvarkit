# VcfStatsJfx

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

VCF statistics


## Usage

```
Usage: vcfstatsjfx [options] Files
  Options:
    --altering, --damaging
      For Prediction, just display children of SO:0001818 ( 
      protein_altering_variant )
      Default: false
    -fgt, --fgt
      Ignore filtered **GENOTYPES**
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --max-concordance
      Max number of concordance to display. disable if <=0
      Default: 100
    -ncl, --norm-contig-length
      For the 'contig' Panel, normalize on contig length.
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    --predictions-per-sample, -pps
      Show Predictions per sample.
      Default: false
    --prefix
      Title Prefix
      Default: <empty string>
    -s, --seconds
      Save Rscript screen every 's' seconds, if output was defined.
      Default: 15
    --stdin
      if there is no file argument. Read vcf from stdin instead of opening a 
      FileOpen dialog
      Default: false
    --trancheAffected
      tranches for the number of affected. A 'range of integers' is a list of 
      integers in ascending order separated with semicolons.
      Default: [[-Inf/0[, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, [10/20[, [20/50[, [50/100[, [100/200[, [200/300[, [300/400[, [400/500[, [500/1000[, [1000/Inf[]
    --trancheIndelSize
      tranches for the Indel size A 'range of integers' is a list of integers 
      in ascending order separated with semicolons.
      Default: [[-Inf/-1000[, [-1000/-100[, [-100/-50[, [-50/-20[, [-20/-15[, [-15/-10[, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, [10/15[, [15/20[, [20/50[, [50/100[, [100/1000[, [1000/Inf[]
    --version
      print version and exit

```


## Keywords

 * vcf
 * stats



## See also in Biostars

 * [https://www.biostars.org/p/308310](https://www.biostars.org/p/308310)
 * [https://www.biostars.org/p/353051](https://www.biostars.org/p/353051)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfstatsjfx
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfstats/VcfStatsJfx.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfstats/VcfStatsJfx.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfstatsjfx** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Examples

```
 java -jar dist/vcfstatsjfx.jar input.vcf.gz  | Rscript -
 
```


## Screenshot

  *  https://twitter.com/yokofakun/status/983280288238317568


![https://video.twimg.com/tweet_video/DaVQGvXXkAAMSBw.mp4](https://video.twimg.com/tweet_video/DaVQGvXXkAAMSBw.mp4 "animation")


## History

   *  removed JFX/gui because openjdk doesn't support jfx :-(

