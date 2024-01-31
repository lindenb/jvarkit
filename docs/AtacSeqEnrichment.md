# AtacSeqEnrichment

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Transcription Start Site (TSS) Enrichment Score calculation


## Usage

```
Usage: java -jar dist/atacseqenrich.jar  [options] Files
Usage: atacseqenrich [options] Files
  Options:
    --bins
      normalize coverage over TSS using 'x' regions
      Default: 100
    --contig-regex
      use contigs matching this regex
      Default: (chr)?[0-9XY]+
    --extend, -x
      extend tss site by 'x' bases.A distance specified as a positive 
      integer.Commas are removed. The following suffixes are interpreted : 
      b,bp,k,kb,m,mb,g,gb 
      Default: 2000
    --gtf, -gtf
      GFF file to be used as a source of TSS start sites
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --mapq
      min mapping quality
      Default: 1
    -o, --output
      Output file. Optional . Default: stdout
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard/gatk CreateSequenceDictionary or samtools dict
    --strand
      use strand TSS information to swap the before/after position of the read 
      according to the TSS strand.
      Default: false
    --threads
      number of parallel jobs
      Default: 1
    --treshold
      comma separated tresholds to set a status concerning,acceptable/ideal 
      see https://www.encodeproject.org/atac-seq/#standards
      Default: 6,10
    --use-transcript
      use 'transcript' type in GTF instead of 'gene'
      Default: false
    --version
      print version and exit

```


## Keywords

 * bam
 * atacseq


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew atacseqenrich
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20240130

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/atacseq/AtacSeqEnrichment.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/atacseq/AtacSeqEnrichment.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **atacseqenrich** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Transcription Start Site (TSS) Enrichment Score

https://www.encodeproject.org/data-standards/terms/#enrichment

> The TSS enrichment calculation is a signal to noise calculation.
> The reads around a reference set of TSSs are collected to form an aggregate distribution of reads centered on the TSSs and extending to 2000 bp in either direction (for a total of 4000bp).
> This distribution is then normalized by taking the average read depth in the 100 bps at each of the end flanks of the distribution (for a total of 200bp of averaged data) and calculating a fold change at each position over that average read depth. This means that the flanks should start at 1, and if there is high read signal at transcription start sites (highly open regions of the genome) there should be an increase in signal up to a peak in the middle. 


## Example

```
find /path -type f -name "*.bam" > bams.list
java -jar dist/jvarkit.jar atacseqenrich -R ref.fa --gtf jeter.gtf bams.list > output.txt
```

## Output

output is a multipart text file. Each part can be isolated using, for example

```
awk '($1=="RAW")' output.txt | cut -f 2-
```


