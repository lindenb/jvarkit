# ScanRetroCopy

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Scan BAM for retrocopies


## Usage

```
Usage: scanretrocopy [options] Files
  Options:
    --bam
      Optional: save matching read in this bam file
    --bedpe, -P, -J
      Optional. Save possible sites of insertion in this Bed-PE file.
    --both
      Force the constraint that both sides of a deleted intron should have at 
      least '--min-depth' reads
      Default: false
    --coding
      ignore non-coding transcripts.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -k, -K, --kg, -kg
      UCSC knownGene File/URL. The knowGene format is a compact alternative to 
      GFF/GTF because one transcript is described using only one line.	Beware 
      chromosome names are formatted the same as your REFERENCE. A typical 
      KnownGene file is 
      http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz 
      .If you only have a gff file, you can try to generate a knownGene file 
      with [http://lindenb.github.io/jvarkit/Gff2KnownGene.html](http://lindenb.github.io/jvarkit/Gff2KnownGene.html)
      Default: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/wgEncodeGencodeBasicV19.txt.gz
    -n, --min-cigar-size
      Minimal cigar element length.
      Default: 6
    --min-depth, -D
      In a transcript one must found at least 'D' reads with a clip-length> 
      'min-cigar-size'. 
      Default: 1
    -o, --output
      Output file. Optional . Default: stdout
    --partition
      Data partitioning using the SAM Read Group (see 
      https://gatkforums.broadinstitute.org/gatk/discussion/6472/ ) . It can 
      be any combination of sample, library....
      Default: sample
      Possible Values: [readgroup, sample, library, platform, center, sample_by_platform, sample_by_center, sample_by_platform_by_center, any]
  * -r, -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    --version
      print version and exit
    --bai, -bai, --with-bai
      Use random access BAM using the bai and using the knownGene data. May be 
      slow at startup
      Default: false

```


## Keywords

 * sam
 * bam
 * cigar
 * clip
 * sv
 * retrocopy


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew scanretrocopy
```

The java jar file will be installed in the `dist` directory.


## Creation Date

2019-01-25

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/retrocopy/ScanRetroCopy.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/retrocopy/ScanRetroCopy.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/retrocopy/ScanRetroCopyTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/retrocopy/ScanRetroCopyTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **scanretrocopy** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example ##

```
$  java -jar  dist/scanretrocopy.jar --bai -R human_g1k_v37.fasta \
	"http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/high_coverage_alignment/HG00096.wgs.ILLUMINA.bwa.GBR.high_cov_pcr_free.20140203.bam"
	
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG00096
2	179296982	.	C	<RETROCOPY>	108	.	AC=1;AF=0.500;AN=2;DP=108;END=179300872;MAXLEN=121;RCP=ENST00000325748.4|-|Exon1|CTCAGTTCAT|CTATATCCAA|Exon2,ENST00000438687.3|-|Exon1|CTCAGTTCAT|CTATATCCAA|Exon2,ENST00000487082.1|-|Exon1|CTCAGTTCAT|CTATATCCAA|Exon2,ENST00000432031.2|-|Exon1|CTCAGTTCAT|CTATATCCAA|Exon2;SVLEN=3891;SVTYPE=DEL	GT:DP:MAXLEN	0/1:108:121
2	179301047	.	C	<RETROCOPY>	48	.	AC=1;AF=0.500;AN=2;DP=48;END=179306337;MAXLEN=60;RCP=ENST00000325748.4|-|Exon2|CTACATTTGT|TAAAGAAATG|Exon3,ENST00000432031.2|-|Exon2|CTACATTTGT|TAAAGAAATG|Exon3,ENST00000438687.3|-|Exon2|CTACATTTGT|TAAAGAAATG|Exon3,ENST00000487082.1|-|Exon2|CTACATTTGT|TAAAGAAATG|Exon3;SVLEN=5291;SVTYPE=DEL	GT:DP:MAXLEN	0/1:48:60
	
```

## Note to self

get a report per gene:

```
find dir -type f  -name "*.tsv" -exec cut -f 4,13 '{}' ';' | awk '{G[$1]=1;S[$2]=1;H[sprintf("%s~%s",$1,$2)]=1;} END{for(g in G) {printf("%s\t",g);n=0;for(s in S) {k=sprintf("%s~%s",g,s);if(k in H){n++;printf("%s;",s);}} printf("\t%d\n",n);}}' | sort -t $'\t' -k3,3n 
```

more flexibility with a `jjs` script

```
var br = new java.io.BufferedReader( new java.io.InputStreamReader(java.lang.System.in));
var line;
var samples={};
var genes={};
var map={};
while((line=br.readLine())!=null) {
	if(line.startsWith("#")) continue;
	//if(line.contains("LowQual")) continue;
	if(line.contains("NON_CODING")) continue;
	var columns = line.split(/\t/);
	var qual=parseInt(columns[5]);


	var infos=columns[7].split(/\;/);
	var sample="";
	for(var i in infos)
		{
		var info=infos[i];
		if(info.startsWith("SAMPLES="))
			{
			var eq=info.indexOf("=");
			sample=info.substr(eq+1);
			samples[sample]=1;
			break;
			}
		}
	for(var i in infos)
		{
		var info=infos[i];
		if(info.startsWith("RCP="))
			{
			var eq=info.indexOf("=");
			var rcps=info.substring(eq+1).split(/[,]/);
			for(var j in rcps)
				{
				var gene =rcps[j].split(/\|/)[0];
				if(gene in genes)
					{
					genes[gene].score = Math.max(genes[gene].score,qual);
					}
				else
					{
					genes[gene]={
						"score":qual,
						"samples":{}
						};
					}
				genes[gene].samples[sample]=1;
				}
			}
		}
	}

var out=java.lang.System.out;
for(var gene in genes)
 {
 var n=0;
 var array=[];
 for(var sample in genes[gene].samples) {
    array.push(sample);
    n++;
    }
 if(n>1) continue;
 out.print(gene);
 out.print("\t");
 out.print(genes[gene].score);
 out.print("\t");
  out.print(array.join(";"));
 out.print("\t");
 out.print(n);
 out.println();
 }

```



