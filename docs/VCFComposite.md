# VCFComposite

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

(in developpement) Finds Variants involved in a Het Compound Disease


## Usage

```
Usage: vcfcomposite [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    -e, -E, --extractors
      Gene Extractors Name. Space/semicolon/Comma separated
      Default: ANN/GeneId VEP/GeneId
    --filter
      [20180718] set FILTER for the variants that are not part of a composite 
      mutation. Blank = filter out non-composites
      Default: NOT_COMPOSITE
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -g, --genes
      Optional tabular text report for genes
    -gf, --genotype-filter
      A Java EXpression Language (JEXL) expressions to filter a genotye in a 
      VCF. Empty string will accept all genotypes. Expression returning a TRUE 
      will accept the genotypes. See 
      https://gatkforums.broadinstitute.org/gatk/discussion/1255 
      Default: <empty string> (ACCEPT ALL)
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -l, --list
      list all available gene extractors
    -max, --max, --max-variants
      [20180718] Max variants per gene. If different than -1, used to set a 
      maximum number of variants per gene; The idea is to filter out the gene 
      having a large number of variants.
      Default: -1
    --maxRecordsInRam
      When writing  files that need to be sorted, this will specify the number 
      of records stored in RAM before spilling to disk. Increasing this number 
      reduces the number of file  handles needed to sort a file, and increases 
      the amount of RAM needed
      Default: 50000
    -u, --one-unaffected-het
      There must be at least one **Unaffected** sample with a HET Genotype for 
      a given variant.
      Default: false
    -o, --out
      Output file. Optional . Default: stdout
  * -p, -ped, --pedigree
      A pedigree file. tab delimited. Columns: family,id,father,mother, 
      sex:(0:unknown;1:male;2:female), phenotype 
      (-9|?|.:unknown;1|affected|case:affected;0|unaffected|control:unaffected) 
    -r, --report
      Optional tabular text report for pairs of variants
    -a, --select-pair
      How to select a pair of variants : any: at least one affected sample 
      must be carried by variant1 and variant2 all: all affected samples be 
      carried by variant1 and variant2
      Default: any
      Possible Values: [any, all]
    -s, --select-variant
      How to select affected samples for *one* variant. This variant will be 
      later challenged with another variant of the same gene. any: at least 
      one affected sample must be HET for the variant all: all affected 
      samples must be HET for the variant.
      Default: any
      Possible Values: [any, all]
    --tmpDir
      tmp working directory. Default: java.io.tmpDir
      Default: []
    -vf, --variant-filter
      A Java EXpression Language (JEXL) expressions to filter the variants 
      from a VCF. Empty string will accept all variants. Expression returning 
      a TRUE will accept the variant. See 
      https://gatkforums.broadinstitute.org/gatk/discussion/1255 
      Default: <empty string> (ACCEPT ALL)
    --version
      print version and exit

```


## Keywords

 * vcf
 * disease
 * annotation
 * pedigree
 * haplotype


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfcomposite
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20170331

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfcomposite/VCFComposite.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfcomposite/VCFComposite.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfcomposite/VCFCompositeTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfcomposite/VCFCompositeTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfcomposite** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Input

input is a VCF file annotated with SNPEff or VEP.

## Example

```

$ cat jeter.ped 
X	S1	S2	S3	0	1
X	S2	0	0	0	0
X	S3	0	0	0	0

$ java -jar dist/vcfcomposite.jar -r report.txt -g gene.txt -p jeter.ped src/test/resources/rotavirus_rf.ann.vcf.gz --filter "" | grep -v "##"
[WARN][VepPredictionParser]NO INFO[CSQ] found in header. This VCF was probably NOT annotated with VEP. But it's not a problem if this tool doesn't need to access VEP Annotations.
[INFO][VCFComposite]reading variants and genes
[INFO][VCFComposite]compile per gene
[INFO][VCFComposite]write variants
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5
RF02	877	.	T	A	3.45	PASS	AC=1;AN=10;ANN=A|missense_variant|MODERATE|UniProtKB/Swiss-Prot:P12472|UniProtKB/Swiss-Prot:P12472|transcript|CAA32213.1|protein_coding|1/1|c.861T>A|p.Asn287Lys|861/2643|861/2643|287/880||,A|upstream_gene_variant|MODIFIER|Gene_1621_1636|Gene_1621_1636|transcript|CAA32214.1|protein_coding||c.-745T>A|||||745|WARNING_TRANSCRIPT_INCOMPLETE,A|upstream_gene_variant|MODIFIER|UniProtKB/Swiss-Prot:P12472|UniProtKB/Swiss-Prot:P12472|transcript|CAA32215.1|protein_coding||c.-1132T>A|||||1132|WARNING_TRANSCRIPT_NO_START_CODON;BQB=1;COMPOSITE=gene|Gene_1621_1636|source|ANN_GeneId|pos|1962|ref|TACA|sample|S1,gene|UniProtKB/Swiss-Prot:P12472|source|ANN_GeneId|pos|1962|ref|TACA|sample|S1;DP=46;DP4=19,21,4,2;HOB=0.02;ICB=0.0439024;MQ=60;MQ0F=0;MQB=1;MQSB=1;RPB=0.841693;SGB=-7.90536;VDB=0.479322	GT:PL	0/1:37,0,50	0/0:0,22,116	0/0:0,22,116	0/0:0,21,94	0/0:0,12,62
RF02	1962	.	TACA	TA	33.43	PASS	AC=1;AN=10;ANN=TA|frameshift_variant|HIGH|UniProtKB/Swiss-Prot:P12472|UniProtKB/Swiss-Prot:P12472|transcript|CAA32213.1|protein_coding|1/1|c.1948_1949delCA|p.His650fs|1948/2643|1948/2643|650/880||,TA|upstream_gene_variant|MODIFIER|UniProtKB/Swiss-Prot:P12472|UniProtKB/Swiss-Prot:P12472|transcript|CAA32215.1|protein_coding||c.-45_-44delCA|||||45|WARNING_TRANSCRIPT_NO_START_CODON,TA|downstream_gene_variant|MODIFIER|Gene_1621_1636|Gene_1621_1636|transcript|CAA32214.1|protein_coding||c.*327_*328delCA|||||327|WARNING_TRANSCRIPT_INCOMPLETE;COMPOSITE=gene|Gene_1621_1636|source|ANN_GeneId|pos|877|ref|T|sample|S1,gene|UniProtKB/Swiss-Prot:P12472|source|ANN_GeneId|pos|877|ref|T|sample|S1;DP=43;DP4=22,11,2,0;HOB=0.02;ICB=0.0439024;IDV=3;IMF=0.3;INDEL;LOF=(UniProtKB/Swiss-Prot:P12472|UniProtKB/Swiss-Prot:P12472|2|0.50);MQ=60;MQ0F=0;MQSB=1;SGB=0.810227;VDB=0.373246	GT:PL	0/1:70,0,159	0/0:0,15,225	0/0:0,15,225	0/0:0,27,231	0/0:0,27,168
RF04	887	.	A	G	5.31	PASS	AC=1;AN=10;ANN=G|missense_variant|MODERATE|Gene_9_2339|Gene_9_2339|transcript|AAB07453.1|protein_coding|1/1|c.878A>G|p.Glu293Gly|878/2331|878/2331|293/776||;BQB=1;COMPOSITE=gene|Gene_9_2339|source|ANN_GeneId|pos|1857|ref|CAGA|sample|S1;DP=48;DP4=16,28,3,1;HOB=0.02;ICB=0.0439024;MQ=60;MQ0F=0;MQB=1;MQSB=1;RPB=0.90467;SGB=3.91248;VDB=0.811811	GT:PL	0/1:40,0,28	0/0:0,24,98	0/0:0,24,98	0/0:0,33,120	0/0:0,42,134
RF04	1857	.	CAGA	CA	39.47	PASS	AC=1;AN=10;ANN=CA|frameshift_variant|HIGH|Gene_9_2339|Gene_9_2339|transcript|AAB07453.1|protein_coding|1/1|c.1850_1851delGA|p.Arg617fs|1850/2331|1850/2331|617/776||;COMPOSITE=gene|Gene_9_2339|source|ANN_GeneId|pos|887|ref|A|sample|S1;DP=45;DP4=12,21,1,1;HOB=0.02;ICB=0.0439024;IDV=2;IMF=0.166667;INDEL;LOF=(Gene_9_2339|Gene_9_2339|1|1.00);MQ=60;MQ0F=0;MQSB=1;SGB=0.810227;VDB=0.969947	GT:PL	0/1:76,0,152	0/0:0,18,194	0/0:0,18,194	0/0:0,15,166	0/0:0,33,255


$ column -t gene.txt
#CHROM  bed.start  bed.end  gene.key                     gene.label                   gene.source  affected.counts  affected.total  affected.samples
RF02    876        1965     Gene_1621_1636               Gene_1621_1636               ANN_GeneId   1                1               S1
RF04    886        1860     Gene_9_2339                  Gene_9_2339                  ANN_GeneId   1                1               S1
RF02    876        1965     UniProtKB/Swiss-Prot:P12472  UniProtKB/Swiss-Prot:P12472  ANN_GeneId   1                1               S1


$ verticalize report.txt 

>>> 2
$1                       #CHROM : RF02
$2                    bed.start : 876
$3                      bed.end : 1965
$4       count.variants.in.gene : 2
$5                     gene.key : Gene_1621_1636
$6                   gene.label : Gene_1621_1636
$7                  gene.source : ANN_GeneId
$8               variant1.start : 877
$9                 variant1.end : 877
$10                variant1.ref : T
$11                variant1.alt : A
$12               variant1.info : ANN=A|missense_variant|MODERATE|UniProtKB/Swiss-Prot:P12472|UniProtKB/Swiss-Prot:P12472|transcript|CAA32213.1|protein_coding|1/1|c.861T>A|p.Asn287Lys|861/2643|861/2643|287/880||,A|upstream_gene_variant|MODIFIER|Gene_1621_1636|Gene_1621_1636|transcript|CAA32214.1|protein_coding||c.-745T>A|||||745|WARNING_TRANSCRIPT_INCOMPLETE,A|upstream_gene_variant|MODIFIER|UniProtKB/Swiss-Prot:P12472|UniProtKB/Swiss-Prot:P12472|transcript|CAA32215.1|protein_coding||c.-1132T>A|||||1132|WARNING_TRANSCRIPT_NO_START_CODON
$13    variant1.gt[S1].affected : HET
$14  variant1.gt[S2].unaffected : HOM_REF
$15  variant1.gt[S3].unaffected : HOM_REF
$16              variant2.start : 1962
$17                variant2.end : 1965
$18                variant2.ref : TACA
$19                variant2.alt : TA
$20               variant2.info : ANN=TA|frameshift_variant|HIGH|UniProtKB/Swiss-Prot:P12472|UniProtKB/Swiss-Prot:P12472|transcript|CAA32213.1|protein_coding|1/1|c.1948_1949delCA|p.His650fs|1948/2643|1948/2643|650/880||,TA|upstream_gene_variant|MODIFIER|UniProtKB/Swiss-Prot:P12472|UniProtKB/Swiss-Prot:P12472|transcript|CAA32215.1|protein_coding||c.-45_-44delCA|||||45|WARNING_TRANSCRIPT_NO_START_CODON,TA|downstream_gene_variant|MODIFIER|Gene_1621_1636|Gene_1621_1636|transcript|CAA32214.1|protein_coding||c.*327_*328delCA|||||327|WARNING_TRANSCRIPT_INCOMPLETE
$21    variant2.gt[S1].affected : HET
$22  variant2.gt[S2].unaffected : HOM_REF
$23  variant2.gt[S3].unaffected : HOM_REF
<<< 2

>>> 3
$1                       #CHROM : RF04
$2                    bed.start : 886
$3                      bed.end : 1860
$4       count.variants.in.gene : 2
$5                     gene.key : Gene_9_2339
$6                   gene.label : Gene_9_2339
$7                  gene.source : ANN_GeneId
$8               variant1.start : 887
$9                 variant1.end : 887
$10                variant1.ref : A
$11                variant1.alt : G
$12               variant1.info : ANN=G|missense_variant|MODERATE|Gene_9_2339|Gene_9_2339|transcript|AAB07453.1|protein_coding|1/1|c.878A>G|p.Glu293Gly|878/2331|878/2331|293/776||
$13    variant1.gt[S1].affected : HET
$14  variant1.gt[S2].unaffected : HOM_REF
$15  variant1.gt[S3].unaffected : HOM_REF
$16              variant2.start : 1857
$17                variant2.end : 1860
$18                variant2.ref : CAGA
$19                variant2.alt : CA
$20               variant2.info : ANN=CA|frameshift_variant|HIGH|Gene_9_2339|Gene_9_2339|transcript|AAB07453.1|protein_coding|1/1|c.1850_1851delGA|p.Arg617fs|1850/2331|1850/2331|617/776||
$21    variant2.gt[S1].affected : HET
$22  variant2.gt[S2].unaffected : HOM_REF
$23  variant2.gt[S3].unaffected : HOM_REF
<<< 3

>>> 4
$1                       #CHROM : RF02
$2                    bed.start : 876
$3                      bed.end : 1965
$4       count.variants.in.gene : 2
$5                     gene.key : UniProtKB/Swiss-Prot:P12472
$6                   gene.label : UniProtKB/Swiss-Prot:P12472
$7                  gene.source : ANN_GeneId
$8               variant1.start : 877
$9                 variant1.end : 877
$10                variant1.ref : T
$11                variant1.alt : A
$12               variant1.info : ANN=A|missense_variant|MODERATE|UniProtKB/Swiss-Prot:P12472|UniProtKB/Swiss-Prot:P12472|transcript|CAA32213.1|protein_coding|1/1|c.861T>A|p.Asn287Lys|861/2643|861/2643|287/880||,A|upstream_gene_variant|MODIFIER|Gene_1621_1636|Gene_1621_1636|transcript|CAA32214.1|protein_coding||c.-745T>A|||||745|WARNING_TRANSCRIPT_INCOMPLETE,A|upstream_gene_variant|MODIFIER|UniProtKB/Swiss-Prot:P12472|UniProtKB/Swiss-Prot:P12472|transcript|CAA32215.1|protein_coding||c.-1132T>A|||||1132|WARNING_TRANSCRIPT_NO_START_CODON
$13    variant1.gt[S1].affected : HET
$14  variant1.gt[S2].unaffected : HOM_REF
$15  variant1.gt[S3].unaffected : HOM_REF
$16              variant2.start : 1962
$17                variant2.end : 1965
$18                variant2.ref : TACA
$19                variant2.alt : TA
$20               variant2.info : ANN=TA|frameshift_variant|HIGH|UniProtKB/Swiss-Prot:P12472|UniProtKB/Swiss-Prot:P12472|transcript|CAA32213.1|protein_coding|1/1|c.1948_1949delCA|p.His650fs|1948/2643|1948/2643|650/880||,TA|upstream_gene_variant|MODIFIER|UniProtKB/Swiss-Prot:P12472|UniProtKB/Swiss-Prot:P12472|transcript|CAA32215.1|protein_coding||c.-45_-44delCA|||||45|WARNING_TRANSCRIPT_NO_START_CODON,TA|downstream_gene_variant|MODIFIER|Gene_1621_1636|Gene_1621_1636|transcript|CAA32214.1|protein_coding||c.*327_*328delCA|||||327|WARNING_TRANSCRIPT_INCOMPLETE
$21    variant2.gt[S1].affected : HET
$22  variant2.gt[S2].unaffected : HOM_REF
$23  variant2.gt[S3].unaffected : HOM_REF
<<< 4

```



