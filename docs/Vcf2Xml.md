# Vcf2Xml

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Convert VCF to XML


## Usage

```
Usage: vcf2xml [options] Files
  Options:
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

 * vcf
 * xml


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcf2xml
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcf2xml/Vcf2Xml.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcf2xml/Vcf2Xml.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcf2xml** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



## Example

```bash
$  curl "https://raw.github.com/arq5x/gemini/master/test/test.region.vep.vcf" |\
   java -jar dist/vcf2xml.jar   |\
   xmllint --format -
```

### Result

```xml
<?xml version="1.0"?>
<vcf>
  <header>
    <infos>
      <info key="AC" countType="A">Allele count in genotypes, for each ALT allele, in the same order as listed</info>
      <info key="AF" countType="A">Allele Frequency, for each ALT allele, in the same order as listed</info>
      <info key="AN" countType="INTEGER" count="1">Total number of alleles in called genotypes</info>
      <info key="BaseQRankSum" countType="INTEGER" count="1">Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities</info>
      <info key="DP" countType="INTEGER" count="1">Approximate read depth; some reads may have been filtered</info>
      <info key="DS" countType="INTEGER" count="0">Were any of the samples downsampled?</info>
      <info key="Dels" countType="INTEGER" count="1">Fraction of Reads Containing Spanning Deletions</info>
      <info key="FS" countType="INTEGER" count="1">Phred-scaled p-value using Fisher's exact test to detect strand bias</info>
      <info key="HRun" countType="INTEGER" count="1">Largest Contiguous Homopolymer Run of Variant Allele In Either Direction</info>
      <info key="HaplotypeScore" countType="INTEGER" count="1">Consistency of the site with at most two segregating haplotypes</info>
      <info key="InbreedingCoeff" countType="INTEGER" count="1">Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation</info>
      <info key="MQ" countType="INTEGER" count="1">RMS Mapping Quality</info>
      <info key="MQ0" countType="INTEGER" count="1">Total Mapping Quality Zero Reads</info>
      <info key="MQRankSum" countType="INTEGER" count="1">Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities</info>
      <info key="QD" countType="INTEGER" count="1">Variant Confidence/Quality by Depth</info>
      <info key="ReadPosRankSum" countType="INTEGER" count="1">Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias</info>
      <info key="CSQ" countType="UNBOUNDED">Consequence type as predicted by VEP. Format: Consequence|Codons|Amino_acids|Gene|HGNC|Feature|EXON|PolyPhen|SIFT</info>
    </infos>
    <formats>
      <format key="AD" countType="UNBOUNDED">Allelic depths for the ref and alt alleles in the order listed</format>
      <format key="DP" countType="INTEGER" count="1">Approximate read depth (reads with MQ=255 or with bad mates are filtered)</format>
      <format key="GQ" countType="INTEGER" count="1">Genotype Quality</format>
      <format key="GT" countType="INTEGER" count="1">Genotype</format>
      <format key="PL" countType="G" value="">Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification</format>
    </formats>

    <filters/>
    <contigs>
      <contig tid="0">
        <key>chr1</key>
      </contig>
      <contig tid="1">
        <key>chr10</key>
      </contig>
      <!-- (...) -->
      <contig tid="92">
        <key>chrY</key>
      </contig>
    </contigs>
    <samples>
      <sample index="0">M10475</sample>
      <sample index="1">M10478</sample>
      <sample index="2">M10500</sample>
      <sample index="3">M128215</sample>
    </samples>
    <metas>
      <meta key="UnifiedGenotyper">"analysis_type=UnifiedGenotyper input_file=[bam/M10478.conc.on.pos.realigned.bam, bam/M10475.conc.on.pos.realigned.bam, bam/M10500.conc.on.pos.realigned.bam, bam/M128215.conc.on.pos.realigned.bam] read_buffer_size=null phone_home=STANDARD read_filter=[] intervals=null excludeIntervals=null interval_set_rule=UNION interval_merging=ALL reference_sequence=/m/cphg-quinlan/cphg-quinlan/shared/genomes/hg19/bwa/gatk/hg19_gatk.fa rodBind=[] nonDeterministicRandomSeed=false downsampling_type=BY_SAMPLE downsample_to_fraction=null downsample_to_coverage=250 baq=OFF baqGapOpenPenalty=40.0 performanceLog=null useOriginalQualities=false defaultBaseQualities=-1 validation_strictness=SILENT unsafe=null num_threads=10 num_cpu_threads=null num_io_threads=null num_bam_file_handles=null read_group_black_list=null pedigree=[] pedigreeString=[] pedigreeValidationType=STRICT allow_intervals_with_unindexed_bam=false logging_level=INFO log_to_file=null help=false genotype_likelihoods_model=BOTH p_nonref_model=EXACT heterozygosity=0.0010 pcr_error_rate=1.0E-4 genotyping_mode=DISCOVERY output_mode=EMIT_VARIANTS_ONLY standard_min_confidence_threshold_for_calling=30.0 standard_min_confidence_threshold_for_emitting=30.0 computeSLOD=false alleles=(RodBinding name= source=UNBOUND) min_base_quality_score=17 max_deletion_fraction=0.05 multiallelic=false max_alternate_alleles=5 min_indel_count_for_genotyping=5 indel_heterozygosity=1.25E-4 indelGapContinuationPenalty=10.0 indelGapOpenPenalty=45.0 indelHaplotypeSize=80 bandedIndel=false indelDebug=false ignoreSNPAlleles=false dbsnp=(RodBinding name= source=UNBOUND) out=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub NO_HEADER=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub sites_only=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub debug_file=null metrics_file=null annotation=[] excludeAnnotation=[] filter_mismatching_base_and_quals=false"</meta>
      <meta key="reference">file:///m/cphg-quinlan/cphg-quinlan/shared/genomes/hg19/bwa/gatk/hg19_gatk.fa</meta>
    </metas>
  </header>
  <variations>
    <variation>
      <chrom>chr1</chrom>
      <start>10001</start>
      <end>10001</end>
      <ref>T</ref>
      <alt>TC</alt>
      <qual>175.91000000000003</qual>
      <infos>
        <BaseQRankSum>4.975</BaseQRankSum>
        <HaplotypeScore>218.6157</HaplotypeScore>
        <CSQ>
          <Consequence>upstream_gene_variant</Consequence>
          <Gene>ENSG00000223972</Gene>
          <HGNC>DDX11L1</HGNC>
          <Feature>ENST00000456328</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000227232</Gene>
          <HGNC>WASH7P</HGNC>
          <Feature>ENST00000488147</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000227232</Gene>
          <HGNC>WASH7P</HGNC>
          <Feature>ENST00000541675</Feature>
        </CSQ>
        <CSQ>
          <Consequence>upstream_gene_variant</Consequence>
          <Gene>ENSG00000223972</Gene>
          <HGNC>DDX11L1</HGNC>
          <Feature>ENST00000450305</Feature>
        </CSQ>
        <CSQ>
          <Consequence>upstream_gene_variant</Consequence>
          <Gene>ENSG00000223972</Gene>
          <HGNC>DDX11L1</HGNC>
          <Feature>ENST00000515242</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000227232</Gene>
          <HGNC>WASH7P</HGNC>
          <Feature>ENST00000538476</Feature>
        </CSQ>
        <CSQ>
          <Consequence>upstream_gene_variant</Consequence>
          <Gene>ENSG00000223972</Gene>
          <HGNC>DDX11L1</HGNC>
          <Feature>ENST00000518655</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000227232</Gene>
          <HGNC>WASH7P</HGNC>
          <Feature>ENST00000438504</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000227232</Gene>
          <HGNC>WASH7P</HGNC>
          <Feature>ENST00000423562</Feature>
        </CSQ>
        <QD>2.31</QD>
        <MQ>35.31</MQ>
        <AC>4</AC>
        <FS>12.516</FS>
        <HRun>0</HRun>
        <MQRankSum>-0.238</MQRankSum>
        <ReadPosRankSum>2.910</ReadPosRankSum>
        <DP>76</DP>
        <AF>0.50</AF>
        <MQ0>0</MQ0>
        <AN>8</AN>
      </infos>
      <genotypes>
        <genotype available="true" called="true" het="true" hom="false" homRef="false" homVar="false" mixed="false" noCall="false" nonInformative="false" filtered="false" phased="false" sample="M128215">
          <AD>
            <value>15</value>
            <value>5</value>
          </AD>
          <DP>24</DP>
          <GQ>0</GQ>
          <PL>
            <value index="1">49</value>
            <value index="2">0</value>
            <value index="3">0</value>
          </PL>
          <alleles>
            <allele ref="true">T</allele>
            <allele>TC</allele>
          </alleles>
        </genotype>
        <genotype available="true" called="true" het="true" hom="false" homRef="false" homVar="false" mixed="false" noCall="false" nonInformative="false" filtered="false" phased="false" sample="M10475">
          <AD>
            <value>10</value>
            <value>2</value>
          </AD>
          <DP>15</DP>
          <GQ>10</GQ>
          <PL>
            <value index="1">25</value>
            <value index="2">0</value>
            <value index="3">10</value>
          </PL>
          <alleles>
            <allele ref="true">T</allele>
            <allele>TC</allele>
          </alleles>
        </genotype>
        <genotype available="true" called="true" het="true" hom="false" homRef="false" homVar="false" mixed="false" noCall="false" nonInformative="false" filtered="false" phased="false" sample="M10500">
          <AD>
            <value>10</value>
            <value>10</value>
          </AD>
          <DP>21</DP>
          <GQ>7</GQ>
          <PL>
            <value index="1">111</value>
            <value index="2">0</value>
            <value index="3">7</value>
          </PL>
          <alleles>
            <allele ref="true">T</allele>
            <allele>TC</allele>
          </alleles>
        </genotype>
        <genotype available="true" called="true" het="true" hom="false" homRef="false" homVar="false" mixed="false" noCall="false" nonInformative="false" filtered="false" phased="false" sample="M10478">
          <AD>
            <value>10</value>
            <value>4</value>
          </AD>
          <DP>16</DP>
          <GQ>5</GQ>
          <PL>
            <value index="1">40</value>
            <value index="2">0</value>
            <value index="3">5</value>
          </PL>
          <alleles>
            <allele ref="true">T</allele>
            <allele>TC</allele>
          </alleles>
        </genotype>
      </genotypes>
    </variation>
    <variation>
      <chrom>chr1</chrom>
      <start>10056</start>
      <end>10056</end>
      <ref>A</ref>
      <alt>C</alt>
      <qual>47.27</qual>
      <infos>
        <BaseQRankSum>-1.865</BaseQRankSum>
        <HaplotypeScore>10.0031</HaplotypeScore>
        <CSQ>
          <Consequence>upstream_gene_variant</Consequence>
          <Gene>ENSG00000223972</Gene>
          <HGNC>DDX11L1</HGNC>
          <Feature>ENST00000456328</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000227232</Gene>
          <HGNC>WASH7P</HGNC>
          <Feature>ENST00000488147</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000227232</Gene>
          <HGNC>WASH7P</HGNC>
          <Feature>ENST00000541675</Feature>
        </CSQ>
        <CSQ>
          <Consequence>upstream_gene_variant</Consequence>
          <Gene>ENSG00000223972</Gene>
          <HGNC>DDX11L1</HGNC>
          <Feature>ENST00000450305</Feature>
        </CSQ>
        <CSQ>
          <Consequence>upstream_gene_variant</Consequence>
          <Gene>ENSG00000223972</Gene>
          <HGNC>DDX11L1</HGNC>
          <Feature>ENST00000515242</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000227232</Gene>
          <HGNC>WASH7P</HGNC>
          <Feature>ENST00000538476</Feature>
        </CSQ>
        <CSQ>
          <Consequence>upstream_gene_variant</Consequence>
          <Gene>ENSG00000223972</Gene>
          <HGNC>DDX11L1</HGNC>
          <Feature>ENST00000518655</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000227232</Gene>
          <HGNC>WASH7P</HGNC>
          <Feature>ENST00000438504</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000227232</Gene>
          <HGNC>WASH7P</HGNC>
          <Feature>ENST00000423562</Feature>
        </CSQ>
        <QD>0.65</QD>
        <MQ>30.59</MQ>
        <AC>1</AC>
        <FS>7.175</FS>
        <HRun>0</HRun>
        <MQRankSum>-1.940</MQRankSum>
        <ReadPosRankSum>1.402</ReadPosRankSum>
        <Dels>0.00</Dels>
        <DP>199</DP>
        <AF>0.13</AF>
        <MQ0>0</MQ0>
        <AN>8</AN>
      </infos>
      <genotypes>
        <genotype available="true" called="true" het="false" hom="true" homRef="true" homVar="false" mixed="false" noCall="false" nonInformative="false" filtered="false" phased="false" sample="M128215">
          <AD>
            <value>56</value>
            <value>0</value>
          </AD>
          <DP>56</DP>
          <GQ>99</GQ>
          <PL>
            <value index="1">0</value>
            <value index="2">144</value>
            <value index="3">1518</value>
          </PL>
          <alleles>
            <allele ref="true">A</allele>
            <allele ref="true">A</allele>
          </alleles>
        </genotype>
        <genotype available="true" called="true" het="false" hom="true" homRef="true" homVar="false" mixed="false" noCall="false" nonInformative="false" filtered="false" phased="false" sample="M10475">
          <AD>
            <value>35</value>
            <value>2</value>
          </AD>
          <DP>39</DP>
          <GQ>47</GQ>
          <PL>
            <value index="1">0</value>
            <value index="2">47</value>
            <value index="3">739</value>
          </PL>
          <alleles>
            <allele ref="true">A</allele>
            <allele ref="true">A</allele>
          </alleles>
        </genotype>
        <genotype available="true" called="true" het="false" hom="true" homRef="true" homVar="false" mixed="false" noCall="false" nonInformative="false" filtered="false" phased="false" sample="M10500">
          <AD>
            <value>30</value>
            <value>0</value>
          </AD>
          <DP>31</DP>
          <GQ>90</GQ>
          <PL>
            <value index="1">0</value>
            <value index="2">90</value>
            <value index="3">983</value>
          </PL>
          <alleles>
            <allele ref="true">A</allele>
            <allele ref="true">A</allele>
          </alleles>
        </genotype>
        <genotype available="true" called="true" het="true" hom="false" homRef="false" homVar="false" mixed="false" noCall="false" nonInformative="false" filtered="false" phased="false" sample="M10478">
          <AD>
            <value>62</value>
            <value>11</value>
          </AD>
          <DP>73</DP>
          <GQ>83</GQ>
          <PL>
            <value index="1">83</value>
            <value index="2">0</value>
            <value index="3">1414</value>
          </PL>
          <alleles>
            <allele ref="true">A</allele>
            <allele>C</allele>
          </alleles>
        </genotype>
      </genotypes>
    </variation>
    <variation>
      <chrom>chr1</chrom>
      <start>10121</start>
      <end>10121</end>
      <ref>A</ref>
      <alt>C</alt>
      <qual>111.08000000000001</qual>
      <infos>
        <BaseQRankSum>0.234</BaseQRankSum>
        <HaplotypeScore>9.1907</HaplotypeScore>
        <CSQ>
          <Consequence>upstream_gene_variant</Consequence>
          <Gene>ENSG00000223972</Gene>
          <HGNC>DDX11L1</HGNC>
          <Feature>ENST00000456328</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000227232</Gene>
          <HGNC>WASH7P</HGNC>
          <Feature>ENST00000488147</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000227232</Gene>
          <HGNC>WASH7P</HGNC>
          <Feature>ENST00000541675</Feature>
        </CSQ>
        <CSQ>
          <Consequence>upstream_gene_variant</Consequence>
          <Gene>ENSG00000223972</Gene>
          <HGNC>DDX11L1</HGNC>
          <Feature>ENST00000450305</Feature>
        </CSQ>
        <CSQ>
          <Consequence>upstream_gene_variant</Consequence>
          <Gene>ENSG00000223972</Gene>
          <HGNC>DDX11L1</HGNC>
          <Feature>ENST00000515242</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000227232</Gene>
          <HGNC>WASH7P</HGNC>
          <Feature>ENST00000538476</Feature>
        </CSQ>
        <CSQ>
          <Consequence>upstream_gene_variant</Consequence>
          <Gene>ENSG00000223972</Gene>
          <HGNC>DDX11L1</HGNC>
          <Feature>ENST00000518655</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000227232</Gene>
          <HGNC>WASH7P</HGNC>
          <Feature>ENST00000438504</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000227232</Gene>
          <HGNC>WASH7P</HGNC>
          <Feature>ENST00000423562</Feature>
        </CSQ>
        <QD>4.63</QD>
        <MQ>29.08</MQ>
        <AC>1</AC>
        <FS>16.063</FS>
        <HRun>0</HRun>
        <MQRankSum>0.521</MQRankSum>
        <ReadPosRankSum>-1.401</ReadPosRankSum>
        <Dels>0.00</Dels>
        <DP>102</DP>
        <AF>0.13</AF>
        <MQ0>0</MQ0>
        <AN>8</AN>
      </infos>
      <genotypes>
        <genotype available="true" called="true" het="true" hom="false" homRef="false" homVar="false" mixed="false" noCall="false" nonInformative="false" filtered="false" phased="false" sample="M128215">
          <AD>
            <value>15</value>
            <value>9</value>
          </AD>
          <DP>24</DP>
          <GQ>99</GQ>
          <PL>
            <value index="1">147</value>
            <value index="2">0</value>
            <value index="3">228</value>
          </PL>
          <alleles>
            <allele ref="true">A</allele>
            <allele>C</allele>
          </alleles>
        </genotype>
        <genotype available="true" called="true" het="false" hom="true" homRef="true" homVar="false" mixed="false" noCall="false" nonInformative="false" filtered="false" phased="false" sample="M10475">
          <AD>
            <value>20</value>
            <value>0</value>
          </AD>
          <DP>20</DP>
          <GQ>39</GQ>
          <PL>
            <value index="1">0</value>
            <value index="2">39</value>
            <value index="3">379</value>
          </PL>
          <alleles>
            <allele ref="true">A</allele>
            <allele ref="true">A</allele>
          </alleles>
        </genotype>
        <genotype available="true" called="true" het="false" hom="true" homRef="true" homVar="false" mixed="false" noCall="false" nonInformative="false" filtered="false" phased="false" sample="M10500">
          <AD>
            <value>10</value>
            <value>0</value>
          </AD>
          <DP>10</DP>
          <GQ>30</GQ>
          <PL>
            <value index="1">0</value>
            <value index="2">30</value>
            <value index="3">305</value>
          </PL>
          <alleles>
            <allele ref="true">A</allele>
            <allele ref="true">A</allele>
          </alleles>
        </genotype>
        <genotype available="true" called="true" het="false" hom="true" homRef="true" homVar="false" mixed="false" noCall="false" nonInformative="false" filtered="false" phased="false" sample="M10478">
          <AD>
            <value>42</value>
            <value>1</value>
          </AD>
          <DP>48</DP>
          <GQ>99</GQ>
          <PL>
            <value index="1">0</value>
            <value index="2">102</value>
            <value index="3">1011</value>
          </PL>
          <alleles>
            <allele ref="true">A</allele>
            <allele ref="true">A</allele>
          </alleles>
        </genotype>
      </genotypes>
    </variation>
    <variation>
      <chrom>chr1</chrom>
      <start>10177</start>
      <end>10177</end>
      <ref>A</ref>
      <alt>C</alt>
      <qual>163.46</qual>
      <infos>
        <BaseQRankSum>1.016</BaseQRankSum>
        <HaplotypeScore>4.4949</HaplotypeScore>
        <CSQ>
          <Consequence>upstream_gene_variant</Consequence>
          <Gene>ENSG00000223972</Gene>
          <HGNC>DDX11L1</HGNC>
          <Feature>ENST00000456328</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000227232</Gene>
          <HGNC>WASH7P</HGNC>
          <Feature>ENST00000488147</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000227232</Gene>
          <HGNC>WASH7P</HGNC>
          <Feature>ENST00000541675</Feature>
        </CSQ>
        <CSQ>
          <Consequence>upstream_gene_variant</Consequence>
          <Gene>ENSG00000223972</Gene>
          <HGNC>DDX11L1</HGNC>
          <Feature>ENST00000450305</Feature>
        </CSQ>
        <CSQ>
          <Consequence>upstream_gene_variant</Consequence>
          <Gene>ENSG00000223972</Gene>
          <HGNC>DDX11L1</HGNC>
          <Feature>ENST00000515242</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000227232</Gene>
          <HGNC>WASH7P</HGNC>
          <Feature>ENST00000538476</Feature>
        </CSQ>
        <CSQ>
          <Consequence>upstream_gene_variant</Consequence>
          <Gene>ENSG00000223972</Gene>
          <HGNC>DDX11L1</HGNC>
          <Feature>ENST00000518655</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000227232</Gene>
          <HGNC>WASH7P</HGNC>
          <Feature>ENST00000438504</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000227232</Gene>
          <HGNC>WASH7P</HGNC>
          <Feature>ENST00000423562</Feature>
        </CSQ>
        <QD>3.21</QD>
        <MQ>28.39</MQ>
        <AC>3</AC>
        <FS>6.680</FS>
        <HRun>2</HRun>
        <MQRankSum>2.046</MQRankSum>
        <ReadPosRankSum>-2.074</ReadPosRankSum>
        <Dels>0.00</Dels>
        <DP>57</DP>
        <AF>0.38</AF>
        <MQ0>0</MQ0>
        <AN>8</AN>
      </infos>
      <genotypes>
        <genotype available="true" called="true" het="true" hom="false" homRef="false" homVar="false" mixed="false" noCall="false" nonInformative="false" filtered="false" phased="false" sample="M128215">
          <AD>
            <value>4</value>
            <value>3</value>
          </AD>
          <DP>7</DP>
          <GQ>18</GQ>
          <PL>
            <value index="1">81</value>
            <value index="2">0</value>
            <value index="3">18</value>
          </PL>
          <alleles>
            <allele ref="true">A</allele>
            <allele>C</allele>
          </alleles>
        </genotype>
        <genotype available="true" called="true" het="true" hom="false" homRef="false" homVar="false" mixed="false" noCall="false" nonInformative="false" filtered="false" phased="false" sample="M10475">
          <AD>
            <value>10</value>
            <value>4</value>
          </AD>
          <DP>14</DP>
          <GQ>56</GQ>
          <PL>
            <value index="1">56</value>
            <value index="2">0</value>
            <value index="3">143</value>
          </PL>
          <alleles>
            <allele ref="true">A</allele>
            <allele>C</allele>
          </alleles>
        </genotype>
        <genotype available="true" called="true" het="false" hom="true" homRef="true" homVar="false" mixed="false" noCall="false" nonInformative="false" filtered="false" phased="false" sample="M10500">
          <AD>
            <value>6</value>
            <value>0</value>
          </AD>
          <DP>6</DP>
          <GQ>15</GQ>
          <PL>
            <value index="1">0</value>
            <value index="2">15</value>
            <value index="3">125</value>
          </PL>
          <alleles>
            <allele ref="true">A</allele>
            <allele ref="true">A</allele>
          </alleles>
        </genotype>
        <genotype available="true" called="true" het="true" hom="false" homRef="false" homVar="false" mixed="false" noCall="false" nonInformative="false" filtered="false" phased="false" sample="M10478">
          <AD>
            <value>24</value>
            <value>6</value>
          </AD>
          <DP>30</DP>
          <GQ>69</GQ>
          <PL>
            <value index="1">69</value>
            <value index="2">0</value>
            <value index="3">329</value>
          </PL>
          <alleles>
            <allele ref="true">A</allele>
            <allele>C</allele>
          </alleles>
        </genotype>
      </genotypes>
    </variation>
    <variation>
      <chrom>chr1</chrom>
      <start>10180</start>
      <end>10180</end>
      <ref>T</ref>
      <alt>C</alt>
      <qual>79.53</qual>
      <infos>
        <BaseQRankSum>-0.295</BaseQRankSum>
        <HaplotypeScore>4.7635</HaplotypeScore>
        <CSQ>
          <Consequence>upstream_gene_variant</Consequence>
          <Gene>ENSG00000223972</Gene>
          <HGNC>DDX11L1</HGNC>
          <Feature>ENST00000456328</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000227232</Gene>
          <HGNC>WASH7P</HGNC>
          <Feature>ENST00000488147</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000227232</Gene>
          <HGNC>WASH7P</HGNC>
          <Feature>ENST00000541675</Feature>
        </CSQ>
        <CSQ>
          <Consequence>upstream_gene_variant</Consequence>
          <Gene>ENSG00000223972</Gene>
          <HGNC>DDX11L1</HGNC>
          <Feature>ENST00000450305</Feature>
        </CSQ>
        <CSQ>
          <Consequence>upstream_gene_variant</Consequence>
          <Gene>ENSG00000223972</Gene>
          <HGNC>DDX11L1</HGNC>
          <Feature>ENST00000515242</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000227232</Gene>
          <HGNC>WASH7P</HGNC>
          <Feature>ENST00000538476</Feature>
        </CSQ>
        <CSQ>
          <Consequence>upstream_gene_variant</Consequence>
          <Gene>ENSG00000223972</Gene>
          <HGNC>DDX11L1</HGNC>
          <Feature>ENST00000518655</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000227232</Gene>
          <HGNC>WASH7P</HGNC>
          <Feature>ENST00000438504</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000227232</Gene>
          <HGNC>WASH7P</HGNC>
          <Feature>ENST00000423562</Feature>
        </CSQ>
        <QD>2.41</QD>
        <MQ>28.08</MQ>
        <AC>2</AC>
        <FS>1.926</FS>
        <HRun>2</HRun>
        <MQRankSum>0.710</MQRankSum>
        <ReadPosRankSum>-0.814</ReadPosRankSum>
        <Dels>0.00</Dels>
        <DP>54</DP>
        <AF>0.25</AF>
        <MQ0>0</MQ0>
        <AN>8</AN>
      </infos>
      <genotypes>
        <genotype available="true" called="true" het="false" hom="true" homRef="true" homVar="false" mixed="false" noCall="false" nonInformative="false" filtered="false" phased="false" sample="M128215">
          <AD>
            <value>7</value>
            <value>0</value>
          </AD>
          <DP>7</DP>
          <GQ>12</GQ>
          <PL>
            <value index="1">0</value>
            <value index="2">12</value>
            <value index="3">121</value>
          </PL>
          <alleles>
            <allele ref="true">T</allele>
            <allele ref="true">T</allele>
          </alleles>
        </genotype>
        <genotype available="true" called="true" het="false" hom="true" homRef="true" homVar="false" mixed="false" noCall="false" nonInformative="false" filtered="false" phased="false" sample="M10475">
          <AD>
            <value>12</value>
            <value>2</value>
          </AD>
          <DP>14</DP>
          <GQ>5</GQ>
          <PL>
            <value index="1">0</value>
            <value index="2">5</value>
            <value index="3">215</value>
          </PL>
          <alleles>
            <allele ref="true">T</allele>
            <allele ref="true">T</allele>
          </alleles>
        </genotype>
        <genotype available="true" called="true" het="true" hom="false" homRef="false" homVar="false" mixed="false" noCall="false" nonInformative="false" filtered="false" phased="false" sample="M10500">
          <AD>
            <value>6</value>
            <value>1</value>
          </AD>
          <DP>7</DP>
          <GQ>13</GQ>
          <PL>
            <value index="1">13</value>
            <value index="2">0</value>
            <value index="3">95</value>
          </PL>
          <alleles>
            <allele ref="true">T</allele>
            <allele>C</allele>
          </alleles>
        </genotype>
        <genotype available="true" called="true" het="true" hom="false" homRef="false" homVar="false" mixed="false" noCall="false" nonInformative="false" filtered="false" phased="false" sample="M10478">
          <AD>
            <value>17</value>
            <value>9</value>
          </AD>
          <DP>26</DP>
          <GQ>99</GQ>
          <PL>
            <value index="1">107</value>
            <value index="2">0</value>
            <value index="3">218</value>
          </PL>
          <alleles>
            <allele ref="true">T</allele>
            <allele>C</allele>
          </alleles>
        </genotype>
      </genotypes>
    </variation>
    <variation>
      <chrom>chr1</chrom>
      <start>10234</start>
      <end>10234</end>
      <ref>C</ref>
      <alt>T</alt>
      <qual>49.1</qual>
      <infos>
        <BaseQRankSum>-2.564</BaseQRankSum>
        <HaplotypeScore>7.4331</HaplotypeScore>
        <CSQ>
          <Consequence>upstream_gene_variant</Consequence>
          <Gene>ENSG00000223972</Gene>
          <HGNC>DDX11L1</HGNC>
          <Feature>ENST00000456328</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000227232</Gene>
          <HGNC>WASH7P</HGNC>
          <Feature>ENST00000488147</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000227232</Gene>
          <HGNC>WASH7P</HGNC>
          <Feature>ENST00000541675</Feature>
        </CSQ>
        <CSQ>
          <Consequence>upstream_gene_variant</Consequence>
          <Gene>ENSG00000223972</Gene>
          <HGNC>DDX11L1</HGNC>
          <Feature>ENST00000450305</Feature>
        </CSQ>
        <CSQ>
          <Consequence>upstream_gene_variant</Consequence>
          <Gene>ENSG00000223972</Gene>
          <HGNC>DDX11L1</HGNC>
          <Feature>ENST00000515242</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000227232</Gene>
          <HGNC>WASH7P</HGNC>
          <Feature>ENST00000538476</Feature>
        </CSQ>
        <CSQ>
          <Consequence>upstream_gene_variant</Consequence>
          <Gene>ENSG00000223972</Gene>
          <HGNC>DDX11L1</HGNC>
          <Feature>ENST00000518655</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000227232</Gene>
          <HGNC>WASH7P</HGNC>
          <Feature>ENST00000438504</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000227232</Gene>
          <HGNC>WASH7P</HGNC>
          <Feature>ENST00000423562</Feature>
        </CSQ>
        <QD>1.36</QD>
        <MQ>31.63</MQ>
        <AC>2</AC>
        <FS>5.371</FS>
        <HRun>1</HRun>
        <MQRankSum>1.310</MQRankSum>
        <ReadPosRankSum>3.118</ReadPosRankSum>
        <Dels>0.00</Dels>
        <DP>45</DP>
        <AF>0.25</AF>
        <MQ0>0</MQ0>
        <AN>8</AN>
      </infos>
      <genotypes>
        <genotype available="true" called="true" het="false" hom="true" homRef="true" homVar="false" mixed="false" noCall="false" nonInformative="false" filtered="false" phased="false" sample="M128215">
          <AD>
            <value>3</value>
            <value>2</value>
          </AD>
          <DP>5</DP>
          <GQ>6</GQ>
          <PL>
            <value index="1">0</value>
            <value index="2">6</value>
            <value index="3">68</value>
          </PL>
          <alleles>
            <allele ref="true">C</allele>
            <allele ref="true">C</allele>
          </alleles>
        </genotype>
        <genotype available="true" called="true" het="true" hom="false" homRef="false" homVar="false" mixed="false" noCall="false" nonInformative="false" filtered="false" phased="false" sample="M10475">
          <AD>
            <value>6</value>
            <value>1</value>
          </AD>
          <DP>7</DP>
          <GQ>19</GQ>
          <PL>
            <value index="1">19</value>
            <value index="2">0</value>
            <value index="3">80</value>
          </PL>
          <alleles>
            <allele ref="true">C</allele>
            <allele>T</allele>
          </alleles>
        </genotype>
        <genotype available="true" called="true" het="false" hom="true" homRef="true" homVar="false" mixed="false" noCall="false" nonInformative="false" filtered="false" phased="false" sample="M10500">
          <AD>
            <value>4</value>
            <value>0</value>
          </AD>
          <DP>4</DP>
          <GQ>12</GQ>
          <PL>
            <value index="1">0</value>
            <value index="2">12</value>
            <value index="3">137</value>
          </PL>
          <alleles>
            <allele ref="true">C</allele>
            <allele ref="true">C</allele>
          </alleles>
        </genotype>
        <genotype available="true" called="true" het="true" hom="false" homRef="false" homVar="false" mixed="false" noCall="false" nonInformative="false" filtered="false" phased="false" sample="M10478">
          <AD>
            <value>22</value>
            <value>6</value>
          </AD>
          <DP>29</DP>
          <GQ>71</GQ>
          <PL>
            <value index="1">71</value>
            <value index="2">0</value>
            <value index="3">463</value>
          </PL>
          <alleles>
            <allele ref="true">C</allele>
            <allele>T</allele>
          </alleles>
        </genotype>
      </genotypes>
    </variation>
    <variation>
      <chrom>chr16</chrom>
      <start>72057282</start>
      <end>72057282</end>
      <ref>A</ref>
      <alt>G</alt>
      <qual>1098.89</qual>
      <infos>
        <BaseQRankSum>-0.026</BaseQRankSum>
        <HaplotypeScore>1.9993</HaplotypeScore>
        <CSQ>
          <Consequence>intron_variant</Consequence>
          <Gene>ENSG00000102967</Gene>
          <HGNC>DHODH</HGNC>
          <Feature>ENST00000219240</Feature>
        </CSQ>
        <CSQ>
          <Consequence>intron_variant</Consequence>
          <Consequence>nc_transcript_variant</Consequence>
          <Gene>ENSG00000102967</Gene>
          <HGNC>DHODH</HGNC>
          <Feature>ENST00000571392</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000102967</Gene>
          <HGNC>DHODH</HGNC>
          <Feature>ENST00000572003</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000102967</Gene>
          <HGNC>DHODH</HGNC>
          <Feature>ENST00000573843</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000102967</Gene>
          <HGNC>DHODH</HGNC>
          <Feature>ENST00000573922</Feature>
        </CSQ>
        <CSQ>
          <Consequence>intron_variant</Consequence>
          <Gene>ENSG00000102967</Gene>
          <HGNC>DHODH</HGNC>
          <Feature>ENST00000574309</Feature>
        </CSQ>
        <CSQ>
          <Consequence>intron_variant</Consequence>
          <Gene>ENSG00000102967</Gene>
          <HGNC>DHODH</HGNC>
          <Feature>ENST00000572887</Feature>
        </CSQ>
        <QD>13.74</QD>
        <MQ>36.24</MQ>
        <AC>5</AC>
        <FS>7.890</FS>
        <HRun>3</HRun>
        <MQRankSum>0.529</MQRankSum>
        <ReadPosRankSum>1.005</ReadPosRankSum>
        <Dels>0.00</Dels>
        <DP>80</DP>
        <AF>0.63</AF>
        <MQ0>0</MQ0>
        <AN>8</AN>
      </infos>
      <genotypes>
        <genotype available="true" called="true" het="true" hom="false" homRef="false" homVar="false" mixed="false" noCall="false" nonInformative="false" filtered="false" phased="false" sample="M128215">
          <AD>
            <value>15</value>
            <value>23</value>
          </AD>
          <DP>40</DP>
          <GQ>99</GQ>
          <PL>
            <value index="1">501</value>
            <value index="2">0</value>
            <value index="3">332</value>
          </PL>
          <alleles>
            <allele ref="true">A</allele>
            <allele>G</allele>
          </alleles>
        </genotype>
        <genotype available="true" called="true" het="false" hom="true" homRef="false" homVar="true" mixed="false" noCall="false" nonInformative="false" filtered="false" phased="false" sample="M10475">
          <AD>
            <value>0</value>
            <value>18</value>
          </AD>
          <DP>18</DP>
          <GQ>39</GQ>
          <PL>
            <value index="1">422</value>
            <value index="2">39</value>
            <value index="3">0</value>
          </PL>
          <alleles>
            <allele>G</allele>
            <allele>G</allele>
          </alleles>
        </genotype>
        <genotype available="true" called="true" het="true" hom="false" homRef="false" homVar="false" mixed="false" noCall="false" nonInformative="false" filtered="false" phased="false" sample="M10500">
          <AD>
            <value>5</value>
            <value>1</value>
          </AD>
          <DP>6</DP>
          <GQ>20</GQ>
          <PL>
            <value index="1">20</value>
            <value index="2">0</value>
            <value index="3">148</value>
          </PL>
          <alleles>
            <allele ref="true">A</allele>
            <allele>G</allele>
          </alleles>
        </genotype>
        <genotype available="true" called="true" het="true" hom="false" homRef="false" homVar="false" mixed="false" noCall="false" nonInformative="false" filtered="false" phased="false" sample="M10478">
          <AD>
            <value>8</value>
            <value>7</value>
          </AD>
          <DP>16</DP>
          <GQ>99</GQ>
          <PL>
            <value index="1">201</value>
            <value index="2">0</value>
            <value index="3">173</value>
          </PL>
          <alleles>
            <allele ref="true">A</allele>
            <allele>G</allele>
          </alleles>
        </genotype>
      </genotypes>
    </variation>
    <variation>
      <chrom>chr16</chrom>
      <start>72057435</start>
      <end>72057435</end>
      <ref>C</ref>
      <alt>T</alt>
      <qual>572.98</qual>
      <infos>
        <BaseQRankSum>-2.270</BaseQRankSum>
        <HaplotypeScore>4.5319</HaplotypeScore>
        <CSQ>
          <Consequence>missense_variant</Consequence>
          <Codons>Cgg/Tgg</Codons>
          <Amino_acids>R/W</Amino_acids>
          <Gene>ENSG00000102967</Gene>
          <HGNC>DHODH</HGNC>
          <Feature>ENST00000219240</Feature>
          <EXON>8/9</EXON>
          <PolyPhen>probably_damaging(0.984)</PolyPhen>
          <SIFT>deleterious(0)</SIFT>
        </CSQ>
        <CSQ>
          <Consequence>non_coding_exon_variant</Consequence>
          <Consequence>nc_transcript_variant</Consequence>
          <Gene>ENSG00000102967</Gene>
          <HGNC>DHODH</HGNC>
          <Feature>ENST00000571392</Feature>
          <EXON>3/4</EXON>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000102967</Gene>
          <HGNC>DHODH</HGNC>
          <Feature>ENST00000572003</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000102967</Gene>
          <HGNC>DHODH</HGNC>
          <Feature>ENST00000573843</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000102967</Gene>
          <HGNC>DHODH</HGNC>
          <Feature>ENST00000573922</Feature>
        </CSQ>
        <CSQ>
          <Consequence>intron_variant</Consequence>
          <Gene>ENSG00000102967</Gene>
          <HGNC>DHODH</HGNC>
          <Feature>ENST00000574309</Feature>
        </CSQ>
        <CSQ>
          <Consequence>missense_variant</Consequence>
          <Codons>Cgg/Tgg</Codons>
          <Amino_acids>R/W</Amino_acids>
          <Gene>ENSG00000102967</Gene>
          <HGNC>DHODH</HGNC>
          <Feature>ENST00000572887</Feature>
          <EXON>8/9</EXON>
          <PolyPhen>probably_damaging(0.926)</PolyPhen>
          <SIFT>deleterious(0)</SIFT>
        </CSQ>
        <QD>8.07</QD>
        <MQ>36.53</MQ>
        <AC>1</AC>
        <FS>0.000</FS>
        <HRun>0</HRun>
        <MQRankSum>0.596</MQRankSum>
        <ReadPosRankSum>0.927</ReadPosRankSum>
        <Dels>0.00</Dels>
        <DP>260</DP>
        <AF>0.13</AF>
        <MQ0>0</MQ0>
        <AN>8</AN>
      </infos>
      <genotypes>
        <genotype available="true" called="true" het="false" hom="true" homRef="true" homVar="false" mixed="false" noCall="false" nonInformative="false" filtered="false" phased="false" sample="M128215">
          <AD>
            <value>65</value>
            <value>0</value>
          </AD>
          <DP>65</DP>
          <GQ>99</GQ>
          <PL>
            <value index="1">0</value>
            <value index="2">180</value>
            <value index="3">1982</value>
          </PL>
          <alleles>
            <allele ref="true">C</allele>
            <allele ref="true">C</allele>
          </alleles>
        </genotype>
        <genotype available="true" called="true" het="true" hom="false" homRef="false" homVar="false" mixed="false" noCall="false" nonInformative="false" filtered="false" phased="false" sample="M10475">
          <AD>
            <value>37</value>
            <value>33</value>
          </AD>
          <DP>71</DP>
          <GQ>99</GQ>
          <PL>
            <value index="1">609</value>
            <value index="2">0</value>
            <value index="3">666</value>
          </PL>
          <alleles>
            <allele ref="true">C</allele>
            <allele>T</allele>
          </alleles>
        </genotype>
        <genotype available="true" called="true" het="false" hom="true" homRef="true" homVar="false" mixed="false" noCall="false" nonInformative="false" filtered="false" phased="false" sample="M10500">
          <AD>
            <value>67</value>
            <value>1</value>
          </AD>
          <DP>68</DP>
          <GQ>99</GQ>
          <PL>
            <value index="1">0</value>
            <value index="2">120</value>
            <value index="3">1460</value>
          </PL>
          <alleles>
            <allele ref="true">C</allele>
            <allele ref="true">C</allele>
          </alleles>
        </genotype>
        <genotype available="true" called="true" het="false" hom="true" homRef="true" homVar="false" mixed="false" noCall="false" nonInformative="false" filtered="false" phased="false" sample="M10478">
          <AD>
            <value>56</value>
            <value>0</value>
          </AD>
          <DP>56</DP>
          <GQ>99</GQ>
          <PL>
            <value index="1">0</value>
            <value index="2">135</value>
            <value index="3">1543</value>
          </PL>
          <alleles>
            <allele ref="true">C</allele>
            <allele ref="true">C</allele>
          </alleles>
        </genotype>
      </genotypes>
    </variation>
    <variation>
      <chrom>chr16</chrom>
      <start>72059269</start>
      <end>72059269</end>
      <ref>T</ref>
      <alt>C</alt>
      <qual>39.18</qual>
      <infos>
        <FS>0.000</FS>
        <AC>2</AC>
        <HaplotypeScore>0.0000</HaplotypeScore>
        <HRun>0</HRun>
        <Dels>0.00</Dels>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000102967</Gene>
          <HGNC>DHODH</HGNC>
          <Feature>ENST00000219240</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000102967</Gene>
          <HGNC>DHODH</HGNC>
          <Feature>ENST00000571392</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000102967</Gene>
          <HGNC>DHODH</HGNC>
          <Feature>ENST00000572003</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000102967</Gene>
          <HGNC>DHODH</HGNC>
          <Feature>ENST00000573843</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000102967</Gene>
          <HGNC>DHODH</HGNC>
          <Feature>ENST00000573922</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000102967</Gene>
          <HGNC>DHODH</HGNC>
          <Feature>ENST00000574309</Feature>
        </CSQ>
        <CSQ>
          <Consequence>downstream_gene_variant</Consequence>
          <Gene>ENSG00000102967</Gene>
          <HGNC>DHODH</HGNC>
          <Feature>ENST00000572887</Feature>
        </CSQ>
        <DP>2</DP>
        <QD>19.59</QD>
        <AF>1.00</AF>
        <MQ0>0</MQ0>
        <MQ>37.00</MQ>
        <AN>2</AN>
      </infos>
      <genotypes>
        <genotype available="true" called="false" het="false" hom="false" homRef="false" homVar="false" mixed="false" noCall="true" nonInformative="true" filtered="false" phased="false" sample="M128215">
          <alleles/>
        </genotype>
        <genotype available="true" called="false" het="false" hom="false" homRef="false" homVar="false" mixed="false" noCall="true" nonInformative="true" filtered="false" phased="false" sample="M10475">
          <alleles/>
        </genotype>
        <genotype available="true" called="true" het="false" hom="true" homRef="false" homVar="true" mixed="false" noCall="false" nonInformative="false" filtered="false" phased="false" sample="M10500">
          <AD>
            <value>0</value>
            <value>2</value>
          </AD>
          <DP>2</DP>
          <GQ>6</GQ>
          <PL>
            <value index="1">70</value>
            <value index="2">6</value>
            <value index="3">0</value>
          </PL>
          <alleles>
            <allele>C</allele>
            <allele>C</allele>
          </alleles>
        </genotype>
        <genotype available="true" called="false" het="false" hom="false" homRef="false" homVar="false" mixed="false" noCall="true" nonInformative="true" filtered="false" phased="false" sample="M10478">
          <alleles/>
        </genotype>
      </genotypes>
    </variation>
  </variations>
</vcf>
```
## Example 2 inserting into MongoDB:

Here is a simple XSLT stylesheet that transforms the above XML to MongoDB statements:

```XML
<?xml version='1.0'  encoding="ISO-8859-1"?>
<xsl:stylesheet xmlns:xsl='http://www.w3.org/1999/XSL/Transform'
	exclude-result-prefixes="x"
	xmlns:x="http://exslt.org/strings"
	extension-element-prefixes="x"
	version='1.0' > 
<xsl:output method="text"/>

<xsl:template match="/">
<xsl:apply-templates select="vcf"/>
</xsl:template>


<xsl:template match="vcf">
<xsl:text>use vcf
</xsl:text>
<xsl:apply-templates select="variations"/>
</xsl:template>


<xsl:template match="variations">
<xsl:apply-templates select="variation"/>
</xsl:template>

<xsl:template match="variation">
<xsl:text> db.variants.insert({</xsl:text>
<xsl:apply-templates select="chrom" mode="txt"/>
<xsl:text>,</xsl:text>
<xsl:apply-templates select="start" mode="num"/>
<xsl:text>,</xsl:text>
<xsl:apply-templates select="end" mode="num"/>

<xsl:if test="id">
<xsl:text>,</xsl:text>
<xsl:apply-templates select="id" mode="txt"/>
</xsl:if>

<xsl:if test="ref">
<xsl:text>,</xsl:text>
<xsl:apply-templates select="ref" mode="txt"/>
</xsl:if>

<xsl:text>,"alt":[</xsl:text>
<xsl:for-each select="alt">
<xsl:if test="position()&gt;1">,</xsl:if>
<xsl:text>"</xsl:text>
<xsl:value-of select="."/>
<xsl:text>"</xsl:text>
</xsl:for-each>
<xsl:text>]</xsl:text>

<xsl:if test="qual">
<xsl:text>,</xsl:text>
<xsl:apply-templates select="qual" mode="num"/>
</xsl:if>

<xsl:text>,"genotypes":[</xsl:text>
<xsl:for-each select="genotypes/genotype">
<xsl:if test="position()&gt;1">,</xsl:if>
<xsl:apply-templates select="."/>
</xsl:for-each>
<xsl:text>]</xsl:text>


<xsl:text>})
</xsl:text>
</xsl:template>


<xsl:template match="genotype">
<xsl:text>{</xsl:text>
<xsl:apply-templates select="@sample" mode="txt"/>


<xsl:if test="count(alleles/allele)&gt;0">
<xsl:text>,"alleles":[</xsl:text>
<xsl:for-each select="alleles/allele">
<xsl:if test="position()&gt;1">,</xsl:if>
<xsl:text>"</xsl:text>
<xsl:value-of select="."/>
<xsl:text>"</xsl:text>
</xsl:for-each>
<xsl:text>]</xsl:text>
</xsl:if>

<xsl:text>}</xsl:text>
</xsl:template>

<xsl:template match="*|@*" mode="txt">
<xsl:text>"</xsl:text>
<xsl:value-of select="local-name(.)"/>
<xsl:text>":"</xsl:text>
<xsl:value-of select="."/><!-- should escape characters here -->
<xsl:text>"</xsl:text>
</xsl:template>

<xsl:template match="*|@*" mode="num">
<xsl:text>"</xsl:text>
<xsl:value-of select="local-name(.)"/>
<xsl:text>":</xsl:text>
<xsl:value-of select="."/>
<xsl:text></xsl:text>
</xsl:template>
</xsl:stylesheet>
```
### generated script

```bash
$ xsltproc vcf2mongo.xsl vcf.xml 
use vcf
 db.variants.insert({"chrom":"chr1","start":10001,"end":10001,"ref":"T","alt":["TC"],"qual":175.91000000000003,"genotypes":[{"sample":"M128215","alleles":["T","TC"]},{"sample":"M10475","alleles":["T","TC"]},{"sample":"M10500","alleles":["T","TC"]},{"sample":"M10478","alleles":["T","TC"]}]})
 db.variants.insert({"chrom":"chr1","start":10056,"end":10056,"ref":"A","alt":["C"],"qual":47.27,"genotypes":[{"sample":"M128215","alleles":["A","A"]},{"sample":"M10475","alleles":["A","A"]},{"sample":"M10500","alleles":["A","A"]},{"sample":"M10478","alleles":["A","C"]}]})
 db.variants.insert({"chrom":"chr1","start":10121,"end":10121,"ref":"A","alt":["C"],"qual":111.08000000000001,"genotypes":[{"sample":"M128215","alleles":["A","C"]},{"sample":"M10475","alleles":["A","A"]},{"sample":"M10500","alleles":["A","A"]},{"sample":"M10478","alleles":["A","A"]}]})
 db.variants.insert({"chrom":"chr1","start":10177,"end":10177,"ref":"A","alt":["C"],"qual":163.46,"genotypes":[{"sample":"M128215","alleles":["A","C"]},{"sample":"M10475","alleles":["A","C"]},{"sample":"M10500","alleles":["A","A"]},{"sample":"M10478","alleles":["A","C"]}]})
 db.variants.insert({"chrom":"chr1","start":10180,"end":10180,"ref":"T","alt":["C"],"qual":79.53,"genotypes":[{"sample":"M128215","alleles":["T","T"]},{"sample":"M10475","alleles":["T","T"]},{"sample":"M10500","alleles":["T","C"]},{"sample":"M10478","alleles":["T","C"]}]})
 db.variants.insert({"chrom":"chr1","start":10234,"end":10234,"ref":"C","alt":["T"],"qual":49.1,"genotypes":[{"sample":"M128215","alleles":["C","C"]},{"sample":"M10475","alleles":["C","T"]},{"sample":"M10500","alleles":["C","C"]},{"sample":"M10478","alleles":["C","T"]}]})
 db.variants.insert({"chrom":"chr16","start":72057282,"end":72057282,"ref":"A","alt":["G"],"qual":1098.89,"genotypes":[{"sample":"M128215","alleles":["A","G"]},{"sample":"M10475","alleles":["G","G"]},{"sample":"M10500","alleles":["A","G"]},{"sample":"M10478","alleles":["A","G"]}]})
 db.variants.insert({"chrom":"chr16","start":72057435,"end":72057435,"ref":"C","alt":["T"],"qual":572.98,"genotypes":[{"sample":"M128215","alleles":["C","C"]},{"sample":"M10475","alleles":["C","T"]},{"sample":"M10500","alleles":["C","C"]},{"sample":"M10478","alleles":["C","C"]}]})
 db.variants.insert({"chrom":"chr16","start":72059269,"end":72059269,"ref":"T","alt":["C"],"qual":39.18,"genotypes":[{"sample":"M128215"},{"sample":"M10475"},{"sample":"M10500","alleles":["C","C"]},{"sample":"M10478"}]})
```

### insert into **mongodb**
```bash
$ xsltproc vcf2mongo.xsl vcf.xml | mongo
```
### select data from mongodb
```
> use vcf;
## select variants on chromosome 1 from 10057 to 10234
> db.variants.find({"chrom":"chr1","start" :{ $gt: 10057, $lte: 10234, } })
{ "_id" : ObjectId("5267e19a7bc3eca84c83784d"), "chrom" : "chr1", "start" : 10121, "end" : 10121, "ref" : "A", "alt" : [  "C" ], "qual" : 111.08000000000001, "genotypes" : [ 	{ 	"sample" : "M128215", 	"alleles" : [ 	"A", 	"C" ] }, 	{ 	"sample" : "M10475", 	"alleles" : [ 	"A", 	"A" ] }, 	{ 	"sample" : "M10500", 	"alleles" : [ 	"A", 	"A" ] }, 	{ 	"sample" : "M10478", 	"alleles" : [ 	"A", 	"A" ] } ] }
{ "_id" : ObjectId("5267e19a7bc3eca84c83784e"), "chrom" : "chr1", "start" : 10177, "end" : 10177, "ref" : "A", "alt" : [  "C" ], "qual" : 163.46, "genotypes" : [ 	{ 	"sample" : "M128215", 	"alleles" : [ 	"A", 	"C" ] }, 	{ 	"sample" : "M10475", "alleles" : [ 	"A", 	"C" ] }, 	{ 	"sample" : "M10500", 	"alleles" : [ 	"A", 	"A" ] }, 	{ 	"sample" : "M10478", 	"alleles" : [ 	"A", 	"C" ] } ] }
{ "_id" : ObjectId("5267e19a7bc3eca84c83784f"), "chrom" : "chr1", "start" : 10180, "end" : 10180, "ref" : "T", "alt" : [  "C" ], "qual" : 79.53, "genotypes" : [ 	{ 	"sample" : "M128215", 	"alleles" : [ 	"T", 	"T" ] }, 	{ 	"sample" : "M10475", "alleles" : [ 	"T", 	"T" ] }, 	{ 	"sample" : "M10500", 	"alleles" : [ 	"T", 	"C" ] }, 	{ 	"sample" : "M10478", 	"alleles" : [ 	"T", 	"C" ] } ] }
{ "_id" : ObjectId("5267e19a7bc3eca84c837850"), "chrom" : "chr1", "start" : 10234, "end" : 10234, "ref" : "C", "alt" : [  "T" ], "qual" : 49.1, "genotypes" : [ 	{ 	"sample" : "M128215", 	"alleles" : [ 	"C", 	"C" ] }, 	{ 	"sample" : "M10475", "alleles" : [ 	"C", 	"T" ] }, 	{ 	"sample" : "M10500", 	"alleles" : [ 	"C", 	"C" ] }, 	{ 	"sample" : "M10478", 	"alleles" : [ 	"C", 	"T" ] } ] }

### select chrom,start,end for  variants on chromosome 1 from 10057 to 10234
> db.variants.find({"chrom":"chr1","start" :{ $gt: 10057, $lte: 10234, } },{"chrom":1,"start":1,"end":1})
{ "_id" : ObjectId("5267e19a7bc3eca84c83784d"), "chrom" : "chr1", "start" : 10121, "end" : 10121 }
{ "_id" : ObjectId("5267e19a7bc3eca84c83784e"), "chrom" : "chr1", "start" : 10177, "end" : 10177 }
{ "_id" : ObjectId("5267e19a7bc3eca84c83784f"), "chrom" : "chr1", "start" : 10180, "end" : 10180 }
{ "_id" : ObjectId("5267e19a7bc3eca84c837850"), "chrom" : "chr1", "start" : 10234, "end" : 10234 }
```

