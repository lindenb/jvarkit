/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


*/
package com.github.lindenb.jvarkit.tools.vcf2xml;

import java.io.OutputStream;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.dict.DictionaryXmlSerializer;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.net.Hyperlink;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParserFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParserFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;
/**
BEGIN_DOC


## Example

```bash
$  curl "https://raw.github.com/arq5x/gemini/master/test/test.region.vep.vcf" |\
   java -jar dist/jvarkit.jar vcf2xml   |\
   xmllint --format -
```

### Result

This is an old exampe: the format/schema may have changed since I created the tool

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

END_DOC

 */
@Program(
	name="vcf2xml",
	description="Convert VCF to XML",
	keywords={"vcf","xml"},
	modificationDate = "20251124",
	jvarkit_amalgamion = true
	)
public class Vcf2Xml extends Launcher
	{
	private static final Logger LOG=Logger.of(Vcf2Xml.class);

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"--hide"},description="features to hide '(genotypes|gt|het|nocall|homref|mixed|homvar|info|filter|id|header)'")
	private String hide_str = "";
	@Parameter(names={"--contig-regex"},description="keep chromosomes matching that regular expression")
	private String contig_regex=null;
	@Parameter(names={"--min-contig-length"},description="keep chromosomes which length is greater than 'x'")
	private int min_contig_length=0;
	@Parameter(names={"--omit-xml-desclaration"},description="Don't print XM declaration.")
	private boolean omit_xml_decl = false;

		
		
	private static void element(final XMLStreamWriter w,String name,Object o) throws XMLStreamException {
		if(o==null) return;
		final String s=String.valueOf(o);
		if(StringUtils.isBlank(s)) return;
		w.writeStartElement(name);
		w.writeCharacters(s);
		w.writeEndElement();
		}
	   
		
	@SuppressWarnings("rawtypes")
	private static void writeAttributes(final XMLStreamWriter w, final Map<String, Object> hash) throws XMLStreamException {
		if (hash == null)
			return;
		for (String key : hash.keySet()) {
			boolean is_array = false;
			final Object o = hash.get(key);
			final Collection col;
			if (o == null)
				continue;
			if (o.getClass().isArray()) {
				final Object array[] = (Object[]) o;
				col = Arrays.asList(array);
				is_array = true;
			} else if (o instanceof Collection) {
				final Collection array = (Collection) o;
				col = array;
				is_array = true;
			} else {
				col = Collections.singletonList(o);
			}
			if (col.isEmpty())
				continue;
			w.writeStartElement("attribute");
			w.writeAttribute("id", key);
			int idx=1;
			for (Object v: col) {
				w.writeStartElement("value");
				if (is_array) {
					w.writeAttribute("index", String.valueOf(idx++));
					}
				w.writeCharacters(String.valueOf(v));
				w.writeEndElement();
				}
			w.writeEndElement();

		}
	}
		

	

	@Override
	public int doWork(final List<String> args) {
		try {
			final String input = super.oneFileOrNull(args);
			
			final Set<String> hide= Arrays.stream(this.hide_str.split("[  ,;|]"))
					.map(S->S.toLowerCase())
					.filter(S->!StringUtils.isBlank(S))
					.collect(Collectors.toSet());
			
			if(hide.contains("genotype")) hide.add("gt");
			if(hide.contains("genotypes")) hide.add("gt");
			if(hide.contains("infos")) hide.add("info");
			if(hide.contains("href")) hide.add("url");
			if(hide.contains("filters")) hide.add("filter");
			if(hide.contains("no_call")) hide.add("nocall");
			if(hide.contains("hom_ref")) hide.add("homref");
			if(hide.contains("hom_var")) hide.add("homvar");
			if(hide.contains("ann")) hide.add("snpeff");
			
			try(VCFIterator iter = input==null || input.equals("-")?new VCFIteratorBuilder().open(stdin()):new VCFIteratorBuilder().open(input)) {
				final VCFHeader header = iter.getHeader();
				final SAMSequenceDictionary dict0 = header.getSequenceDictionary();
				final Hyperlink hyperlink = dict0==null?Hyperlink.empty():Hyperlink.compile(dict0);
				
				final Pattern contigRegex = (this.contig_regex==null?null:Pattern.compile(this.contig_regex));
				final SAMSequenceDictionary dict = dict0==null?null:new SAMSequenceDictionary(
						dict0.getSequences().stream().
						filter(SSR->min_contig_length<=0 || SSR.getSequenceLength()>=min_contig_length).
						filter(SSR->contigRegex==null?true:contigRegex.matcher(SSR.getSequenceName()).matches()).
						collect(Collectors.toList())
						);
				if(dict==null || dict.isEmpty()) {
					LOG.warn("empty dictionary");
					}
				
				final Map<String,Long> contig2genomicindex = new HashMap<>(dict==null?10:dict.size());
				final long genome_length = dict==null?0L:dict.getReferenceLength();
				final VepPredictionParser vepParser = new VepPredictionParserFactory(header).get();
				final AnnPredictionParser annParser = new AnnPredictionParserFactory(header).get();
				
				try(OutputStream out= super.openPathOrStdoutAsStream(this.outputFile)) {
					final XMLOutputFactory xmlfactory= XMLOutputFactory.newInstance();
					final String encoding = "UTF-8";
					final XMLStreamWriter w= xmlfactory.createXMLStreamWriter(out, encoding);
					if(!omit_xml_decl) w.writeStartDocument(encoding, "1.0");
					w.writeStartElement("vcf");
					if(input!=null) w.writeAttribute("src", input);
					
					/* VCF HEADER *****************************************/
					if(!hide.contains("header")) {
						w.writeStartElement("header");
						if(header.getInfoHeaderLines()!=null)
							{
							w.writeStartElement("infos");
							for(final VCFInfoHeaderLine h:header.getInfoHeaderLines())
								{
								w.writeStartElement("info");
								w.writeAttribute("id",h.getID());
								w.writeAttribute("countType", String.valueOf( h.getCountType()));
								if(h.getCountType()==VCFHeaderLineCount.INTEGER)
									{
									w.writeAttribute("count", String.valueOf( h.getCount()));
									}
								if(!StringUtils.isBlank(h.getValue())) {
									w.writeAttribute("value",h.getValue());
									}
								w.writeCharacters(h.getDescription());
								w.writeEndElement();
								}
							w.writeEndElement();
							}
						
						if(header.getFormatHeaderLines()!=null && !hide.contains("genotype"))
							{
							w.writeStartElement("formats");
							for(final VCFFormatHeaderLine h: header.getFormatHeaderLines())
								{
								w.writeStartElement("format");
								w.writeAttribute("id",h.getID());
								w.writeAttribute("countType",h.getCountType().name());
								if(h.getCountType()==VCFHeaderLineCount.INTEGER)
									{
									w.writeAttribute("count",String.valueOf(h.getCount()));
									}
								if(!StringUtils.isBlank(h.getValue())) {
									w.writeAttribute("value",h.getValue());
									}
								w.writeCharacters(h.getDescription());
								w.writeEndElement();
								}
							w.writeEndElement();
							}
						
						if(header.getFilterLines()!=null && !hide.contains("filter"))
							{
							w.writeStartElement("filters");
							for(VCFFilterHeaderLine h: header.getFilterLines())
								{
								w.writeStartElement("filter");
								w.writeAttribute("id",h.getID());
								if(!StringUtils.isBlank(h.getDescription())) {
									w.writeCharacters(h.getDescription());
									}
								w.writeEndElement();
								}
							w.writeEndElement();
							}
						
						//write dict
						if(dict!=null)
							{
							new DictionaryXmlSerializer().writeDictionary(w, dict);
							long x=0;
							for(SAMSequenceRecord ssr: dict.getSequences()) {
								contig2genomicindex.put(ssr.getContig(), x);
								x+= ssr.getLengthOnReference();
								}
							
							} // end dict
						
						if(header.getSampleNamesInOrder()!=null && !hide.contains("gt"))
							{
							w.writeStartElement("samples");
							for(final String name:header.getSampleNamesInOrder())
								{
								w.writeStartElement("sample");
								w.writeAttribute("index", String.valueOf(header.getSampleNameToOffset().get(name)));
								w.writeCharacters(name);
								w.writeEndElement();
								}
							w.writeEndElement();
							}
						if(header.getMetaDataInInputOrder()!=null)
							{
							w.writeStartElement("metas");
							for(final VCFHeaderLine meta:header.getMetaDataInInputOrder())
								{
								if(meta.getKey().equals( "INFO"))continue;
								if(meta.getKey().equals("FORMAT"))continue;
								if(meta.getKey().equals("contig"))continue;
								if(meta.getKey().equals("FILTER"))continue;
								if(meta.getKey().equals("fileformat"))continue;
								w.writeStartElement("meta");
								w.writeAttribute("key", meta.getKey());
								if(!StringUtils.isBlank(meta.getValue())) {
									w.writeCharacters(meta.getValue());
									}
								w.writeEndElement();
								}
							w.writeEndElement();
							}
						w.writeEndElement();//header
						}
					/* variants */
					w.writeStartElement("variants");
					while(iter.hasNext()) {
						final VariantContext ctx = iter.next();
						final SAMSequenceRecord ssr;
						if(dict!=null) {
							ssr = dict.getSequence(ctx.getContig());
							if(ssr==null) continue;
							}
						else
							{
							if(contigRegex!=null && !contigRegex.matcher(ctx.getContig()).matches()) continue; 
							ssr=null;
							}
						w.writeStartElement("variant");
						if(dict!=null && contig2genomicindex!=null) {
							final long x1 =  contig2genomicindex.get(ctx.getContig()).longValue() + ctx.getStart();
							final long x2 =  contig2genomicindex.get(ctx.getContig()).longValue() + ctx.getEnd();
							w.writeAttribute("x1", String.valueOf(x1));
							w.writeAttribute("x2", String.valueOf(x2));
							w.writeAttribute("f1", String.valueOf(x1/(double)genome_length));
							w.writeAttribute("f2", String.valueOf(x2/(double)genome_length));
							}
						int len = ctx.getLengthOnReference();
						if(ctx.hasAttribute("SVLEN") ) {
							len = Math.max(len,ctx.getAttributeAsIntList("SVLEN",0).stream().mapToInt(n->Math.abs(n)).max().orElse(0));
							}
						w.writeAttribute("length", String.valueOf(len));
						
						
						if(!hide.contains("url")) {
							final Optional<String> url= hyperlink.apply(ctx);
							if(url.isPresent()) {
								w.writeEmptyElement("url");
								w.writeAttribute("href", url.get());
								}
							else if(ctx.getID()!=null && ctx.getID().matches("rs[0-9]+"))
								{
								w.writeEmptyElement("url");
								w.writeAttribute("href",  "https://www.ncbi.nlm.nih.gov/snp/"+ctx.getID());
								}
							}
						
						
						
						w.writeStartElement("chrom");
						if(ssr!=null) {
							w.writeAttribute("tid", String.valueOf(ssr.getSequenceIndex()));
							}
						w.writeCharacters(ctx.getContig());
						w.writeEndElement();
						
						w.writeStartElement("start");
						w.writeCharacters(String.valueOf(ctx.getStart()));
						w.writeEndElement();
						
						w.writeStartElement("end");
						w.writeCharacters(String.valueOf(ctx.getEnd()));
						w.writeEndElement();
						

						
						if(ctx.hasID() && !hide.contains("id"))
							{
							for(String id: CharSplitter.SEMICOLON.split(ctx.getID())) {
								if(id.equals(".") || StringUtils.isBlank(id)) continue;
								w.writeStartElement("id");
								w.writeCharacters(id);
								w.writeEndElement();
								}
							}
						
						
						
						w.writeStartElement("alleles");
						w.writeAttribute("count", String.valueOf(ctx.getNAlleles()));
						for(Allele a:ctx.getAlleles())
							{
							w.writeStartElement("allele");
							if(a.isReference()) {
								w.writeAttribute("ref", "true");
								}
							w.writeAttribute("idx", String.valueOf(ctx.getAlleleIndex(a)));
							w.writeCharacters(a.getDisplayString());
							w.writeEndElement();
							}
						w.writeEndElement();
						
						if(ctx.hasLog10PError() && !hide.contains("qual"))
							{
							w.writeStartElement("qual");
							w.writeCharacters(String.valueOf(ctx.getPhredScaledQual()));
							w.writeEndElement();
							}
						
						if((ctx.isFiltered() || ctx.filtersWereApplied()) && !hide.contains("filter"))
							{
							w.writeStartElement("filters");
							if(ctx.isFiltered())
								{
								for(String s: ctx.getFilters())
									{
									w.writeStartElement("filter");
									w.writeCharacters(s);
									w.writeEndElement();
									}
								}
							else if(ctx.filtersWereApplied())
								{
								w.writeStartElement("filter");
								w.writeCharacters(VCFConstants.PASSES_FILTERS_v4);
								w.writeEndElement();
								}
							w.writeEndElement();
							}
						/* INFO COLUMN********************************************* */
						if(ctx.getAttributes()!=null  && !hide.contains("info"))
							{
							w.writeStartElement("infos");
							writeAttributes(w,ctx.getAttributes());

							// VEP
							if(!hide.contains("vep") && vepParser.isValid() && ctx.hasAttribute(vepParser.getTag())) {
								final List<VepPredictionParser.VepPrediction> preds= vepParser.getPredictions(ctx);
								if(!preds.isEmpty()) {
									for(VepPredictionParser.VepPrediction pred:preds) {
										w.writeStartElement("vep");
										for(final String cat:vepParser.getCategories()) {
											final String v = pred.getByCol(cat);
											if(StringUtils.isBlank(v)) continue;
											List<String> cols = Collections.singletonList(v);
											if(cat.equals("Consequence")) cols= CharSplitter.AMP.splitAsStringList(v);
											for(String s: cols) {
												if(StringUtils.isBlank(s)) continue;
												w.writeStartElement("property");
												w.writeAttribute("name", cat);
												w.writeCharacters(s);
												w.writeEndElement();
												}
											}
										w.writeEndElement();
										}
									}
								}
							// SNPEFF
							if(!hide.contains("snpeff") && annParser.isValid() && ctx.hasAttribute(annParser.getTag())) {
								final List<AnnPredictionParser.AnnPrediction> preds= annParser.getPredictions(ctx);
								if(!preds.isEmpty()) {
									for(AnnPredictionParser.AnnPrediction pred:preds) {
										w.writeStartElement("ann");
										for(int x=0;x< AnnPredictionParser.getColumnCount();++x) {
											final String v = pred.at(x);
											if(StringUtils.isBlank(v)) continue;
											List<String> cols = Collections.singletonList(v);
											if( AnnPredictionParser.getLabel(x).equals("Annotation")) cols= CharSplitter.AMP.splitAsStringList(v);
											for(String s: cols) {
												if(StringUtils.isBlank(s)) continue;
												w.writeStartElement("property");
												w.writeAttribute("name", AnnPredictionParser.getLabel(x));
												w.writeCharacters(s);
												w.writeEndElement();
												}
											}
										w.writeEndElement();
										}
									}
								}
							w.writeEndElement();
							} // INFO
						
						
						
						if(ctx.hasGenotypes() && !hide.contains("gt"))
							{
							w.writeStartElement("genotypes");
							for(String sample:ctx.getSampleNames())
								{
								final Genotype g=ctx.getGenotype(sample);
								if(g==null) continue;
								if(hide.contains("het") && g.isHet()) continue;
								if(hide.contains("homvar") && g.isHomVar()) continue;
								if(hide.contains("homref") && g.isHomRef()) continue;
								if(hide.contains("nocall") && g.isNoCall()) continue;
								if(hide.contains("mixed") && g.isMixed()) continue;
								
								
								w.writeStartElement("genotype");
								w.writeAttribute("type",g.getType().name());
								if(g.isFiltered()) w.writeAttribute("filtered",String.valueOf(g.isFiltered()));
								if(g.isPhased()) w.writeAttribute("phased",String.valueOf(g.isPhased()));
								w.writeAttribute("sample",g.getSampleName());
								if(g.hasAD() && !hide.contains("ad"))
									{
									w.writeStartElement("AD");
									for(int ad:g.getAD())
										{
										element(w,"value", ad);
										}
									w.writeEndElement();
									}
								if(g.hasDP())
									{
									}
								if(g.hasGQ() && !hide.contains("gq"))
									{
									w.writeStartElement("GQ");
									w.writeCharacters(String.valueOf(g.getGQ()));
									w.writeEndElement();
									}
								if(g.hasPL() && !hide.contains("pl"))
									{
									w.writeStartElement("PL");
									int index=0;
									for(int v:g.getPL())
										{
										w.writeStartElement("value");
										w.writeAttribute("index",String.valueOf(index++));
										w.writeCharacters(String.valueOf(v));
										w.writeEndElement();
										}
									w.writeEndElement();
									}
								
								
								
								w.writeStartElement("alleles");
								for(Allele a:g.getAlleles())
									{
									w.writeStartElement("allele");
									if(ctx.getAlleleIndex(a)>=0) w.writeAttribute("index", String.valueOf(ctx.getAlleleIndex(a)));
									if(a.isReference()) w.writeAttribute("ref", String.valueOf(a.isReference()));
									if(a.isSymbolic()) w.writeAttribute("symbolic","true");
									w.writeCharacters(a.getDisplayString());
									w.writeEndElement();
									}
								w.writeEndElement();
								
								
								writeAttributes(w,g.getExtendedAttributes());
								w.writeEndElement();//genotype
								}
							w.writeEndElement();// genotypes
							}
						
						w.writeEndElement();//variant
						}
					
					w.writeEndElement();// variations
					w.writeEndElement();//vcf
					w.writeEndDocument();
					w.flush();
					w.close();
					out.flush();
					}
				}
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	public static void main(final String[] args)
		{
		new Vcf2Xml().instanceMainWithExit(args);
		}
}
