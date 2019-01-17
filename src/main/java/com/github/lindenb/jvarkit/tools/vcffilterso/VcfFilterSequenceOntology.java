/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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


History:
* 2014 creation
* 2015 knime interface

*/
package com.github.lindenb.jvarkit.tools.vcffilterso;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;
import javax.xml.bind.annotation.XmlType;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;
import com.github.lindenb.jvarkit.util.vcf.DelegateVariantContextWriter;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VariantAttributesRecalculator;
import com.github.lindenb.jvarkit.util.vcf.VariantContextWriterFactory;
import htsjdk.variant.vcf.VCFIterator;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParserFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.MyPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParserFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParserFactory;


/**
 * VCF filter on Sequence Ontology
BEGIN_DOC

## Examples

### Example 

list the variants having a "*feature_elongation*"  ( SO:0001907  ) http://www.sequenceontology.org/browser/current_release/term/SO:0001907

here the variant is selected because *stop_lost* is a children of *SO:0001907* .

```
$ curl -s -k "https://raw.github.com/arq5x/gemini/master/test/test5.vep.snpeff.vcf" | \
java -jar dist/vcffilterso.jar  -A "SO:0001907"  |\
grep -v "##"

#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	M10475	M10478	M10500M128215
chr10	1142208	.	T	C	3404.30	.	AC=8;AF=1.00;AN=8;CSQ=intron_variant|||ENSG00000047056|WDR37|ENST00000263150|||,downstream_gene_variant|||ENSG00000047056|WDR37|ENST00000436154|||,intron_variant|||ENSG00000047056|WDR37|ENST00000358220|||,stop_lost|Tga/Cga|* / R|ENSG00000047056|WDR37|ENST00000381329|9/9||;DP=122;Dels=0.00;EFF=DOWNSTREAM(MODIFIER||||208|WDR37|protein_coding|CODING|ENST00000436154|),INTRON(MODIFIER||||494|WDR37|protein_coding|CODING|ENST00000263150|9),INTRON(MODIFIER||||494|WDR37|protein_coding|CODING|ENST00000358220|9),STOP_LOST(HIGH|MISSENSE|Tga/Cga|*250R|249|WDR37|protein_coding|CODING|ENST00000381329|);FS=0.000;HRun=0;HaplotypeScore=2.6747;MQ=36.00;MQ0=0;QD=27.90	GT:AD:DP:GQ:PL	1/1:1,37:39:87.16:940,87,0	1/1:0,29:29:78.20:899,78,0	1/1:0,24:24:66.14:729,66,0	1/1:0,30:30:75.18:836,75,0
```

### Example 

invert the query:

 ```
$ curl -s -k "https://raw.github.com/arq5x/gemini/master/test/test5.vep.snpeff.vcf" | \
java -jar dist/vcffilterso.jar -v -A "SO:0001907"  |\
grep -v "##"

#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	M10475	M10478	M10500M128215
chr1	145273345	.	T	C	289.85	.	AC=3;AF=0.38;AN=8;BaseQRankSum=1.062;CSQ=missense_variant|Tct/Cct|S/P|ENSG00000213240|NOTCH2NL|ENST00000369340|4/6|benign(0.238)|tolerated(0.45),missense_variant|Tct/Cct|S/P|ENSG00000213240|NOTCH2NL|ENST00000362074|3/5|benign(0.238)|tolerated(0.45),missense_variant&amp;NMD_transcript_variant|Tct/Cct|S/P|ENSG00000255168||ENST00000468030|3/23|benign(0.416)|tolerated(0.55),missense_variant|Tct/Cct|S/P|ENSG00000213240|NOTCH2NL|ENST00000344859|3/6|possibly_damaging(0.545)|tolerated(0.44);DP=1000;DS;Dels=0.00;EFF=EXON(MODIFIER|||||RP11-458D21.5|nonsense_mediated_decay|NON_CODING|ENST00000468030|3),NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|Tct/Cct|S67P|230|NOTCH2NL|protein_coding|CODING|ENST00000344859|),NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|Tct/Cct|S67P|236|NOTCH2NL|protein_coding|CODING|ENST00000362074|),NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|Tct/Cct|S67P|236|NOTCH2NL|protein_coding|CODING|ENST00000369340|);FS=3.974;HRun=1;HaplotypeScore=17.4275;MQ=29.25;MQ0=0;MQRankSum=-1.370;QD=0.39;ReadPosRankSum=-1.117	GT:AD:DP:GQ:PL	0/0:226,22:250:99:0,158,4259	0/1:224,24:250:5.77:6,0,5314	0/1:219,28:249:57.30:57,0,5027	0/1:215,34:250:99:269,0,3796
chr1	156011444	.	T	C	2523.46	.	AC=4;AF=0.50;AN=8;BaseQRankSum=-0.490;CSQ=missense_variant|atA/atG|I/M|ENSG00000160803|UBQLN4|ENST00000368309|10/11|benign(0.012)|tolerated(0.3),downstream_gene_variant|||ENSG00000160803|UBQLN4|ENST00000459954|||,missense_variant|Atc/Gtc|I/V|ENSG00000160803|UBQLN4|ENST00000368307|6/7|unknown(0)|tolerated(0.88);DP=204;Dels=0.00;EFF=DOWNSTREAM(MODIFIER|||||UBQLN4|processed_transcript|CODING|ENST00000459954|),NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|Atc/Gtc|I148V|226|UBQLN4|protein_coding|CODING|ENST00000368307|),NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|atA/atG|I495M|601|UBQLN4|protein_coding|CODING|ENST00000368309|);FS=4.328;HRun=0;HaplotypeScore=4.3777;MQ=35.24;MQ0=0;MQRankSum=-0.101;QD=14.93;ReadPosRankSum=1.575	GT:AD:DP:GQ:PL	0/1:24,15:40:99:214,0,443	0/1:32,36:68:99:702,0,794	1/1:1,59:61:99:1656,132,0	0/0:34,1:35:69.10:0,69,717
chr5	64982321	.	T	C	61.12	.	AC=4;AF=1.00;AN=4;CSQ=intron_variant|||ENSG00000197860|SGTB|ENST00000381007|||,intron_variant|||ENSG00000197860|SGTB|ENST00000506816|||;DP=4;Dels=0.00;EFF=INTRON(MODIFIER||||194|SGTB|protein_coding|CODING|ENST00000506816|5),INTRON(MODIFIER||||304|SGTB|protein_coding|CODING|ENST00000381007|5);FS=0.000;HRun=0;HaplotypeScore=0.0000;MQ=37.00;MQ0=0;QD=20.37	GT:AD:DP:GQ:PL	1/1:0,2:2:6:58,6,0	1/1:0,1:1:3.01:37,3,0	./.	./.
chr10	126678092	.	G	A	89.08	.	AC=1;AF=0.13;AN=8;BaseQRankSum=-3.120;CSQ=stop_gained|Caa/Taa|Q/*|ENSG00000175029|CTBP2|ENST00000531469|11/11||,stop_gained|Caa/Taa|Q/*|ENSG00000175029|CTBP2|ENST00000309035|9/9||,downstream_gene_variant|||ENSG00000019995|ZRANB1|ENST00000359653|||,stop_gained|Caa/Taa|Q/*|ENSG00000175029|CTBP2|ENST00000494626|11/11||,stop_gained|Caa/Taa|Q/*|ENSG00000175029|CTBP2|ENST00000337195|11/11||,stop_gained|Caa/Taa|Q/*|ENSG00000175029|CTBP2|ENST00000334808|9/9||,stop_gained|Caa/Taa|Q/*|ENSG00000175029|CTBP2|ENST00000411419|11/11||,downstream_gene_variant|||ENSG00000175029|CTBP2|ENST00000395705|||;DP=185;Dels=0.00;EFF=DOWNSTREAM(MODIFIER||||708|ZRANB1|protein_coding|CODING|ENST00000359653|),DOWNSTREAM(MODIFIER|||||CTBP2|processed_transcript|CODING|ENST00000395705|),STOP_GAINED(HIGH|NONSENSE|Caa/Taa|Q445*|445|CTBP2|protein_coding|CODING|ENST00000337195|),STOP_GAINED(HIGH|NONSENSE|Caa/Taa|Q445*|445|CTBP2|protein_coding|CODING|ENST00000411419|),STOP_GAINED(HIGH|NONSENSE|Caa/Taa|Q445*|445|CTBP2|protein_coding|CODING|ENST00000494626|),STOP_GAINED(HIGH|NONSENSE|Caa/Taa|Q445*|445|CTBP2|protein_coding|CODING|ENST00000531469|),STOP_GAINED(HIGH|NONSENSE|Caa/Taa|Q513*|513|CTBP2|protein_coding|CODING|ENST00000334808|),STOP_GAINED(HIGH|NONSENSE|Caa/Taa|Q985*|985|CTBP2|protein_coding|CODING|ENST00000309035|);FS=3.490;HRun=0;HaplotypeScore=3.3843;MQ=25.32;MQ0=0;MQRankSum=6.568;QD=2.02;ReadPosRankSum=-5.871	GT:AD:DP:GQ:PL	0/0:64,3:67:99:0,165,1505	0/0:11,1:12:7.31:0,7,240	0/0:52,10:62:54.97:0,55,1263	0/1:35,9:44:99:125,0,693
chr10	135210791	.	T	C	65.41	.	AC=4;AF=0.50;AN=8;BaseQRankSum=2.054;CSQ=intron_variant|||ENSG00000148824||ENST00000317502|||,upstream_gene_variant|||ENSG00000148824||ENST00000492266|||,intron_variant&amp;nc_transcript_variant|||ENSG00000148824||ENST00000460848|||,intron_variant|||ENSG00000254536|MTG1|ENST00000468317|||,intron_variant&amp;nc_transcript_variant|||ENSG00000148824||ENST00000477902|||,intron_&amp;nc_transcript_variant|||ENSG00000148824||ENST00000473735|||,intron_variant&amp;nc_transcript_variant|||ENSG00000148824||ENST00000498790|||,intron_variant|||ENSG00000148824||ENST00000432508|||,intron_variant&amp;nc_transcript_variant|||ENSG00000148824||ENST00000498334|||,intron_variant&amp;nc_transcript_variant|||ENSG00000148824||ENST00000495014|||;DP=11;Dels=0.00;EFF=INTRON(MODIFIER||||283|MTG1|protein_coding|CODING|ENST00000432508|3),INTRON(MODIFIER||||334|MTG1|protein_coding|CODING|ENST00000317502|3),INTRON(MODIFIER||||339|MTG1|protein_coding|CODING|ENST00000468317|4),INTRON(MODIFIER|||||MTG1|processed_transcript|CODING|ENST00000460848|3),INTRON(MODIFIER|||||MTG1|processed_transcript|CODING|ENST00000473735|3),INTRON(MODIFIER|||||MTG1|processed_transcript|CODING|ENST00000477902|3),INTRON(MODIFIER|||||MTG1|processed_transcript|CODING|ENST00000495014|3),INTRON(MODIFIER|||||MTG1|processed_transcript|CODING|ENST00000498334|3),INTRON(MODIFIER|||||MTG1|processed_transcript|CODING|ENST00000498790|3),UPSTREAM(MODIFIER|||||MTG1|processed_transcript|CODING|ENST00000492266|);FS=0.000;HRun=0;HaplotypeScore=0.2489;MQ=35.12;MQ0=0;MQRankSum=0.248;QD=16.35;ReadPosRankSum=-1.001	GT:AD:DP:GQ:PL	0/0:4,0:4:9:0,9,84	1/1:0,3:3:6.02:74,6,0	1/1:0,1:1:3.01:37,3,0	0/0:3,0:3:9.02:0,9,100
chr13	48873835	.	G	A	58.95	.	AC=4;AF=1.00;AN=4;CSQ=upstream_gene_variant|||ENSG00000139687|RB1|ENST00000467505|||,upstream_gene_variant|||ENSG00000139687|RB1|ENST00000525036|||,downstream_gene_variant|||ENSG00000231473|LINC00441|ENST00000436963|||,upstream_gene_variant|||ENSG00000139687|RB1|ENST00000267163|||,intron_variant&amp;nc_transcript_variant|||ENSG00000231473|LINC00441|ENST00000433480|||;DP=3;Dels=0.00;EFF=DOWNSTREAM(MODIFIER|||||LINC00441|processed_transcript|NON_CODING|ENST00000436963|),INTRON(MODIFIER|||||LINC00441|processed_transcript|NON_CODING|ENST00000433480|2),UPSTREAM(MODIFIER||||928|RB1|protein_coding|CODING|ENST00000267163|),UPSTREAM(MODIFIER|||||RB1|nonsense_mediated_decay|CODING|ENST00000467505|),UPSTREAM(MODIFIER|||||RB1|retained_intron|CODING|ENST00000525036|);FS=0.000;HRun=1;HaplotypeScore=0.0000;MQ=37.00;MQ0=0;QD=19.65	GT:AD:DP:GQ:PL	./.	./.	1/1:0,2:2:6.01:62,6,0	1/1:0,1:1:3.01:31,3,0
chr20	36779424	.	G	A	128.76	.	AC=1;AF=0.13;AN=8;BaseQRankSum=0.610;CSQ=non_coding_exon_variant&amp;nc_transcript_variant|||ENSG00000198959|TGM2|ENST00000468262|4/10||,non_coding_exon_variant&amp;nc_transcript_variant|||ENSG00000198959|TGM2|ENST00000485572|4/13||,non_coding_exon_variant&amp;nc_transcript_variant|||ENSG00000198959|TGM2|ENST00000474777|4/6||,stop_gained|Cag/Tag|Q/*|ENSG00000198959|TGM2|ENST00000373403|5/7||,stop_gained|Cag/Tag|Q/*|ENSG00000198959|TGM2|ENST00000361475|4/13||,stop_gained|Cag/Tag|Q/*|ENSG00000198959|TGM2|ENST00000453095|5/5||,stop_gained|Cag/Tag|Q/*|ENSG00000198959|TGM2|ENST00000536701|3/12||,stop_gained|Cag/Tag|Q/*|ENSG00000198959|TGM2|ENST00000536724|3/12||;DP=196;Dels=0.00;EFF=EXON(MODIFIER|||||TGM2|processed_transcript|CODING|ENST00000485572|4),EXON(MODIFIER|||||TGM2|retained_intron|CODING|ENST00000468262|4),EXON(MODIFIER|||||TGM2|retained_intron|CODING|ENST00000474777|4),STOP_GAINED(HIGH|NONSENSE|Cag/Tag|Q157*|158|TGM2|protein_coding|CODING|ENST00000453095|),STOP_GAINED(HIGH|NONSENSE|Cag/Tag|Q157*|279|TGM2|protein_coding|CODING|ENST00000373403|),STOP_GAINED(HIGH|NONSENSE|Cag/Tag|Q157*|687|TGM2|protein_coding|CODING|ENST00000361475|),STOP_GAINED(HIGH|NONSENSE|Cag/Tag|Q76*|606|TGM2|protein_coding|CODING|ENST00000536701|),STOP_GAINED(HIGH|NONSENSE|Cag/Tag|Q97*|627|TGM2|protein_coding|CODING|ENST00000536724|);FS=1.447;HRun=0;HaplotypeScore=4.5749;MQ=36.22;MQ0=0;MQRankSum=-0.814;QD=3.90;ReadPosRankSum=-0.570	GT:AD:DP:GQ:PL	0/0:49,1:52:63.68:0,64,969	0/0:17,0:17:30.05:0,30,320	0/0:93,0:94:99:0,216,2384	0/1:24,9:33:99:165,0,505
chrX	17819377	.	T	C	7515.25	.	AC=8;AF=1.00;AN=8;CSQ=downstream_gene_variant|||ENSG00000248906||ENST00000509491|||,missense_variant|Atg/Gtg|M/V|ENSG00000131831|RAI2|ENST00000360011|3/3|unknown(0)|tolerated(1),missense_variant|Atg/Gtg|M/V|ENSG00000131831|RAI2|ENST00000331511|3/3|unknown(0)|tolerated(1),missense_variant|Atg/Gtg|M/V|ENSG00000131831|RAI2|ENST00000415486|3/3|unknown(0)|tolerated(1),missense_variant|Atg/Gtg|M/V|ENSG00000131831|RAI2|ENST00000451717|2/2|unknown(0)|tolerated(1),missense_variant|Atg/Gtg|M/V|ENSG00000131831|RAI2|ENST00000545871|3/3|unknown(0)|tolerated(1);DP=319;Dels=0.00;EFF=DOWNSTREAM(MODIFIER|||||RP3-389A20.4|processed_transcript|NON_CODING|ENST00000509491|),NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|Atg/Gtg|M202V|480|RAI2|protein_coding|CODING|ENST00000415486|),NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|Atg/Gtg|M252V|530|RAI2|protein_coding|CODING|ENST00000331511|),NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|Atg/Gtg|M252V|530|RAI2|protein_coding|CODING|ENST00000360011|),NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|Atg/Gtg|M252V|530|RAI2|protein_coding|CODING|ENST00000451717|),NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|Atg/Gtg|M252V|530|RAI2|protein_coding|CODING|ENST00000545871|);FS=0.000;HRun=1;HaplotypeScore=7.7850;MQ=36.33;MQ0=0;QD=23.56	GT:AD:DP:GQ:PL	1/1:0,125:126:99:2343,237,0	1/1:0,26:26:78.14:837,78,0	1/1:0,90:92:99:2640,244,0	1/1:0,74:75:99:1695,171,0
```

### Example 

list the available SO:terms:

```
$  java -jar dist/vcffilterso.jar -S

SO:0001650	inframe_variant
SO:0001791	exon_variant
SO:0001792	non_coding_exon_variant
SO:0001593	minus_2_frameshift_variant
SO:0001594	plus_1_frameshift_variant
SO:0001595	plus_2_frameshift_variant
SO:0001596	transcript_secondary_structure_variant
SO:0001597	compensatory_transcript_secondary_structure_variant
SO:0001598	translational_product_structure_variant
SO:0001599	3D_polypeptide_structure_variant
(...)
```

## See also

 * GroupByGene
 * VCFPredictions
 * http://www.sequenceontology.org/browser/obob.cgi

## History

 * 2018-02-07 refactored a large part of the code
 * 2017 moved to jcommander


 * END_DOC
 *
 */
@Program(name="vcffilterso",
	keywords={"vcf","filter","sequenceontology","prediction","so"},
	description="Filter a VCF file annotated with SNPEff or VEP with terms from Sequence-Ontology. Reasoning : Children of user's SO-terms will be also used.<"
	)
public class VcfFilterSequenceOntology
	extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfFilterSequenceOntology.class).make();
	private static final String GT_FILTER_RESET_TO_NOCALL="NO_CALL";
	
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@ParametersDelegate
	private CtxWriterFactory component = new CtxWriterFactory();
	@Parameter(names={"-S","--showacn"},description="list the available SO accession and exit.")
	private boolean showList = false;


	@XmlType(name="vcffilterso")
	@XmlRootElement(name="vcffilterso")
	@XmlAccessorType(XmlAccessType.FIELD)
	public static class CtxWriterFactory 
		implements VariantContextWriterFactory
			{
			
			@Parameter(names={"-i","--invert"},description="invert SO:Term selection")
			private boolean invertSoTerms = false;
		
			@Parameter(names={"-d","--noreasoning"},description="disable reasoning: do not use SO term's children.")
			private boolean disableReasoning = false;
		
			@Parameter(names={"-fi","--filterin"},description="Do not discard variant but add this FILTER its' prediction is found in the database")
			private String filterIn = "";
		
			@Parameter(names={"-fo","--filterout"},description="Do not discard variant but add this FILTER its' prediction is NOT found in the database")
			private String filterOut = "";
		
			@Parameter(names={"-r","--rmatt","--remove-attribute"},description="Do not remove the variant itself, just remove the mismatching Prediction in the INFO column: e.g: CSQ=OK,OK,NO,OK -> CSQ=OK,OK,OK.")
			private boolean removeUnusedAttribute = false;
		
			@Parameter(names={"-R","--rmnoatt","--remove-variant-if-no-INFO"},description="remove the variant if option -r was used and the is no more attribute")
			private boolean removeVariantIfNoMoreAttribute = false;
		
			@XmlElement(name="accession")
			@Parameter(names={"-A","--acn","--accession"},description="add this SO:ACN. e.g.: 'SO:0001818' Protein altering variant [http://www.sequenceontology.org/miso/current_svn/term/SO:0001818](http://www.sequenceontology.org/miso/current_svn/term/SO:0001818)")
			private List<String> userTermsAsString = new ArrayList<>();
		
			@XmlElement(name="accession-file")
			@Parameter(names={"-f","--acnfile"},description="file of SO accession numbers, one per line")
			private File userAcnFile = null;
		
			@XmlElement(name="filter-genotypes")
			@Parameter(hidden=true,names={"-fg","--filter-genotype"},description="[20180205] Experimental. Filter genotypes having NO ALT allele carrying a matching prediction."
					+ "Only works when I can extract a valid Allele from a prediction. Use with care. "
					+ "Idea is to FILTER out genotype '0/2' of multialleleic variant where only '0/1' is of interest."
					+ "Special FILTER named '"+GT_FILTER_RESET_TO_NOCALL+"' will set the Genotype to NO_CALL. "
					+ "Use with care, it's not always possible to find the allele corresponding to an annotation."
					)
			private String filterGenotypesStr=null;

			@ParametersDelegate
			private VariantAttributesRecalculator recalculator = new VariantAttributesRecalculator();

			
			@XmlTransient
			@Parameter(names={"-owluri","--owluri"},description="Experimental. If not empty, don't use the internal SO ontology but load a OWL description of the ontology. Tested with https://github.com/The-Sequence-Ontology/SO-Ontologies/raw/master/releases/so-xp.owl/so-xp-simple.owl")
			private String owluri = "";
			
			@XmlTransient
			private SequenceOntologyTree sequenceOntologyTree = SequenceOntologyTree.createDefault();

			
			/* all sequence terms */
			@XmlTransient
			private final Set<SequenceOntologyTree.Term> user_terms=new HashSet<SequenceOntologyTree.Term>();

			
			private class CtxWriter extends DelegateVariantContextWriter
				{
				private final List<AbstractPredictionHandler> predictionHandlers = new ArrayList<>();				
				private final boolean removeVariantIfNoMoreAttribute = CtxWriterFactory.this.removeVariantIfNoMoreAttribute;
				private final boolean removeUnusedAttribute = CtxWriterFactory.this.removeUnusedAttribute;
				private final String filterGenotypesStr = CtxWriterFactory.this.filterGenotypesStr;
				private final String filterIn = CtxWriterFactory.this.filterIn;
				private final String filterOut = CtxWriterFactory.this.filterOut;

				private abstract class AbstractPredictionHandler
					{
					final List<Object> predStrings = new ArrayList<>();
					final Set<Allele> matching_alleles = new HashSet<>();
					boolean keepFlag = false;
					abstract boolean isValid();
					abstract String getTag();
					
					AbstractPredictionHandler reset(final VariantContext ctx) {
						this.predStrings.clear();
						this.keepFlag = false;
						this.matching_alleles.clear();
						this.matching_alleles.add(ctx.getReference());
						return this;
						}
					boolean supportFilterAlleles() { return false;}
					abstract AbstractPredictionHandler visit(final VariantContext ctx,final VariantContextBuilder vcb);
					void updateInfo(final VariantContext ctx,final VariantContextBuilder vcb)
						{
						if(CtxWriter.this.removeUnusedAttribute) {
							vcb.rmAttribute(this.getTag());
							if(!this.predStrings.isEmpty()) {
								vcb.attribute(this.getTag(), this.predStrings);
								}
							}
						}
					
					}
				/** handler for VEP */
				private class VepPredictionHandler extends AbstractPredictionHandler
					{
					private final VepPredictionParser parser;
					VepPredictionHandler(final VCFHeader header)
						{
						this.parser =new VepPredictionParserFactory().
								header(header).get().
								sequenceOntologyTree(CtxWriterFactory.this.sequenceOntologyTree);
						}
					@Override
					String getTag() {
						return this.parser.getTag();
						}
					@Override
					boolean isValid() { return this.parser.isValid();}
					@Override
					boolean supportFilterAlleles() {
						return true;
						}
					AbstractPredictionHandler visit(final VariantContext ctx,final VariantContextBuilder vcb)
						{
						if(!ctx.hasAttribute(this.getTag())) return this;
						for(final VepPredictionParser.VepPrediction pred : this.parser.getPredictions(ctx))
							{

							if(pred==null) continue;
							if(hasUserTemLabel(pred.getSOTerms()))
								{
								if(isRecodingGenotypes()) {
									if(pred.getAllele()!=null) this.matching_alleles.add(pred.getAllele());
									}
								this.keepFlag=true;
								if(CtxWriter.this.removeUnusedAttribute) {
									this.predStrings.add(pred.getOriginalAttributeAsString());
									}
								}
							else
								{
								//nothing
								}
							}
						updateInfo(ctx,vcb);
						return this;
						}
					}
				
				/** handler for SNPEFF */
				private class SnpEffPredictionHandler extends AbstractPredictionHandler
					{
					private final SnpEffPredictionParser parser;
					SnpEffPredictionHandler(final VCFHeader header)
						{
						this.parser =new SnpEffPredictionParserFactory().
								header(header).get().
								sequenceOntologyTree(CtxWriterFactory.this.sequenceOntologyTree);
						}
					@Override
					boolean supportFilterAlleles() {
						return true;
						}
					@Override
					String getTag() {
						return this.parser.getTag();
						}
					@Override
					boolean isValid() { return this.parser.isValid();}
					@Override
					AbstractPredictionHandler visit(final VariantContext ctx,final VariantContextBuilder vcb)
						{
						if(!ctx.hasAttribute(this.getTag())) return this;
						for(final SnpEffPredictionParser.SnpEffPrediction pred : this.parser.getPredictions(ctx))
							{
							if(pred==null) continue;
							if(hasUserTemLabel(pred.getSOTerms()))
								{
								if(isRecodingGenotypes()) {
									final Allele alt = pred.getAllele();
									if(alt!=null) this.matching_alleles.add(alt);
									}
								this.keepFlag=true;
								if(CtxWriter.this.removeUnusedAttribute) {
									this.predStrings.add(pred.getOriginalAttributeAsString());
									}
								}
							}
						updateInfo(ctx,vcb);
						return this;
						}
					}
				/** handler for my Prediction */
				private class MyPredictionHandler extends AbstractPredictionHandler
					{
					private final MyPredictionParser parser;
					private final boolean _is_valid;
					MyPredictionHandler(final VCFHeader header)
						{
						this.parser = new MyPredictionParser(header).
								sequenceOntologyTree(CtxWriterFactory.this.sequenceOntologyTree);
						this._is_valid = header.getInfoHeaderLine(this.parser.getTag())!=null;
						}
					@Override
					String getTag() {
						return this.parser.getTag();
						}
					@Override
					boolean isValid() { return _is_valid;}
					@Override
					AbstractPredictionHandler visit(final VariantContext ctx,final VariantContextBuilder vcb)
						{
						if(!ctx.hasAttribute(this.getTag())) return this;
						for(final MyPredictionParser.MyPrediction pred : this.parser.getPredictions(ctx))
							{
							if(pred==null) continue;
							if(hasUserTemLabel(pred.getSOTerms()))
								{
								this.keepFlag=true;
								if(CtxWriter.this.removeUnusedAttribute) {
									this.predStrings.add(pred.getOriginalAttributeAsString());
									}
								}
							}
						updateInfo(ctx,vcb);
						return this;
						}
					}
				
				/** prediction for ANN person */
				private class AnnPredictionHandler extends AbstractPredictionHandler
					{
					private final AnnPredictionParser parser;
					AnnPredictionHandler(final VCFHeader header)
						{
						this.parser = new AnnPredictionParserFactory().header(header).get().
								sequenceOntologyTree(CtxWriterFactory.this.sequenceOntologyTree);
						}
					@Override
					String getTag() {
						return this.parser.getTag();
						}
					@Override
					boolean isValid() { return this.parser.isValid();}
					@Override
					boolean supportFilterAlleles() {
						return true;
						}
					
					AbstractPredictionHandler visit(final VariantContext ctx,final VariantContextBuilder vcb)
						{
						if(!ctx.hasAttribute(this.getTag())) return this;
						for(final AnnPredictionParser.AnnPrediction pred : this.parser.getPredictions(ctx))
							{
							if(pred==null) continue;
							if(hasUserTemLabel(pred.getSOTerms()))
								{
								if(isRecodingGenotypes() && !StringUtil.isBlank(pred.getAllele()))
									{
									this.matching_alleles.add(Allele.create(pred.getAllele(),false));
									}
								this.keepFlag=true;
								if(CtxWriter.this.removeUnusedAttribute) {
									this.predStrings.add(pred.getOriginalAttributeAsString());
									}
								}
							}
						updateInfo(ctx,vcb);
						return this;
						}
					}
				
				
				CtxWriter(final VariantContextWriter delegate) {
					super(delegate);
					}
				
				@Override
				public void writeHeader(final VCFHeader header) {
					final VCFHeader header2= new VCFHeader(header);
					
					final String termlist = String.join(", ",CtxWriterFactory.this.user_terms.stream().
								map(S->S.getAcn()+"("+S.getLabel()+")").
								collect(Collectors.toSet()))
								;
					if(!CtxWriterFactory.this.filterIn.isEmpty()) {
						header2.addMetaDataLine(new VCFFilterHeaderLine(CtxWriterFactory.this.filterIn,
								"Variant having SO terms:"+ termlist));
						}
					if(!CtxWriterFactory.this.filterOut.isEmpty()) {
						header2.addMetaDataLine(new VCFFilterHeaderLine(CtxWriterFactory.this.filterOut,
								"Variant non having SO terms :" + termlist));
						}
					
					if(!StringUtil.isBlank(this.filterGenotypesStr) &&
						!GT_FILTER_RESET_TO_NOCALL.equals(this.filterGenotypesStr))
						{
						header2.addMetaDataLine(new VCFFormatHeaderLine(
								VCFConstants.GENOTYPE_FILTER_KEY,
								1,VCFHeaderLineType.String,
								"Genotype was filterered by vcffilterso : "+this.filterGenotypesStr
								));
						}
					
					AbstractPredictionHandler ph = new VepPredictionHandler(header);
					if(ph.isValid()) this.predictionHandlers.add(ph);
					ph = new SnpEffPredictionHandler(header);
					if(ph.isValid()) this.predictionHandlers.add(ph);
					ph = new MyPredictionHandler(header);
					if(ph.isValid()) this.predictionHandlers.add(ph);
					ph = new AnnPredictionHandler(header);
					if(ph.isValid()) this.predictionHandlers.add(ph);
					
					CtxWriterFactory.this.recalculator.setHeader(header2);
					super.writeHeader(header2);
					}
				
				@Override
				public void add(final VariantContext ctx) {
					boolean variant_has_one_matching_pred = true;
					if(this.predictionHandlers.isEmpty() && !this.removeUnusedAttribute) {
						variant_has_one_matching_pred = false;
						}
					
					
					final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
					
					/* loop over each handler to detect the matching predictions */
					for(final AbstractPredictionHandler handler:this.predictionHandlers) {
						handler.reset(ctx).visit(ctx,vcb);
					}
					
					/* FILTER genotypes having NO causal ATL allele */
					if(isRecodingGenotypes() && 
							!this.predictionHandlers.isEmpty()) {
						/* samples to be FILTERED */
						final Set<String> invalidSamples = new HashSet<>(ctx.getNSamples());
						/* loop over each genotype */
						for(final String sample: ctx.getSampleNames())
							{
							boolean sample_is_ok = false;
							for(final AbstractPredictionHandler handler:this.predictionHandlers) {							
								if(!handler.supportFilterAlleles()) continue;
								final Genotype gt = ctx.getGenotype(sample);
								if(gt==null || gt.isNoCall() || gt.isHomRef() || gt.isFiltered())
									{
									sample_is_ok = true;
									break;
									}
								if(gt.getAlleles().stream().
									filter(A->A.isCalled() && !A.isReference()).
									anyMatch(A->handler.matching_alleles.contains(A)))
									{
									sample_is_ok = true;
									break;
									}
								}
							if(!sample_is_ok) {
								invalidSamples.add(sample);
								}
							}
						final Function<Genotype,Genotype> convertGt = G->{
							/* sample is not invalid */
							if(!invalidSamples.contains(G.getSampleName()))
								{
								return G;
								}
							/* sample is invalid and we reset to NO_CALL (./.) */
							else if(GT_FILTER_RESET_TO_NOCALL.equals(this.filterGenotypesStr))
								{
								return GenotypeBuilder.createMissing(G.getSampleName(), G.getPloidy());
								}
							/* sample is invalid and we set the GT FILTER */
							else
								{
								
								return new GenotypeBuilder(G).filter(this.filterGenotypesStr).make();
								}	
							};
						
						/* update genotypes */
						vcb.genotypes(ctx.getGenotypes().stream().
								map(convertGt).
								collect(Collectors.toList()));
						}
					
					/* all attributes have been removed ? should we keep this variant ?*/
					if( this.removeUnusedAttribute &&
						this.predictionHandlers.stream().allMatch(P->P.predStrings.isEmpty())
						)
						{
						if(this.removeVariantIfNoMoreAttribute) return;
						variant_has_one_matching_pred=false;
						}
					
					/* all handlers failed */
					if(variant_has_one_matching_pred && 
						!this.removeUnusedAttribute &&
						this.predictionHandlers.stream().allMatch(P->!P.keepFlag))
						{
						variant_has_one_matching_pred=false;
						}
					
					
					if(!StringUtil.isBlank(this.filterIn))
						{
						if(variant_has_one_matching_pred){
							vcb.filter(this.filterIn);
							}
						else if( !ctx.filtersWereApplied()) {
							vcb.passFilters();
							}
						}
					else  if(!StringUtil.isBlank(this.filterOut)) {
						if(variant_has_one_matching_pred && !ctx.filtersWereApplied()) {
							vcb.passFilters();
							}
						else if(!variant_has_one_matching_pred) {
							vcb.filter(CtxWriterFactory.this.filterOut);
							}
						}
					else if(!variant_has_one_matching_pred)
						{
						return ;
						}
					
					final VariantContext ctx3;
					if(GT_FILTER_RESET_TO_NOCALL.equals(this.filterGenotypesStr))
					 	{
						ctx3 = CtxWriterFactory.this.recalculator.apply(vcb.make());
					 	}
					 else
					 	{
						ctx3= vcb.make(); 
					 	}
					super.add(ctx3);
					}
				
				private boolean hasUserTemLabel(final Collection<SequenceOntologyTree.Term> ctxTerms)
					{
					if(ctxTerms==null || ctxTerms.isEmpty()) return false;
					return ctxTerms.
							stream().
							anyMatch(S->user_terms.contains(S));
					}				
				}
			
			@Override
			public int initialize() {
				
				try {
					this.sequenceOntologyTree = readSequenceOntologyTree(this.owluri);
				} catch(final IOException err)
					{
					LOG.error(err);
					return -1;
					}
				
				
				final boolean reasoning = !this.disableReasoning;
				if(this.userAcnFile!=null)
					{
					try {
						this.parseAccessionsFile(this.userAcnFile);
						} 
					catch(final IOException err)
						{
						LOG.error(err);
						return -1;
						}
					}
				final Function<String,SequenceOntologyTree.Term> acn2term = (ACN)->{
					SequenceOntologyTree.Term T = this.sequenceOntologyTree.getTermByAcn(ACN);
					if(T==null)
						{
						T = this.sequenceOntologyTree.getTermByLabel(ACN);
						}
					if(T==null)
						{
						throw new JvarkitException.UserError("Unknown SO:Accession/label \""+ACN+"\"");
						}
					return T;
					};
				final Set<SequenceOntologyTree.Term> tmpSet1 = new HashSet<>();
					this.userTermsAsString.stream().
					map(S->S.trim()).
					filter(S->!S.isEmpty()).
					flatMap(S->Arrays.stream(S.split("[,;& \t]"))).
					filter(S->!S.isEmpty()).
					map(acn2term).
					forEach(t->{
						tmpSet1.add(t);
						if(reasoning) tmpSet1.addAll(t.getAllDescendants());					
					});
				
				
				if(this.invertSoTerms)
					{
					final Set<SequenceOntologyTree.Term> tmpSet2 = new HashSet<>(this.sequenceOntologyTree.getTerms());
					tmpSet2.removeAll(tmpSet1);
					this.user_terms.addAll(tmpSet2);
					
					}
				else
					{
					this.user_terms.addAll(tmpSet1);
					}
				
				if(this.user_terms.isEmpty())
					{
					LOG.warn("No SO: term found ");
					}
				LOG.info("Will be using :"+this.user_terms.stream().
						map(T->T.getAcn()+"("+T.getLabel()+")").
						collect(Collectors.joining(" ")));
				
				this.userTermsAsString.clear();//we don't need this anymore

						
				
				if(!StringUtil.isBlank(this.filterIn) && !StringUtil.isBlank(this.filterOut)) {
					LOG.error("Option filterIn && filterOut both defined.");
					return -1;
				}
				
				return 0;
				}
			
			@Override
			public VariantContextWriter open(final VariantContextWriter delegate) {
				return new CtxWriter(delegate);
				}
			
			@Override
			public void close() throws IOException {
				}
			
			private boolean isRecodingGenotypes() {
				return !StringUtil.isBlank(this.filterGenotypesStr);
			}
			
			private void parseAccessionsFile(final File f) throws IOException
				{
				final BufferedReader in=IOUtils.openFileForBufferedReading(f);
				in.lines().
					filter(L->!L.startsWith("#")).
					map(L->L.trim()).
					filter(L->!L.isEmpty()).
					forEach(L->{
						this.userTermsAsString.add(L);
					});
				in.close();
				}

			
			}
	
	public VcfFilterSequenceOntology()
		{
		}
	
	@Override
	protected int doVcfToVcf(final String inputName,final VCFIterator iter, final VariantContextWriter delegate) {	
		final VariantContextWriter out = this.component.open(delegate);
		final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(iter.getHeader()).logger(LOG);
		out.writeHeader(iter.getHeader());
		while(iter.hasNext())
			{
			out.add(progress.watch(iter.next()));
			}
		out.close();
		progress.finish();
		return 0;
		}
			
	/** return SequenceOntologyTree
	 * 
	 * @param owluri if null, use getInstance
	 * @return a SequenceOntologyTree
	 * @throws IOException
	 */
	private static SequenceOntologyTree readSequenceOntologyTree(final String owluri) throws IOException {
		if(!StringUtil.isBlank(owluri)) {
			LOG.info("loading 'SO' tree from "+owluri);
			SequenceOntologyTree sequenceOntologyTree = SequenceOntologyTree.fromUri(owluri.trim());
			LOG.info("Done loading SO.");
			return sequenceOntologyTree;
			}
		else
			{
			return SequenceOntologyTree.getInstance();
			}
		}
	
	@Override
	public int doWork(final List<String> args) {
		try {
			if(this.showList)
				{
				final SequenceOntologyTree tree = readSequenceOntologyTree(this.component.owluri);
				final PrintWriter pw=super.openFileOrStdoutAsPrintWriter(this.outputFile);
				for(final SequenceOntologyTree.Term t: tree.getTerms())
					{
					pw.println(t.getAcn()+"\t"+t.getLabel());
					}
				pw.flush();
				pw.close();
				return 0;
				}
			
			if(this.component.initialize()!=0) return -1;
			
			return doVcfToVcf(args, this.outputFile);
			}
		catch(final Exception err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.component);
			}
		}
		
	
	public static void main(final String[] args)
		{
		new VcfFilterSequenceOntology().instanceMainWithExit(args);
		}
	}
