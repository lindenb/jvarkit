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
*/
package com.github.lindenb.jvarkit.tools.ensemblreg;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.seekablestream.SeekableStreamFactory;

import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BedFeature;
import org.broad.igv.bbfile.BigBedIterator;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;

import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.igv.SeekableStreamAdaptor;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import htsjdk.variant.vcf.VCFIterator;
/**
BEGIN_DOC

## Example

original vcf:

```
$ curl -kLs "https://raw.githubusercontent.com/jamescasbon/PyVCF/master/vcf/test/gatk.vcf" | cut -f 1-8
##fileformat=VCFv4.1
(...)
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr22	42522392	rs28371738	G	A	2951.95	.	AC=2;AF=0.143;AN=14;BaseQRankSum=0.375;DB;DP=1506;DS;Dels=0.00;FS=0.000;HRun=0;HaplotypeScore=123.5516;MQ=253.92;MQ0=0;MQRankSum=0.685;QD=5.90;ReadPosRankSum=0.590
chr22	42522613	rs1135840	G	C	11611.03	.	AC=6;AF=0.429;AN=14;BaseQRankSum=16.289;DB;DP=1518;DS;Dels=0.03;FS=0.000;HRun=0;HaplotypeScore=142.5716;MQ=242.46;MQ0=0;MQRankSum=2.010;QD=9.16;ReadPosRankSum=-1.731
chr22	42522755	.	C	G	36.98	.	AC=1;AF=0.071;AN=14;BaseQRankSum=-14.866;DP=1527;DS;Dels=0.01;FS=0.000;HRun=0;HaplotypeScore=253.4254;MQ=197.36;MQ0=2;MQRankSum=-10.810;QD=0.15;ReadPosRankSum=-17.244
chr22	42523003	rs116917064	A	G	7113.55	.	AC=8;AF=0.571;AN=14;BaseQRankSum=6.026;DB;DP=1433;DS;Dels=0.00;FS=0.000;HRun=1;HaplotypeScore=101.7894;MQ=182.04;MQ0=0;MQRankSum=-2.501;QD=4.96;ReadPosRankSum=8.294
chr22	42523077	.	A	G	54.31	.	AC=1;AF=0.071;AN=14;BaseQRankSum=-0.563;DP=1521;DS;Dels=0.00;FS=0.000;HRun=0;HaplotypeScore=54.8434;MQ=164.04;MQ0=1;MQRankSum=-2.419;QD=2.59;ReadPosRankSum=-1.229
(...)
```

execute:

```bash
$  curl -kLs "https://raw.githubusercontent.com/jamescasbon/PyVCF/master/vcf/test/gatk.vcf" |\
   java -jar dist/vcfensemblreg.jar > out.vcf

```


out.vcf:

```
##fileformat=VCFv4.1
(...)
##INFO=<ID=ATF3,Number=1,Type=Float,Description="Overlap summary of ATF3 ChipSeq binding peaks across available datasets http://ngs.sanger.ac.uk/production/ensembl/regulation//hg19/tf
bs/ATF3.bw">
##INFO=<ID=BAF155,Number=1,Type=Float,Description="Overlap summary of BAF155 ChipSeq binding peaks across available datasets http://ngs.sanger.ac.uk/production/ensembl/regulation//hg1
9/tfbs/BAF155.bw">
##INFO=<ID=BAF170,Number=1,Type=Float,Description="Overlap summary of BAF170 ChipSeq binding peaks across available datasets http://ngs.sanger.ac.uk/production/ensembl/regulation//hg1
9/tfbs/BAF170.bw">
##INFO=<ID=BuildOverview,Number=1,Type=String,Description="Ensembl Regulatory annotation of regional function http://ngs.sanger.ac.uk/production/ensembl/regulation//hg19/overview/RegB
uild.bb">
(...)
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr22	42522392	rs28371738	G	A	2951.95	.	AC=2;AF=0.143;AN=14;BaseQRankSum=0.375;BuildOverview=open_646005|UnannotatedOpenChromatinRegions;DB;DP=1506;DS;
Dels=0.00;FS=0.000;HRun=0;HaplotypeScore=123.5516;MQ=253.92;MQ0=0;MQRankSum=0.685;QD=5.90;ReadPosRankSum=0.590;Segway_17_18=1.0;Segway_17_20=2.0;Segway_17_21=3.0;Segway_17_22=1.0;Segw
ay_17_5=9.0;Segway_17_7=1.0;Segway_17_A549_projected=open_646005|InactiveRegions;Segway_17_GM12878_projected=open_646005|InactiveRegions;Segway_17_H1HESC_projected=open_646005|Inactiv
eRegions;Segway_17_HELAS3_projected=open_646005|InactiveRegions;Segway_17_HEPG2_projected=open_646005|InactiveRegions;Segway_17_HMEC_projected=open_646005|InactiveRegions;Segway_17_HS
MM_projected=open_646005|InactiveRegions;Segway_17_HUVEC_projected=open_646005|InactiveRegions;Segway_17_K562_projected=open_646005|InactiveRegions;Segway_17_NHA_projected=open_646005
|InactiveRegions;Segway_17_NHDFAD_projected=open_646005|InactiveRegions;Segway_17_NHEK_projected=open_646005|UnannotatedActiveOpenChromatinRegions;Segway_17_NHLF_projected=open_646005
|InactiveRegions
chr22	42522613	rs1135840	G	C	11611.03	.	AC=6;AF=0.429;AN=14;BaseQRankSum=16.289;DB;DP=1518;DS;Dels=0.03;FS=0.000;HRun=0;HaplotypeScore=142.5716
;MQ=242.46;MQ0=0;MQRankSum=2.010;QD=9.16;ReadPosRankSum=-1.731;Segway_17_18=1.0;Segway_17_20=4.0;Segway_17_21=3.0;Segway_17_22=1.0;Segway_17_5=7.0;Segway_17_7=1.0
chr22	42522755	.	C	G	36.98	.	AC=1;AF=0.071;AN=14;BaseQRankSum=-14.866;DP=1527;DS;Dels=0.01;FS=0.000;HRun=0;HaplotypeScore=253.4254;MQ=197.36;MQ0=2;M
QRankSum=-10.810;QD=0.15;ReadPosRankSum=-17.244;Segway_17_15=1.0;Segway_17_20=3.0;Segway_17_21=1.0;Segway_17_22=1.0;Segway_17_3=1.0;Segway_17_5=10.0
chr22	42523003	rs116917064	A	G	7113.55	.	AC=8;AF=0.571;AN=14;BaseQRankSum=6.026;DB;DP=1433;DS;Dels=0.00;FS=0.000;HRun=1;HaplotypeScore=101.7894;MQ=182.0
4;MQ0=0;MQRankSum=-2.501;QD=4.96;ReadPosRankSum=8.294;Segway_17_15=1.0;Segway_17_18=1.0;Segway_17_20=2.0;Segway_17_22=3.0;Segway_17_3=1.0;Segway_17_5=9.0
chr22	42523077	.	A	G	54.31	.	AC=1;AF=0.071;AN=14;BaseQRankSum=-0.563;DP=1521;DS;Dels=0.00;FS=0.000;HRun=0;HaplotypeScore=54.8434;MQ=164.04;MQ0=1;MQR
ankSum=-2.419;QD=2.59;ReadPosRankSum=-1.229;Segway_17_18=1.0;Segway_17_20=1.0;Segway_17_22=4.0;Segway_17_24=1.0;Segway_17_3=1.0;Segway_17_5=9.0;Segway_17_DND41_segments=24_gene_86550|
TranscriptionAssociated
chr22	42523209	rs28371730	T	C	15556.89	.	AC=8;AF=0.571;AN=14;BaseQRankSum=3.458;DB;DP=1509;DS;Dels=0.01;FS=0.000;HRun=0;HaplotypeScore=120.8206;
MQ=221.07;MQ0=0;MQRankSum=-4.945;QD=10.31;ReadPosRankSum=0.639;Segway_17_18=1.0;Segway_17_20=1.0;Segway_17_22=4.0;Segway_17_24=1.0;Segway_17_3=1.0;Segway_17_5=9.0;Segway_17_DND41_segm
ents=24_gene_86550|TranscriptionAssociated
chr22	42523211	rs2004511	T	C	2445.52	.	AC=2;AF=0.143;AN=14;BaseQRankSum=10.587;DB;DP=1509;DS;Dels=0.00;FS=0.000;HRun=0;HaplotypeScore=102.7564;MQ=221.
50;MQ0=0;MQRankSum=-6.926;QD=4.89;ReadPosRankSum=2.057;Segway_17_18=1.0;Segway_17_20=1.0;Segway_17_22=4.0;Segway_17_24=1.0;Segway_17_3=1.0;Segway_17_5=9.0;Segway_17_DND41_segments=24_
gene_86550|TranscriptionAssociated
chr22	42523409	rs1985842	G	T	6801.90	.	AC=6;AF=0.429;AN=14;BaseQRankSum=20.509;BuildOverview=open_646006|UnannotatedOpenChromatinRegions;DB;DP=1454;DS
;Dels=0.00;FS=0.000;HRun=0;HaplotypeScore=150.8967;MQ=200.12;MQ0=0;MQRankSum=4.472;QD=5.65;ReadPosRankSum=9.396;Segway_17_18=1.0;Segway_17_20=3.0;Segway_17_21=1.0;Segway_17_22=4.0;Seg
way_17_24=2.0;Segway_17_5=6.0;Segway_17_A549_projected=open_646006|InactiveRegions;Segway_17_DND41_segments=24_gene_86550|TranscriptionAssociated;Segway_17_GM12878_projected=open_6460
06|InactiveRegions;Segway_17_H1HESC_projected=open_646006|InactiveRegions;Segway_17_HELAS3_projected=open_646006|InactiveRegions;Segway_17_HEPG2_projected=open_646006|InactiveRegions;
Segway_17_HMEC_projected=open_646006|InactiveRegions;Segway_17_HSMMT_segments=24_gene_74042|TranscriptionAssociated;Segway_17_HSMM_projected=open_646006|InactiveRegions;Segway_17_HUVE
C_projected=open_646006|InactiveRegions;Segway_17_K562_projected=open_646006|InactiveRegions;Segway_17_NHA_projected=open_646006|InactiveRegions;Segway_17_NHDFAD_projected=open_646006
|InactiveRegions;Segway_17_NHEK_projected=open_646006|InactiveRegions;Segway_17_NHLF_projected=open_646006|InactiveRegions
chr22	42523805	rs28371725	C	T	1637.33	.	AC=1;AF=0.071;AN=14;BaseQRankSum=-0.379;BuildOverview=open_646007|UnannotatedOpenChromatinRegions;DB;DP=1516;DS
;Dels=0.00;FS=0.000;HRun=2;HaplotypeScore=77.2321;MQ=226.05;MQ0=0;MQRankSum=2.862;QD=6.55;ReadPosRankSum=0.064;Segway_17_11=1.0;Segway_17_13=1.0;Segway_17_15=1.0;Segway_17_20=5.0;Segw
ay_17_21=3.0;Segway_17_22=2.0;Segway_17_3=2.0;Segway_17_5=2.0;Segway_17_A549_projected=open_646007|InactiveRegions;Segway_17_GM12878_projected=open_646007|InactiveRegions;Segway_17_H1
HESC_projected=open_646007|InactiveRegions;Segway_17_HELAS3_projected=open_646007|UnannotatedActiveOpenChromatinRegions;Segway_17_HEPG2_projected=open_646007|InactiveRegions;Segway_17
_HMEC_projected=open_646007|UnannotatedActiveOpenChromatinRegions;Segway_17_HSMM_projected=open_646007|InactiveRegions;Segway_17_HUVEC_projected=open_646007|InactiveRegions;Segway_17_
K562_projected=open_646007|InactiveRegions;Segway_17_NHA_projected=open_646007|InactiveRegions;Segway_17_NHDFAD_projected=open_646007|InactiveRegions;Segway_17_NHEK_projected=open_646
007|UnannotatedActiveOpenChromatinRegions;Segway_17_NHLF_projected=open_646007|InactiveRegions
chr22	42523943	rs16947	A	G	23661.10	.	AC=8;AF=0.571;AN=14;BaseQRankSum=4.602;BuildOverview=open_646007|UnannotatedOpenChromatinRegions;DB;DP=1514;DS;
Dels=0.00;FS=0.000;HRun=0;HaplotypeScore=38.3217;MQ=238.64;MQ0=0;MQRankSum=2.485;QD=15.63;ReadPosRankSum=3.749;Segway_17_11=2.0;Segway_17_13=1.0;Segway_17_15=1.0;Segway_17_18=1.0;Segw
ay_17_20=6.0;Segway_17_21=1.0;Segway_17_22=2.0;Segway_17_3=1.0;Segway_17_5=2.0;Segway_17_A549_projected=open_646007|InactiveRegions;Segway_17_GM12878_projected=open_646007|InactiveReg
ions;Segway_17_H1HESC_projected=open_646007|InactiveRegions;Segway_17_HELAS3_projected=open_646007|UnannotatedActiveOpenChromatinRegions;Segway_17_HEPG2_projected=open_646007|Inactive
Regions;Segway_17_HMEC_projected=open_646007|UnannotatedActiveOpenChromatinRegions;Segway_17_HSMM_projected=open_646007|InactiveRegions;Segway_17_HUVEC_projected=open_646007|InactiveR
egions;Segway_17_K562_projected=open_646007|InactiveRegions;Segway_17_NHA_projected=open_646007|InactiveRegions;Segway_17_NHDFAD_projected=open_646007|InactiveRegions;Segway_17_NHEK_p
rojected=open_646007|UnannotatedActiveOpenChromatinRegions;Segway_17_NHLF_projected=open_646007|InactiveRegions
chr22	42524150	.	C	G	3758.65	.	AC=8;AF=0.571;AN=14;BaseQRankSum=24.314;BuildOverview=open_646008|UnannotatedOpenChromatinRegions;DP=1506;DS;Dels=0.00;
FS=0.000;HRun=1;HaplotypeScore=172.5901;MQ=242.92;MQ0=0;MQRankSum=11.537;QD=2.50;ReadPosRankSum=-9.185;Segway_17_11=3.0;Segway_17_13=1.0;Segway_17_18=2.0;Segway_17_20=10.0;Segway_17_5
=1.0;Segway_17_A549_projected=open_646008|InactiveRegions;Segway_17_GM12878_projected=open_646008|InactiveRegions;Segway_17_H1HESC_projected=open_646008|UnannotatedActiveOpenChromatin
Regions;Segway_17_HELAS3_projected=open_646008|InactiveRegions;Segway_17_HEPG2_projected=open_646008|InactiveRegions;Segway_17_HMEC_projected=open_646008|InactiveRegions;Segway_17_HSM
M_projected=open_646008|InactiveRegions;Segway_17_HUVEC_projected=open_646008|InactiveRegions;Segway_17_K562_projected=open_646008|InactiveRegions;Segway_17_NHA_projected=open_646008|
InactiveRegions;Segway_17_NHDFAD_projected=open_646008|InactiveRegions;Segway_17_NHEK_projected=open_646008|InactiveRegions;Segway_17_NHLF_projected=open_646008|InactiveRegions
chr22	42524435	rs1807313	T	A	5252.25	.	AC=3;AF=0.214;AN=14;BaseQRankSum=-0.192;DB;DP=1526;DS;Dels=0.01;FS=0.000;HRun=1;HaplotypeScore=152.3866;MQ=242.
06;MQ0=0;MQRankSum=1.923;QD=9.99;ReadPosRankSum=3.008;Segway_17_11=4.0;Segway_17_13=1.0;Segway_17_15=1.0;Segway_17_18=1.0;Segway_17_20=5.0;Segway_17_21=1.0;Segway_17_24=1.0;Segway_17_
3=1.0;Segway_17_5=1.0;Segway_17_7=1.0;Segway_17_HEPG2_segments=24_gene_85862|TranscriptionAssociated
chr22	42524696	rs58440431	T	C	6423.61	.	AC=2;AF=0.143;AN=14;BaseQRankSum=3.119;BuildOverview=open_646009|UnannotatedOpenChromatinRegions;DB;DP=1509;DS;
Dels=0.00;FS=0.000;HRun=0;HaplotypeScore=53.0005;MQ=230.78;MQ0=0;MQRankSum=2.825;QD=12.85;ReadPosRankSum=2.051;Segway_17_11=3.0;Segway_17_13=1.0;Segway_17_15=1.0;Segway_17_18=1.0;Segw
ay_17_20=4.0;Segway_17_21=2.0;Segway_17_24=1.0;Segway_17_5=3.0;Segway_17_7=1.0;Segway_17_A549_projected=open_646009|InactiveRegions;Segway_17_GM12878_projected=open_646009|InactiveReg
ions;Segway_17_H1HESC_projected=open_646009|InactiveRegions;Segway_17_HELAS3_projected=open_646009|UnannotatedActiveOpenChromatinRegions;Segway_17_HEPG2_projected=open_646009|Inactive
Regions;Segway_17_HEPG2_segments=24_gene_85862|TranscriptionAssociated;Segway_17_HMEC_projected=open_646009|UnannotatedActiveOpenChromatinRegions;Segway_17_HSMM_projected=open_646009|
UnannotatedActiveOpenChromatinRegions;Segway_17_HUVEC_projected=open_646009|UnannotatedActiveOpenChromatinRegions;Segway_17_K562_projected=open_646009|UnannotatedActiveOpenChromatinRe
gions;Segway_17_NHA_projected=open_646009|InactiveRegions;Segway_17_NHDFAD_projected=open_646009|InactiveRegions;Segway_17_NHEK_projected=open_646009|InactiveRegions;Segway_17_NHLF_pr
ojected=open_646009|InactiveRegions
```


END_DOC
 */
@Program(
	name="vcfensemblreg",
	description="Annotate a VCF with the UCSC genome hub tracks for Ensembl Regulation.",
	keywords={"vcf","ensembl","regulation"}
	)
public class VcfEnsemblReg extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfEnsemblReg.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	
	private static class Track
		{
		String id;
		String type=null;
		String parent=null;
		String shortLabel=null;
		String longLabel=null;
		URL url=null;
		}
	private List<Track> tracks=new ArrayList<Track>();
	
	@Parameter(names="-d",description="trackDB.txt URL.")
	private String trackDBUrl="http://ngs.sanger.ac.uk/production/ensembl/regulation/hg19/trackDb.txt";

	private static final String BrightRed="187,0,0";
	private static final String LightRed="255,105,105";
	private static final String Orange="250,202,0";
	private static final String Blue="10,190,254";
	private static final String Gold="209,157,0";
	private static final String Yellow="255,252,4";
	private static final String LightGray="225,225,225";
	private static final String Gray="127,127,127";
	private static final String DarkGreen="0,176,80";
	

	/** http://ngs.sanger.ac.uk/production/ensembl/regulation/hg19/segmentation_summaries/ */
	private String segway_17SegmentationSummaries(String color)
	{
	if(color==null) return null;
		switch(color)
		{
		case BrightRed: return  "ActivePromoters";
		case LightRed :return "ProximalEnhancer";
		case Orange : return "DistalEnhancer";
		case Blue : return "DistalCTF";
		case DarkGreen: return " TranscriptionAssociated";
		case Gray: return "PolycombRepressed";
		case LightGray: return null;//  Weak signal  | Heterochromatin/Repetitive/Copy Number Variation  
		default: System.err.println("UNKNOWN Seg COLOR "+color);return null;
		}
	}
	/** http://ngs.sanger.ac.uk/production/ensembl/regulation/hg19/projected_segmentations/ */
	private String projectedSegments(String color)
	{
	if(color==null) return null;
		switch(color)
		{
		case BrightRed: return  "PredictedActivePromoters";
		case LightRed :return "PredicteActivedPromoterFlankingRegions";
		case Orange : return "PredictedActiveEnhancers";
		case Blue : return "ActiveCTCFBindingSite";
		case Gold : return "UnannotatedActiveTFBS";
		case Yellow: return "UnannotatedActiveOpenChromatinRegions" ;
		case LightGray: return "InactiveRegions"; 
		default: LOG.warning("projectedSegments: undefined color: "+color); return null;
		}
	}
	/** http://ngs.sanger.ac.uk/production/ensembl/regulation/hg19/overview/ */
	private String regBuildOverview(String color)
	{
	if(color==null) return null;
		switch(color)
		{
		case BrightRed: return  "PredictedPromoters";
		case LightRed :return "PredictedPromoterFlankingRegions";
		case Orange : return "PredictedEnhancers";
		case Blue : return "CTCFBindingSite";
		case Gold : return "UnannotatedTFBS";
		case Yellow: return "UnannotatedOpenChromatinRegions" ;
		default: LOG.warning("regBuildOverview: undefined color: "+color);return null;
		}
	}
	
	/** http://ngs.sanger.ac.uk/production/ensembl/regulation/hg19/segmentations/ */
	private String segway_17CellSegments(String color)
		{
		if(color==null) return null;
			switch(color)
			{
			case BrightRed: return  "ActivePromoters";
			case LightRed :return "ProximalEnhancer";
			case Orange : return "DistalEnhancer";
			case Blue : return "DistalCTF";
			case DarkGreen: return "TranscriptionAssociated";
			case Gray: return "PolycombRepressed";
			case LightGray: return null;//  Weak signal  | Heterochromatin/Repetitive/Copy Number Variation  
			default: System.err.println("UNKNOWN Seg COLOR "+color);return null;
			}
		}

	private void annotate(Track track,File inf,File outf) throws IOException
		{
		boolean contained=false;
		LOG.info("Processing "+track.id+" ("+track.shortLabel+") "+track.url);
		VCFIterator in=VCFUtils.createVCFIteratorFromFile(inf);
		VCFHeader header=in.getHeader();
		VCFInfoHeaderLine info=null;
		
		
		SeekableStream sstream=SeekableStreamFactory.getInstance().getStreamFor(track.url);
		BBFileReader bigFile = new BBFileReader(track.url.toString(), new SeekableStreamAdaptor(sstream));
		
		VariantContextWriter w1=VCFUtils.createVariantContextWriter(outf);
		
		if(bigFile.isBigWigFile())
			{
			info = new VCFInfoHeaderLine(track.id, 1, VCFHeaderLineType.Float,String.valueOf(track.longLabel)+" "+track.url);
			}
		else
			{
			info = new VCFInfoHeaderLine(track.id, 1, VCFHeaderLineType.String,String.valueOf(track.longLabel)+" "+track.url);
			}
		
		header.addMetaDataLine(info);
		w1.writeHeader(in.getHeader());
		
		while(in.hasNext())
			{
			VariantContext ctx=in.next();
			String chrom=ctx.getContig();
			if(!chrom.startsWith("chr")) chrom="chr"+chrom;
			if(!chrom.matches("(chrX|chrY|chr[0-9]|chr1[0-9]|chr2[12])"))
				{
				w1.add(ctx);
				}
			else if(bigFile.isBigWigFile())
				{
				BigWigIterator iter=bigFile.getBigWigIterator(
						chrom,
						ctx.getStart()-1,
						chrom,
						ctx.getStart(),
						contained
						);
				Float wigValue=null;
				while(iter!=null && iter.hasNext() && wigValue==null)
					{
					WigItem item=iter.next();
					wigValue=item.getWigValue();
					}
				if(wigValue==null)
					{
					w1.add(ctx);
					continue;
					}
				VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				vcb.attribute(track.id, wigValue);
				w1.add(vcb.make());
				}
			else
				{
				BigBedIterator iter = bigFile.getBigBedIterator(
						chrom,
						ctx.getStart()-1,
						chrom,
						ctx.getStart(),
						contained
						);
				Set<String> bedValues=new HashSet<String>();
				while(iter!=null && iter.hasNext())
					{
					BedFeature item=iter.next();
					String rest[]=item.getRestOfFields();
					if(rest==null || rest.length!=6)
						{
						System.err.println(track.id+" "+Arrays.toString(item.getRestOfFields()));
						continue;
						}
					String color=null;
					if(track.parent!=null)
						{
						if(track.parent.startsWith("Segway_17SegmentationSummaries"))
							{
							color = segway_17SegmentationSummaries(rest[5]);
							}
						else if(track.parent.startsWith("ProjectedSegments") )
							{
							color = projectedSegments(rest[5]);
							}
						else if(track.parent.startsWith("RegBuildOverview") )
							{
							color = regBuildOverview(rest[5]);
							}
						else if(track.parent.startsWith("Segway_17CellSegments") )
							{
							color = segway_17CellSegments(rest[5]);
							}
						else
							{
							System.err.println("Unknown parent:"+track.parent);
							}
						}
					
					
					if(color==null) continue;
					bedValues.add(rest[0]+"|"+color);
					} 
				if(bedValues.isEmpty())
					{
					w1.add(ctx);
					continue;
					}
				StringBuilder sb=new StringBuilder();
				for(String s:bedValues) { if(sb.length()!=0) sb.append(",");sb.append(s);}
				VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				vcb.attribute(track.id, sb.toString());
				w1.add(vcb.make());
				}
			}
		
		sstream.close();
		in.close();
		w1.close();
		}
	@Override
	protected int doVcfToVcf(String inputName, VCFIterator in, VariantContextWriter out) {
		try {
		File tmpDir= new File(System.getProperty("java.io.tmpdir"));
		File tmp1=File.createTempFile("tmp_",".vcf", tmpDir);
		File tmp2=File.createTempFile("tmp_",".vcf", tmpDir);
		tmp1.deleteOnExit();
		tmp2.deleteOnExit();
		
		VariantContextWriter w1=VCFUtils.createVariantContextWriter(tmp1);
		VCFHeader header=in.getHeader();
		header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
		JVarkitVersion.getInstance().addMetaData(getClass().getSimpleName(), header);
		w1.writeHeader(header);
		while(in.hasNext())
			{
			w1.add(in.next());
			}
		in.close();
		w1.close();
		
	
		for(Track track:this.tracks)
			{
			if(track.url==null) continue;
			annotate(track, tmp1,tmp2);
			tmp1.delete();

			tmp2.renameTo(tmp1);
			}
		in=VCFUtils.createVCFIteratorFromFile(tmp1);
		out.writeHeader(in.getHeader());
		while(in.hasNext())
			{
			out.add(in.next());
			}
		in.close();
		tmp1.delete();
		tmp2.delete();
		return 0;
		} catch(Exception err ) {
			LOG.error(err);
			return -1;
		}
		}
	
	private void parseTrackDB(final URL url) throws IOException
		{
		LOG.info("Parsing "+url);
		BufferedReader r=new BufferedReader(new InputStreamReader(url.openStream(),"UTF-8"));
		String line;
		Track track=null;
		for(;;)
			{
			line=r.readLine();
			
			while(line!=null &&
					!line.isEmpty() &&
					Character.isWhitespace(line.charAt(0)))
					line=line.substring(1);
				
			if(line!=null && line.isEmpty()) continue;
			
			if(	line==null ||
				line.trim().startsWith("track "))
				{
				if(track!=null &&
						track.url!=null &&
						track.type!=null
						)
					{
					this.tracks.add(track);
					}
				if(line==null) break;
				track=null;
				}
			int w=line.indexOf(' ');
			if(w==-1)
				{
				LOG.info("No whitespace in "+line+" ?");
				continue;
				}
			if(track==null) track=new Track();
			final String key=line.substring(0,w).trim();
			final String val=line.substring(w+1).trim();
			switch(key)
				{
				case "track":track.id=val;break;
				case "shortLabel":track.shortLabel=val;break;
				case "longLabel":track.longLabel=val;break;
				case "bigDataUrl":track.url=new URL(url,val);break;
				case "type":track.type=val;break;
				case "parent":track.parent=val.split("[ \t]")[0];break;
				default:break;
				}

			}
		
		CloserUtil.close(r);
		}
	
	@Override
	public int doWork(final List<String> args) {
		try
			{
			parseTrackDB(new URL(trackDBUrl));
			
			return doVcfToVcf(args,outputFile);
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		
		}

	public static void main(String[] args)
		{
		new VcfEnsemblReg().instanceMainWithExit(args);
		}
}
