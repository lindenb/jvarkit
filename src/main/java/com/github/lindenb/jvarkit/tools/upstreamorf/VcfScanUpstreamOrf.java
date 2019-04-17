/*

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
The MIT License (MIT)
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/
package com.github.lindenb.jvarkit.tools.upstreamorf;

import java.io.BufferedReader;
import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.AbstractCharSequence;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.Paranoid;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.Algorithms;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.GeneticCode;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;
/**
BEGIN_DOC

## inspiration

part of this code was inspired from: https://github.com/ImperialCardioGenetics/uORFs/blob/master/5primeUTRannotator/five_prime_UTR_annotator.pm

## Examples

### Example 1

```
 wget -q -O - "https://storage.googleapis.com/gnomad-public/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.chr1.vcf.bgz" |\
 bcftools annotate -x "INFO,FILTER" |\
 java -jar /home/lindenb/src/jvarkit-git/dist/vcfscanupstreamorf.jar \
 	-R human_g1k_v37.fasta  --uorf-only  --canonical  

#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	89333	rs1008713359	A	G	283.15	.	UORF_DEL_ATG=transcript:ENST00000466430.1|strand:-|utr-start:89295|utr-end:120932|alt:C|atg-pos:89334|cap-atg:2708|atg-cds:40|atg-frame:atg-out-of-cds-frame|kozak-seq:GAAATGC|kozak-strength:Moderate|stop-frame:not-in-frame-stop|stop-pos:89318|atg-stop:16|pep:.
1	89359	rs1327179626	C	T	3839.47	.	UORF_ADD_ATG=transcript:ENST00000466430.1|strand:-|utr-start:89295|utr-end:120932|alt:A|atg-pos:89359|cap-atg:2683|atg-cds:65|atg-frame:atg-out-of-cds-frame|kozak-seq:TGCATGT|kozak-strength:Weak|stop-frame:in-frame-stop|stop-pos:89356|atg-stop:3|pep:M
1	89391	rs1332733110	T	C	2045.80	.	UORF_DEL_ATG=transcript:ENST00000466430.1|strand:-|utr-start:89295|utr-end:120932|alt:G|atg-pos:89391|cap-atg:2651|atg-cds:97|atg-frame:atg-out-of-cds-frame|kozak-seq:TGAATGA|kozak-strength:Weak|stop-frame:in-frame-stop|stop-pos:89382|atg-stop:9|pep:VNK
1	89555	rs1200434471	A	C	283.62	.	UORF_ADD_ATG=transcript:ENST00000466430.1|strand:-|utr-start:89295|utr-end:120932|alt:G|atg-pos:89557|cap-atg:2485|atg-cds:263|atg-frame:atg-out-of-cds-frame|kozak-seq:GAAATGA|kozak-strength:Moderate|stop-frame:in-frame-stop|stop-pos:89452|atg-stop:105|pep:MKSQNVSQKIIYNVCVRKRQYPSNFESLHQKENSK,transcript:ENST00000495576.1|strand:-|utr-start:89551|utr-end:91105|alt:G|atg-pos:89557|cap-atg:1312|atg-cds:7|atg-frame:atg-out-of-cds-frame|kozak-seq:GAAATGA|kozak-strength:Moderate|stop-frame:not-in-frame-stop|stop-pos:89556|atg-stop:1|pep:.
1	89560	rs1234719556	C	A	448.62	.	UORF_DEL_ATG=transcript:ENST00000466430.1|strand:-|utr-start:89295|utr-end:120932|alt:T|atg-pos:89562|cap-atg:2480|atg-cds:268|atg-frame:atg-out-of-cds-frame|kozak-seq:TGAATGA|kozak-strength:Weak|stop-frame:in-frame-stop|stop-pos:89541|atg-stop:21|pep:IKLKVKM,transcript:ENST00000495576.1|strand:-|utr-start:89551|utr-end:91105|alt:T|atg-pos:89562|cap-atg:1307|atg-cds:12|atg-frame:atg-in-cds-frame|kozak-seq:TGAATGA|kozak-strength:Weak|stop-frame:not-in-frame-stop|stop-pos:89555|atg-stop:7|pep:.
1	89624	rs1166058274	T	C	606.05	.	UORF_DEL_ATG=transcript:ENST00000466430.1|strand:-|utr-start:89295|utr-end:120932|alt:G|atg-pos:89624|cap-atg:2418|atg-cds:330|atg-frame:atg-in-cds-frame|kozak-seq:ACAATGA|kozak-strength:Moderate|stop-frame:in-frame-stop|stop-pos:89609|atg-stop:15|pep:VKELF,transcript:ENST00000495576.1|strand:-|utr-start:89551|utr-end:91105|alt:G|atg-pos:89624|cap-atg:1245|atg-cds:74|atg-frame:atg-out-of-cds-frame|kozak-seq:ACAATGA|kozak-strength:Moderate|stop-frame:in-frame-stop|stop-pos:89609|atg-stop:15|pep:VKELF
1	89718	rs865856422	A	G	1466.11	.	UORF_DEL_ATG=transcript:ENST00000466430.1|strand:-|utr-start:89295|utr-end:120932|alt:C|atg-pos:89719|cap-atg:2323|atg-cds:425|atg-frame:atg-out-of-cds-frame|kozak-seq:AAAATGA|kozak-strength:Moderate|stop-frame:in-frame-stop|stop-pos:89710|atg-stop:9|pep:TKL,transcript:ENST00000495576.1|strand:-|utr-start:89551|utr-end:91105|alt:C|atg-pos:89719|cap-atg:1150|atg-cds:169|atg-frame:atg-out-of-cds-frame|kozak-seq:AAAATGA|kozak-strength:Moderate|stop-frame:in-frame-stop|stop-pos:89710|atg-stop:9|pep:TKL
1	89831	rs1209426147	A	G	372.62	.	UORF_DEL_ATG=transcript:ENST00000466430.1|strand:-|utr-start:89295|utr-end:120932|alt:C|atg-pos:89832|cap-atg:2210|atg-cds:538|atg-frame:atg-out-of-cds-frame|kozak-seq:CTTATGT|kozak-strength:Weak|stop-frame:in-frame-stop|stop-pos:89811|atg-stop:21|pep:TFAIYHT,transcript:ENST00000495576.1|strand:-|utr-start:89551|utr-end:91105|alt:C|atg-pos:89832|cap-atg:1037|atg-cds:282|atg-frame:atg-in-cds-frame|kozak-seq:CTTATGT|kozak-strength:Weak|stop-frame:in-frame-stop|stop-pos:89811|atg-stop:21|pep:TFAIYHT
1	89945	rs1376722481	G	C	297.51	.	UORF_ADD_ATG=transcript:ENST00000466430.1|strand:-|utr-start:89295|utr-end:120932|alt:G|atg-pos:89947|cap-atg:2095|atg-cds:653|atg-frame:atg-out-of-cds-frame|kozak-seq:AATATGC|kozak-strength:Moderate|stop-frame:in-frame-stop|stop-pos:89803|atg-stop:144|pep:MPLASVSHLAKPRLRSGKMEAISSWERRQRRWEYYVATYVCNLPYLAL,transcript:ENST00000495576.1|strand:-|utr-start:89551|utr-end:91105|alt:G|atg-pos:89947|cap-atg:922|atg-cds:397|atg-frame:atg-out-of-cds-frame|kozak-seq:AATATGC|kozak-strength:Moderate|stop-frame:in-frame-stop|stop-pos:89803|atg-stop:144|pep:MPLASVSHLAKPRLRSGKMEAISSWERRQRRWEYYVATYVCNLPYLAL
1	90032	rs866094671	C	T	14378.50	.	UORF_ADD_ATG=transcript:ENST00000466430.1|strand:-|utr-start:89295|utr-end:120932|alt:A|atg-pos:90032|cap-atg:2010|atg-cds:738|atg-frame:atg-in-cds-frame|kozak-seq:TTCATGG|kozak-strength:Moderate|stop-frame:in-frame-stop|stop-pos:89951|atg-stop:81|pep:MGQLVSRAARETKPQCTFYSLCAHQTC,transcript:ENST00000495576.1|strand:-|utr-start:89551|utr-end:91105|alt:A|atg-pos:90032|cap-atg:837|atg-cds:482|atg-frame:atg-out-of-cds-frame|kozak-seq:TTCATGG|kozak-strength:Moderate|stop-frame:in-frame-stop|stop-pos:89951|atg-stop:81|pep:MGQLVSRAARETKPQCTFYSLCAHQTC
```

### extract to bed format

```
track name="uORF" description="uORF for http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/wgEncodeGencodeBasicV19.txt.gz"
#chrom	chromStart	chromEnd	name	score	strand	thickStart	thickEnd	itemRgb	blockCount	blockSizes	blockStarts
chr1	34553	36081	ENST00000417324.1.uorf	100	-	35140	35736	0,0,255	1	1528	0
chr1	89294	120932	ENST00000466430.1.uorf	500	-	91254	91491	0,255,0	1	31638	0
chr1	89550	91105	ENST00000495576.1.uorf	100	-	90431	90590	0,0,255	1	1555	0
chr1	139789	140339	ENST00000493797.1.uorf	500	-	139816	140223	0,255,0	1	550	0
chr1	141473	149707	ENST00000484859.1.uorf	1000	-	146708	146978	255,0,0	1	8234	0
chr1	142807	146831	ENST00000490997.1.uorf	1000	-	142988	146482	255,0,0	1	4024	0
chr1	157783	157887	ENST00000410691.1.uorf	0	-	157848	157887	0,0,0	1	104	0
chr1	236111	267253	ENST00000424587.2.uorf	100	-	236759	236918	0,0,255	1	31142	0
chr1	453632	460480	ENST00000450983.1.uorf	1000	-	453980	454166	255,0,0	1	6848	0
chr1	521368	523833	ENST00000417636.1.uorf	500	-	522285	523620	0,255,0	1	2465	0
chr1	529838	532878	ENST00000357876.5.uorf	500	-	530001	532684	0,255,0	1	3040	0
chr1	562756	564390	ENST00000452176.1.uorf	500	-	562878	562995	0,255,0	1	1634	0
chr1	646721	655580	ENST00000414688.1.uorf	500	-	647189	655553	0,255,0	1	8859	0
chr1	677192	685396	ENST00000416385.1.uorf	100	-	682910	683180	0,0,255	1	8204	0
chr1	693612	693716	ENST00000411249.1.uorf	0	-	693689	693716	0,0,0	1	104	0
chr1	694411	700305	ENST00000417659.1.uorf	100	-	700133	700208	0,0,255	1	5894	0
chr1	700236	714006	ENST00000428504.1.uorf	500	-	705034	709660	0,255,0	1	13770	0
chr1	736258	745541	ENST00000447500.1.uorf	500	-	741231	745515	0,255,0	1	9283	0
chr1	745488	753092	ENST00000435300.1.uorf	500	-	752900	753047	0,255,0	1	7604	0
chr1	761585	762902	ENST00000473798.1.uorf	500	-	762082	762571	0,255,0	1	1317	0
chr1	803450	812283	ENST00000446136.1.uorf	100	-	810390	812268	0,0,255	1	8833	0
chr1	852249	855072	ENST00000417705.1.uorf	100	-	852976	854794	0,0,255	1	2823	0
chr1	889805	894689	ENST00000487214.1.uorf	1000	-	889839	894620	255,0,0	1	4884	0
chr1	916546	917473	ENST00000341290.2.uorf	1000	-	916549	917473	255,0,0	1	927	0
chr1	931345	933431	ENST00000606034.1.uorf	100	-	931510	932137	0,0,255	1	2086	0
chr1	935353	935552	ENST00000428771.2.uorf	500	-	935487	935544	0,255,0	1	199	0
chr1	947376	948573	ENST00000458555.1.uorf	100	-	947459	947507	0,0,255	1	1197	0
chr1	997587	998668	ENST00000442292.2.uorf	500	-	997810	998119	0,255,0	1	1081	0
chr1	1019305	1051623	ENST00000482816.1.uorf	1000	-	1019401	1026923	255,0,0	1	32318	0
chr1	1026923	1027554	ENST00000379320.1.uorf	100	-	1027028	1027400	0,0,255	1	631	0
chr1	1026923	1041507	ENST00000379319.1.uorf	100	-	1041338	1041410	0,0,255	1	14584	0
(...)
```

END_DOC

*/
@Program(name="vcfscanupstreamorf",
description="Scan BAM for upstream-ORF. Inspired from https://github.com/ImperialCardioGenetics/uORFs ",
keywords={"vcf","uorf"},
creationDate="2019-02-18",
modificationDate="2019-02-25"
)
public class VcfScanUpstreamOrf extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfScanUpstreamOrf.class).make();
	private static final Paranoid PARANOID = Paranoid.createThrowingInstance();
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT+
			". If there is no argument and 'output' is a directory or it ends with '.zip' an archive containing the fasta+bed of the uORF is created and the program exits.")
	private File outputFile = null;
	@Parameter(names={"-r","-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private File faidx = null;
	@Parameter(names={"-k","-K","--kg","-kg"},description=KnownGene.OPT_KNOWNGENE_DESC)
	private String knownGeneUri = KnownGene.getDefaultUri();
	@Parameter(names={"--canonical"},description="reduce the number of transcripts. Keep one if some share the same UTR")
	private boolean canonical_utr = false;
	@Parameter(names={"--exclude-cds"},description="remove a uORF it if enterely overlaps a coding region of the exon of an alternative transcript.")
	private boolean exclude_cds_overlaping_alternative = false;
	@Parameter(names={"--uorf-only"},description="only print variants having something to say about an uorf")
	private boolean print_uorf_only = false;
	@Parameter(names={"--plus"},description="only transcripts on '+' strand. For debugging only",hidden=true)
	private boolean plus_strand_only = false;
	@Parameter(names={"--dac"},description="disable scan for ATG creation")
	private boolean disable_atg_create = false;
	@Parameter(names={"--dad"},description="disable scan for ATG deletion")
	private boolean disable_atg_delete = false;
	@Parameter(names={"--dsc"},description="disable scan for STOP creation")
	private boolean disable_stop_create = false;
	@Parameter(names={"--dsd"},description="disable scan for STOP deletion")
	private boolean disable_stop_delete = false;
	@Parameter(names={"--kal"},description="disable scan for Kozak change")
	private boolean disable_kozak_alteration= false;
	@Parameter(names={"--strong"},description="only accept events that are related to 'Strong' Kozak pattern.")
	private boolean only_strong_kozak= false;

	
	
	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;
	private ContigNameConverter refCtgNameConverter =null;
	private GenomicSequence genomicSequence=null;
	private final IntervalTreeMap<List<KnownGene>> knownGeneMap = new IntervalTreeMap<>();
	private final GeneticCode geneticCode = GeneticCode.getStandard();
	

	private final VCFInfoHeaderLine infoAddATG = new VCFInfoHeaderLine(
			"UORF_ADD_ATG",
			VCFHeaderLineCount.UNBOUNDED,
			VCFHeaderLineType.String,
			"Information about de-novo ATG creation"
			);
	private final VCFInfoHeaderLine infoDelATG = new VCFInfoHeaderLine(
			"UORF_DEL_ATG",
			VCFHeaderLineCount.UNBOUNDED,
			VCFHeaderLineType.String,
			"Information about ATG deletion"
			);
	private final VCFInfoHeaderLine infoAddSTOP = new VCFInfoHeaderLine(
			"UORF_ADD_STOP",
			VCFHeaderLineCount.UNBOUNDED,
			VCFHeaderLineType.String,
			"Information about de-novo STOP creation"
			);
	private final VCFInfoHeaderLine infoDelSTOP = new VCFInfoHeaderLine(
			"UORF_DEL_STOP",
			VCFHeaderLineCount.UNBOUNDED,
			VCFHeaderLineType.String,
			"Information about STOP deletion"
			);
	
	private final VCFInfoHeaderLine infoKozakAlteration = new VCFInfoHeaderLine(
			"UORF_KOZAK_ALTERATION",
			VCFHeaderLineCount.UNBOUNDED,
			VCFHeaderLineType.String,
			"Information about kozak sequences alteration"
			);

	
	enum KozakStrength {
		Strong,Moderate,Weak,nil;
		
		String getColor() {
			switch(this)
			{
			case Strong: return "255,0,0";
			case Moderate: return "0,255,0";
			case Weak: return "0,0,255";
			default:  return "0,0,0";
			}
		}
		
	}
	
	/* kozak consensus length */
	private static final int KOZAK_LENGTH=7;
	/* position of ATG in the kozak sequence */
	private static final int KOZAK_ATG=3;
	
	/* a Kozak consensus */
	private class KozakSequence
		extends AbstractCharSequence
		{
		final CharSequence delegate;
		final int offset_beg;
		KozakSequence(final  CharSequence delegate,final int offset_atg) {
			this.delegate = delegate;
			this.offset_beg = offset_atg-KOZAK_ATG;
			}
		
		@Override
		public char charAt(int index) {
			final int idx= this.offset_beg + index;
			return idx >=0 && idx <delegate.length()?delegate.charAt(idx):'N';
			}
		
		private boolean hasATG() {
			if(charAt(KOZAK_ATG  )!='A') return false;
			if(charAt(KOZAK_ATG+1)!='T') return false;
			if(charAt(KOZAK_ATG+2)!='G') return false;
			return true;
			}

		
		KozakStrength getStrength() {
			char c0 = charAt(0);
			if(c0=='N') return KozakStrength.nil;
			char c6 = charAt(6);
			if(c6=='N') return KozakStrength.nil;
			if(!hasATG()) return KozakStrength.nil;
			
			
			boolean ok_chart_at_0 = c0 == 'A' || c0 == 'G' ;
			
  			if ( ok_chart_at_0 && c6 == 'G'){
  				return  KozakStrength.Strong;
  				}
  			else if (ok_chart_at_0 || c6 == 'G'){
  				return  KozakStrength.Moderate;
  				}
  			return KozakStrength.Weak;
			}
		
		@Override
		public int length() {
			return KOZAK_LENGTH;
			}
		
		void annotate(final Map<String,String> sb) {
			sb.put("kozak-seq",this.toString());
			sb.put("kozak-strength",this.getStrength().name());
			}
		
		}
	private final static char NO_BASE='\0';
	
	private interface Locatable0
		{
		/** associated gene */
		public KnownGene getGene();
		/** chrom */
		public String getContig();
		/** 0 based chrom start */
		public int getChromStart();
		/** 0 based past chrom end */
		public int getChromEnd();
		
		}
	
	/** wrapper about the information about a RNA sequence, it can be
	 * the RNA itself or the uORF
	 */
	private abstract class AbstractRNASequence 
		extends AbstractCharSequence
		implements Locatable0
		{
		protected AbstractRNASequence() {
			}
		
		protected void throwIfNotInRange(final int idx0) {
			if(idx0 >=0 && idx0 < this.length()) return;
			throw new IndexOutOfBoundsException("Bad index" +idx0+" in RNA length="+this.length());
		}
		
		/** return the associated UCSC gene */
		public abstract KnownGene getGene();
	
		/** return true if the is an AT at this position */
		boolean isATG(final int index) {
			return index>=0 && 
					index+2 < this.length() &&
					charAt(index+0)=='A' &&
					charAt(index+1)=='T' && 
					charAt(index+2)=='G'
					;
			}
		/** return true if translation at index is stop */
		boolean isStop(final int index) {
			return index>=0 && 
					index+2 < this.length() && 
					VcfScanUpstreamOrf.this.geneticCode.isStopCodon(
							charAt(index),
							charAt(index+1),
							charAt(index+2)
							)
					;
			}
		
		/** find best ATG upstream of stop 'best' is first in frame with atg */
		int findBestATG(final int stop_pos) {
			PARANOID.assertTrue(isStop(stop_pos));
			int best=-1;
			for(int idx=stop_pos-3;idx>=0;idx-=3) {
				if(isStop(idx))
					{
					//there is another stop at this position in frame, not a good candidate
					return -1;
					}
				else if(!isATG(idx)) 
					{
					//nothing
					}
				else if(best==-1) {
					best= idx;
					}
				else  {
					final KozakSequence kozak1= getKozakContext(best);
					final KozakSequence kozak2= getKozakContext(idx);
					if(kozak2.getStrength().compareTo(kozak1.getStrength())<0) {
						best= idx;
						}
					}
				}
			return best;
			}
		
		/** find best stop starting from ATG. 'best' is first in frame with atg */
		int findBestStop(final int ATG_pos) {
			if(isStop(ATG_pos)) throw new IllegalStateException("cannot be a stop !");
			int stop=-1;
			for(int idx=ATG_pos+3; idx+2 <  this.length();idx+=3) {
				if(isStop(idx)) return stop;
				}
			return -1;
			}
		/** get the strand of the gene */
		protected final boolean isPositiveStrand() {
			return this.getGene().isPositiveStrand();
			}
		/** get the contig of the gene */
		@Override
		public final String getContig() {
			return getGene().getContig();
			}
		/** return kozak context at ATG position . Return null if cannot build context **/
		KozakSequence getKozakContext(int atg0) {
			return new KozakSequence(this,atg0);
			}
		
		/** return peptide between pos_atg and post_stop or '.' if they're not in frame */
		String translate(final int pos_ATG,final int pos_stop) {
			if(pos_ATG>=pos_stop) throw new IllegalStateException("atg"+pos_ATG+">stop="+pos_stop);
			if(pos_ATG%3 != pos_stop%3) return ".";
			final StringBuilder sb=new StringBuilder();
			for(int i=pos_ATG;i<pos_stop;i+=3) {
				sb.append(VcfScanUpstreamOrf.this.geneticCode.translate(
						charAt(i+0),
						charAt(i+1),
						charAt(i+2))
						);
				}
			return sb.toString();
			}
		/** convert 0-based position of orf to 0-based position of genome */
		protected abstract int convertOrfToGenomicIndex0(final int index);
		
		/** convert genomic index to mRNA index */
		protected abstract int convertGenomicIndex0ToOrf(int genomic0);
		
		/** return true if genomic position is transcripted in this RNA */
		abstract boolean containsGenomicPos0(int genomic0);
		}
	
	
	/** wraps a UCSC gene into a CharSequence */
	private abstract class AbstractKnownGeneRNA 
		extends AbstractRNASequence
		{
		protected AbstractKnownGeneRNA() {
			}

		@Override
		public abstract KnownGene getGene();
		
		@Override
		public final int getChromStart() {
			return getGene().getTxStart();
			}
		@Override
		public final int getChromEnd() {
			return getGene().getTxEnd();
			}
		}
	
	/** wraps a UCSC gene into a CharSequence */
	private class KnownGeneRNA 
		extends AbstractKnownGeneRNA
		{
		private final KnownGene knownGene; 
		/** maps RNA index to genomic index in ascending order */
		private final int genomic_indexes0[];
		/** maps RNA index to base, in RNA order */
		private final char base_buffer[];
		
		KnownGeneRNA(final KnownGene knownGene) {
			this.knownGene = knownGene;
			final List<Integer> L=new ArrayList<>(3_000);
			for(int exon_index=0;exon_index< knownGene.getExonCount();++exon_index) {
				int beg = knownGene.getExonStart(exon_index);
				final int end = knownGene.getExonEnd(exon_index);
				while(beg<end)
					{
					L.add(beg);
					++beg;
					}
				}
				
			this.genomic_indexes0 = new int[L.size()];
			for(int x=0;x< L.size();++x) this.genomic_indexes0[x] = L.get(x);
			this.base_buffer = new char[this.genomic_indexes0.length];
			Arrays.fill(this.base_buffer, NO_BASE);
			}
		
		@Override
		public KnownGene getGene() {
			return this.knownGene;
			}
		
		
		@Override
		public int length() {
			return this.base_buffer.length;
			}
		
		@Override
		boolean containsGenomicPos0(int genomic0) {
			final int idx =  Arrays.binarySearch(this.genomic_indexes0,genomic0);
			return idx>= 0 &&  idx < this.length();
			}
		
		@Override
		protected int convertGenomicIndex0ToOrf(final int g0) {
			final int x= Arrays.binarySearch(this.genomic_indexes0, g0);
			throwIfNotInRange(x);
			PARANOID.assertEq(this.genomic_indexes0[x],g0);
			return isPositiveStrand() ? x: (this.length()-1)-x;
			}
		
		@Override
		protected int convertOrfToGenomicIndex0(final int index) {
			throwIfNotInRange(index);
			return this.genomic_indexes0[
			        isPositiveStrand()?
					index:
					(length()-1)-index
					];
			}
		
		@Override
		public char charAt(final int index_in_rna) {
			throwIfNotInRange(index_in_rna);
			if(this.base_buffer[index_in_rna] != NO_BASE) {
				return this.base_buffer[index_in_rna];
				}
			final int g0  = this.convertOrfToGenomicIndex0(index_in_rna);
			char c = Character.toUpperCase(VcfScanUpstreamOrf.this.genomicSequence.charAt(g0));
			if(getGene().isNegativeStrand()) {
				c=AcidNucleics.complement(c);
				}
			this.base_buffer[index_in_rna]=c;
			
			return c;
			}	
		}
	
	private class MutatedKnownGeneRNA extends AbstractKnownGeneRNA
		{
		private final MutatedUTR delegate;
		MutatedKnownGeneRNA(final MutatedUTR delegate)  {
			this.delegate = delegate;
			}

		@Override
		public KnownGene getGene() {
			return getWildRNA().getGene();
			}
		private KnownGeneRNA getWildRNA() {
			return delegate.delegate.getRNA();
		}
		
		@Override
		public char charAt(int index) {
			throwIfNotInRange(index);
			if(index==delegate.alt_position0) return delegate.alt_base;
			return getWildRNA().charAt(index);
			}
		@Override
		boolean containsGenomicPos0(int genomic0) {
			return getWildRNA().containsGenomicPos0(genomic0);
			}
		@Override
		protected int convertGenomicIndex0ToOrf(int genomic0) {
			return getWildRNA().convertGenomicIndex0ToOrf(genomic0);
			}
		@Override
		protected int convertOrfToGenomicIndex0(int index) {
			return getWildRNA().convertOrfToGenomicIndex0(index);
			}
		@Override
		public int length() {
			return getWildRNA().length();
			}
		}
	
	/** base class for UORF or mutated-UROF */
	private abstract class AbstractUTRSequence extends AbstractRNASequence
		{
		protected abstract AbstractKnownGeneRNA getRNA();
		
		@Override
		public final KnownGene getGene() {
			return getRNA().getGene();
			}
		
		public final int getChromStart() {
			return this.isPositiveStrand()?
					this.getGene().getTxStart():
					this.getGene().getCdsEnd();
			}
		public final int getChromEnd() {
			return this.isPositiveStrand()?
					this.getGene().getCdsStart():
					this.getGene().getTxEnd();
			}
		
		
		@Override
		final boolean containsGenomicPos0(final int genomic0) {
			if(genomic0<this.getChromStart()) return false;
			if(genomic0>=this.getChromEnd()) return false;
			PARANOID.assertTrue(getRNA().containsGenomicPos0(genomic0));
			final int idx =  getRNA().convertGenomicIndex0ToOrf(genomic0);
			return idx>= 0 &&  idx < this.length();
			}
	
		@Override
		protected final int convertOrfToGenomicIndex0(final int idx) {
			throwIfNotInRange(idx);
			// uORF is always in 5' of the RNA
			final int g0 = this.getRNA().convertOrfToGenomicIndex0(idx);
			PARANOID.assertTrue(containsGenomicPos0(g0));
			return g0;
			}
		
		@Override
		protected final int convertGenomicIndex0ToOrf(int genomic0) {
			final int idx = getRNA().convertGenomicIndex0ToOrf(genomic0);
			throwIfNotInRange(idx);
			return idx;
			}
		
		/** add annotation about stop */
		protected final void annotStop(final Map<String,String> sb,final int pos_ATG,final int pos_stop) {
			PARANOID.assertGe(pos_ATG, 0);
			if(pos_stop<0) return;
			PARANOID.assertEq(pos_ATG%3,pos_stop%3);
			
			sb.put("stop-pos",String.valueOf(1+this.getRNA().convertOrfToGenomicIndex0(pos_stop)));
			sb.put("atg~stop",String.valueOf(pos_stop-pos_ATG));
			sb.put("stop",pos_stop>=this.length()?"in-cds":"in-uorf");
			sb.put("pep",this.getRNA().translate(pos_ATG,pos_stop));
				
			}
		
		
		/** find best OenReading frame in the region */
		public OpenReadingFrame getBestOpenReadingFrame() {
			OpenReadingFrame best=null;
			for(int i=0;i +2< this.length();i++) {
				if(!isATG(i)) continue;
				int stop=-1;
				for(int j=i+3;j+2<this.length();j+=3) {
					if(isStop(j)) {
						stop=j;
						break;
					}
				}
				if(stop==-1) continue;
				final OpenReadingFrame orf = new OpenReadingFrame(this,i,stop);
				if(best==null || orf.isBetterThan(best)) {
					best=orf;
					}
				}
			return best;
			}
		}
	
	
	/** an ORF in an UpstreamORF, a subsectiobn of an mRNA */
	private class OpenReadingFrame extends AbstractRNASequence 
		{
		private final AbstractUTRSequence uorf;
		private final int atg0;
		private final int stop0;
		//private final int chromStart0;
		OpenReadingFrame(final AbstractUTRSequence uorf,final int atg0,final int stop0) {
			this.uorf = uorf;
			this.atg0 = atg0;
			this.stop0 = stop0;
			PARANOID.assertLt(atg0,stop0).
				assertTrue((stop0-atg0)%3==0).
				assertTrue(uorf.isATG(atg0)).
				assertTrue(uorf.isStop(stop0));
			}
		public AbstractUTRSequence getUORF() {
			return this.uorf;
			}
		
		public AbstractKnownGeneRNA getRNA() {
			return this.getUORF().getRNA();
			}
		@Override
		public final KnownGene getGene() {
			return this.getRNA().getGene();
			}
		
		@Override
		public final int getChromStart() {
			return this.getRNA().getChromStart();
			}
		
		@Override
		public final int getChromEnd() {
			return this.getRNA().getChromEnd();
			}
		
		@Override
		protected final int convertGenomicIndex0ToOrf(int genomic0) {
			return this.getRNA().convertGenomicIndex0ToOrf(genomic0);
			}
		
		@Override
		boolean containsGenomicPos0(int genomic0) {
			return this.getRNA().containsGenomicPos0(genomic0);
			}
		@Override
		protected int convertOrfToGenomicIndex0(int index) {
			return this.getRNA().convertOrfToGenomicIndex0(index);
			}
		
		@Override
		public int length() {
			return this.stop0 - this.atg0;
			}
		
		@Override
		public char charAt(final int index) {
			if(index>=this.length()) return 'N';// if stop is beyond uORF, in whole mRNA
			throwIfNotInRange(index);
			return this.getRNA().charAt(this.atg0+index);
			}
		/** get kozak context for the ATG of this sequence */
		KozakSequence getKozakSequence() {
			return this.getKozakContext(this.atg0);
			}
		KozakStrength getKozakStrength() {
			return getKozakSequence().getStrength();
			}
		public boolean isBetterThan(final OpenReadingFrame o) {
			int i= Integer.compare(this.length(), o.length());
			if(i!=0) return i>0;
			i= getKozakStrength().compareTo(o.getKozakStrength());
			if(i!=0) return i<0;
			return false;
			}
	
		// private int[] getCdsStartEnd0() {}
		
		
		private String getLabel() {
			return new StringBuilder(this.getGene().getName()).
					append("|strand:").
					append(this.isPositiveStrand()?"+":"-").
					toString();
					
			}

		String translate() {
			final StringBuilder sb= new StringBuilder();
			for(int x=0;x+2< this.length();x+=3) {
				sb.append(geneticCode.translate(charAt(x), charAt(x+1), charAt(x+2)));
				}
			return sb.toString();
			}
		
		
		void printFastaDNA(final PrintWriter pw) {
			pw.println(">"+getLabel()+"|pep:"+translate());
			pw.println(this.toString());
			}
		void printFastaPep(final PrintWriter pw) {
			pw.println(">"+getLabel());
			pw.println(this.translate());
			}
		void printBed(final PrintWriter pw) {
			pw.print(this.getContig()); pw.print('\t');//chrom
			pw.print(this.getChromStart()); pw.print('\t');//chromStart of uORF
			pw.print(this.getChromEnd()); pw.print('\t');// chromEnd of uORF
			pw.print(this.getGene().getName()+".uorf"); pw.print('\t');//name
			switch(getKozakStrength())
				{
				case Strong: pw.print(1000);break;
				case Moderate: pw.print(500);break;
				case Weak: pw.print(100);break;
				default: pw.print(0);break;
				}
			pw.print('\t');//score
			
			final char strand;
			final int thickStart;
			final int thickEnd;
			
			if(this.isPositiveStrand()) {
				strand = '+';
				thickStart= this.convertOrfToGenomicIndex0(this.atg0);
				thickEnd = this.convertOrfToGenomicIndex0(this.stop0);
				PARANOID.assertLt(thickStart,thickEnd);
				PARANOID.assertTrue(thickStart< getGene().getCdsStart());
				}
			else
				{
				strand='-';
				thickStart = this.convertOrfToGenomicIndex0(this.stop0)+1;
				thickEnd = this.convertOrfToGenomicIndex0(this.atg0) +1;
				PARANOID.assertLt(this.atg0,this.stop0);
				PARANOID.assertLt(thickStart,thickEnd);
				PARANOID.assertLe(getGene().getCdsEnd(),thickEnd);
				}
			PARANOID.assertTrue(thickStart<thickEnd);
			pw.print(strand);//strand
			pw.print('\t');
			pw.print(thickStart);//thickStart
			pw.print('\t');
			pw.print(thickEnd);//thickStart
		
			pw.print('\t');
			pw.print(this.getKozakStrength().getColor()); pw.print('\t');//itemRgb
			pw.print('\t');// score 
			pw.print(getGene().getExonCount());pw.print("\t");// blockCount
			for(int x=0;x< getGene().getExonCount();++x) {
				if(x>0) pw.print(",");
				 pw.print(getGene().getExonEnd(x)-getGene().getExonStart(x));
				}
			pw.print('\t');// blockSizes
			for(int x=0;x< getGene().getExonCount();++x) {
				if(x>0) pw.print(",");
				pw.print(getGene().getExonStart(x)-this.getChromStart());
				} // blockStarts
			pw.println();
			}
		}
	
	/** the 5' UTR of an mRNA */
	private class UpstreamORF
		extends AbstractUTRSequence
		{
		final KnownGeneRNA mRNA; 
		final int _length;
		UpstreamORF(final KnownGeneRNA mRNA) {
			this.mRNA=mRNA;
			if(mRNA.isPositiveStrand()) {
				final int cdsStart = mRNA.getGene().getCdsStart();
				this._length =  Algorithms.lower_bound(mRNA.genomic_indexes0, cdsStart);
				}
			else
				{
				final int cdsEnd = mRNA.getGene().getCdsEnd();
				final int x1 =  mRNA.length();
				final int x0 = Algorithms.lower_bound(mRNA.genomic_indexes0, cdsEnd);
				this._length = x1-x0;
				
				final int z[]=new int[]{0,0,0,1,1,2,2,3};
				PARANOID.assertEq(Algorithms.lower_bound(z, 0),0);
				PARANOID.assertEq(Algorithms.upper_bound(z, 0),3);
				//PARANOID.assertEq(Algorithms.upper_bound(z, 0),2);
				
				
				int x3=this.mRNA.genomic_indexes0.length;
				int x4=x3;
				while(x3-1>=0) {
					if(mRNA.genomic_indexes0[x3-1]>=cdsEnd) {
						x4=x3-1;
						}
					x3--;
					}
				PARANOID.assertEq(x4,x0);
				}
			}		
		
		
		@Override
		protected final KnownGeneRNA getRNA() {
			return this.mRNA;
			}
		
		@Override
		public final int length() {
			return this._length;
			}
		

		@Override
		public final char charAt(final int index) {
			if( index>=0 &&   index< this.length())
				{
				return mRNA.charAt(index);
				}
			throw new IndexOutOfBoundsException(String.valueOf(index));
			}
		
		
		}
	
	/** mutated version of UpstreamORF */
	private class MutatedUTR extends AbstractUTRSequence
		{
		final MutatedKnownGeneRNA mutatedRNA;
		final UpstreamORF delegate;
		final int genomic_index0_of_variant;
		final char alt_base;
		final int alt_position0;
		final Set<String> denovo_atg_set = new LinkedHashSet<>();
		final Set<String> remove_atg_set = new LinkedHashSet<>();
		final Set<String> remove_stop_set = new LinkedHashSet<>();
		final Set<String> denovo_stop_set = new LinkedHashSet<>();
		final Set<String> kozak_alterations_set = new LinkedHashSet<>();
		
		MutatedUTR(final UpstreamORF delegate,final VariantContext ctx,final int alt_index) {
			this.delegate=delegate;
			this.genomic_index0_of_variant = ctx.getStart()-1;
			final char c=Character.toUpperCase((char)ctx.getAlleles().get(alt_index).getBases()[0]);
			this.alt_base= delegate.isPositiveStrand() ? c: AcidNucleics.complement(c);
			this.alt_position0 = delegate.convertGenomicIndex0ToOrf(this.genomic_index0_of_variant);
			if(delegate.convertOrfToGenomicIndex0(this.alt_position0)!=this.genomic_index0_of_variant) throw new IllegalStateException();
			this.mutatedRNA = new MutatedKnownGeneRNA(this);
			}
		
		@Override
		protected AbstractKnownGeneRNA getRNA() {
			return this.mutatedRNA;
			}
		
		@Override
		public char charAt(final int index) {
			return  index==this.alt_position0?
					this.alt_base:
					this.delegate.charAt(index);
			}
		
		@Override
		public int length() {
			return delegate.length();
			}
		
		boolean isChanging() {
			return !this.denovo_atg_set.isEmpty() ||
				!this.denovo_stop_set.isEmpty() ||
				!this.remove_atg_set.isEmpty() ||
				!this.remove_stop_set.isEmpty() ||
				!this.kozak_alterations_set.isEmpty()
				;
			}
		
		/* return atg shift position if alt base creates a new ATG: 0,1,2 or -1 if there is no new ATG */
		private Map<String,String> getAnnotPrefix(final int pos_ATG) {
			final int dist_from_cap = pos_ATG;
			final int dist_to_cds_start = this.length() - pos_ATG ;
			
			final Map<String,String> hash = new LinkedHashMap<>();
			hash.put("transcript", this.getGene().getName());
			hash.put("strand",this.isPositiveStrand()?"+":"-");
			hash.put("utr-start",String.valueOf(this.getChromStart()+1));
			hash.put("utr-end",String.valueOf(this.getChromEnd()));
			hash.put("alt",String.valueOf(this.alt_base));
			hash.put("alt-pos",String.valueOf(1+this.delegate.convertOrfToGenomicIndex0(pos_ATG)));
			hash.put("cap~atg",String.valueOf(dist_from_cap));
			hash.put("atg~cds",String.valueOf(dist_to_cds_start));
			hash.put("atg-frame", dist_to_cds_start % 3 !=0 ? "atg-out-of-cds-frame" :  "atg-in-cds-frame");
			return hash;
			}

		
		void invoke() {
			if(!disable_atg_delete) removeATG();
			if(!disable_stop_delete) removeStop();
			if(!disable_atg_create) findDeNovoATG();
			if(!disable_stop_create) deNovoStop();
			if(!disable_kozak_alteration) kozakAlteration();
			}
		private String toString(final Map<String,String> map) {
			return map.entrySet().stream().map(K->K.getKey()+":"+K.getValue()).collect(Collectors.joining("|"));
		}
		
		void removeStop() {
			for(int i=0;i<3;i++) {
				// stop removed
				if(!this.isStop(this.alt_position0 - i /* yes minus */)) continue;
				if(this.delegate.isStop(this.alt_position0 - i /* yes minus */)) continue;
				
				final int pos_stop = this.alt_position0 - i;
				
				final int pos_ATG = this.findBestATG(pos_stop);
				if(pos_ATG==-1) continue;
				PARANOID.assertLt(pos_ATG,pos_stop);
				
				final KozakSequence kozakSequence = this.getKozakContext(pos_ATG);
				if(!acceptKozak(kozakSequence)) {
					continue;
					}
				final Map<String,String> sb =getAnnotPrefix(pos_ATG);
				kozakSequence.annotate(sb);
				this.delegate.annotStop(sb, pos_ATG, pos_stop);
				
				this.remove_stop_set.add(toString(sb));
				}
			}
		/** find de NOVO stop added to uORF */
		void deNovoStop() {
			for(int i=0;i<3;i++) {
				// stop added
				if(!this.isStop(this.alt_position0 - i /* yes minus */)) continue;
			    if(this.delegate.isStop(this.alt_position0 - i /* yes minus */)) continue;
			    
				final int pos_stop = this.alt_position0 - i;
				
				
				final int pos_ATG = this.findBestATG(pos_stop);
				if(pos_ATG==-1) continue;
				PARANOID.assertLt(pos_ATG,pos_stop);
				
				final KozakSequence kozakSequence = this.getKozakContext(pos_ATG);
				if(!acceptKozak(kozakSequence)) continue;
				
				final Map<String,String> sb =getAnnotPrefix(pos_ATG);
				
				kozakSequence.annotate(sb);
				
				this.annotStop(sb, pos_ATG, pos_stop);
				
				this.denovo_stop_set.add(toString(sb));
				}
			}

		/** ATG in uORF removed */
		void removeATG() {
			for(int i=0;i<3;i++) {
				if(!this.delegate.isATG(this.alt_position0 - i /* yes minus */)) continue;
				if(this.isATG(this.alt_position0 - i)) continue;
				
				
				final int pos_ATG = this.alt_position0 - i;
				final KozakSequence kozakSequence = this.delegate.getKozakContext(pos_ATG);
				if(!acceptKozak(kozakSequence)) continue;
				
				// yes we use mRNA.findBest because stop can be beyond CDS
				final int pos_stop = this.delegate.getRNA().findBestStop(pos_ATG);
				if(pos_stop==-1) continue;
				PARANOID.assertLt(pos_ATG,pos_stop);
				
				final Map<String,String> sb = getAnnotPrefix(pos_ATG);
				kozakSequence.annotate(sb);
				annotStop(sb, pos_ATG, pos_stop);
				
				this.remove_atg_set.add(toString(sb));
				}
			}
		
		/** ATG in uORF added */
		void findDeNovoATG() {
			for(int i=0;i<3;i++) {
				//must be atg denovo
				if(!this.isATG(this.alt_position0 - i /* yes minus */)) continue;					
				if(this.delegate.isATG(this.alt_position0 - i /* yes minus */)) continue;					
				final int pos_A = this.alt_position0 - i;
				
				final KozakSequence kozakSequence = this.getKozakContext(pos_A);
				if(!acceptKozak(kozakSequence)) continue;
				// yes we use mRNA.findBest because stop can be beyond CDS
				final int stop = this.getRNA().findBestStop(pos_A);
				if(stop==-1) continue;
				
				PARANOID.assertLt(pos_A,stop);
				
				final Map<String,String> sb =getAnnotPrefix(pos_A);
				kozakSequence.annotate(sb);
				this.annotStop(sb,pos_A,stop);
				this.denovo_atg_set.add(toString(sb));
				}
			}
		
		void kozakAlteration() {
			for(int i=0;i< KOZAK_LENGTH;i++) {
				final int kozak_start = this.alt_position0 - i;
				if(kozak_start<0 || kozak_start+KOZAK_LENGTH>=this.length()) continue;
				final int pos_A = kozak_start+ KOZAK_ATG;
				if(!isATG(pos_A)) continue;
				if(!this.delegate.isATG(pos_A)) continue;
				final KozakSequence kozak1 = this.delegate.getKozakContext(pos_A);
				if(kozak1.getStrength()==KozakStrength.nil) continue;
				
				final KozakSequence kozak2 = this.getKozakContext(pos_A);
				if(kozak2.getStrength()==KozakStrength.nil) continue;
				
				if(!acceptKozak(kozak1) && !acceptKozak(kozak2)) {
					continue;
					}
				
				final int stop = this.getRNA().findBestStop(pos_A);
				if(stop<0) continue;

				
				final Map<String,String> sb = getAnnotPrefix(pos_A);
				
				sb.put("kozak-ref-seq",kozak1.toString());
				sb.put("kozak-ref-strength",kozak1.getStrength().name());
				sb.put("kozak-alt-seq",kozak2.toString());
				sb.put("kozak-alt-strength",kozak2.getStrength().name());
				
				annotStop(sb, pos_A, stop);
				
				this.kozak_alterations_set.add(toString(sb));
				}
			}
		}
	
	private boolean acceptKozak(final KozakSequence k) {
		if(this.only_strong_kozak && !KozakStrength.Strong.equals(k.getStrength())) return false;
		return true;
	}
	
	/** rconvert transcript to interval return null if there is no UTR */
	private Interval kgToUTRInterval(final KnownGene kg) {
		if(kg.isNonCoding()) {
			throw new IllegalStateException();
			}
		else if(kg.isPositiveStrand()) {
			if(kg.getTxStart()>=kg.getCdsStart()) return null;
			return new Interval(kg.getContig(),kg.getTxStart()+1,kg.getCdsStart(),false,kg.getName());
			}
		else  if(kg.isNegativeStrand())
			{
			if(kg.getCdsEnd()>=kg.getEnd()) return null;
			return new Interval(kg.getContig(),kg.getCdsEnd()/* no +1 */,kg.getTxEnd(),true,kg.getName());
			}
		else
			{
			throw new IllegalStateException("bad strand");
			}
		}
	
	/** return true of all bases in uORF candidate overlap a CDS in kg */
	private boolean overlapCDS(final KnownGene candidate,final KnownGene kg) {
		final KnownGeneRNA mrna1 = new KnownGeneRNA(candidate);
		final UpstreamORF uorf1 = new UpstreamORF(mrna1);
		final KnownGeneRNA mrna2 = new KnownGeneRNA(kg);
		
		for(int i=0;i< uorf1.length();i++) {
			final int g1 = uorf1.convertOrfToGenomicIndex0(i);
			if(g1 < kg.getCdsStart()) return false;
			if(g1 >= kg.getCdsEnd()) return false;
			if(!mrna2.containsGenomicPos0(g1)) return false;
			}
		return true;
		}
	
	
	@Override
	protected int doVcfToVcf(
		final String inputName,
		final VCFIterator iter,
		final VariantContextWriter out
		) {
		try {
			/** build vcf header */						
			final VCFHeader header= iter.getHeader();
			header.addMetaDataLine(this.infoAddATG);
			header.addMetaDataLine(this.infoAddSTOP);
			header.addMetaDataLine(this.infoDelATG);
			header.addMetaDataLine(this.infoDelSTOP);
			header.addMetaDataLine(this.infoKozakAlteration);
			JVarkitVersion.getInstance().addMetaData(this, header);
			final ProgressFactory.Watcher<VariantContext> progress = ProgressFactory.newInstance().dictionary(header).logger(LOG).build();
			
			out.writeHeader(header);
			
			while(iter.hasNext()) {
				final VariantContext ctx = progress.apply(iter.next());
				if(!ctx.isVariant()) {
					if(!this.print_uorf_only) out.add(ctx);
					continue;
					}
				if(!ctx.isSNP()) {
					if(!this.print_uorf_only) out.add(ctx);
					continue;
					}
				final String refContig = this.refCtgNameConverter.apply(ctx.getContig());
				if(StringUtils.isBlank(refContig)) {
					if(!this.print_uorf_only) out.add(ctx);
					continue;
					}
				
				/* new reference sequence */
				if(this.genomicSequence==null || !this.genomicSequence.getChrom().equals(refContig)) {
					this.genomicSequence = new GenomicSequence(this.indexedFastaSequenceFile, refContig);
					}

				
				final Interval interval = new Interval(refContig,ctx.getStart(),ctx.getEnd()); 
				final List<KnownGene> kgGenes = this.knownGeneMap.getOverlapping(interval).
						stream().
						flatMap(C->C.stream()).
						collect(Collectors.toList());
				
				if(kgGenes.isEmpty()) {
					if(!this.print_uorf_only) out.add(ctx);
					continue;
					}
				
				final List<UpstreamORF> uorfs = kgGenes.
						stream().
						map(KG->new KnownGeneRNA(KG)).
						filter(KG->KG.containsGenomicPos0(ctx.getStart()-1)).
						map(KG->new UpstreamORF(KG)).
						filter(KG->KG.containsGenomicPos0(ctx.getStart()-1)).
						sorted((A,B)->Integer.compare(A.getChromStart(),B.getChromStart())).
						collect(Collectors.toList());
				
				if(uorfs.isEmpty()) {
					if(!this.print_uorf_only) out.add(ctx);
					continue;
					}
								
				final List<MutatedUTR> mutatedList = new ArrayList<>();
				for(final UpstreamORF uorf: uorfs) {
					for(int alt_idx=1 /* 0==REF */;alt_idx< ctx.getAlleles().size();++alt_idx) {
						final Allele alt_allele = ctx.getAlleles().get(alt_idx);
						if(alt_allele.isSymbolic() || !alt_allele.isCalled() || alt_allele.length()!=1) continue;
						if(!AcidNucleics.isATGC(alt_allele.getDisplayBases()[0])) continue;
						final MutatedUTR mutated = new MutatedUTR(uorf, ctx, alt_idx);
						mutated.invoke();
						mutatedList.add(mutated);
						}
					}
				
				if(mutatedList.isEmpty() || mutatedList.stream().noneMatch(M->M.isChanging())) {
					if(!this.print_uorf_only) out.add(ctx);
					continue;
					}
				final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				
				List<String> ann = new ArrayList<>(
							mutatedList.
							stream().
							flatMap(M->M.remove_atg_set.stream()).
							collect(Collectors.toCollection(LinkedHashSet::new)));
				if(!ann.isEmpty()) {
					vcb.attribute(this.infoDelATG.getID(), ann);
					}
				
				ann = new ArrayList<>(
						mutatedList.
						stream().
						flatMap(M->M.remove_stop_set.stream()).
						collect(Collectors.toCollection(LinkedHashSet::new)));
				if(!ann.isEmpty()) {
					vcb.attribute(this.infoDelSTOP.getID(), ann);
					}
				
				ann = new ArrayList<>(
						mutatedList.
						stream().
						flatMap(M->M.denovo_atg_set.stream()).
						collect(Collectors.toCollection(LinkedHashSet::new)));
				if(!ann.isEmpty()) {
					vcb.attribute(this.infoAddATG.getID(), ann);
					}
				
				ann = new ArrayList<>(
						mutatedList.
						stream().
						flatMap(M->M.denovo_stop_set.stream()).
						collect(Collectors.toCollection(LinkedHashSet::new)));
				if(!ann.isEmpty()) {
					vcb.attribute(this.infoAddSTOP.getID(), ann);
					}
				
				ann = new ArrayList<>(
						mutatedList.
						stream().
						flatMap(M->M.kozak_alterations_set.stream()).
						collect(Collectors.toCollection(LinkedHashSet::new)));
				if(!ann.isEmpty()) {
					vcb.attribute(this.infoKozakAlteration.getID(), ann);
					}
				
				out.add(vcb.make());
				}
			
			progress.close();
			
			return 0;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			}
		}
	
	private int saveArchive(final File outFile,final SAMSequenceDictionary dict) {
		try {
			final ArchiveFactory archive=ArchiveFactory.open(outFile);
			final ContigDictComparator ctgCmp=  new ContigDictComparator(dict);
			final Comparator<Locatable> locCmp1=(A,B)->{
				int i= ctgCmp.compare(A.getContig(), B.getContig());
				if(i!=0) return i;
				i =  Integer.compare(A.getStart(), B.getStart());
				if(i!=0) return i;
				return  Integer.compare(A.getEnd(), B.getEnd());
				};
			final Comparator<KnownGene> locCmp2=(A,B)->{
				return  locCmp1.compare(kgToUTRInterval(A),kgToUTRInterval(B));
				};
			
			final PrintWriter pw1=archive.openWriter("uorf-dna.fa");
			final PrintWriter pw2=archive.openWriter("uorf-pep.fa");
			final PrintWriter pw3=archive.openWriter("uorf.bed");
			pw3.println("track name=\"uORF\" description=\"uORF for "+this.knownGeneUri+"\"");
			pw3.println("#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts");
			
			final ProgressFactory.Watcher<Locatable> progress=ProgressFactory.newInstance().dictionary(dict).logger(LOG).build();
			for(final KnownGene kg:this.knownGeneMap.values().
					stream().
					flatMap(L->L.stream()).
					sorted(locCmp2).
					collect(Collectors.toList())) {
				/* new reference sequence */
				if(this.genomicSequence==null || !this.genomicSequence.getChrom().equals(kg.getContig())) {
					this.genomicSequence = new GenomicSequence(this.indexedFastaSequenceFile, kg.getContig());
					}
				
				final KnownGeneRNA mRNA=new KnownGeneRNA(kg);
				final UpstreamORF uorf=new UpstreamORF(mRNA);
				progress.apply(new Interval(uorf.getContig(),uorf.getChromStart()+1,uorf.getChromEnd()));
				final OpenReadingFrame orf = uorf.getBestOpenReadingFrame();
				if(orf==null) continue;
				orf.printFastaDNA(pw1);
				orf.printFastaPep(pw2);
				orf.printBed(pw3);
				};
			pw1.flush(); pw1.close();
			pw2.flush(); pw2.close();
			pw3.flush(); pw3.close();
			progress.close();
			
			archive.close();
			return 0;
		} catch(Exception err) {
			LOG.error(err);
			return -1;
		}
	}
	
	@Override
	public int doWork(final List<String> args) {
		try {
			this.indexedFastaSequenceFile = new IndexedFastaSequenceFile(this.faidx);
			final SAMSequenceDictionary refDict = SequenceDictionaryUtils.extractRequired(this.indexedFastaSequenceFile);
			this.refCtgNameConverter= ContigNameConverter.fromOneDictionary(refDict);
						
			LOG.info("Loading "+this.knownGeneUri);
			try(BufferedReader br= IOUtils.openURIForBufferedReading(this.knownGeneUri)) {
				String line;
				final Set<String> unknownContigs = new HashSet<>();
				final CharSplitter tab=CharSplitter.TAB;
				// tmp IntervalTreeMap for gene, will be used to remove uORF overlapping alternate transcript with CDS */
				final IntervalTreeMap<List<KnownGene>> tmpTreeMap = new IntervalTreeMap<>();
				while((line=br.readLine())!=null)
					{
					if(StringUtils.isBlank(line) || line.startsWith("#"))continue;
					
					final String tokens[]=tab.split(line);
					final KnownGene kg=new KnownGene(tokens);
					if(kg.isNonCoding()) continue;
					
					if(this.plus_strand_only && !kg.isPositiveStrand()) continue;
					
					final String refContig=this.refCtgNameConverter.apply(kg.getContig());
					if(StringUtils.isBlank(refContig)) {
						if(unknownContigs.add(kg.getContig())) {
							LOG.warn("unknown contig: "+kg.getContig());
							}
						continue;
						}
					
					kg.setChrom(refContig);
					final Interval interval = new Interval(
						refContig,
						kg.getTxStart()+1,
						kg.getTxEnd()
						);
					
					List<KnownGene> L =  tmpTreeMap.get(interval);
					if(L==null) {
						L=new ArrayList<>();
						tmpTreeMap.put(interval,L);
						}
					L.add(kg);
					}
				
				for(final KnownGene kg:tmpTreeMap.values().
						stream().
						flatMap(L->L.stream()).
						collect(Collectors.toList())) {
					final Interval interval = kgToUTRInterval(kg);
					if(interval==null || (1+interval.getEnd()-interval.getStart())<1) continue;
					
					if(this.exclude_cds_overlaping_alternative) {
						if(tmpTreeMap.getOverlapping(interval).
							stream().
							flatMap(L->L.stream()).
							filter(G->G!=kg).//same object in memory
							anyMatch(K->overlapCDS(kg,K))
							) {
							//LOG.debug("overlap "+kg.getName());
							continue;
							}
						}
					
					List<KnownGene> L =  this.knownGeneMap.get(interval);
					if(L==null) {
						L=new ArrayList<>();
						this.knownGeneMap.put(interval,L);
						}
					if(this.canonical_utr) {
						if(L.stream().
							map(K->kgToUTRInterval(K)).
							filter(K->K!=null).
							anyMatch(R->R.isNegativeStrand()==interval.isNegativeStrand() && R.equals(interval))
							) continue;
						}
					L.add(kg);
					}
				
				LOG.info("number of transcripts :"+this.knownGeneMap.values().stream().flatMap(L->L.stream()).count());
				}
  
			if(this.knownGeneMap.isEmpty()) {
				LOG.error("no transcripts found in "+this.knownGeneUri);
				return -1;
				}
			
			if(this.outputFile!=null && args.isEmpty() && (this.outputFile.isDirectory() || this.outputFile.getName().endsWith(".zip"))) {
				return this.saveArchive(this.outputFile,refDict);
				}
			else
				{
				return super.doVcfToVcf(args, this.outputFile);
				}
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.indexedFastaSequenceFile);
			}
		
		}
	
	public static void main(final String[] args) {
		new VcfScanUpstreamOrf().instanceMainWithExit(args);
	}
}