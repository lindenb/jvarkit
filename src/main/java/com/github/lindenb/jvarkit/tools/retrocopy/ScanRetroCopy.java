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
package com.github.lindenb.jvarkit.tools.retrocopy;

import java.io.BufferedReader;
import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;


import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.lang.AbstractCharSequence;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
/**
BEGIN_DOC

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



END_DOC

*/
@Program(name="scanretrocopy",
description="Scan BAM for retrocopies",
keywords={"sam","bam","cigar","clip","sv","retrocopy"},
creationDate="2019-01-25",
modificationDate="2019-02-14"
)
public class ScanRetroCopy extends Launcher
	{
	private static final Logger LOG = Logger.build(ScanRetroCopy.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	
	@Parameter(names={"-r","-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private File faidx = null;
	@Parameter(names={"-k","-K","--kg","-kg"},description=KnownGene.OPT_KNOWNGENE_DESC)
	private String knownGeneUri = KnownGene.getDefaultUri();
	@Parameter(names={"-n","--min-cigar-size"},description="Minimal cigar element length.")
	private int minCigarSize = 6;
	@Parameter(names={"-m"},description="Ignore reads having a clip lower than this value",hidden=true)
	private int _priv_ignoreCigarSize = 1;
	@Parameter(names={"--bai","-bai","--with-bai"},description="Use random access BAM using the bai and using the knownGene data. May be slow at startup")
	private boolean use_bai;
	@Parameter(names={"--partition"},description=SAMRecordPartition.OPT_DESC)
	private SAMRecordPartition partiton=SAMRecordPartition.sample;
	@Parameter(names={"--coding"},description="ignore non-coding transcripts.")
	private boolean onlyCodingTranscript=false;
	@Parameter(names={"--malus"},description="use malus value in score. bad idea. Due to alernative splicing, there is often a cigar 'M' containing the next exon.",hidden=true)
	private boolean use_malus=false;
	@Parameter(names={"--min-depth","-D"},description="In a transcript one must found at least 'D' reads with a clip-length> 'min-cigar-size'.")
	private int min_depth=1;
	@Parameter(names={"--low-depth","-d"},description="Min number of reads to set FILTER=PASS.",hidden=true)
	private int low_depth_threshold=10;
	@Parameter(names={"--bedpe","-P","-J"},description="Optional. Save possible sites of insertion in this Bed-PE file.")
	private File saveBedPeTo=null;
	@Parameter(names={"--insertion-distance","-id"},description="for insertion,s merge sites with a distance lower than this value. "+
		DistanceParser.OPT_DESCRIPTION,
		converter=DistanceParser.StringConverter.class,
		splitter=NoSplitter.class,
		hidden=true
		)
	private int merge_distance = 1_000;
	@Parameter(names={"--bam"},description="Optional: save matching read in this bam file")
	private File saveBamTo = null;
	@Parameter(names={"--both"},description="Force the constraint that both sides of a deleted intron should have at least '--min-depth' reads ")
	private boolean force_both_side=false;

	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;
	private ContigNameConverter refCtgNameConverter =null;
	private GenomicSequence genomicSequence=null;
	private PrintWriter saveInsertionsPw = null;
	private final IntervalTreeMap<List<KnownGene>> knownGenesMap = new IntervalTreeMap<>();
	private final List<Match> intronBuffer=new ArrayList<>(100_000);
	

	private static final String ATT_BEST_MATCHING_LENGTH="MAXLEN";
	private static final String ATT_FILTER_NONDOCODING="NON_CODING";
	private static final String ATT_SAMPLES="SAMPLES";
	private static final String ATT_LOW_DEPTH_FILTER="LowQual";
	private static final String ATT_INSERTION="INS";
	private final static String ATT_KG_STRAND= "STRAND";
	private final static String ATT_INTRONS_INFO="SPLICED";
	private final static String ATT_INTRONS_COUNT="ITC";
	private final static String ATT_INTRONS_CANDIDATE_COUNT="ICC";
	private final static String ATT_INTRONS_CANDIDATE_FRACTION="ICF";
	private final static String ATT_NOT_ALL_INTRONS="NOT_ALL_INTRONS";
	private final static String ATT_GT_INTRON="INTRONS";
	
	
	
	/* one-based sequence for a cigarelement */
	private class CigarLocatable 
		extends AbstractCharSequence
		implements Locatable
		{	
		final SAMRecord record;
		final int cigar_index;
		final int chromStart1;
		final int readStart0;
		CigarLocatable(final String refContig,final SAMRecord record,final int cigar_index) {
			this.record=record;
			this.cigar_index=cigar_index;
			
			int ref1 = this.record.getUnclippedStart();
			int read0 = 0;
			for(int i=0;i< this.cigar_index;i++) {
				final CigarElement ce = this.record.getCigar().getCigarElement(i);
				final CigarOperator op =ce.getOperator();
				if(op.consumesReferenceBases() || op.equals(CigarOperator.S)) {
					ref1 += ce.getLength();
					}
				if(op.consumesReadBases()) {
					read0 += ce.getLength();
					}
 				}
			this.chromStart1 =ref1;
			this.readStart0 =read0;
			}
		CigarElement getCigarElement() {
			return this.record.getCigar().getCigarElement(this.cigar_index);
			}
		@Override
		public String getContig() {
			return genomicSequence.getChrom();//may be not the same as record.getContig
			}
		@Override
		public int getStart() {
			return this.chromStart1;
			}
		@Override
		public int getEnd() {
			return getStart()+this.size()-1;
			}
		@Override
		public final int length() {
			return this.size();
			}
		public int size() {
			return this.getCigarElement().getLength();
			}
		@Override
		public final char charAt(int index) {
			return readBaseAt0(index);
			}
		public char readBaseAt0(int readPos0) {
			if(readPos0<0 || readPos0>=this.size()) throw new IndexOutOfBoundsException(String.valueOf(readPos0+"/size="+size()));
			final byte bases[]=this.record.getReadBases();
			return Character.toUpperCase((char)bases[this.readStart0+readPos0]);
			}
		}
	
	/** exon with one based coordinate */
	private class ExonOne extends AbstractCharSequence implements Locatable
		{
		private final KnownGene.Exon delegate;
		ExonOne(final KnownGene.Exon delegate) {
			this.delegate = delegate;
			}
		@Override
		public String getContig() {
			return delegate.getGene().getContig();
			}
		@Override
		public int getStart() {
			return this.delegate.getStart()+1;
			}
		@Override
		public int getEnd() {
			return this.delegate.getEnd();
			}
		
		/** implements charAt in CHROMOSOME space */
		public char charAt1(int gpos1) {
			if(gpos1<getStart()) {
				LOG.error("charAt1 out of bound ??"+gpos1+" <"+getStart());
				return 'N';
				}
			if(gpos1>getEnd()) {
				LOG.error("charAt1 out of bound ??"+gpos1+" >"+getEnd());
				return 'N';
				}
			if(gpos1<1 || gpos1>genomicSequence.length()) return 'N';
			return genomicSequence.charAt(gpos1-1);
			}
		@Override
		public int length() {
			return CoordMath.getLength(getStart(), getEnd());
			}
		@Override
		/** WARNING: implements charAt in exon space, NOT chromosome space */
		public char charAt(int index) {
			return charAt1(getStart()+index);
			}
		}	
	
	private class Match implements Locatable
		{
		final KnownGene knownGene;
		final String sampleName;
		final int intron_index;
		Interval junction = null;
		final byte side;/* 5 or 3 */
		final int clip_length;
		
		Match(final KnownGene.Intron intron,final String sampleName,final SAMRecord rec,final byte side,final int clip_length) {
			this.knownGene = intron.getGene();
			this.sampleName = sampleName;
			this.intron_index = intron.getIndex();
			this.clip_length = clip_length;
			this.side = side;
			
			if(rec.getReadPairedFlag()&& !rec.getProperPairFlag() && !rec.getMateUnmappedFlag()) {
				final String mateCtg = ScanRetroCopy.this.refCtgNameConverter.apply(rec.getMateReferenceName());
				if(StringUtils.isBlank(mateCtg)) return;
				
				final int mateEnd= SAMUtils.getMateCigar(rec)!=null?
						SAMUtils.getMateAlignmentEnd(rec):
						rec.getMateAlignmentStart()
						;
				
				final Interval mateInterval = new Interval(mateCtg,rec.getMateAlignmentStart(),mateEnd);		
				
				if(mateCtg.equals(this.knownGene.getContig()) && ScanRetroCopy.this.knownGenesMap.
					getOverlapping(this.knownGene).
					stream().
					flatMap(L->L.stream()).
					anyMatch(R->R.overlaps(mateInterval))) {			
					return;
					}	
				
				this.junction = mateInterval;
				}
			}
		

		@Override
		public String getContig() {
			return this.knownGene.getContig();
			}
		
		@Override
		public int getStart() {
			return this.knownGene.getIntronStart(this.intron_index)+1;
			}
		
		@Override
		public int getEnd() {
			return this.knownGene.getIntronEnd(this.intron_index);
			}
		/*
		boolean isInterval(final Locatable loc) {
			
			
			return loc.getContig().equals(this.getContig()) &&
					loc.getStart() == this.getStart() &&
					loc.getEnd() == this.getEnd()
					;
			}
		
		Interval toInterval() {
			return new Interval(this);
			}*/
		
		@Override
		public String toString() {
			return  this.knownGene.getName()+" "+this.knownGene.getContig()+":"+ this.knownGene.getIntronStart(this.intron_index)+"-"+this.knownGene.getIntronEnd(this.intron_index)+" "+this.sampleName+" "+this.side+" "+this.clip_length;
			}
		
		}
	
		
	private class JunctionInfo
		{
		final List<Integer> clip_lengths = new ArrayList<>();
		void visit(final Match match) {
			this.clip_lengths.add(match.clip_length);
			}
		// a exon exon junction is valid if there is at least 'min_depth' cigar with length >= minCigarSize
		boolean hasValidDepth() {
			return bestDepth() >= (long)ScanRetroCopy.this.min_depth;
			}
		
		IntStream all() {
			return this.clip_lengths.stream().mapToInt(X->X.intValue());
		}
		
		int bestDepth() {
			return (int)all().
					filter(X->X>=ScanRetroCopy.this.minCigarSize).
					count();
			}
		
		int longestClip() {
			return all().max().orElse(0);
			}
		double average() {
			return all().average().orElse(0.0);
			}
		}
	
	private class GeneInfo
		{
		final List<JunctionInfo> intron_5_side;
		final List<JunctionInfo> intron_3_side;
		GeneInfo(final KnownGene kg) {
			//we must find one exon-exon junction with a least 'min_depth'
			this.intron_5_side = new ArrayList<>(kg.getIntronCount());
			this.intron_3_side = new ArrayList<>(kg.getIntronCount());
			for(int i=0;i< kg.getIntronCount();i++) {
				this.intron_5_side.add(new JunctionInfo());
				this.intron_3_side.add(new JunctionInfo());
				}
			}
		Stream<JunctionInfo> all() {
			return Stream.concat(intron_5_side.stream(),intron_3_side.stream());
		}
		
		Stream<JunctionInfo> intron(int intron_index) {
			return Arrays.asList(this.intron_5_side.get(intron_index),this.intron_3_side.get(intron_index)).stream();
		}
		
		void visit(final Match match) {
			switch((int)match.side) {
				case 5: this.intron_5_side.get(match.intron_index).visit(match);break;
				case 3: this.intron_3_side.get(match.intron_index).visit(match);break;
				default: throw new IllegalStateException();
				}
			}
		
		
		
		boolean hasValidDepth(int intron_index) {
			return intron(intron_index).anyMatch(M->M.hasValidDepth());
			}
		
		boolean hasValidDepth() {
			if(ScanRetroCopy.this.force_both_side)
				{
				return this.intron_5_side.stream().anyMatch(X->X.hasValidDepth()) &&
					   this.intron_3_side.stream().anyMatch(X->X.hasValidDepth());
				}
			else
				{
				return all().anyMatch(X->X.hasValidDepth());
				}
			}
		int bestDepth() {
			return all().filter(V->V.hasValidDepth()).mapToInt(M->M.bestDepth()).max().orElse(0);
			}
		int longestClip() {
			return all().filter(V->V.hasValidDepth()).mapToInt(M->M.longestClip()).max().orElse(0);
			}
		Genotype makeGenotype(final KnownGene kg,final String sample,List<Allele> refalt) {
			final GenotypeBuilder gb ;
			if(!this.hasValidDepth())
				{
				gb=new GenotypeBuilder(sample, Arrays.asList(refalt.get(0),refalt.get(0)));
				}
			else
				{
				gb=new GenotypeBuilder(sample,refalt);
				}
			gb.DP(bestDepth());
			gb.attribute(ATT_BEST_MATCHING_LENGTH, longestClip());
			final StringBuilder sb=new StringBuilder();
			final String formatDbl="%.2f";
			
			for(int intron_index=0;intron_index < kg.getIntronCount();++intron_index) {
				if(intron_index>0) sb.append("|");
				sb.append("i").append(intron_index).append(",");
				sb.append(hasValidDepth(intron_index)?"Y":"N").append(",");
				
				sb.append(this.intron_5_side.get(intron_index).bestDepth());
				sb.append(",");
				sb.append(this.intron_3_side.get(intron_index).bestDepth());
				sb.append(",");
				sb.append(this.intron_5_side.get(intron_index).longestClip());
				sb.append(",");
				sb.append(this.intron_3_side.get(intron_index).longestClip());
				sb.append(",");
				sb.append(String.format(formatDbl,this.intron_5_side.get(intron_index).average()));
				sb.append(",");
				sb.append(String.format(formatDbl,this.intron_3_side.get(intron_index).average()));

				
				}
			gb.attribute(ATT_GT_INTRON,sb.toString());
			/*
			gb.AD(new int[] {
				this.intron_5_side.stream().filter(X->hasValidDepth()).mapToInt(X->X.bestDepth()).max().orElse(0)
				,	
				this.intron_3_side.stream().filter(X->hasValidDepth()).mapToInt(X->X.bestDepth()).max().orElse(0)
				});*/
			return gb.make();
			}
		}
	
	private void dump(final VariantContextWriter vcw,final Locatable before) {
		final Allele alt= Allele.create("<RETROCOPY>", false);
		
		/* get a list of overlapping gene as string + coding state*/
		final Function<Locatable,String> findGenes = R->{
			final String s1 = ScanRetroCopy.this.knownGenesMap.getOverlapping(R).
					stream().
					flatMap(G->G.stream()).
					map(G->G.getName()).
					sorted().
					collect(Collectors.joining(";"));
			final boolean coding = knownGenesMap.getOverlapping(R).
					stream().
					flatMap(G->G.stream()).
					anyMatch(G->!G.isNonCoding());

			return (s1.isEmpty()?".":s1)+"\t"+(coding?".":ATT_FILTER_NONDOCODING);
			};
		
		// genes to be considered for this dump
		final Set<KnownGene> candidateGenes = this.intronBuffer.
				stream().
				map(K->K.knownGene).
				filter(K->before==null || K.getEnd() < before.getStart()).
				collect(Collectors.toCollection(()->new TreeSet<KnownGene>((A,B)-> {
					final int i= Integer.compare(A.getStart(), B.getStart());
					if(i!=0) return i;
					return A.getName().compareTo(B.getName());
					})));
		
		final Set<String> candidateSamples = this.intronBuffer.
				stream().
				filter(M->candidateGenes.contains(M.knownGene)).
				map(K->K.sampleName).
				collect(Collectors.toSet());
		
		// loop over genes
		for(final KnownGene kg: candidateGenes) {
			boolean filter_set=false;
			
			final Map<String,GeneInfo> sample2info = new HashMap<>(candidateSamples.size());
			for(final String sn:candidateSamples) sample2info.put(sn, new GeneInfo(kg));
			// visit all matches for this gene
			this.intronBuffer.
				stream().
				filter(M->M.knownGene.getName().equals(kg.getName())).
				forEach(M->sample2info.get(M.sampleName).visit(M));
			
			// we need at least one junction with a min depth
			if(sample2info.values().stream().noneMatch(GI->GI.hasValidDepth())) {
				continue;
				}
			
			// ok good candidate
			final VariantContextBuilder vcb = new VariantContextBuilder();
			vcb.chr(kg.getContig());
			vcb.start(kg.getStart());
			vcb.stop(kg.getEnd());
			vcb.id(kg.getName());
			final Allele ref= Allele.create((byte)genomicSequence.charAt(kg.getTxStart()), true);
			final List<Allele> alleles = Arrays.asList(ref,alt);

			final int max_depth = sample2info.values().stream().mapToInt(X->X.bestDepth()).max().orElse(0);
			vcb.attribute(VCFConstants.DEPTH_KEY,max_depth);
			vcb.log10PError(max_depth/-10.0);

			if(max_depth < ScanRetroCopy.this.low_depth_threshold)
				{
				vcb.filter(ATT_LOW_DEPTH_FILTER+ScanRetroCopy.this.low_depth_threshold);
				filter_set = true;
				}
			
			final int AC=(int)sample2info.values().stream().filter(X->X.hasValidDepth()).count();
			final int AN=2*sample2info.size();
			vcb.attribute(VCFConstants.ALLELE_NUMBER_KEY,AN);
			vcb.attribute(VCFConstants.ALLELE_COUNT_KEY,AC);
			if(AN>0) vcb.attribute(VCFConstants.ALLELE_FREQUENCY_KEY,AC/(double)AN);
			vcb.attribute(VCFConstants.SVTYPE,"DEL");
			vcb.attribute(VCFConstants.END_KEY,kg.getEnd());
			vcb.attribute(ATT_KG_STRAND,kg.isNegativeStrand()?"minus":"plus");
			vcb.attribute(ATT_BEST_MATCHING_LENGTH,sample2info.values().stream().filter(X->X.hasValidDepth()).mapToInt(M->M.longestClip()).max().orElse(0));
			vcb.attribute("SVLEN",kg.getLengthOnReference());
			
			
			vcb.alleles(alleles);
			
			vcb.attribute(ATT_SAMPLES, new ArrayList<>(
					sample2info.entrySet().stream().
					filter(KV->KV.getValue().hasValidDepth()).
					map(KV->KV.getKey()).
					collect(Collectors.toCollection(TreeSet::new))));

			// introns sequences
			final List<String> intronInfos=new ArrayList<>(kg.getIntronCount());
			for(int intron_idx=0;intron_idx < kg.getIntronCount();++intron_idx) {
				final KnownGene.Intron the_intron = kg.getIntron(intron_idx);
				final int tmp_idx = intron_idx;
				if(sample2info.values().stream().noneMatch(X->X.hasValidDepth(tmp_idx))) continue;
				
				final CharSequence intronSequence = ScanRetroCopy.this.genomicSequence.subSequence(the_intron.getStart(),the_intron.getEnd());
				final StringBuilder sb=new StringBuilder(the_intron.getName().replaceAll("[ ]","_"));
				sb.append("|");
				sb.append(String.valueOf(the_intron.getStart()+1));
				sb.append("|");
				sb.append(String.valueOf(the_intron.getEnd()+1));
				sb.append("|");
				sb.append(StringUtils.left(intronSequence,ScanRetroCopy.this.minCigarSize));
				sb.append("|");
				sb.append(StringUtils.right(intronSequence,ScanRetroCopy.this.minCigarSize));
				
				
				
				intronInfos.add(sb.toString());
				}
			vcb.attribute(ATT_INTRONS_INFO,intronInfos);
			vcb.attribute(ATT_INTRONS_CANDIDATE_COUNT,intronInfos.size());
			vcb.attribute(ATT_INTRONS_COUNT,kg.getIntronCount());
			vcb.attribute(ATT_INTRONS_CANDIDATE_FRACTION,intronInfos.size()/(float)kg.getIntronCount());
			if(kg.getIntronCount()!=intronInfos.size()) {
				vcb.filter(ATT_NOT_ALL_INTRONS);
				filter_set=true;
				}
			
			
			/* build genotypes */
			final List<Genotype> genotypes = new ArrayList<>(sample2info.size());
			for(final String sample: sample2info.keySet()) {
				final GeneInfo geneInfo = sample2info.get(sample);
				genotypes.add(geneInfo.makeGenotype(kg,sample,alleles));
				}
			
			
			/* insertions */
			final List<Interval> insertions = this.intronBuffer.
					stream().
					filter(M->M.junction!=null).
					filter(M->M.knownGene.getName().equals(kg.getName())).
					filter(M->sample2info.get(M.sampleName).hasValidDepth(M.intron_index)).
					map(M->M.junction).
					sorted().
					collect(Collectors.toCollection(ArrayList::new));
			
			if(!insertions.isEmpty())
				{
				final Set<String> jset = new HashSet<>();
				int i=0;
				while(i<insertions.size()) {
					Interval insertion = insertions.get(i);
					int j=i+1;
					int count_evidence =1;
					while(j<insertions.size() ) {
						final Interval m = insertions.get(j);
						if(!insertion.withinDistanceOf(m, ScanRetroCopy.this.merge_distance)) break;
						insertion= new Interval(insertion.getContig(),Math.min(insertion.getStart(), m.getStart()),Math.max(insertion.getEnd(), m.getEnd()));
						insertions.remove(j);
						++count_evidence;
						}
					i=j;
					final List<KnownGene> mateGenes = ScanRetroCopy.this.knownGenesMap.getOverlapping(insertion).
								stream().
								flatMap(G->G.stream()).
								sorted((A,B)->A.getName().compareTo(B.getName())).
								collect(Collectors.toList());
					final StringBuilder sb=new StringBuilder(insertion.getContig()+":"+insertion.getStart()+"-"+insertion.getEnd());
					sb.append("|");
					sb.append(count_evidence);
					sb.append("|");
					sb.append(mateGenes.isEmpty()?".":mateGenes.stream().map(G->G.getName()).collect(Collectors.joining("&")));
					sb.append("|");
					sb.append(!mateGenes.isEmpty() && mateGenes.stream().allMatch(G->G.isNonCoding())?ATT_FILTER_NONDOCODING:".");
					sb.append("|");
					
					final String distanceStr;
					
					if(insertion.overlaps(kg)) {
						distanceStr = "0";
						}
					else  if(insertion.contigsMatch(kg)) {
						distanceStr = String.valueOf(Math.abs(Math.min(kg.getStart()-insertion.getEnd(),insertion.getStart()-kg.getEnd())));
						}
					else
						{
						distanceStr = "NOT_SAME_CONTIG";
						}
					sb.append(distanceStr);
					jset.add(sb.toString());
					
					// save insertions
					saveInsertionsPw.print(kg.getContig());
					saveInsertionsPw.print("\t");
					saveInsertionsPw.print(kg.getStart()-1);
					saveInsertionsPw.print("\t");
					saveInsertionsPw.print(kg.getEnd());
					saveInsertionsPw.print("\t");
					saveInsertionsPw.print(insertion.getContig());
					saveInsertionsPw.print("\t");
					saveInsertionsPw.print(insertion.getStart()-1);
					saveInsertionsPw.print("\t");
					saveInsertionsPw.print(insertion.getEnd());
					saveInsertionsPw.print("\t");
					saveInsertionsPw.print(".");//name
					saveInsertionsPw.print("\t");
					saveInsertionsPw.print(count_evidence);//score
					saveInsertionsPw.print("\t");
					saveInsertionsPw.print(".");//strand 1
					saveInsertionsPw.print("\t");
					saveInsertionsPw.print(".");//strand 2
					saveInsertionsPw.print("\t");

					// "Any number of additional, user-defined fields ..."
					saveInsertionsPw.print(kg.getName());
					saveInsertionsPw.print("\t");
					saveInsertionsPw.print(findGenes.apply(insertion));
					saveInsertionsPw.print("\t");
					if(insertion.overlaps(kg)) 
						{
						saveInsertionsPw.print("0");
						}
					else if(insertion.contigsMatch(kg)) {
						saveInsertionsPw.print(Math.abs(Math.min(kg.getStart()-insertion.getEnd(),insertion.getStart()-kg.getEnd())));
						}
					else
						{
						saveInsertionsPw.print("NOT_SAME_CONTIG");
						}
					/* TODO
					saveInsertionsPw.print("\t");
					saveInsertionsPw.print( sampleMap.values().stream().mapToInt(M->M.countSupportingReads()).sum());
					saveInsertionsPw.print("\t");
					saveInsertionsPw.print(String.join(";",sampleMap.keySet()));
					*/
					saveInsertionsPw.println();
					}
				vcb.attribute(ATT_INSERTION,new ArrayList<>(jset));
				}
			
			
			
			if(kg.isNonCoding()) {
				vcb.filter(ATT_FILTER_NONDOCODING);
				filter_set=true;
				}
			
			if(!filter_set) {
				vcb.passFilters();
				}
			
			vcb.genotypes(genotypes);
			vcw.add(vcb.make());
			
			// cleanup
			this.intronBuffer.removeIf(M->M.knownGene.getName().equals(kg.getName()));
			}		
		
		
		if(before!=null) {
			// remove transcript if there is not enough evidence(s).
			this.intronBuffer.removeIf(M->M.knownGene.getEnd() < before.getStart());
			} 
		else
			{
			this.intronBuffer.clear();
			//this.kgId2knownGenes.clear();
			}

		
		
		}
	
	private boolean isCandidateCigarElement(final CigarElement C) {
		return C.getOperator().equals(CigarOperator.S) && C.getLength()>=this._priv_ignoreCigarSize;
	}

	
	@Override
	public int doWork(final List<String> args) {
		if(this.min_depth <1) {
			LOG.error("Bad min depth");
			return -1;
			}
		
		SamReader sr = null;
		VariantContextWriter vcw0=null;
		CloseableIterator<SAMRecord> iter = null;
		SAMFileWriter sfw = null;
		try {
			this.indexedFastaSequenceFile = new IndexedFastaSequenceFile(this.faidx);
			final SAMSequenceDictionary refDict = SequenceDictionaryUtils.extractRequired(this.indexedFastaSequenceFile);
			this.refCtgNameConverter= ContigNameConverter.fromOneDictionary(refDict);

			/* READ KNOWGENES FILES */
			LOG.info("Loading "+this.knownGeneUri);
			try(BufferedReader br= IOUtils.openURIForBufferedReading(this.knownGeneUri)) {
				String line;
				final CharSplitter tab=CharSplitter.TAB;
				while((line=br.readLine())!=null)
					{
					if(StringUtils.isBlank(line))continue;
					final String tokens[]=tab.split(line);
					final KnownGene kg=new KnownGene(tokens);
					if(kg.getExonCount()<2) continue;
					if(this.onlyCodingTranscript && kg.getCdsStart()==kg.getCdsEnd())continue; 
					final String ctg = this.refCtgNameConverter.apply(kg.getContig());
					if(StringUtils.isBlank(ctg)) continue;
					kg.setChrom(ctg);
					final Interval interval = new Interval(ctg,kg.getTxStart()+1,kg.getTxEnd(),kg.isNegativeStrand(),kg.getName());
					List<KnownGene> L=  this.knownGenesMap.get(interval);
					if(L==null) {
						L=new ArrayList<KnownGene>();
						this.knownGenesMap.put(interval,L);
						}
					L.add(kg);
					}
				
				}

			if(this.knownGenesMap.isEmpty()) {
				LOG.error("no gene found in "+this.knownGeneUri);
				return -1;
				}
			LOG.info("Number of transcripts: "+ this.knownGenesMap.values().stream().flatMap(L->L.stream()).count());
			
			sr = super.openSamReader(oneFileOrNull(args));
			final SAMFileHeader samFileHeader = sr.getFileHeader();
			if(!samFileHeader.getSortOrder().equals(SAMFileHeader.SortOrder.coordinate)) {
				LOG.error("input is not sorted on coordinate but \""+samFileHeader.getSortOrder()+"\"");
				return -1;
			}
			
			if(this.saveBamTo!=null) {
				sfw = new SAMFileWriterFactory().
						makeSAMOrBAMWriter(samFileHeader, true, this.saveBamTo);
				}
			
			if(this.use_bai && !sr.hasIndex()) {
				LOG.warning("Cannot used bai because input is not indexed");
				}
			
			if(this.use_bai && sr.hasIndex())
				{
				LOG.info("building intervals...");
				final SAMSequenceDictionary samdict= SequenceDictionaryUtils.extractRequired(samFileHeader);

				final ContigNameConverter samConvert = ContigNameConverter.fromOneDictionary(samdict);
				final List<QueryInterval> intervalsL = this.knownGenesMap.values().
						stream().
						flatMap(K->K.stream()).
						filter(KG->samConvert.apply(KG.getContig())!=null).
						flatMap(KG->KG.getExons().stream()).
						flatMap(exon->{
							// we need the reads overlapping the exon bounds
							final int tid=samdict.getSequenceIndex(samConvert.apply(exon.getGene().getContig()));
							final QueryInterval q1=new QueryInterval(tid,exon.getStart()+1,exon.getStart()+1);
							final QueryInterval q2=new QueryInterval(tid,exon.getEnd(),exon.getEnd());
							return Arrays.stream(new QueryInterval[]{q1,q2});
						}).
						sorted().
						collect(Collectors.toList());
				
				final QueryInterval intervals[]=QueryInterval.optimizeIntervals(intervalsL.toArray(new QueryInterval[intervalsL.size()]));
				intervalsL.clear();//GC
				LOG.debug("Query bam using "+intervals.length+" random access intervals. Please wait...");
				iter = sr.queryOverlapping(intervals);
				}
			else
				{
				iter= sr.iterator();
				}
			final Set<String> samples = this.partiton.getPartitions(samFileHeader);
			if(samples.isEmpty()) {
				LOG.error("No sample was defined in the read group of the input bam.");
				return -1;
				}
			/** build vcf header */
			final Set<VCFHeaderLine> metaData = new HashSet<>();
			metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_KEY,true));
			metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_QUALITY_KEY,true));
			metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.DEPTH_KEY,true));
			metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_ALLELE_DEPTHS,true));
			metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY,true));
			metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_NUMBER_KEY,true));
			metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY,true));
			metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY,true));
			metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_FREQUENCY_KEY,true));
			metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY,true));
			metaData.add(new VCFInfoHeaderLine(VCFConstants.SVTYPE, 1, VCFHeaderLineType.String,"Variation type"));
			metaData.add(new VCFInfoHeaderLine("SVLEN", 1, VCFHeaderLineType.Integer,"Variation Length"));
			metaData.add(new VCFInfoHeaderLine(ATT_BEST_MATCHING_LENGTH, 1,VCFHeaderLineType.Integer,"Best Matching length"));
			metaData.add(new VCFFormatHeaderLine(ATT_BEST_MATCHING_LENGTH, 1,VCFHeaderLineType.Integer,"Best Matching length"));
			metaData.add(new VCFFormatHeaderLine(ATT_GT_INTRON, 1,VCFHeaderLineType.String,
						"Introns info: (intron-0-idx,valid,dp-5,dp-3,max-len-5,max-len-3,avg-5,avg-3)*"));
			
			
			//metaData.add(new VCFFormatHeaderLine(ATT_COUNT_SUPPORTING_READS, 2,VCFHeaderLineType.Integer,"Count supporting reads [intron-left/intron-right]"));
			//metaData.add(new VCFInfoHeaderLine(ATT_RETRO_DESC, VCFHeaderLineCount.UNBOUNDED,VCFHeaderLineType.String,
			//		"Retrocopy attributes: transcript-id|strand|exon-left|exon-left-bases|exon-right-bases|exon-right"));
			metaData.add(new VCFInfoHeaderLine(ATT_INSERTION, VCFHeaderLineCount.UNBOUNDED,VCFHeaderLineType.String,
					"Possible place of insertion:"+ "chr:start-end|count-evidence|mate-genes|non-coding|distance"));
			metaData.add(new VCFInfoHeaderLine(ATT_KG_STRAND, 1, VCFHeaderLineType.String,"KnownGene strand."));
			metaData.add(new VCFInfoHeaderLine(ATT_INTRONS_COUNT, 1, VCFHeaderLineType.Integer,"Number of introns for the Transcript"));
			metaData.add(new VCFInfoHeaderLine(ATT_INTRONS_CANDIDATE_COUNT, 1, VCFHeaderLineType.Integer,"Number of introns found retrocopied for the transcript"));
			metaData.add(new VCFInfoHeaderLine(ATT_INTRONS_CANDIDATE_FRACTION, 1, VCFHeaderLineType.Float,"Fraction of introns found retrocopied for the transcript"));
			metaData.add(new VCFInfoHeaderLine(ATT_INTRONS_INFO, VCFHeaderLineCount.UNBOUNDED,VCFHeaderLineType.String,
					"Introns found: chr|start|end|seq-left|seq-right"));

			
			
			
			metaData.add(new VCFFilterHeaderLine(ATT_FILTER_NONDOCODING,"Only non-coding transcripts"));
			metaData.add(new VCFInfoHeaderLine(ATT_SAMPLES,VCFHeaderLineCount.UNBOUNDED,VCFHeaderLineType.String,"Samples found. partition:"+this.partiton.name()));
			metaData.add(new VCFFilterHeaderLine(ATT_LOW_DEPTH_FILTER+this.low_depth_threshold,"Number of read is lower than :"+this.min_depth));
			metaData.add(new VCFFilterHeaderLine(ATT_NOT_ALL_INTRONS,"Not all introns were found retrocopied"));

			
			
			
			
			final VCFHeader header=new VCFHeader(metaData, samples);
			JVarkitVersion.getInstance().addMetaData(this, header);
			header.setSequenceDictionary(refDict);
			
			/* open vcf for writing*/
			vcw0=super.openVariantContextWriter(this.outputFile);
			final VariantContextWriter vcw=vcw0;
			vcw.writeHeader(header);
			
			/* save gene writer */
			if(this.saveBedPeTo!=null) {
				this.saveInsertionsPw = super.openFileOrStdoutAsPrintWriter(this.saveBedPeTo);
				}
			else
				{
				this.saveInsertionsPw = new PrintWriter(new NullOuputStream());
				}
			
			final ProgressFactory.Watcher<SAMRecord> progress = ProgressFactory.newInstance().dictionary(samFileHeader).logger(LOG).build();
			
			while(iter.hasNext()) {
				final SAMRecord rec = progress.apply(iter.next());
				if(rec.getReadUnmappedFlag()) continue;
				if(rec.isSecondaryOrSupplementary()) continue;
				if(rec.getDuplicateReadFlag()) continue;
				final byte bases[]=rec.getReadBases();
				if(bases==null || bases==SAMRecord.NULL_SEQUENCE) continue;
				final Cigar cigar = rec.getCigar();
				if(cigar==null || cigar.numCigarElements()<2) continue;
				final String refContig = this.refCtgNameConverter.apply(rec.getContig());
				boolean save_read_to_bam = false;
				
				if(StringUtils.isBlank(refContig)) continue;
				
				/* get sample */
				final String sampleName = this.partiton.getPartion(rec, null);
				if(StringUtils.isBlank(sampleName)) continue;

				/* this is a new reference sequence */
				if(this.genomicSequence==null || !this.genomicSequence.getChrom().equals(refContig)) {
					if(this.genomicSequence!=null) {
						
						/* DUMP things BEFORE changing the reference sequence!!! */						
						/* dump buffer */
						dump(vcw,null);
						}
					/* map transcript-name to their  transcript */
					/*this.kgId2knownGenes.clear();
					this.knownGenesMap.values().
						stream().
						flatMap(L->L.stream()).
						filter(G->refContig.equals(G.getContig())).
						forEach(K->this.kgId2knownGenes.put(K.getName(), K));*/
					/* now, we can change genomicSequence */
					this.genomicSequence = new GenomicSequence(this.indexedFastaSequenceFile, refContig);
					}
			
				
				
				final CigarElement leftCigar = cigar.getCigarElement(0);
				final CigarElement rightCigar = cigar.getCigarElement(cigar.numCigarElements()-1);
				
				/* both ends are not candidate */
				if(!isCandidateCigarElement(leftCigar) && !isCandidateCigarElement(rightCigar) ) {
					continue;
				}
				
				final List<KnownGene> genes = this.knownGenesMap.getOverlapping(
						new Interval(refContig,rec.getUnclippedStart(),rec.getUnclippedEnd())
						).stream().
						flatMap(L->L.stream()).
						collect(Collectors.toList());
				
				
				/* time to time dump the buffer for the transcripts before the current ones*/
				if(!genes.isEmpty()) {
					// get all genes overlapping the current set of genes
					int minTxStart= genes.stream().mapToInt(K->K.getTxStart()).min().getAsInt();
					int maxTxStart= genes.stream().mapToInt(K->K.getTxEnd()).max().getAsInt();
					// update minTxtStart to get the lowest gene overlapping the set of genes
					minTxStart = this.knownGenesMap.getOverlapping(
							new Interval(refContig,minTxStart,maxTxStart)
							).stream().
							flatMap(L->L.stream()).
							mapToInt(K->K.getStart()).
							min().
							getAsInt();
					//not max, because we only need the 5' side
					dump(vcw,new Interval(refContig,minTxStart,minTxStart));
					}
				
				
				/* test each side of the clipped read */
				for(int side=0;side<2 && !genes.isEmpty();++side) {
					final CigarElement ce_side = (side==0?leftCigar:rightCigar);
					if(!isCandidateCigarElement(ce_side)) continue;
					for(final KnownGene knownGene:genes) {
						for(int exonIndex=0;exonIndex< knownGene.getExonCount();exonIndex++) {
							if(side==0) /* looking at cigar string in 5' */
								{
								if(exonIndex==0) continue;
								
								//last 'M' element
								final CigarLocatable cigarM = new CigarLocatable(refContig, rec,1);
								
								//last cigar element
								final CigarLocatable cigarS = new CigarLocatable(refContig, rec,0);
								// current exon
								final ExonOne exonRight = new ExonOne(knownGene.getExon(exonIndex));
								if(!cigarM.overlaps(exonRight)) continue;
								if(!(exonRight.getStart() >= cigarM.getStart())) continue;
								// get next exon
								final ExonOne exonLeft = new ExonOne(knownGene.getExon(exonIndex-1));
								if(exonLeft.getLengthOnReference() < this.minCigarSize) continue;
								
								/* end of cigar 'M' can have same bases than the prev exon. */
								final int malus = exonRight.getStart() - cigarM.getStart();
								
								int genomic1 = exonLeft.getEnd()-malus;
								if(genomic1<exonLeft.getStart() || genomic1>exonLeft.getEnd()) {
									continue;
								}
								
								int matchLength= (this.use_malus?malus:0);
								int readIdx0=cigarS.size()-1;
								// loop over sequence
								while(readIdx0 >=0 && genomic1 >= exonLeft.getStart()) {
									final char read_base = cigarS.readBaseAt0(readIdx0);
									final char genome_base = exonLeft.charAt1(genomic1);
									if(read_base!=genome_base)
										{
										break;
										}
									readIdx0--;
									matchLength++;
									genomic1--;
									}
								
								if(matchLength<this._priv_ignoreCigarSize) continue;
								
								final KnownGene.Intron intron=knownGene.getIntron(exonIndex-1); 
								
								//find match or create new
								final Match match = new Match(intron, sampleName, rec,(byte)3,matchLength);
								this.intronBuffer.add(match);
								save_read_to_bam = true;
								}
							else /* test last cigar */
								{
								if(exonIndex+1>=knownGene.getExonCount()) continue;
								//last 'M' element
								final CigarLocatable cigarM = new CigarLocatable(refContig, rec, cigar.numCigarElements()-2);
								
								//last cigar element
								final CigarLocatable cigarS = new CigarLocatable(refContig, rec, cigar.numCigarElements()-1);
								// current exon
								final ExonOne exonLeft = new ExonOne(knownGene.getExon(exonIndex));
								
								if(!cigarM.overlaps(exonLeft)) continue;
								if(!(exonLeft.getEnd() <= cigarM.getEnd())) continue;
								// get next exon
								final ExonOne exonRight = new ExonOne(knownGene.getExon(exonIndex+1));
								if(exonRight.getLengthOnReference() < this.minCigarSize) continue;
								
								/* end of cigar 'M' can have same bases than the next exon. */
								final int malus = cigarM.getEnd()-exonLeft.getEnd();
								
								int genomic1 = exonRight.getStart()+malus;
								if(genomic1<exonRight.getStart() || genomic1>exonRight.getEnd()) {
									continue;
								}
								
								int matchLength= (this.use_malus?malus:0);
								int readIdx0=0;
								// loop over sequence
								while(readIdx0 <cigarS.size() && genomic1 <= exonRight.getEnd()) {
									final char read_base = cigarS.readBaseAt0(readIdx0);
									final char genome_base = exonRight.charAt1(genomic1);
									if(read_base!=genome_base)
										{
										break;
										}
									readIdx0++;
									matchLength++;
									genomic1++;
									}
								
								if(matchLength<this._priv_ignoreCigarSize) continue;

								//find match or create new
								
								final KnownGene.Intron intron=knownGene.getIntron(exonIndex); 
								
								final Match match = new Match(intron, sampleName, rec,(byte)5,matchLength);
								this.intronBuffer.add(match);
								save_read_to_bam = true;
								}
							}
						}
					} //end for side
				if(save_read_to_bam && sfw!=null) sfw.addAlignment(rec);
				}
			/* dump buffer */
			dump(vcw,null);
			
			
			progress.close();
			vcw.close();
			iter.close();
			iter=null;
			sr.close();
			sr=null;
			this.saveInsertionsPw.flush();
			this.saveInsertionsPw.close();
			this.saveInsertionsPw=null; 
			
			if(sfw!=null) {
				sfw.close();
				sfw=null;
				}
			return 0;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(sr);
			CloserUtil.close(vcw0);
			CloserUtil.close(sfw);
			CloserUtil.close(this.indexedFastaSequenceFile);
			CloserUtil.close(this.saveInsertionsPw); 
			}
		}
	
	public static void main(final String[] args) {
		new ScanRetroCopy().instanceMainWithExit(args);
	}
}
