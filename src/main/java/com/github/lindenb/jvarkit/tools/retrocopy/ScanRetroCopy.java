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
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.lang.AbstractCharSequence;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
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
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
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
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
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
```

END_DOC

*/
@Program(name="scanretrocopy",
description="Scan BAM for retrocopies",
keywords={"sam","bam","cigar","sv","retrocopy"},
creationDate="2019-01-25"
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
	@Parameter(names={"-n","--min-cigar-size"},description="Minimal cigar element size.")
	private int minCigarSize = 10;
	@Parameter(names={"--bai","-bai","--with-bai"},description="Use random access BAM using the bai and using the knownGene data.")
	private boolean use_bai;
	@Parameter(names={"--partition"},description=SAMRecordPartition.OPT_DESC)
	private SAMRecordPartition partiton=SAMRecordPartition.sample;
	@Parameter(names={"--coding"},description="ignore non-coding transcript ")
	private boolean onlyCodingTranscript=false;
	@Parameter(names={"--save-gene","-S"},description="Optional. save per-gene info in this file.")
	private File saveGeneTo=null;
	@Parameter(names={"--malus"},description="use malus value in score. bad idea. Due to alernative splicing, there is often a cigar 'M' containing the next exon.",hidden=true)
	private boolean use_malus=false;

	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;
	private ContigNameConverter refCtgNameConverter =null;
	private GenomicSequence genomicSequence=null;
	private PrintWriter saveGenePw=null;

	private static final String ATT_BEST_MATCHING_LENGTH="MAXLEN";
	private static final String ATT_RETRO_DESC="RCP";
	private final Predicate<CigarElement> isCandidateCigarElement=(C)->C.getOperator().equals(CigarOperator.S) && C.getLength()>=this.minCigarSize;
	private final Map<String,GeneInfo> geneIdToInfo=new HashMap<>();
	
	/*  what was found for each gene */
	private static class GeneInfo
		{
		final KnownGene gene;
		final Set<Interval> intronSet=new HashSet<>();
		GeneInfo(final KnownGene gene) {
			this.gene = gene;
			}
		
		}
	
	/* per-sample-information */
	private static class PerSample {
		int countSupportingReads = 0;
		int bestLength=0;
		}
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
		@Override
		public String toString() {
			return getCigarElement().getOperator()+"["+getContig()+":"+getStart()+"-"+getEnd()+"]";
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
		
		public String getName() {
			return "Exon"+(this.delegate.getIndex()+1);
		}
		@Override
		public String toString() {
			return delegate.getName()+"["+getContig()+":"+getStart()+"-"+getEnd()+"]";
			}
		/** implements charAt in CHROMOSOME space */
		public char charAt1(int gpos1) {
			if(gpos1<getStart() || gpos1> getEnd()) throw new IndexOutOfBoundsException(gpos1+" "+getStart()+" "+getEnd());
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
		
		
	private class Match
		{
		final String contig;
		final int chromStart0;
		final int chromEnd0;
		final Set<String> attributes = new HashSet<>();
		final Map<String,PerSample> sampleMap = new HashMap<>();
		Match(final KnownGene.Intron intron) {
			this.contig=intron.getGene().getContig();
			this.chromStart0 = intron.getStart();
			this.chromEnd0 = intron.getEnd();
			}
		
		PerSample getSample(final String sample) {
			PerSample p = this.sampleMap.get(sample);
			if(p==null) {
				p=new PerSample();
				this.sampleMap.put(sample, p);
				}
			return p;
			}
		
		VariantContext build() {
			final VariantContextBuilder vcb = new VariantContextBuilder();
			vcb.chr(this.contig);
			vcb.start(this.chromStart0+1);
			vcb.stop(this.chromEnd0+1);
			final Allele ref= Allele.create((byte)genomicSequence.charAt(chromStart0), true);
			final Allele alt= Allele.create("<RETROCOPY>", false);
			final List<Allele> alleles = Arrays.asList(ref,alt);
			vcb.alleles(alleles);
			vcb.attribute(ATT_RETRO_DESC,new ArrayList<>(this.attributes));
			final int sum_count = this.sampleMap.values().stream().mapToInt(M->M.countSupportingReads).sum();
			vcb.attribute(VCFConstants.DEPTH_KEY,sum_count);
			vcb.attribute(ATT_BEST_MATCHING_LENGTH,this.sampleMap.values().stream().mapToInt(M->M.bestLength).max().orElse(0));
			vcb.attribute(VCFConstants.ALLELE_NUMBER_KEY,2);
			vcb.attribute(VCFConstants.ALLELE_COUNT_KEY,1);
			vcb.attribute(VCFConstants.ALLELE_FREQUENCY_KEY,0.5);
			vcb.attribute(VCFConstants.SVTYPE,"DEL");
			vcb.attribute(VCFConstants.END_KEY,chromEnd0+1);
			vcb.attribute("SVLEN",(chromEnd0-chromStart0)+1);
			vcb.log10PError(sum_count/-10.0);
			final List<Genotype> genotypes = new ArrayList<>(this.sampleMap.size());
			for(final String sample:this.sampleMap.keySet()) {
				final PerSample perSample = this.sampleMap.get(sample);
				final GenotypeBuilder gb = new GenotypeBuilder(sample, alleles);/** always het */
				gb.DP(perSample.countSupportingReads);
				gb.attribute(ATT_BEST_MATCHING_LENGTH, perSample.bestLength);
				genotypes.add(gb.make());
				}
			vcb.genotypes(genotypes);
			return vcb.make();
			}
		}
	private void saveGeneInfo(final GeneInfo info) {
		saveGenePw.print(info.gene.getContig());
		saveGenePw.print('\t');
		saveGenePw.print(info.gene.getStart());
		saveGenePw.print('\t');
		saveGenePw.print(info.gene.getEnd());
		saveGenePw.print('\t');
		saveGenePw.print(info.gene.getName());
		saveGenePw.print('\t');
		saveGenePw.print(info.gene.getStrand().encodeAsChar());
		saveGenePw.print('\t');
		saveGenePw.print(info.gene.getExonCount());
		saveGenePw.print('\t');
		saveGenePw.print(info.gene.getIntronCount());//intron
		saveGenePw.print('\t');
		saveGenePw.print(info.intronSet.size());//intron
		saveGenePw.print('\t');
		saveGenePw.print(info.intronSet.size()==info.gene.getIntronCount()?"ALL_INTRONS":".");//intron
		saveGenePw.print('\t');
		saveGenePw.print(info.intronSet.stream().sorted().map(I->I.getContig()+":"+(I.getStart()-1)+"-"+I.getEnd()).collect(Collectors.joining(" ")));//intron
		saveGenePw.println();
		}
	private void reportGene(final KnownGene gene,final Match match) {
		GeneInfo g=this.geneIdToInfo.get(gene.getName());
		if(g==null) {
			g=new  GeneInfo(gene);
			this.geneIdToInfo.put(gene.getName(),g);
			}
		g.intronSet.add(new Interval(gene.getContig(),match.chromStart0+1,match.chromEnd0));
		}
	
	@Override
	public int doWork(final List<String> args) {
		SamReader sr = null;
		VariantContextWriter vcw0=null;
		CloseableIterator<SAMRecord> iter = null;
		final IntervalTreeMap<List<KnownGene>> knownGeneMap = new IntervalTreeMap<>();

		try {
			this.indexedFastaSequenceFile = new IndexedFastaSequenceFile(this.faidx);
			final SAMSequenceDictionary refDict = SequenceDictionaryUtils.extractRequired(this.indexedFastaSequenceFile);
			this.refCtgNameConverter= ContigNameConverter.fromOneDictionary(refDict);

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
					List<KnownGene> L=  knownGeneMap.get(interval);
					if(L==null) {
						L=new ArrayList<KnownGene>();
						knownGeneMap.put(interval,L);
						}
					L.add(kg);
					}
				
				}

			if(knownGeneMap.isEmpty()) {
				LOG.error("no gene found in "+this.knownGeneUri);
				return -1;
				}
			LOG.info("Number of transcripts: "+ knownGeneMap.values().stream().flatMap(L->L.stream()).count());
			
			final List<Match> matchBuffer=new ArrayList<>();
			sr = super.openSamReader(oneFileOrNull(args));
			final SAMFileHeader samFileHeader = sr.getFileHeader();
			if(!samFileHeader.getSortOrder().equals(SAMFileHeader.SortOrder.coordinate)) {
				LOG.error("input is not sorted on coordinate");
				return -1;
			}
			
			if(this.use_bai && !sr.hasIndex()) {
				LOG.warning("Cannot used bai because input is not indexed");
				}
			
			if(this.use_bai && sr.hasIndex())
				{
				LOG.info("building intervals...");
				final SAMSequenceDictionary samdict= SequenceDictionaryUtils.extractRequired(samFileHeader);

				final ContigNameConverter samConvert = ContigNameConverter.fromOneDictionary(samdict);
				final List<QueryInterval> intervalsL = knownGeneMap.values().
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
				LOG.debug("Query bam using "+intervals.length+" random access intervals.");
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
			metaData.add(new VCFInfoHeaderLine(ATT_RETRO_DESC, VCFHeaderLineCount.UNBOUNDED,VCFHeaderLineType.String,
					"Retrocopy attributes: transcript-id|strand|exon-1|exon-2"));

			
			final VCFHeader header=new VCFHeader(metaData, samples);
			header.setSequenceDictionary(refDict);
			
			/* open vcf for writing*/
			vcw0=super.openVariantContextWriter(this.outputFile);
			final VariantContextWriter vcw=vcw0;
			vcw.writeHeader(header);
			
			/* save gene writer */
			if(this.saveGeneTo!=null) {
				this.saveGenePw = super.openFileOrStdoutAsPrintWriter(this.saveGeneTo);
				}
			else
				{
				this.saveGenePw = new PrintWriter(new NullOuputStream());
				}
				
			final ProgressFactory.Watcher<SAMRecord> progress = ProgressFactory.newInstance().dictionary(samFileHeader).logger(LOG).build();
			long count_reads = 0L;
		
			while(iter.hasNext()) {
				final SAMRecord rec = progress.apply(iter.next());
				if(rec.getReadUnmappedFlag()) continue;
				if(rec.isSecondaryOrSupplementary()) continue;
				if(rec.getDuplicateReadFlag()) continue;
				final byte bases[]=rec.getReadBases();
				if(bases==null || bases==SAMRecord.NULL_SEQUENCE) continue;
				final Cigar cigar = rec.getCigar();
				if(cigar==null || cigar.numCigarElements()<2 || !cigar.isClipped()) continue;
				final String refContig = this.refCtgNameConverter.apply(rec.getContig());
				
				if(StringUtils.isBlank(refContig)) continue;
				
				
				
				final CigarElement leftCigar = cigar.getCigarElement(0);
				final CigarElement rightCigar = cigar.getCigarElement(cigar.numCigarElements()-1);
				
				/* both ends are not candidate */
				if(!isCandidateCigarElement.test(leftCigar) && !isCandidateCigarElement.test(rightCigar) ) continue;
				/* get sample */
				final String sampleName = this.partiton.getPartion(rec, null);
				if(StringUtils.isBlank(sampleName)) continue;
				
				final List<KnownGene> genes = knownGeneMap.getOverlapping(
						new Interval(refContig,rec.getUnclippedStart(),rec.getUnclippedEnd())
						).stream().
						flatMap(L->L.stream()).
						collect(Collectors.toList());
				if(genes.isEmpty()) continue;
				
				count_reads++;

				if(this.genomicSequence==null || !this.genomicSequence.getChrom().equals(refContig)) {
					this.genomicSequence = new GenomicSequence(this.indexedFastaSequenceFile, refContig);
					/* dump buffer */
					matchBuffer.stream().map(B->B.build()).forEach(V->vcw.add(V));
					matchBuffer.clear();
					/* dump genes */
					this.geneIdToInfo.values().stream().forEach(G->saveGeneInfo(G));
					this.geneIdToInfo.clear();
					}
				else if(count_reads%1_000_000==0)
					{
					final int min_start=  genes.stream().mapToInt(G->G.getStart()-1_000_000/*arbitrary, useful */).min().orElse(0);
					/* dump buffer. DOesn't work */
					int i=0;
					while(i < matchBuffer.size())
						{
						final Match match=matchBuffer.get(i);
						if((match.chromEnd0+1)+this.minCigarSize < min_start)
							{
							vcw.add(match.build());
							matchBuffer.remove(i);
							}
						else
							{
							i++;
							}
						}
					
					}
				
				/* test each side of the clipped read */
				for(int side=0;side<2;++side) {
					final CigarElement ce_side = (side==0?leftCigar:rightCigar);
					if(!ce_side.getOperator().equals(CigarOperator.S)) continue;
					if(ce_side.getLength()< this.minCigarSize) continue;
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
								
								if(matchLength<this.minCigarSize) continue;
								
								final KnownGene.Intron intron=knownGene.getIntron(exonIndex-1); 
								
								//find match or create new
								Match match = matchBuffer.stream().
										filter(B->B.chromStart0==intron.getStart() && B.chromEnd0==intron.getEnd()).
										findFirst().orElse(null);
								if(match==null)
									{
									//LOG.debug("MEW MATCH0 "+knownGene.getName());
									match = new Match(intron);
									matchBuffer.add(match);
									reportGene(knownGene,match);
									}
								final PerSample perSample = match.getSample(sampleName);
															
								
								match.attributes.add(
										knownGene.getName()+"|"+
									    knownGene.getStrand().encodeAsChar()+"|"+
										exonLeft.getName()+"|"+ StringUtils.right(exonLeft,this.minCigarSize)+"|"+
										StringUtils.left(exonRight,this.minCigarSize)+"|"+exonRight.getName()
										);
								perSample.countSupportingReads++;
								perSample.bestLength=Math.max(perSample.bestLength, matchLength);
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
								
								if(matchLength<this.minCigarSize) continue;
								//find match or create new
								
								final KnownGene.Intron intron=knownGene.getIntron(exonIndex); 
								
								Match match = matchBuffer.stream().filter(B->B.chromStart0==intron.getStart() && B.chromEnd0==intron.getEnd()).findFirst().orElse(null);
								if(match==null)
									{
									//LOG.debug("MEW MATCH1 "+knownGene.getName());
									match = new Match(intron);
									matchBuffer.add(match);
									reportGene(knownGene,match);
									}
								final PerSample perSample = match.getSample(sampleName);
								
								match.attributes.add(
										knownGene.getName()+"|"+
									    knownGene.getStrand().encodeAsChar()+"|"+
									    exonLeft.getName()+"|"+ StringUtils.right(exonLeft,this.minCigarSize)+"|"+
										StringUtils.left(exonRight,this.minCigarSize)+"|"+exonRight.getName()
										);
								
								perSample.countSupportingReads++;
								perSample.bestLength=Math.max(perSample.bestLength, matchLength);
								}
							}
						}
					}
				}
			/* dump buffer */
			matchBuffer.stream().map(B->B.build()).forEach(V->vcw.add(V));
			matchBuffer.clear();
			/* dump gene */
			this.geneIdToInfo.values().stream().forEach(G->saveGeneInfo(G));
			this.geneIdToInfo.clear();
			
			progress.close();
			vcw.close();
			iter.close();
			iter=null;
			sr.close();
			sr=null;
			this.saveGenePw.flush();
			this.saveGenePw.close();
			this.saveGenePw=null;
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
			CloserUtil.close(this.indexedFastaSequenceFile);
			CloserUtil.close(this.saveGenePw);
			}
		}
	
	public static void main(final String[] args) {
		new ScanRetroCopy().instanceMainWithExit(args);
	}
}
