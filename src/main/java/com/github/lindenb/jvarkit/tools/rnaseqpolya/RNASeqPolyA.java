/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.rnaseqpolya;

import java.io.InputStream;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.BiFunction;
import java.util.function.ToIntFunction;
import java.util.function.UnaryOperator;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.gff.Gff3Codec;
import htsjdk.tribble.gff.Gff3Feature;
import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

/**
BEGIN_DOC

## Motivation

find evidence of poly-A tails in RNASeq data.

## Input

input is a set of path to indexed BAM/CRAM files or a file with the `.list` suffix containing the path to the BAM/CRAM files (one per line)

## Output

ouput is a VCF file. Each variant is a transcript.

## Example

```bash

```

END_DOC
*/
@Program(name="rnaseqpolya",
	description="find poly-A tail in RNASeq data",
	keywords={"bam","sam","rnaseq","polya"},
	creationDate="20210913",
	modificationDate="20210914"
	)
public class RNASeqPolyA extends Launcher {
	private static final Logger LOG = Logger.build(RNASeqPolyA.class).make();
	
	@Parameter(names={"-out","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-g","--gff","--gff3"},description="GFF3 file containing the exons.",required=true)
	private Path gffPath = null;
	@Parameter(names= {"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path faidx;
	@Parameter(names= {"-p","--primer"},description="Search for poly-A 'A{x}' dandling part of the read. -1 : ignore and search for poly-A just after the exon boundary . If it's found, we count the number of A starting with this pattern.")
	private int polyA_primer_size = -1;
	@Parameter(names= {"-indels","--indels"},description="Ignore reads with indels.")
	private boolean ignore_with_indels = false;
	@Parameter(names= {"--disable-index"},description="Disable use of BAM index")
	private boolean disable_bam_index = false;

	@ParametersDelegate
	private WritingVariantsDelegate writingVariantsDelegate = new WritingVariantsDelegate();
	
	
	private static class ExonCount {
		int n_tested_reads = 0;
		int n_tested_reads_with_A = 0;
		int max_length_polyA=0;
	}
	
	private static class GeneInfo {
		String contig;
		String geneId;
		String geneName = null;
		String biotype = null;
		final Map<String,LastExon> transcripts = new HashMap<>();
		
		public int getmaxPolyA() {
			return transcripts.values().stream().mapToInt(T->T.getMaxPolyA()).max().orElse(0);
			}
	}
	
	private static class LastExon implements Locatable {
		GeneInfo gene;
		int start;
		int end;
		Strand strand = Strand.NONE;
		String transcriptId;
		final Map<String,ExonCount> sample2count = new HashMap<>();
		@Override
		public String getContig() {
			return gene.contig;
			}
		@Override
		public int getStart() {
			return start;
			}
		@Override
		public int getEnd() {
			return end;
			}
		public int getDP() {
			return this.sample2count.values().stream().mapToInt(C->C.n_tested_reads).sum();
			}

		public int getMaxPolyA() {
			return this.sample2count.values().stream().mapToInt(C->C.max_length_polyA).max().orElse(0);
			}
		}
	
	
	
	private RNASeqPolyA()
		{
		}
	
	private void decodeGff3Feature(final Gff3Feature geneFeat,final ContigNameConverter converter,final Map<String,GeneInfo> geneMap) {
		if(geneFeat==null) return;
		final BiFunction<Gff3Feature, String,String> getKey =  (FEAT,KEY) -> {
			List<String> atts = FEAT.getAttribute(KEY);
			return atts==null || atts.size()!=1?null:atts.get(0);
			};
		
		if(!geneFeat.getType().equals("gene"))  {
			for(Gff3Feature other:geneFeat.getChildren()) {
				decodeGff3Feature(other,converter,geneMap);
				}
			return;
			}
			
		final String chrom = converter.apply(geneFeat.getContig());
		if(StringUtils.isBlank(chrom)) return;
		
		final String geneID = geneFeat.getID();
		if(StringUtils.isBlank(geneID)) {
			LOG.warn("not ID for "+geneFeat);
			return;
		}
		if(geneMap.containsKey(geneID)) throw new IllegalArgumentException("duplicate gene ID "+geneID);
		final GeneInfo gene = new GeneInfo();
		geneMap.put(geneID, gene);
		gene.geneId = geneID;
		gene.contig = chrom;
		gene.geneName = getKey.apply(geneFeat, "Name");
		gene.biotype = getKey.apply(geneFeat, "biotype");
		
		
		for(Gff3Feature trFeat:geneFeat.getChildren())
			{
			/* ignore this , there is plenty of type under gene 
			if(!(trFeat.getType().equals("transcript") || trFeat.getType().equals("mRNA"))) {
				continue;
				}
			*/
			final String transcriptId = trFeat.getID();
			if(StringUtils.isBlank(transcriptId))  {
				LOG.warn("not ID for "+transcriptId);
				continue;
			}
			if(gene.transcripts.containsKey(transcriptId)) throw new IllegalArgumentException("duplicate transcriptId ID "+transcriptId);

			
			for(final Gff3Feature exFeat:trFeat.getChildren())
				{
				if(!exFeat.getType().equals("exon")) {
					continue;
				}
				
				LastExon lastExon  = gene.transcripts.get(transcriptId);
				if(lastExon==null) {
					lastExon = new LastExon();
					lastExon.gene = gene;
					lastExon.transcriptId = transcriptId;
					lastExon.strand = exFeat.getStrand();
					lastExon.start = exFeat.getStart();
					lastExon.end = exFeat.getEnd();
					gene.transcripts.put(transcriptId,lastExon);
					}
				else if(!lastExon.strand.equals(exFeat.getStrand())) {
					throw new IllegalArgumentException("conflict strand for "+transcriptId);
					}
				if(
					(lastExon.strand.equals(Strand.POSITIVE) && exFeat.getEnd() > lastExon.getEnd()) ||
					(lastExon.strand.equals(Strand.NEGATIVE) && exFeat.getStart() < lastExon.getStart())) {
					lastExon.start = exFeat.getStart();
					lastExon.end = exFeat.getEnd();
					}
				}
			}
		
			
			
	}
	
	@Override
	public int doWork(final List<String> args) {
		final Map<String,GeneInfo> geneToTrancripts = new HashMap<>();
		try
			{
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(this.faidx);
			final ContigNameConverter converter = ContigNameConverter.fromOneDictionary(dict);
			
			final Gff3Codec gff3 = new Gff3Codec(Gff3Codec.DecodeDepth.DEEP);
			try(InputStream is = IOUtils.openPathForReading(this.gffPath)) {
				final AsciiLineReader asciiLineReader = AsciiLineReader.from(is);
				final LineIterator lr= new LineIteratorImpl(asciiLineReader);
				while(!gff3.isDone(lr)) {
					decodeGff3Feature(gff3.decode(lr),converter,geneToTrancripts);
					}
				gff3.close(lr);
				asciiLineReader.close();
			}
			
			if(geneToTrancripts.isEmpty()) {
				LOG.warn("no transcript was found.");
				//continue , empty VCF must be produced
				}
			
			final IntervalTreeMap<LastExon> exonMap = new IntervalTreeMap<>();
			//fill exonMap
			geneToTrancripts.values().
				stream().
				flatMap(K->K.transcripts.values().stream()).
				forEach(X-> exonMap.put(new Interval(X), X));
			
			final ToIntFunction<String> toTid = C->{
				final SAMSequenceRecord ssr = dict.getSequence(C);
				if(ssr==null) throw new JvarkitException.ContigNotFoundInDictionary(C, dict);
				return ssr.getSequenceIndex();
			};
			
			final QueryInterval[] intervals = this.disable_bam_index?
					null:
					QueryInterval.optimizeIntervals(
					exonMap.keySet().
						stream().
						map(R->new QueryInterval(toTid.applyAsInt(R.getContig()), R.getStart(), R.getEnd())).
						toArray(N->new QueryInterval[N])
					);
			
			final SamReaderFactory srf = super.createSamReaderFactory().
					referenceSequence(this.faidx);
			
			final String primerAAA;
			if( polyA_primer_size>0) {
				primerAAA = StringUtils.repeat(this.polyA_primer_size, 'A');
			} else
			{
				primerAAA = null;
			}
			
			
			final List<Path> inputs = IOUtils.unrollPaths(args);
			final Set<String> samples = new HashSet<>();
			for(final Path bamPath: inputs) {
				try(SamReader sr = srf.open(bamPath)) {
					final SAMFileHeader header0 =  sr.getFileHeader();
					SequenceUtil.assertSequenceDictionariesEqual(dict, SequenceDictionaryUtils.extractRequired(header0));
					final String sample =header0.getReadGroups().stream().map(RG->RG.getSample()).filter(S->!StringUtils.isBlank(S)).findFirst().orElse(IOUtils.getFilenameWithoutCommonSuffixes(bamPath));
					if(samples.contains(sample)) {
						LOG.error("duplicate sample "+sample);
						return -1;
						}
					samples.add(sample);
					final ProgressFactory.Watcher<SAMRecord> progress = ProgressFactory.newInstance().dictionary(dict).build();
					try(CloseableIterator<SAMRecord> iter= (intervals==null?sr.iterator():sr.query(intervals, false))) {
						while(iter.hasNext()) {
							final SAMRecord rec = progress.apply(iter.next());
							if(rec.getReadUnmappedFlag()) continue;
							final Collection<LastExon> lastExons = exonMap.getOverlapping(rec);
							if(lastExons.isEmpty()) continue;
							final Cigar cigar = rec.getCigar();
							if(cigar==null || cigar.isEmpty()) continue;
							final byte[] bases = rec.getReadBases();
							if(bases==null || SAMRecord.NULL_SEQUENCE.equals(bases)) continue;
							for(LastExon exon : lastExons) {
								ExonCount count = exon.sample2count.get(sample);
								if(count==null) {
									count = new ExonCount();
									exon.sample2count.put(sample,count);
									}
								
								StringBuilder sb = new StringBuilder() ;
								boolean indel_flag = false;
								int ref1 = rec.getUnclippedStart();
								int read0 = 0;
								for(CigarElement ce:cigar) {
									if(exon.strand.equals(Strand.NEGATIVE) && ref1> exon.start) break;
									if(this.ignore_with_indels && indel_flag) break;
									final CigarOperator op = ce.getOperator();
									switch(op) {
										case P: break;
										case I: read0+=ce.getLength(); indel_flag  = true; break;
										case D:case N: ref1 +=ce.getLength(); indel_flag = true; break;
										case H:
											for(int i=0;i< ce.getLength();i++) {
												if(exon.strand.equals(Strand.POSITIVE) && ref1> exon.end) {
													sb.append('N');
													}
												else if(exon.strand.equals(Strand.NEGATIVE) && ref1< exon.start) {
													sb.append('N');
													}
												read0++;
												ref1++;
												}
											break;
										case S: case M: case X: case EQ:
											for(int i=0;i< ce.getLength();i++) {
												if(exon.strand.equals(Strand.POSITIVE) && ref1> exon.end) {
													sb.append((char)Character.toUpperCase(bases[read0]));
													}
												else if(exon.strand.equals(Strand.NEGATIVE) && ref1< exon.start) {
													sb.append((char)Character.toUpperCase(bases[read0]));
													}
												read0++;
												ref1++;
												}
											break;
										default:throw new IllegalStateException();
										}
									}
								// premature end or start
								if((exon.strand.equals(Strand.POSITIVE) && ref1 < exon.end) ||
								   (exon.strand.equals(Strand.NEGATIVE) && ref1 < exon.start) ||
								   (this.ignore_with_indels && indel_flag)) {
									continue;
								}
								
								++count.n_tested_reads;
								String polyA;
								if(exon.strand.equals(Strand.NEGATIVE)) {
									polyA = AcidNucleics.reverseComplement(sb);
									}
								else
									{
									polyA = sb.toString();
									}
								
								if(primerAAA!=null) {
									final int pos = polyA.indexOf(primerAAA);
									if(pos>0) polyA = polyA.substring(pos);
								}
								
								int count_polyA = 0;
								for(int i=0;i< polyA.length();i++) {
									if(polyA.charAt(i)!='A') break;
									count_polyA++;
									}
								if(count_polyA>0) count.n_tested_reads_with_A++;
								count.max_length_polyA = Math.max(count.max_length_polyA, count_polyA);
								}
							}
						progress.close();
						}
					}
				}
				
			final Set<VCFHeaderLine> metaData = new HashSet<>();
			final VCFInfoHeaderLine infoGeneId = new VCFInfoHeaderLine("GENE",1,VCFHeaderLineType.String,"Gene ID in "+this.gffPath);
			metaData.add(infoGeneId);
			final VCFInfoHeaderLine infoTranscriptId = new VCFInfoHeaderLine("TRANSCRIPT",1,VCFHeaderLineType.String,"Transcript ID in "+this.gffPath);
			metaData.add(infoTranscriptId);
			final VCFInfoHeaderLine infoStrand = new VCFInfoHeaderLine("STRAND",1,VCFHeaderLineType.String,"Strand");
			metaData.add(infoStrand);
			final VCFInfoHeaderLine infoTranscriptMaxPolyA = new VCFInfoHeaderLine("TRANSCRIPT_MAX",1,VCFHeaderLineType.Integer,"Max poly A in Transcript");
			metaData.add(infoTranscriptMaxPolyA);
			final VCFInfoHeaderLine infoGeneMaxPolyA = new VCFInfoHeaderLine("GENE_MAX",1,VCFHeaderLineType.Integer,"Max poly A in Gene");
			metaData.add(infoGeneMaxPolyA);
			final VCFInfoHeaderLine infoGeneName = new VCFInfoHeaderLine("GENE_NAME",1,VCFHeaderLineType.String,"Gene Name");
			metaData.add(infoGeneName);
			final VCFInfoHeaderLine infoBiotype = new VCFInfoHeaderLine("GENE_BIOTYPE",1,VCFHeaderLineType.String,"Gene Biotype");
			metaData.add(infoBiotype);

			final VCFFormatHeaderLine fmtMaxPolyA = new VCFFormatHeaderLine("MAX",1,VCFHeaderLineType.Integer,"Max poly A");
			metaData.add(fmtMaxPolyA);
			final VCFFormatHeaderLine fmtReadPolyA = new VCFFormatHeaderLine("DPA",1,VCFHeaderLineType.Integer,"Read with at least one A");
			metaData.add(fmtReadPolyA);

			
			VCFStandardHeaderLines.addStandardInfoLines(metaData, true, VCFConstants.DEPTH_KEY,VCFConstants.END_KEY);
			VCFStandardHeaderLines.addStandardFormatLines(metaData, true, VCFConstants.DEPTH_KEY);
			
			final VCFHeader header = new VCFHeader(metaData,samples.stream().sorted().collect(Collectors.toList()));
			header.setSequenceDictionary(dict);
			JVarkitVersion.getInstance().addMetaData(this, header);
			
			final UnaryOperator<String> afterColon= S->{
				int colon = S.indexOf(":");
				return colon==-1?S:S.substring(colon+1);
			};
			
			final List<Allele> ALLELES = Collections.singletonList(Allele.create("N",true));
			try(VariantContextWriter w=writingVariantsDelegate.dictionary(dict).open(this.outputFile)) {
				w.writeHeader(header);
				exonMap.values().stream().sorted(new ContigDictComparator(dict).createLocatableComparator()).forEach(T->{
					if(T.getDP()==0) return;
					
					final VariantContextBuilder vcb = new VariantContextBuilder();
					vcb.chr(T.getContig());
					vcb.start(T.getStart());
					vcb.stop(T.getEnd());
					vcb.id(afterColon.apply(T.transcriptId));
					vcb.attribute(VCFConstants.END_KEY, T.getEnd());
					vcb.attribute(infoGeneId.getID(), afterColon.apply(T.gene.geneId));
					vcb.attribute(infoTranscriptId.getID(),  afterColon.apply(T.transcriptId));
					vcb.attribute(infoStrand.getID(), T.strand.name());
					if(!StringUtils.isBlank(T.gene.geneName)) {
						vcb.attribute(infoGeneName.getID(),T.gene.geneName);
						}
					if(!StringUtils.isBlank(T.gene.biotype)) {
						vcb.attribute(infoBiotype.getID(),T.gene.biotype);
						}
					
					final List<Genotype> genotypes = new ArrayList<>(samples.size());
					for(String sn: samples) {
						final ExonCount count = T.sample2count.get(sn);
						final GenotypeBuilder gb =new GenotypeBuilder(sn);
						gb.attribute(fmtMaxPolyA.getID(),count==null?0:count.max_length_polyA);
						gb.attribute(fmtReadPolyA.getID(),count==null?0:count.n_tested_reads_with_A);
						gb.DP(count==null?0:count.n_tested_reads);
						genotypes.add(gb.make());
						}
					vcb.alleles(ALLELES);
					vcb.genotypes(genotypes);
					
					vcb.attribute(
						VCFConstants.DEPTH_KEY,
						T.getDP()
						);

					final int score = T.getMaxPolyA();
					if(score>0) vcb.log10PError(score/-10.0);
					
					vcb.attribute(
						infoTranscriptMaxPolyA.getID(),
						score
						);
					
					vcb.attribute(
						infoGeneMaxPolyA.getID(),
						T.gene.getmaxPolyA()
						);

					
					w.add(vcb.make());
				});
			}
				
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			}
		}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new RNASeqPolyA().instanceMainWithExit(args);

	}

}
