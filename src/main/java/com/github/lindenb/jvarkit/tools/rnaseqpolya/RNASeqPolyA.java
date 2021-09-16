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
import java.util.Objects;
import java.util.Set;
import java.util.function.BiFunction;
import java.util.function.ToIntFunction;
import java.util.function.UnaryOperator;
import java.util.stream.Collectors;

import com.beust.jcommander.DynamicParameter;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;
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
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CoordMath;
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
import htsjdk.variant.vcf.VCFHeaderLineCount;
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
$ find dir1 -type f -name "*.bam" > in.list
$ java -jar rnaseqpolya.jar -p 5 --reference ref.fa --gff3 in.gff  --out out.vcf.gz  in.list 
```

## See also

https://twitter.com/yokofakun/status/1438137720484900869

![twitter](https://pbs.twimg.com/media/E_VKkVJWEAMA6xy?format=png&name=900x900 "Screenshot")

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
	
	@Parameter(names={"-o","-out","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-g","--gff","--gff3"},description="GFF3 file containing the exons.",required=true)
	private Path gffPath = null;
	@Parameter(names= {"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path faidx;
	@Parameter(names= {"-p","--primer"},description="Search for poly-A 'A{x}' dandling part of the read. -1 : ignore and search for poly-A just after the exon boundary . If it's found, we count the number of A starting with this pattern.")
	private int polyA_primer_size = -1;
	@Parameter(names= {"-indels","--indels"},description="Ignore reads containing indels.")
	private boolean ignore_with_indels = false;
	@Parameter(names= {"--disable-index"},description="Disable use of BAM index")
	private boolean disable_bam_index = false;
	@DynamicParameter(names = "-D", description = "extra parameters. Undocumented.",hidden=true)
	private Map<String, String> dynaParams = new HashMap<>();
	@Parameter(names= {"--filter-reads"},description="remove duplicate, supplementary, malformed reads")
	private boolean default_read_filter = false;
	@Parameter(names= {"-d","--duplicate-ends"},description="keep only one transcript if transcripts share the same 3' coordinate")
	private boolean remove_duplicate_transcripts = false;
	@Parameter(names= {"-x","--overlapping-exon"},description="Ignore exon if an exon from another transcript overlaps the end of the last exon.")
	private boolean remove_if_overlapping_exon = false;
	@Parameter(names= {"-C","--contig"},description="limit to this contig/chromosome")
	private String limit_contig=null;

	
	@ParametersDelegate
	private WritingVariantsDelegate writingVariantsDelegate = new WritingVariantsDelegate();
	
	
	private static class ExonCount {
		int n_tested_reads = 0;
		int n_tested_reads_with_A = 0;
		int max_length_polyA=0;
		long sum_polyA=0L;
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
		Set<String> otherIds = null;
		final Map<String,ExonCount> sample2count = new HashMap<>();
		
		boolean isPlusStrand() { return strand.equals(Strand.POSITIVE);}
		boolean isMinusStrand() { return strand.equals(Strand.NEGATIVE);}
		
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
		/* get last position of exon */
		public int getPosition() {
			return isPlusStrand()?getEnd():getStart();
			}
		
		boolean isAfterExon(int ref1) {
			return (
					(this.isPlusStrand() && ref1> getEnd()) ||
					(this.isMinusStrand() && ref1< getStart())
					);
			}
		}
	
	
	
	private RNASeqPolyA()
		{
		}
	
	private void decodeGff3Feature(final Gff3Feature geneFeat,final ContigNameConverter converter,final Map<String,GeneInfo> geneMap) {
		if(geneFeat==null) return;
		final String debugTranscript  = dynaParams.getOrDefault("debug.transcript","");

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
		if(!StringUtils.isBlank(this.limit_contig) && !chrom.equals(this.limit_contig)) return;
		
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
		if(StringUtils.isBlank(gene.geneName)) gene.geneName = getKey.apply(geneFeat, "gene_name");
		gene.biotype = getKey.apply(geneFeat, "biotype");
		if(StringUtils.isBlank(gene.biotype)) gene.biotype = getKey.apply(geneFeat, "gene_type");

		final List<Interval> all_exons = new ArrayList<>();
		
		for(Gff3Feature trFeat:geneFeat.getChildren()) {
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
			
			if(!StringUtils.isBlank(debugTranscript) && !debugTranscript.equals(transcriptId))  {
				continue;
				}
			if(gene.transcripts.containsKey(transcriptId)) throw new IllegalArgumentException("duplicate transcriptId ID "+transcriptId);

			
			for(final Gff3Feature exFeat:trFeat.getChildren())
				{
				if(!exFeat.getType().equals("exon")) {
					continue;
					}
				
				all_exons.add(new Interval(chrom,exFeat.getStart(),exFeat.getEnd()));
				
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
					(lastExon.isPlusStrand() && exFeat.getEnd() > lastExon.getEnd()) ||
					(lastExon.isMinusStrand() && exFeat.getStart() < lastExon.getStart())) {
					lastExon.start = exFeat.getStart();
					lastExon.end = exFeat.getEnd();
					}
				}
			}
		
		if(this.remove_if_overlapping_exon) {
			final Set<String> toRemove = new HashSet<>();
			for(String id:gene.transcripts.keySet()) {
				final LastExon exi = gene.transcripts.get(id);
				if(all_exons.stream().anyMatch(EX->EX.getStart() < exi.getPosition() && EX.getEnd()>exi.getPosition())) {
					toRemove.add(exi.transcriptId);
					}
				}
			for(String id: toRemove) gene.transcripts.remove(id);
			}
		
		if(this.remove_duplicate_transcripts) {
			final List<String> ids = new ArrayList<>(gene.transcripts.keySet());
			final Set<String> toRemove = new HashSet<>();
			for(int i=0; i+1 < ids.size();i++) {
				final LastExon exi = gene.transcripts.get(ids.get(i));
				if(toRemove.contains(exi.transcriptId)) continue;
				for(int j=i+1; j< ids.size();j++) {
					final LastExon exj = gene.transcripts.get(ids.get(j));
					if( (exi.isPlusStrand() && exj.isPlusStrand() && exi.getEnd()==exj.getEnd()) ||
						(exi.isMinusStrand() && exj.isMinusStrand() && exi.getStart()==exj.getStart())) {
						if(exi.otherIds==null) exi.otherIds=new HashSet<>();
						exi.otherIds.add(exj.transcriptId);
						toRemove.add(exj.transcriptId);
						}
					}
				}
			for(String id: toRemove) gene.transcripts.remove(id);
			}

		}
	
	@Override
	public int doWork(final List<String> args) {
		final Map<String,GeneInfo> geneToTrancripts = new HashMap<>();
		try
			{
			final String debugTranscript  = dynaParams.getOrDefault("debug.transcript","");

			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(this.faidx);
			if(!StringUtils.isBlank(this.limit_contig) && dict.getSequence(this.limit_contig)==null) {
				throw new JvarkitException.ContigNotFoundInDictionary(this.limit_contig, dict);
			}
			
			
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
			
			final List<Path> inputs = IOUtils.unrollPaths(args);

			
			final QueryInterval[] intervals = this.disable_bam_index || inputs.isEmpty()?
					null:
					QueryInterval.optimizeIntervals(
					exonMap.values().
						stream().
						map(R->new QueryInterval(toTid.applyAsInt(R.getContig()), R.getPosition(), R.getPosition())).
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
			
			
			
			final Set<String> samples = new HashSet<>();
			int bam_index=0;
			// loop over the bams
			for(;;) {
				final Path bamFilename = inputs.isEmpty()?null:inputs.get(bam_index);
				try(SamReader sr = inputs.isEmpty()?srf.open(SamInputResource.of(stdin())):srf.open(bamFilename)) {
					final SAMFileHeader header0 =  sr.getFileHeader();
					SequenceUtil.assertSequenceDictionariesEqual(dict, SequenceDictionaryUtils.extractRequired(header0));
					final String sample =header0.getReadGroups().stream().map(RG->RG.getSample()).
							filter(S->!StringUtils.isBlank(S)).findFirst().
							orElse(bamFilename==null?"STDIN":IOUtils.getFilenameWithoutCommonSuffixes(bamFilename));
					if(samples.contains(sample)) {
						LOG.error("duplicate sample "+sample);
						return -1;
						}
					samples.add(sample);
					final ProgressFactory.Watcher<SAMRecord> progress = ProgressFactory.newInstance().dictionary(dict).build();
					try(CloseableIterator<SAMRecord> iter= (intervals==null || inputs.isEmpty()/*stdin*/?sr.iterator():sr.query(intervals, false))) {
						while(iter.hasNext()) {
							final SAMRecord rec = progress.apply(iter.next());
							if(rec.getReadUnmappedFlag()) continue;
							if(!StringUtils.isBlank(this.limit_contig) && !rec.getContig().equals(this.limit_contig)) continue;
							if(this.default_read_filter && !SAMRecordDefaultFilter.accept(rec)) continue;

							final Collection<LastExon> lastExons = exonMap.getOverlapping(rec);
							if(lastExons.isEmpty()) continue;
							final Cigar cigar = rec.getCigar();
							if(cigar==null || cigar.isEmpty()) continue;
							final byte[] bases = rec.getReadBases();
							if(bases==null || SAMRecord.NULL_SEQUENCE.equals(bases)) continue;
							for(LastExon exon : lastExons) {
								if(!StringUtils.isBlank(debugTranscript) && !exon.transcriptId.equals(debugTranscript)) {
									continue; 
								}
								ExonCount count = exon.sample2count.get(sample);
								if(count==null) {
									count = new ExonCount();
									exon.sample2count.put(sample,count);
									}
								
								final StringBuilder sb = new StringBuilder() ;
								boolean indel_flag = false;
								boolean last_exon_in_intron_flag = false;
								boolean match_last_base = false;
								int ref1 = rec.getUnclippedStart();
								int read0 = 0;
								for(CigarElement ce:cigar) {
									if(exon.isMinusStrand() && ref1> exon.start) break;
									if(this.ignore_with_indels && indel_flag) break;
									final CigarOperator op = ce.getOperator();
									switch(op) {
										case P: break;
										case I:
											indel_flag  = true;
											for(int i=0;i< ce.getLength();i++) {
												if(exon.isAfterExon(ref1)) {
													sb.append((char)Character.toUpperCase(bases[read0]));
													}
												read0++;
												}									
											break;
										case D:case N:
											if( (exon.isPlusStrand()  && CoordMath.overlaps(ref1, ref1+ce.getLength()-1,exon.getEnd(), exon.getEnd()+1)) ||
												(exon.isMinusStrand() && CoordMath.overlaps(ref1, ref1+ce.getLength()-1,exon.getStart()-1, exon.getStart()))
												) {
												last_exon_in_intron_flag=true;
												}
											
											ref1 +=ce.getLength();
											indel_flag = true;
											break;
										case H:
											for(int i=0;i< ce.getLength();i++) {
												if(exon.isAfterExon(ref1)) {
													sb.append('N');
													}
												if(ref1==exon.getPosition()) match_last_base = true;
												ref1++;
												}
											break;
										case S: case M: case X: case EQ:
											for(int i=0;i< ce.getLength();i++) {
												if(exon.isAfterExon(ref1)) {
													sb.append((char)Character.toUpperCase(bases[read0]));
													}
												if(ref1==exon.getPosition()) match_last_base = true;
												read0++;
												ref1++;
												}
											break;
										default:throw new IllegalStateException(op.name());
										}
									} //end loop cigar
								// premature end or start
								if(!match_last_base ||
								   (exon.isPlusStrand() && ref1 < exon.getEnd()) ||
								   (exon.isMinusStrand() && ref1 < exon.getStart()) ||
								   (this.ignore_with_indels && indel_flag) || 
								   last_exon_in_intron_flag) {
									continue;
								}
								//if(read0!=bases.length) throw new IllegalStateException("read0:"+read0+" expected "+bases.length+" in "+rec.getReadName());
								//if(ref1!=1+rec.getUnclippedEnd())throw new IllegalStateException("ref1:"+ref1+" expected 1+"+rec.getUnclippedEnd()+" in "+rec.getReadName());

								++count.n_tested_reads;
								String polyA;
								if(exon.isMinusStrand()) {
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
								if(count_polyA>0) {
									count.n_tested_reads_with_A++;
									count.sum_polyA += count_polyA;
									}
								if(count_polyA>count.max_length_polyA) {
									count.max_length_polyA = count_polyA;
									}
								}//end of loop last exon
							}
						progress.close();
						}
					}
				++bam_index;
				if(inputs.isEmpty() || bam_index>=inputs.size()) break;
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
			final VCFInfoHeaderLine infoEndPos = new VCFInfoHeaderLine("POS3",1,VCFHeaderLineType.Integer,"End 3 prime position");
			metaData.add(infoEndPos);
			final VCFInfoHeaderLine infoGeneName = new VCFInfoHeaderLine("GENE_NAME",1,VCFHeaderLineType.String,"Gene Name");
			metaData.add(infoGeneName);
			final VCFInfoHeaderLine infoBiotype = new VCFInfoHeaderLine("GENE_BIOTYPE",1,VCFHeaderLineType.String,"Gene Biotype");
			metaData.add(infoBiotype);
			final int n_last_exon_bases = Math.max(0, Integer.parseInt(dynaParams.getOrDefault("last.n.exons","10")));
			final VCFInfoHeaderLine infoLastExonBases = new VCFInfoHeaderLine("LAST_BASES",1,VCFHeaderLineType.String,"Last exon bases N="+ n_last_exon_bases+". Reverse-complemented for negative strand.");
			metaData.add(infoLastExonBases);
			final int n_after_exon_bases = Math.max(0, Integer.parseInt(dynaParams.getOrDefault("after.n.exons","10")));
			final VCFInfoHeaderLine infoAfterExonBases = new VCFInfoHeaderLine("AFTER_BASES",1,VCFHeaderLineType.String,"Bases after exon N="+ n_after_exon_bases+". Reverse-complemented for negative strand.");
			metaData.add(infoAfterExonBases);
			final VCFInfoHeaderLine infoOtherIdss = new VCFInfoHeaderLine("OTHER_IDS",VCFHeaderLineCount.UNBOUNDED,VCFHeaderLineType.String,"Other transcripts ending at the same coordinate.");
			metaData.add(infoOtherIdss);
			

			final VCFFormatHeaderLine fmtMaxPolyA = new VCFFormatHeaderLine("MAX",1,VCFHeaderLineType.Integer,"Max poly A");
			metaData.add(fmtMaxPolyA);
			final VCFFormatHeaderLine fmtReadPolyA = new VCFFormatHeaderLine("DPA",1,VCFHeaderLineType.Integer,"Read with at least one A");
			metaData.add(fmtReadPolyA);
			final VCFFormatHeaderLine fmtAveragePolyA = new VCFFormatHeaderLine("AVG",1,VCFHeaderLineType.Float,"average length of poly-A for reads carrying at least one A.");
			metaData.add(fmtAveragePolyA);

			

			
			VCFStandardHeaderLines.addStandardInfoLines(metaData, true, VCFConstants.DEPTH_KEY,VCFConstants.END_KEY);
			VCFStandardHeaderLines.addStandardFormatLines(metaData, true, VCFConstants.DEPTH_KEY);
			
			final VCFHeader header = new VCFHeader(metaData,samples.stream().sorted().collect(Collectors.toList()));
			header.setSequenceDictionary(dict);
			JVarkitVersion.getInstance().addMetaData(this, header);
			
			final UnaryOperator<String> afterColon= S->{
				if(!(S.startsWith("gene:") || S.startsWith("transcript:"))) return S;
				int colon = S.indexOf(":");
				return S.substring(colon+1);
			};
			
			final List<Allele> ALLELES = Collections.singletonList(Allele.create("N",true));
			try(VariantContextWriter w=writingVariantsDelegate.dictionary(dict).open(this.outputFile);
				ReferenceSequenceFile fai=ReferenceSequenceFileFactory.getReferenceSequenceFile(this.faidx)	
				) {
				w.writeHeader(header);
				exonMap.values().stream().sorted(new ContigDictComparator(dict).createLocatableComparator()).forEach(T->{
					if(T.getDP()==0) return;
					
					final String lastBases;
					final String afterBases;
					
					if(T.isPlusStrand()) {
						lastBases = fai.getSubsequenceAt(T.getContig(), Math.max(T.getStart(),T.getEnd()-n_last_exon_bases),T.getEnd()).getBaseString();
						final SAMSequenceRecord ssr = Objects.requireNonNull(dict.getSequence(T.getContig()));
						afterBases = fai.getSubsequenceAt(T.getContig(), T.getEnd()+1,Math.min(T.getEnd()+n_after_exon_bases,ssr.getLengthOnReference())).getBaseString();
						}
					else if(T.isMinusStrand()) {
						lastBases  = AcidNucleics.reverseComplement(fai.getSubsequenceAt(T.getContig(), T.getStart(),Math.min(T.getEnd(),T.getStart()+n_last_exon_bases)).getBaseString());
						afterBases = AcidNucleics.reverseComplement(fai.getSubsequenceAt(T.getContig(), Math.max(1,T.getStart()-n_after_exon_bases), T.getStart()-1).getBaseString());
						}
					else
						{
						lastBases = null;
						afterBases = null;
						}
					
					final VariantContextBuilder vcb = new VariantContextBuilder();
					vcb.chr(T.getContig());
					vcb.start(T.getStart());
					vcb.stop(T.getEnd());
					vcb.id(afterColon.apply(T.transcriptId));
					vcb.attribute(VCFConstants.END_KEY, T.getEnd());
					vcb.attribute(infoGeneId.getID(), afterColon.apply(T.gene.geneId));
					vcb.attribute(infoTranscriptId.getID(),  afterColon.apply(T.transcriptId));
					vcb.attribute(infoStrand.getID(), T.strand.name());
					vcb.attribute(infoEndPos.getID(), T.getPosition());
					if(T.otherIds!=null && !T.otherIds.isEmpty()) {
						vcb.attribute(infoOtherIdss.getID(), T.otherIds.stream().map(afterColon).collect(Collectors.toList()));
						}
					
					if(!StringUtils.isBlank(lastBases)) {
						vcb.attribute(infoLastExonBases.getID(), lastBases);
						}
					if(!StringUtils.isBlank(afterBases)) {
						vcb.attribute(infoAfterExonBases.getID(), afterBases);
						}
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
						gb.attribute(fmtAveragePolyA.getID(), count==null || count.n_tested_reads_with_A==0?0f:count.sum_polyA/(float)count.n_tested_reads_with_A);
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
