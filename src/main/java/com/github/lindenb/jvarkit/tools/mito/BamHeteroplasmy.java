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
package com.github.lindenb.jvarkit.tools.mito;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.TreeSet;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.stats.FisherExactTest;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SequenceUtil;
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

Experimental.

Input is one or more indexed BAM/CRAM files or a file with the suffix '.list' containing the path to the BAMs.


END_DOC
*/
@Program(name="bamheteroplasmy",
description="Call a VCF for the Mitochondria (Experimental)",
creationDate="20190910",
modificationDate="20190912",
keywords={"vcf","sam","bam","mitochondria"}
)
public class BamHeteroplasmy extends Launcher {
	private static final Logger LOG = Logger.build(BamHeteroplasmy.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names = { "-R", "--reference" }, description = INDEXED_FASTA_REFERENCE_DESCRIPTION, required = true)
	private Path faidx=null;
	@Parameter(names = { "-partition", "--partition" }, description = SAMRecordPartition.OPT_DESC)
	private SAMRecordPartition partition = SAMRecordPartition.sample;
	@Parameter(names = { "--organelle" }, description = "Organelle name; if empty, the program tries to find the best name in the REF")
	private String organelleName = null;
	@Parameter(names = { "--supplementary" }, description = "accept supplementary alignments.")
	private boolean acceptSupplementary = false;
	@Parameter(names = { "--secondary" }, description = "accept secondary alignments.")
	private boolean acceptSecondary = false;
	@Parameter(names = { "--discordant" }, description = "accept discordant alignments (mate mapping another contig)")
	private boolean acceptDiscordant = false;
	@Parameter(names = { "--mate-unmapped" }, description = "accept mate unmapped")
	private boolean acceptMateUnmapped = false;
	@Parameter(names = { "-sa","--sa" }, description = "accept read having 'SA:' supplementary alignments mapping another contig")
	private boolean acceptSA = false;
	@Parameter(names = { "-Q","--mapq" }, description = "min mapping quality")
	private int mapq = 30;
	@Parameter(names = { "-B","--baq" }, description = "min base quality")
	private int baseq = 0;
	@Parameter(names = { "--all-alleles" }, description = "when calculating the depth, use all alleles, not just the major/minor")
	private boolean use_all_alleles = false;
	@Parameter(names = { "--min-allele-dp" }, description = "min-dp per alt allele")
	private int min_allele_dp = 20;
	@Parameter(names = { "--max-clipping" }, description = "max clip per read")
	private int max_clip = 2;
	@Parameter(names = { "--fisher-treshold" }, description = "set filter if fisher test for strand bias is lower than 'x'.")
	private double fisher_strand_bias = 0.01;

	
	
	private ReferenceSequence chrMSequence = null;
	private static final char ATGCNatgcn[]=new char[] {'A','T','G','C','N','a','t','g','c','n'};
	private static final int ACGT[]=new int[] {'A','C','G','T'};
	

	
	private  static class Pileup {
		final Map<Character,Integer> counter = new HashMap<>(ATGCNatgcn.length);
		
		
		int count(char c) {
			return counter.getOrDefault(c, 0);
		}
		
		int countPlus(char ch) {
			return  count(Character.toUpperCase(ch));
			}
		
		int countMinus(char ch) {
			return  count(Character.toLowerCase(ch));
			}
		int countPlus(Allele ch) {
			return  countPlus(ch.getBaseString().charAt(0));
			}
		
		int countMinus(Allele ch) {
			return  countPlus(ch.getBaseString().charAt(0));
			}
		int countIgnoreCase(char ch) {
			return  countPlus(ch)+ countMinus(ch);
			}
		
		int count(Allele c) {
			return count(c.getBaseString().charAt(0));
			}
		
		int countIgnoreCase(Allele c) {
			return countIgnoreCase(c.getBaseString().charAt(0));
		}
	}
	
	private class SampleHeteroplasmy
		{
		final String sn;
		final Map<Integer,Pileup> pos1ToPileups;
	
		SampleHeteroplasmy(final String sn) {
			this.sn = sn;
			this.pos1ToPileups = new HashMap<>(BamHeteroplasmy.this.chrMSequence.length());
			}
		void visit(int pos1,boolean negativeStrand,byte base) {
			Pileup p = this.pos1ToPileups.get(pos1);
			if(p==null) {
				p = new Pileup();
				this.pos1ToPileups.put(pos1,p);
				}
			final char c =  negativeStrand ?
				(char)Character.toLowerCase(base):
				(char)Character.toUpperCase(base);
			p.counter.put(c, p.counter.getOrDefault(c,0)+1);
			}
		}

	
	@Override
	public int doWork(List<String> args) {
		VariantContextWriter vcw = null;
		ReferenceSequenceFile referenceSequenceFile = null;
		try {
			referenceSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(this.faidx);
			final SamReaderFactory srf = super.createSamReaderFactory();
			if(this.faidx!=null) srf.referenceSequence(this.faidx);
			final List<Path> all_bams = IOUtils.unrollPaths(args);
			
			final SAMSequenceDictionary dict1 = SequenceDictionaryUtils.extractRequired(this.faidx);
			final SAMSequenceRecord ssr1 ;
			
			if(StringUtils.isBlank(this.organelleName)) {
				ssr1 = dict1.getSequences().stream().filter(SSR->SSR.getSequenceName().matches("(chr)?(M|MT)")).findFirst().orElseThrow(()->new JvarkitException.ContigNotFoundInDictionary("mitochondrial chromosome", dict1));
				}
			else
				{
				ssr1 = dict1.getSequence(this.organelleName);
				if(ssr1==null) throw new JvarkitException.ContigNotFoundInDictionary(this.organelleName, dict1);
				}
			
			this.chrMSequence = Objects.requireNonNull(referenceSequenceFile.getSequence(ssr1.getSequenceName()));
			
			
			final Map<String,SampleHeteroplasmy> sample2heteroplasmy= new HashMap<>();
			
			
			for(final Path bam1File:all_bams) {
				final SamReader bam1 = srf.open(bam1File);
				final SAMFileHeader header = bam1.getFileHeader();
				if(header.getReadGroups().
						stream().
						map(RG->this.partition.apply(RG)).
						filter(S->!StringUtils.isBlank(S)).count()==0)
					{
					LOG.warn("no Read Group in "+bam1File+" for partiton:"+this.partition);
					}
				SequenceUtil.assertSequenceDictionariesEqual(SequenceDictionaryUtils.extractRequired(header), dict1);
				
				final SAMRecordIterator iter=bam1.queryOverlapping(ssr1.getSequenceName(), 1, ssr1.getSequenceLength());
				while(iter.hasNext()) {
					final SAMRecord rec = iter.next();
					if(rec.getReadUnmappedFlag()) continue;
					if(rec.getMappingQuality()< this.mapq) continue;
					if(!this.acceptSupplementary && rec.getSupplementaryAlignmentFlag()) continue;
					if(!this.acceptSecondary && rec.getSecondOfPairFlag()) continue;
					if(!this.acceptDiscordant && rec.getReadPairedFlag() && !rec.getMateUnmappedFlag() && !rec.getMateReferenceName().equals(ssr1.getSequenceName())) continue;
					if(!this.acceptMateUnmapped && rec.getReadPairedFlag() && rec.getMateUnmappedFlag() ) continue;
					if(!this.acceptSA && 
						SAMUtils.getOtherCanonicalAlignments(rec).stream().
							anyMatch(R->!R.getContig().equals(ssr1.getSequenceName()))
						) continue;
						
					final String sampleName = this.partition.getPartion(rec,bam1File.getFileName().toString());
					SampleHeteroplasmy sampleData = sample2heteroplasmy.get(sampleName);
					if(sampleData==null) {
						sampleData = new SampleHeteroplasmy(sampleName);
						sample2heteroplasmy.put(sampleName,sampleData);
					}
					
					final Cigar cigar = rec.getCigar();
					if(cigar==null || cigar.isEmpty()) continue;
					if(cigar.numCigarElements()>1 &&
						cigar.getCigarElements().stream().filter(CE->CE.getOperator().isClipping()).mapToInt(CE->CE.getLength()).sum()> this.max_clip) continue;
					final byte bases[]=rec.getReadBases();
					if(bases==SAMRecord.NULL_SEQUENCE) continue;
					
					final byte qualities[]=rec.getBaseQualities();
					
					int refpos = rec.getUnclippedStart();
					int readpos=0;
					for(final CigarElement ce:cigar) {
						final CigarOperator op =ce.getOperator();
						if(op.consumesReferenceBases()) {
							if(op.consumesReadBases())
								{
								for(int i=0;i< ce.getLength();i++) {
									final int pos1 = refpos+i;
									if(pos1<1 || pos1> this.chrMSequence.length()) continue;
									if(qualities!=SAMRecord.NULL_QUALS && qualities[readpos+i] < this.baseq) continue;
									
									sampleData.visit(pos1,rec.getReadNegativeStrandFlag(),bases[readpos+i]);
									}
								}
							refpos+=ce.getLength();
							}
						if(op.consumesReadBases()) {
							readpos+=ce.getLength();
							}
						}
					}
				iter.close();
				bam1.close();
				}
			
			final Set<VCFHeaderLine> metaData = new HashSet<>();
			metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_KEY,true));
			metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_QUALITY_KEY,true));
			metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.DEPTH_KEY,true));
			metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_ALLELE_DEPTHS,true));
			metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_FILTER_KEY,true));
			metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY,true));
			metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_NUMBER_KEY,true));
			metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY,true));
			metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY,true));
			metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_FREQUENCY_KEY,true));
			metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY,true));
			
			final VCFFormatHeaderLine formatHeteroPlasmy = new VCFFormatHeaderLine(
					"HP",
					1,
					VCFHeaderLineType.Float,
					"Heteroplasmy (read alt/depth)"
					);
			metaData.add(formatHeteroPlasmy);
			final VCFFormatHeaderLine formatFisherStrand = new VCFFormatHeaderLine(
					"FS",
					1,
					VCFHeaderLineType.Float,
					"Fisher Strand"
					);
			metaData.add(formatFisherStrand);
			
			final VCFFormatHeaderLine formatDP8 = new VCFFormatHeaderLine(
					"DP8",
					8,
					VCFHeaderLineType.Integer,
					"number of bases ACGT(formard) acgt(reverse) "
					);
			metaData.add(formatDP8);
			
			final VCFInfoHeaderLine infoSamples = new VCFInfoHeaderLine(
					"SAMPLES",
					VCFHeaderLineCount.UNBOUNDED,
					VCFHeaderLineType.String,
					"Samples with genotype"
					);
			metaData.add(infoSamples);
			
			final VCFInfoHeaderLine infoNoRefSamples = new VCFInfoHeaderLine(
					"NOREFSAMPLES",
					VCFHeaderLineCount.UNBOUNDED,
					VCFHeaderLineType.String,
					"Samples without REF base"
					);
			metaData.add(infoNoRefSamples);
			
			final VCFInfoHeaderLine infoNSamples = new VCFInfoHeaderLine(
					"NSAMPLES",
					1,
					VCFHeaderLineType.Integer,
					"Number of Samples with genotype"
					);
			metaData.add(infoNSamples);
			
			metaData.add(new VCFFilterHeaderLine(formatFisherStrand.getID(),"Fails fisher test for Strands : "+this.fisher_strand_bias));
			
			final VCFHeader header=new VCFHeader(metaData,sample2heteroplasmy.keySet()); 
			header.setSequenceDictionary(dict1);
			JVarkitVersion.getInstance().addMetaData(this, header);
			
			vcw = VCFUtils.createVariantContextWriterToPath(this.outputFile);
			vcw.writeHeader(header);
			
			
			for(int pos1=1; pos1 <= this.chrMSequence.length();pos1++) {
				final byte ref_base = this.chrMSequence.getBases()[pos1-1];
				final Allele ref = Allele.create(ref_base, true);
				for(int alt_base: ACGT) {
					if(alt_base==(int)ref_base) continue;
					final Set<String> filters = new HashSet<>();
					final Set<String> samplesWithGt = new TreeSet<>();
					final Set<String> noRefSamples = new TreeSet<>();
					final Allele alt= Allele.create((byte)alt_base, false);
					int max_gq = 0;
					final VariantContextBuilder vcb = new VariantContextBuilder(null, this.chrMSequence.getName(), pos1, pos1,
							Arrays.asList(ref,alt));
					
					final List<Genotype> genotypes = new ArrayList<>(sample2heteroplasmy.size());
					for(SampleHeteroplasmy compound: sample2heteroplasmy.values()) {
						final Pileup pileup = compound.pos1ToPileups.get(pos1);
						
						if( pileup==null ) {
							continue;
							}
						
						if(pileup.count(ref)==0) {
							noRefSamples.add(compound.sn);
							}
						
						if( pileup.count(alt) < min_allele_dp  ) {
							genotypes.add(GenotypeBuilder.createMissing(compound.sn, 1));
							continue;
							}
						
						
						final int dp ;
						if(use_all_alleles) {
							dp= pileup.counter.values().
									stream().
									mapToInt(I->I).
									sum();
							}
						else
							{
							dp = pileup.countIgnoreCase(ref) +  pileup.countIgnoreCase(alt);
							}
						
						if(dp==0) {
							genotypes.add(GenotypeBuilder.createMissing(compound.sn, 1));
							continue;
							}
						
						final GenotypeBuilder gb = new GenotypeBuilder(compound.sn,Collections.singletonList(alt));

						
						gb.DP(dp);
						gb.AD(new int[] {
								pileup.countIgnoreCase(ref),
								pileup.countIgnoreCase(alt)
							});
						
						final double p = pileup.countIgnoreCase(alt)/dp;
						gb.attribute(formatHeteroPlasmy.getID(), p);
						
					
						double tmp = 1.96 * Math.sqrt( ( p * ( 1.0 - p ) ) / dp );
						double min_confidence = Math.max(p-tmp, 0.0);
						double max_confidence = Math.min(p+tmp, 1.0);
						double distance = max_confidence - min_confidence;
						int gq = (int)(1.0-distance)*99;
						max_gq = Math.max(max_gq, gq);
						gb.GQ(gq);
						
						final int dp8[]=new int[ACGT.length*2];
						for(int y=0;y< ACGT.length;++y) {
							char b = (char)ACGT[y];
							dp8[y] = pileup.countPlus(b);
							dp8[y + ACGT.length ] = pileup.countMinus(b);
							}
						
						gb.attribute(formatDP8.getID(),dp8);
						
							
						final FisherExactTest fisher = FisherExactTest.compute(
								pileup.countPlus(ref), 
								pileup.countMinus(ref), 
								pileup.countPlus(alt), 
								pileup.countMinus(alt)
								);
						gb.attribute(formatFisherStrand.getID(), fisher.getAsDouble());
						
						if(fisher.getAsDouble()< this.fisher_strand_bias) {
							gb.attribute(VCFConstants.GENOTYPE_FILTER_KEY, formatFisherStrand.getID());
							filters.add(formatFisherStrand.getID());
							}
						
						genotypes.add(gb.make());
						samplesWithGt.add(compound.sn);
						} // end of loop compound
					
					//no genotype found here, don't create a variant
					if(samplesWithGt.isEmpty()) continue;
					if(filters.isEmpty()) {
						vcb.passFilters();
					} else
						{
						vcb.filters(filters);
						}
					if(max_gq>0) vcb.log10PError(max_gq/-10.0);
					vcb.attribute(infoNSamples.getID(), samplesWithGt.size());
					vcb.attribute(infoSamples.getID(),new ArrayList<>(samplesWithGt));
					vcb.genotypes(genotypes);
					vcw.add(vcb.make());
					
					
					} //end loop over alt base
				}
			vcw.close();
			vcw=null;
			referenceSequenceFile.close();
			referenceSequenceFile=null;
			return 0;
		} catch (final Throwable e) {
			LOG.error(e);
			return -1;
		} finally {
			CloserUtil.close(vcw);
			CloserUtil.close(referenceSequenceFile);
			}
		}
	public static void main(final String[] args) {
		new BamHeteroplasmy().instanceMainWithExit(args);
	}
}
