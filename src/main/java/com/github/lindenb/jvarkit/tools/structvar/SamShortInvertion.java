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
package com.github.lindenb.jvarkit.tools.structvar;

import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Consumer;
import java.util.function.ToIntBiFunction;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.tools.misc.ConcatSam;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.StringUtil;
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

finds the regions having some short inversions.

input is a set of BAM files. One file ending with '.list' is interpreted as a file containing some path to the bams.

output is a VCF file

## Example:

```
```

END_DOC

**/

@Program(name="samshortinvert",
	description="Scan short inversions in SAM",
	keywords={"sam","bam","sv","inversion"}
	)
public class SamShortInvertion extends Launcher
	{
	private static final Logger LOG = Logger.build(SamShortInvertion.class).make();
	private static final byte SUPPORTING_LEFT=(byte)1;
	private static final byte SUPPORTING_RIGHT=(byte)2;
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-m","--maxsize"},description="max size of inversion")
	private int max_size_inversion = 10_000 ;
	@Parameter(names={"-r","--rgn"},description="limit to that region CHROM:START-END")
	private String intervalStr = null;
	@Parameter(names={"-B","--bed"},description="limit to that bed file")
	private File intervalBed = null;
	@Parameter(names={"-partition","--partition"},description=SAMRecordPartition.OPT_DESC)
	private SAMRecordPartition partition = SAMRecordPartition.sample;
	@Parameter(names={"-x","--extend"},description="extends interval by 'x' pb before merging.")
	private int extend=50;
	@Parameter(names={"-s","-supporting"},description="Don't print the variant if INFO/DP <= 's'")
	private int min_supporting_reads = 1;
	@Parameter(names={"--debug"},description="Debug",hidden=true)
	private boolean debug = false;

	private static class Arc
		{
		String sample;
		int tid;
		int chromStart;
		int chromEnd;
		boolean consummed=false;
		byte type;
		
		int length() {
			return chromEnd-chromStart+1;
		}
		
		@Override
		public String toString() {
			return "("+tid+"):"+chromStart+"-"+chromEnd +" type="+(int)type;
			}
		}
	
	

	private void dump(
			final SAMSequenceDictionary dict,
			final IntervalTreeMap<List<Arc>> database,
			final VariantContextWriter vcw,
			final Set<String> samples,
			final Integer before
			) {
		if(this.debug) LOG.debug("dump");
		final Allele REF = Allele.create("N", true);
		final Allele SPLIT = Allele.create("<INV>", false);
		final ContigDictComparator ctgCmp  = new ContigDictComparator(dict);
		final List<Interval> intervals  = database.keySet().stream().
				map(R-> new Interval(
					R.getContig(),
					Math.max(1,R.getStart() - this.extend),
					R.getEnd() + this.extend
					)).
				filter(R->(before==null?true:R.getEnd() < before.intValue())).
				sorted((A,B)->{
					int i = ctgCmp.compare(A.getContig(), B.getContig());
					if(i!=0) return i;
					i = A.getStart() - B.getStart();
					if(i!=0) return i;
					return A.getEnd() - B.getEnd();
					}).
				collect(Collectors.toList());
		
		for(final Interval interval0:intervals) {
			
			final List<Arc> arcs = database.getOverlapping(interval0).
					stream().
					flatMap(L->L.stream()).
					filter(A->!A.consummed).
					filter(A->
						Math.abs(interval0.getStart()-A.chromStart) <= this.extend &&
						Math.abs(interval0.getEnd()-A.chromEnd) <= this.extend).
					collect(Collectors.toList());
			
			if(arcs.isEmpty()) continue;
			arcs.forEach(A->A.consummed=true);
			
			int maxdp = 0;
			final VariantContextBuilder vcb = new VariantContextBuilder();
			final Set<Allele> alleles = new HashSet<>();
			alleles.add(REF);
			final List<Genotype> genotypes = new ArrayList<>(samples.size());
			
			
			vcb.chr(dict.getSequence(arcs.get(0).tid).getSequenceName());
			final int chromStart = arcs.stream().mapToInt(A->A.chromStart).min().getAsInt();
			vcb.start(chromStart);
			final int chromEnd = arcs.stream().mapToInt(A->A.chromEnd).max().getAsInt();
			vcb.stop(chromEnd);
			
			vcb.attribute(VCFConstants.END_KEY, chromEnd);
			vcb.attribute("SVLEN", 1+chromEnd-chromStart);
			
			int depth = 0;
			int nsamples = 0;
			for(final String sample : samples) {
			
				final List<Arc> sampleArcs = arcs.stream().filter(A->A.sample.equals(sample)).collect(Collectors.toList());
				if(sampleArcs.isEmpty())
					{
					genotypes.add(GenotypeBuilder.createMissing(sample, 2));
					}
				else
					{
					final GenotypeBuilder gb = new GenotypeBuilder(sample);
					alleles.add(SPLIT);
					gb.alleles(Arrays.asList(REF,SPLIT));
					final int countCat1= (int)sampleArcs.stream().filter(A->A.type==SUPPORTING_LEFT).count();
					final int countCat2= (int)sampleArcs.stream().filter(A->A.type==SUPPORTING_RIGHT).count();
					gb.DP(countCat1+countCat2);
					gb.attribute("N5", countCat1);
					gb.attribute("N3", countCat2);
					
					maxdp = Math.max(maxdp, countCat1+countCat2);
					
					depth+=countCat1+countCat2;
					genotypes.add(gb.make());
					++nsamples;
					}
				}
			if(depth<=this.min_supporting_reads) continue;
			
			vcb.genotypes(genotypes);
			vcb.alleles(alleles);
			vcb.attribute(VCFConstants.DEPTH_KEY, depth);
			vcb.attribute("NSAMPLES", nsamples);
			vcb.attribute(VCFConstants.SVTYPE, "INV");
			vcb.attribute("DPMAX", maxdp);
			vcw.add(vcb.make());
			}
		
		}

	
	@Override
	public int doWork(final List<String> args) {
		if(this.max_size_inversion<100) {
			LOG.error("max size insersion must be >=100");
			return -1;
		}
	 	ConcatSam.ConcatSamIterator iter = null;
	 	VariantContextWriter vcw= null;
	 	final  IntervalTreeMap<List<Arc>> database = new IntervalTreeMap<>(); 
		try {
			final ConcatSam.Factory concatFactory = new ConcatSam.Factory();
			concatFactory.setEnableUnrollList(true);
			concatFactory.addInterval(this.intervalStr);
			if(this.intervalBed!=null)
				{
				final BedLineCodec codec = new BedLineCodec();
				final BufferedReader br = IOUtils.openFileForBufferedReading(this.intervalBed);
				br.lines().
					filter(L->!StringUtil.isBlank(L)).
					map(L->codec.decode(L)).
					filter(L->L!=null).
					forEach(B->concatFactory.addInterval(B.getContig()+":"+B.getStart()+"-"+B.getEnd()));
				br.close();
				}
			
			
			iter = concatFactory.open(args);
			
			final SAMSequenceDictionary dict = iter.getFileHeader().getSequenceDictionary();
			
			final Set<String> samples = iter.getFileHeader().getReadGroups().stream().map(RG->RG.getSample()).filter(S->!StringUtil.isBlank(S)).collect(Collectors.toSet());
			if(samples.isEmpty())
				{
				iter.close();
				LOG.error("No samples defined");
				return -1;
				}
			
			final ToIntBiFunction<Locatable, Locatable> distance = (A,B) ->{
				if(CoordMath.overlaps(A.getStart(), A.getEnd(), B.getStart(), B.getEnd())) return 0;
				if(A.getEnd()<B.getStart()) {
					return B.getStart() - A.getEnd();
					}
				else
					{
					return A.getStart() - B.getEnd();
					}
				};
		
			
			final Set<VCFHeaderLine> meta=new HashSet<>();
			VCFStandardHeaderLines.addStandardFormatLines(meta,true,
					VCFConstants.GENOTYPE_KEY,
					VCFConstants.DEPTH_KEY
					);
			VCFStandardHeaderLines.addStandardInfoLines(meta,true,
					VCFConstants.DEPTH_KEY,
					VCFConstants.END_KEY
					
					);
			meta.add(new VCFFormatHeaderLine("N5", 1, VCFHeaderLineType.Integer,"Number of validating clipped reads in 5'"));
			meta.add(new VCFFormatHeaderLine("N3", 1, VCFHeaderLineType.Integer,"Number of validating clipped reads in 3'"));
			meta.add(new VCFInfoHeaderLine("SVLEN", 1, VCFHeaderLineType.Integer,"SV length"));
			meta.add(new VCFInfoHeaderLine("NSAMPLES", 1, VCFHeaderLineType.Integer,"Number of sample having some split reads"));
			meta.add(new VCFInfoHeaderLine("DPMAX", 1, VCFHeaderLineType.Integer,"MAX DP among samples"));
			meta.add(new VCFInfoHeaderLine("SVTYPE", 1, VCFHeaderLineType.String,"Structural variant type"));
			
			
			
			final VCFHeader header=new VCFHeader(meta,samples);
			JVarkitVersion.getInstance().addMetaData(this, header);
			header.setSequenceDictionary(dict);
			vcw = super.openVariantContextWriter(outputFile);
			vcw.writeHeader(header);
			
			
			final ProgressFactory.Watcher<SAMRecord> progress= ProgressFactory.
					newInstance().
					dictionary(iter.getFileHeader().getSequenceDictionary()).
					logger(LOG).
					build();
			final short SA_TAG = SAMTag.SA.getBinaryTag();
			
			String prevContig=null;
			while(iter.hasNext())
				{
				final SAMRecord rec = progress.apply(iter.next());
				
				if(rec.getReadUnmappedFlag()) continue;
				if(rec.getReadFailsVendorQualityCheckFlag()) continue;
				if(rec.isSecondaryOrSupplementary()) continue;
				if(rec.getDuplicateReadFlag()) continue;
				if(rec.getAttribute(SA_TAG) == null) continue;
				
				final Cigar cigar = rec.getCigar();
				if(cigar==null || cigar.isEmpty() || !cigar.isClipped()) continue;
				
				final String sample= this.partition.getPartion(rec, null);
				if(StringUtil.isBlank(sample))continue;
				
				final List<SAMRecord> others = SAMUtils.getOtherCanonicalAlignments(rec).
						stream().
						filter(R->rec.getContig().equals(R.getContig())).
						filter(R->rec.getReadNegativeStrandFlag()!=R.getReadNegativeStrandFlag()).
						filter(R->distance.applyAsInt(rec,R)< this.max_size_inversion).
						collect(Collectors.toList());
				
				if(others.isEmpty()) continue;
				
				if(!rec.getContig().equals(prevContig)) {
					dump(dict,database,vcw,samples,null);
					database.clear();
					prevContig = rec.getContig();
					}
				else
					{
					final int before = rec.getUnclippedStart() - this.max_size_inversion*2;
					dump(dict,database,vcw,samples,before);
					database.entrySet().removeIf(entries->entries.getKey().getEnd()< before);
					}

				final Consumer<Arc> registerArc = (A)->{
					if(A.chromEnd<=A.chromStart) throw new IllegalArgumentException(A.toString());
					final Interval rgn = new Interval(rec.getContig(), A.chromStart,A.chromEnd);
					List<Arc> list = database.get(rgn);
					if(list==null) {
						list = new ArrayList<>();
						database.put(rgn,list);
						}
					list.add(A);					
					};
				
				if(cigar.isLeftClipped())
					{
					for(final SAMRecord rec2:others) {
						// NON if(rec.getEnd()>= rec2.getStart()) continue;
						final Arc arc = new Arc();
						arc.sample  = sample;
						arc.tid = rec.getReferenceIndex();
						arc.chromStart = Math.min(rec.getStart(),rec2.getStart());
						arc.chromEnd = Math.max(rec.getEnd(),rec2.getEnd());
						if(arc.length()> this.max_size_inversion) continue;

						
						arc.type = SUPPORTING_LEFT;
						registerArc.accept(arc);
						}
					}
			
				
				if(cigar.isRightClipped())
					{
					for(final SAMRecord rec2:others) {
						final Arc arc = new Arc();
						arc.sample  = sample;
						arc.tid = rec.getReferenceIndex();
						arc.chromStart = Math.min(rec.getStart(),rec2.getStart());
						arc.chromEnd = Math.max(rec.getEnd(),rec2.getEnd());
						if(arc.length()> this.max_size_inversion) continue;
						
						arc.type = SUPPORTING_RIGHT;
						registerArc.accept(arc);
						}
					}				
				}
			dump(dict,database,vcw,samples,null);
			iter.close();iter=null;
			progress.close();
			vcw.close();vcw=null;
			return 0;
			} 
		catch (final Exception e) {
			LOG.error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);	
			CloserUtil.close(vcw);	
			}
		}

	
	public static void main(final String[] args)
		{
		new SamShortInvertion().instanceMainWithExit(args);
		}
	}
