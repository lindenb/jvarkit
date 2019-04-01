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

import java.io.Closeable;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.OptionalDouble;
import java.util.Set;

import org.apache.commons.math3.stat.descriptive.moment.Mean;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.StructuralVariantType;
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
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFStandardHeaderLines;


/**
BEGIN_DOC

## Motivation

only SVTYPE=DEL or SVTYPE=INS are considered

## Example

```
find DIR -type f -name "*.bam" > bam.list
$ java -jar ${JVARKIT_HOME}/dist/validatecnv.jar \
	-B bam.list --min-read-support-del 1 --median-adjust 0.5 --max-variance 15 --extend 0.33 20190320.MANTA.vcf 
```

END_DOC

 */
@Program(name="validatecnv",
	description="experimental CNV detection. Look variange of depths before/after putative known CNV.",
	keywords= {"cnv","bam","sam","vcf"},
	modificationDate="20190328",
	generate_doc=false
	)
public class ValidateCnv extends Launcher
	{
	private static final Logger LOG = Logger.build(ValidateCnv.class).make();
	
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile=null;
	@Parameter(names={"-B","--bams","--bam"},description="Path to bam. File with suffix .list is interpretted as a file containing a list of paths to bams.")
	private List<Path> bamFiles = new ArrayList<>();
	@Parameter(names={"-x","--extend"},description="Search the boundaries in a region that is 'x'*(CNV-length)")
	private double extendFactor=0.1;	
	@Parameter(names={"-md","--min-dp"},description="At least one of the bounds must have a median-depth greater than this value.")
	private int min_depth = 20;
	@Parameter(names={"-cd","--critical-dp"},description="The middle section shouldn't have any point with a DP lower than this value (prevent HOM_VAR deletions)")
	private int critical_middle_depth = 1;
	@Parameter(names={"--min"},description="Min abs(SV) size to consider." +DistanceParser.OPT_DESCRIPTION ,converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int min_abs_sv_size = 50;
	@Parameter(names={"--max"},description="Max abs(SV) size to consider." +DistanceParser.OPT_DESCRIPTION ,converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int max_abs_sv_size = 1_000_000;
	@Parameter(names={"--dicard"},description="Dicard variant where all Called Genotypes have been filtered.")
	private boolean discard_all_filtered_variant = false;
	@Parameter(names={"--max-variance"},description="Maximum tolerated variance in one genomic segment.")
	private double max_variance = 10.0;
	@Parameter(names={"--median-adjust"},description="TODO.")
	private double median_factor = 0.33;
	@Parameter(names={"--min-read-support-del"},description="min number of read supporting deletion.")
	private int min_read_supporting_del=3;
	
	private final Mean medianCalc =new Mean();

	
	private class Coverage
		{
		final OptionalDouble median;
		final OptionalDouble variance;
		
		Coverage(final double array[],int array_start,int array_end) {
			final int len = array_end-array_start;
			if(len>0) {
				double m = medianCalc.evaluate(array, array_start, len);
				median = OptionalDouble.of(m);
				double v = 0;
				int n=0;
				for(int i=array_start;i< array_end;i++) {
					v+=Math.abs(array[i]-m);
					n++;
					}
				variance = OptionalDouble.of(v/n);
				}
			else
				{
				median = OptionalDouble.empty();
				variance = OptionalDouble.empty();
				}
			}
		@Override
		public String toString() {
			return "median:"+ (median.isPresent()?".":String.valueOf(median.getAsDouble()))+" "+
					"variance:"+(variance.isPresent()?".":String.valueOf(variance.getAsDouble()))
					;
			}
		}
	
	
	/** BAM input */
	private class Input implements Closeable
		{
		final SamReader samReader;
		final SAMFileHeader header;
		final SAMSequenceDictionary dict;
		final ContigNameConverter ctgNameConverter;
		final String sampleName;
		Input(final Path uri) throws IOException {
			final SamInputResource sri = SamInputResource.of(uri);
			this.samReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT).open(sri);
			this.header  = samReader.getFileHeader();
			this.dict = this.header.getSequenceDictionary();
			this.ctgNameConverter = ContigNameConverter.fromOneDictionary(this.dict);
			this.sampleName = this.header.getReadGroups().stream().
						map(R->R.getSample()).
						filter(S->!StringUtil.isBlank(S)).
						findFirst().
						orElseThrow(()->new IllegalArgumentException("no ReadGroup/SM defined in "+uri));
			}
		@Override
		public void close() {
			CloserUtil.close(this.samReader);
			}
		}
	

	private boolean isSameMedianDepth(final OptionalDouble dp1,final OptionalDouble dp2) {
		if(!dp1.isPresent()) return true;
		if(!dp2.isPresent()) return true;
		double v1=dp1.getAsDouble();
		double v1a  = v1*this.median_factor;
		double v2=dp2.getAsDouble();
		double v2a  = v2*this.median_factor;
		if(v2 < v1 - v1a) return false;
		if(v2 > v1 + v1a) return false;
		if(v1 < v2 - v2a) return false;
		if(v1 > v2 + v2a) return false;
		return true;
	}
	private boolean isHightVariance(OptionalDouble va) {
		if(!va.isPresent()) return false;
		double v1=va.getAsDouble();
		return v1> this.max_variance;
	}
	
	@Override
	public int doWork(final List<String> args) {		
		if(this.extendFactor<=0)
			{
			LOG.error("bad extend factor "+this.extendFactor);
			return -1;
			}
		
		
		final Map<String,Input> sample2bam = new HashMap<>();
		final List<Input> inputs = new ArrayList<>();
		VariantContextWriter out = null;
		VCFIterator iterIn = null;
		try
			{
			
			
			iterIn = super.openVCFIterator(oneFileOrNull(args));
			final VCFHeader header0 = iterIn.getHeader();
			if(!header0.hasGenotypingData()) {
				LOG.error("No genotype in input");
				return -1;
				}
			
			for(final Path p: this.bamFiles) {
				for(final Path p2:IOUtils.unrollPath(p)) {
					final Input input = new Input(p2);
					if(sample2bam.containsKey(input.sampleName)) {
						LOG.error("sample "+input.sampleName+" specified twice.");
						input.close();
						return -1;
						}
					if(!header0.getSampleNameToOffset().containsKey(input.sampleName)) {
						LOG.error("sample "+input.sampleName+" is not present in vcf.");
						input.close();
						continue;
					}
					sample2bam.put(input.sampleName, input);
				}
			}
			if(sample2bam.isEmpty()) {
				LOG.error("no bam was defined");
				return -1;
			}
			
			final SAMSequenceDictionary dict/* may be null*/ = header0.getSequenceDictionary();
			
			final Set<VCFHeaderLine> metadata = new HashSet<>();
			
			final VCFInfoHeaderLine infoSVSamples = 
					new VCFInfoHeaderLine("PASSING_SAMPLES",VCFHeaderLineCount.UNBOUNDED,
							VCFHeaderLineType.String,
							"Samples that could carry a SV");
			metadata.add(infoSVSamples);
			
			final VCFFormatHeaderLine leftMedianDepth = 
					new VCFFormatHeaderLine("LDP",1,VCFHeaderLineType.Integer,"Left median depth or -1");
			metadata.add(leftMedianDepth);
			final VCFFormatHeaderLine rightMedianDepth = 
					new VCFFormatHeaderLine("RDP",1,VCFHeaderLineType.Integer,"Right median depth or -1");
			metadata.add(rightMedianDepth);
			final VCFFormatHeaderLine midMedianDepth = 
					new VCFFormatHeaderLine("MDP",1,VCFHeaderLineType.Integer,"Middle median depth or -1");
			metadata.add(midMedianDepth);
			final VCFFormatHeaderLine leftVarianceDepth = 
					new VCFFormatHeaderLine("LVA",1,VCFHeaderLineType.Integer,"Left variance depth or -1");
			metadata.add(leftVarianceDepth);
			final VCFFormatHeaderLine rightVarianceDepth = 
					new VCFFormatHeaderLine("RVA",1,VCFHeaderLineType.Integer,"Right variance depth or -1");
			metadata.add(rightVarianceDepth);
			final VCFFormatHeaderLine midVarianceDepth = 
					new VCFFormatHeaderLine("MVA",1,VCFHeaderLineType.Integer,"Middle variance depth or -1");
			metadata.add(midVarianceDepth);
			final VCFFormatHeaderLine nReadsSupportingDel = 
					new VCFFormatHeaderLine("RSD",1,VCFHeaderLineType.Integer,"Read supporting deletion");
			metadata.add(nReadsSupportingDel);
			
			VCFStandardHeaderLines.addStandardFormatLines(metadata, true,
					VCFConstants.DEPTH_KEY,
					VCFConstants.GENOTYPE_KEY,
					VCFConstants.GENOTYPE_FILTER_KEY
					);
			VCFStandardHeaderLines.addStandardInfoLines(metadata, true,
					VCFConstants.DEPTH_KEY,
					VCFConstants.END_KEY
					);

			
			final VCFHeader header = new VCFHeader(header0);
			JVarkitVersion.getInstance().addMetaData(this, header);
			metadata.stream().forEach(M->header.addMetaDataLine(M));
			
			final ProgressFactory.Watcher<VariantContext> progress = ProgressFactory.newInstance().dictionary(header0).logger(LOG).build();
			out =  VCFUtils.createVariantContextWriterToPath(this.outputFile);
			out.writeHeader(header);
			
			while(iterIn.hasNext())
				{
				final VariantContext ctx = iterIn.next();
				final StructuralVariantType svType = ctx.getStructuralVariantType();
				if(!(svType==StructuralVariantType.DEL || svType==StructuralVariantType.INS) ) {
					continue;
				}
				final int svLen;
				if(ctx.hasAttribute("SVLEN")) {
					svLen = Math.abs(ctx.getAttributeAsInt("SVLEN", -1));
				} else {
					svLen = ctx.getLengthOnReference();
				}
				if(svLen< this.min_abs_sv_size) continue;
				if(svLen> this.max_abs_sv_size) continue;
				
				final int extend = 1+(int)(svLen*this.extendFactor);
				
				
				
				final int leftPos =  Math.max(1, ctx.getStart()-extend);
				final int array_mid_start = ctx.getStart()-leftPos;
				final int array_mid_end = ctx.getEnd()-leftPos;
				int rightPos =  ctx.getEnd()+extend;
				
				if(dict!=null) {
					final SAMSequenceRecord ssr = dict.getSequence(ctx.getContig());
					if(ssr!=null) {
						rightPos = Math.min(rightPos, ssr.getSequenceLength());
						}
					}
				
				//System.err.println(""+leftPos+" "+ctx.getStart()+" "+ctx.getEnd()+" "+rightPos);
				
				final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
				
			    final Set<String> sv_samples = new HashSet<>(header.getSampleNamesInOrder());
				
				final List<Genotype> genotypes = new ArrayList<>(ctx.getNSamples());
				for(final Genotype gt : ctx.getGenotypes())
					{					
					final Input input = sample2bam.get(gt.getSampleName());
					if(input==null) {
						sv_samples.remove(gt.getSampleName());
						genotypes.add(gt);
						continue;
						}
					final String ctg2 =  input.ctgNameConverter.apply(ctx.getContig());
					if(StringUtils.isBlank(ctg2)) {
						sv_samples.remove(gt.getSampleName());
						genotypes.add(gt);
						continue;
						}
					
					final double coverage[] = new double[rightPos-leftPos+1];
					Arrays.fill(coverage, 0.0);
					int n_reads_supporting_deletions = 0;
					try(CloseableIterator<SAMRecord> iter2 = input.samReader.queryOverlapping(
							ctg2,
							leftPos,
							rightPos
							)) {
						while(iter2.hasNext()) {
							final SAMRecord rec = iter2.next();
							if(rec.getReadUnmappedFlag()) continue;
							if(rec.getReadUnmappedFlag()) continue;
							if(rec.getReadFailsVendorQualityCheckFlag()) continue;
							if(rec.getDuplicateReadFlag()) continue;
							if(rec.isSecondaryOrSupplementary()) continue;
							final Cigar cigar = rec.getCigar();
							if(cigar==null || cigar.isEmpty()) continue;
							
							// any read supporting deletion ?
							if(svType==StructuralVariantType.DEL &&
								rec.getReadPairedFlag() && 
								!rec.getProperPairFlag() &&
								rec.getReadLength() < svLen
								)
								{
								final List<SAMRecord> suppl = SAMUtils.getOtherCanonicalAlignments(rec);
								if(rec.getStart() <= ctx.getStart() && 
										suppl.stream().
											filter(R->R.getContig().equals(ctg2)).
											anyMatch(R->R.getEnd()>=ctx.getEnd())
										) {
										n_reads_supporting_deletions++; 
										}
								else if(rec.getEnd() >= ctx.getEnd() &&
										suppl.stream().
										filter(R->R.getContig().equals(ctg2)).
										anyMatch(R->R.getStart()<=ctx.getStart())
										) {
										n_reads_supporting_deletions++; 
										}
								}
							
							int ref=rec.getStart();
							for(final CigarElement ce:cigar) {
								final CigarOperator op = ce.getOperator();
								if(op.consumesReferenceBases())
									{
									if(op.consumesReadBases()) {
										for(int x=0;x < ce.getLength() && ref + x - leftPos < coverage.length ;++x)
											{
											final int p = ref + x - leftPos;
											if(p<0 || p>=coverage.length) continue;
											coverage[p]++;
											}
										}
									ref+=ce.getLength();
									}
								}
							}// end while iter record
						}//end try
					
					final Coverage covL = new Coverage(coverage,0, array_mid_start);
					final Coverage covM = new Coverage(coverage,array_mid_start, array_mid_end);
					final Coverage covR = new Coverage(coverage,array_mid_end, coverage.length);
										
					final GenotypeBuilder gb = new GenotypeBuilder(gt);
					gb.attribute(nReadsSupportingDel.getID(), n_reads_supporting_deletions);
					
					
					
					gb.attribute(leftMedianDepth.getID(),covL.median.isPresent()? (int)covL.median.getAsDouble():-1);
					gb.attribute(midMedianDepth.getID(),covM.median.isPresent()? (int)covM.median.getAsDouble():-1);
					gb.attribute(rightMedianDepth.getID(),covR.median.isPresent()? (int)covR.median.getAsDouble():-1);
					gb.attribute(leftVarianceDepth.getID(),covL.variance.isPresent()? (int)covL.variance.getAsDouble():-1);
					gb.attribute(midVarianceDepth.getID(),covM.variance.isPresent()? (int)covM.variance.getAsDouble():-1);
					gb.attribute(rightVarianceDepth.getID(),covR.variance.isPresent()? (int)covR.variance.getAsDouble():-1);
					
					final Set<String> gtFilters = new HashSet<>();
					
					//prevent HOM_VAR
					for(int j=array_mid_start;j< array_mid_end;j++)
						{
						if(coverage[j]<this.critical_middle_depth)
							{
							gtFilters.add("HOMVAR");
							break;
							}
						}
					// there are some split read that could support the SV in non-called
					if(svType==StructuralVariantType.DEL &&
						n_reads_supporting_deletions < this.min_read_supporting_del ) {
						gtFilters.add("SPLITREAD");
						}
					
					// high variance left
					if(isHightVariance(covL.variance)) {
						gtFilters.add("HIGHVARL");
						}
					
					// high variance left
					if(isHightVariance(covR.variance)) {
						gtFilters.add("HIGHVARR");
						}
					//low dp left
					if(covL.median.isPresent() && covL.median.getAsDouble() < this.min_depth) {
						gtFilters.add("LOWDPL");
						}
					//low dp right
					if(covR.median.isPresent() && covR.median.getAsDouble() < this.min_depth) {
						gtFilters.add("LOWDPR");
						}
					
					// big difference of depth left/right
					if(!isSameMedianDepth(covL.median, covR.median)) {
						gtFilters.add("DIFFLR");
						}
					// no difference between mid and left
					if(isSameMedianDepth(covL.median, covM.median))
						{
						gtFilters.add("DIFFLM");
						}
					// no difference between mid and right
					if(isSameMedianDepth(covM.median, covR.median))
						{
						gtFilters.add("DIFFMR");
						}
					// no diff beween L and M
					if(isSameMedianDepth(covL.median, covM.median)) {
						gtFilters.add("NELM");
						}
					
					// no diff beween M and R
					if(isSameMedianDepth(covM.median, covR.median)) {
						gtFilters.add("NEMR");
						}
					
					if(!gtFilters.isEmpty()) {
						sv_samples.remove(gt.getSampleName());
						gb.filter(String.join("~",gtFilters));
						}
					genotypes.add(gb.make());
					}
				
				if(sv_samples.isEmpty() && this.discard_all_filtered_variant) {
					continue;
					}
				
				if(!sv_samples.isEmpty() ) {
					vcb.attribute(infoSVSamples.getID(),new ArrayList<>(sv_samples));
					if(!ctx.isFiltered()) vcb.passFilters();
				}
				
				vcb.genotypes(genotypes);
				
				out.add(vcb.make());
				}
			progress.close();
			out.close();
			iterIn.close();
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(iterIn);
			CloserUtil.close(inputs);
			CloserUtil.close(out);
			sample2bam.values().forEach(F->F.close());
			}
		}
	
	
	public static void main(final String[] args) {
		new ValidateCnv().instanceMainWithExit(args);
		}
	}
