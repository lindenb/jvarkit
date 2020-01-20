/*
The MIT License (MIT)

Copyright (c) 2020 Pierre Lindenbaum

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
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.BiFunction;
import java.util.function.IntPredicate;
import java.util.stream.IntStream;


import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.math.stats.Percentile;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.samtools.util.SimplePosition;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
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
	description="Experimental CNV Genotyping. Look variance of depths before/after putative known CNV.",
	keywords= {"cnv","bam","sam","vcf","depth"},
	modificationDate="20200120",
	generate_doc=false
	)
public class ValidateCnv extends Launcher
	{
	private static final Logger LOG = Logger.build(ValidateCnv.class).make();
	
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile=null;
	@Parameter(names={"-B","--bams","--bam"},description="Path to bam. File with suffix .list is interpretted as a file containing a list of paths to bams.")
	private List<String> bamFiles = new ArrayList<>();
	@Parameter(names={"-R","--reference"},description="For reading CRAM. "+INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private Path rererencePath = null;
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
	@ParametersDelegate
	private WritingVariantsDelegate writingVariantsDelegate = new WritingVariantsDelegate();

	

	
	
	
	
	private static class Coverage
		{
		final Double median;
		final Double variance;
		Coverage(final Double median,final Double variance) {
			this.median= median;
			this.variance = variance;
		}
		
		@Override
		public String toString() {
			return "median:"+ (median==null?".":String.valueOf(median))+" "+
					"variance:"+(variance==null?".":String.valueOf(variance))
					;
			}
		}
	
	private Coverage calculateCoverage(final double array[],IntPredicate acceptPos) {
		final double cov[] = IntStream.range(0, array.length).
				filter(acceptPos).
				mapToDouble(POS->array[POS]).
				toArray();
		return calculateCoverage(cov,0,cov.length);
	}
	
	
	private Coverage calculateCoverage(final double array[],int array_start,int array_end) {
		final Double median;
		final Double variance;
		final int len = array_end-array_start;
		if(len>0) {
			final double m = Percentile.median().evaluate(array, array_start, len);
			median = m;
			double v = 0;
			int n=0;
			for(int i=array_start;i< array_end;i++) {
				v+=Math.abs(array[i]-m);
				n++;
				}
			variance = (v/n);
			}
		else
			{
			median = null;
			variance = null;
			}
		return new Coverage(median, variance);
	}
	

	private boolean isSameMedianDepth(final Double dp1,final Double dp2) {
		if(dp1==null) return true;
		if(dp2==null) return true;
		double v1=dp1.doubleValue();
		double v1a  = v1*this.median_factor;
		double v2=dp2.doubleValue();
		double v2a  = v2*this.median_factor;
		if(v2 < v1 - v1a) return false;
		if(v2 > v1 + v1a) return false;
		if(v1 < v2 - v2a) return false;
		if(v1 > v2 + v2a) return false;
		return true;
	}
	private boolean isHightVariance(final Double va) {
		if(va==null) return false;
		double v1=va.doubleValue();
		return v1> this.max_variance;
	}
	
	private enum CnvType {UNDEFINED,WILD,DEL1,DEL2,DUP1,DUP2};

	private CnvType getCNVTypeFromNormDepth(double depth) {
		if(Double.isNaN(depth)) return CnvType.UNDEFINED;
		final double dy1=0.3;
		final double dy2=0.9;
		if(depth<1.0 - dy2 ) return CnvType.DEL2;
		if(depth<1.0 - dy1 ) return CnvType.DEL1;
		if(depth>1.0 + dy2 ) return CnvType.DUP2;
		if(depth>1.0 + dy1 ) return CnvType.DUP1;
		return CnvType.WILD;
		
	}
	
	@Override
	public int doWork(final List<String> args) {		
		if(this.extendFactor<=0)
			{
			LOG.error("bad extend factor "+this.extendFactor);
			return -1;
			}
		
		final Map<String,SamReader> sample2bam = new HashMap<>();
		VariantContextWriter out = null;
		VCFIterator iterIn = null;
		final Allele delAllele  = Allele.create("<DEL>", false);
		final Allele dupAllele  = Allele.create("<DUP>", false);
		final Allele refAllele  = Allele.create("N", true);
		try
			{
			final SamReaderFactory samReaderFactory =  SamReaderFactory.makeDefault().
					referenceSequence(this.rererencePath).
					validationStringency(ValidationStringency.LENIENT)
					;
			final List<Path> bamPaths =IOUtils.unrollPaths(this.bamFiles);

			iterIn = super.openVCFIterator(oneFileOrNull(args));
			final VCFHeader header0 = iterIn.getHeader();
			if(!header0.hasGenotypingData()) {
				LOG.error("No genotype in input");
				return -1;
				}
			/* register each bam */
			for(final Path p2: bamPaths) {				
				final SamReader samReader = samReaderFactory.open(p2);
				if(!samReader.hasIndex()) {
					LOG.error("Bam is not indexed : "+p2);
					return-1;
					}
				
				final String sampleName =  samReader.getFileHeader().
							getReadGroups().
							stream().
							map(R->R.getSample()).
							filter(S->!StringUtil.isBlank(S)).
							findFirst().
							orElseThrow(()->new IllegalArgumentException("no ReadGroup/SM defined in "+p2));
				
				
				if(sample2bam.containsKey(sampleName)) {
					LOG.error("sample "+ sampleName +" specified twice.");
					samReader.close();
					return -1;
					}
				if(!header0.getSampleNameToOffset().containsKey(sampleName)) {
					LOG.error("sample "+ sampleName +" is not present in vcf. Ignoring this bam.");
					samReader.close();
					continue;
					}
				sample2bam.put(sampleName, samReader);
				}
			
			if(sample2bam.isEmpty()) {
				LOG.error("no bam was defined");
				return -1;
			}
			
			final SAMSequenceDictionary dict/* may be null*/ = header0.getSequenceDictionary();
			
			final Set<VCFHeaderLine> metadata = new HashSet<>();
			
			final VCFInfoHeaderLine infoSVSamples = 
					new VCFInfoHeaderLine("N_SAMPLES",
							1,
							VCFHeaderLineType.Integer,
							"Number of Samples that could carry a SV");
			metadata.add(infoSVSamples);
			
			final BiFunction<String, String, VCFFormatHeaderLine> makeFmt = (TAG,DESC)-> new VCFFormatHeaderLine(TAG,1,VCFHeaderLineType.Integer,DESC);
			
			/*
			final VCFFormatHeaderLine leftMedianDepth = makeFmt.apply("LDP","Left median depth or -1");
			metadata.add(leftMedianDepth);
			final VCFFormatHeaderLine rightMedianDepth = makeFmt.apply("RDP","Right median depth or -1");
			metadata.add(rightMedianDepth);
			final VCFFormatHeaderLine midMedianDepth = makeFmt.apply("MDP","Middle median depth or -1");
			metadata.add(midMedianDepth);
			final VCFFormatHeaderLine leftVarianceDepth = makeFmt.apply("LVA","Left variance depth or -1");
			metadata.add(leftVarianceDepth);
			final VCFFormatHeaderLine rightVarianceDepth = makeFmt.apply("RVA","Right variance depth or -1");
			metadata.add(rightVarianceDepth);
			final VCFFormatHeaderLine midVarianceDepth = makeFmt.apply("MVA","Middle variance depth or -1");
			metadata.add(midVarianceDepth);
			*/
			final VCFFormatHeaderLine midVarianceDepth = makeFmt.apply("CBN","Middle variance depth or -1");
			metadata.add(midVarianceDepth);
			final VCFFormatHeaderLine formatCN =  new VCFFormatHeaderLine("CN",1,VCFHeaderLineType.Float,"normalized copy-number");
			metadata.add(formatCN);
			final VCFFormatHeaderLine nReadsSupportingDel = makeFmt.apply("RSD","number of split read supporting SV.");
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
			
			final ProgressFactory.Watcher<VariantContext> progress = ProgressFactory.newInstance().
					dictionary(header0).
					logger(LOG).
					build();
			out =  this.writingVariantsDelegate.
					dictionary(header.getSequenceDictionary()).
					open(this.outputFile);
			out.writeHeader(header);
			
			
			
			while(iterIn.hasNext())
				{
				final VariantContext ctx = iterIn.next();
				final StructuralVariantType svType = ctx.getStructuralVariantType();
				if(!(svType==StructuralVariantType.DEL || svType==StructuralVariantType.DUP || svType==StructuralVariantType.INS) ) {
					continue;
					}
				
				if(ctx.getNAlleles()!=2) continue;
				final int svLen;
				if(ctx.hasAttribute("SVLEN")) {
					svLen = Math.abs(ctx.getAttributeAsInt("SVLEN", -1));
				} else {
					svLen = ctx.getLengthOnReference();
					}
				if(svLen< this.min_abs_sv_size) continue;
				if(svLen> this.max_abs_sv_size) continue;
				
				int n_samples_with_cnv  = 0;
				final SimplePosition breakPointLeft = new SimplePosition(ctx.getContig(), ctx.getStart());
				final SimplePosition breakPointRight = new SimplePosition(ctx.getContig(), ctx.getEnd());
				
				final int extend = 1+(int)(svLen*this.extendFactor);
				
				
				
				final int leftPos =  Math.max(1, breakPointLeft.getPosition()-extend);
				final int array_mid_start = breakPointLeft.getPosition()-leftPos;
				final int array_mid_end = breakPointRight.getPosition()-leftPos;
				int rightPos =  breakPointRight.getPosition()+extend;
				
				if(dict!=null) {
					final SAMSequenceRecord ssr = dict.getSequence(ctx.getContig());
					if(ssr!=null) {
						rightPos = Math.min(rightPos, ssr.getSequenceLength());
						}
					}
				
				//System.err.println(""+leftPos+" "+ctx.getStart()+" "+ctx.getEnd()+" "+rightPos);
				
				final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
				vcb.attributes(new HashMap<>());
				vcb.attribute(VCFConstants.END_KEY, ctx.getEnd());
				
			    final Set<Allele> alleles = new HashSet<>();
			    alleles.add(refAllele);
				
				final List<Genotype> genotypes = new ArrayList<>(ctx.getNSamples());
				
				for(final Genotype gt : ctx.getGenotypes())
					{
					final SamReader samReader =sample2bam.get(gt.getSampleName());
					if(samReader==null) {
						continue;
						}

					final double coverage[] = new double[CoordMath.getLength(leftPos,rightPos)];
					Arrays.fill(coverage, 0.0);
					int n_reads_supporting_deletions = 0;
					
					try(CloseableIterator<SAMRecord> iter2 = samReader.queryOverlapping(
							ctx.getContig(),
							leftPos,
							rightPos
							)) {
						while(iter2.hasNext()) {
							final SAMRecord rec = iter2.next();
							if(rec.getReadUnmappedFlag()) continue;
							if(rec.getReadFailsVendorQualityCheckFlag()) continue;
							if(rec.getDuplicateReadFlag()) continue;
							if(rec.isSecondaryOrSupplementary()) continue;
							final Cigar cigar = rec.getCigar();
							if(cigar==null || cigar.isEmpty()) continue;
							// any clip supporting deletion ?
							boolean read_supports_cnv = false;
							final int breakpoint_distance= 10;

							// any clip on left ?
							if(cigar.isLeftClipped() &&
								rec.getUnclippedStart() < rec.getAlignmentStart() && 
								new SimpleInterval(ctx.getContig(), rec.getUnclippedStart(), rec.getAlignmentStart()).
								withinDistanceOf(breakPointLeft,breakpoint_distance)) {
									read_supports_cnv  = true;
								}
							// any clip on right ?
							if(	!read_supports_cnv &&
								cigar.isRightClipped() &&
								rec.getAlignmentEnd() < rec.getUnclippedEnd() && 
								new SimpleInterval(ctx.getContig(), rec.getAlignmentEnd(), rec.getUnclippedEnd()).
								withinDistanceOf(breakPointRight,breakpoint_distance)) {
									read_supports_cnv  = true;
								}
							
							if(read_supports_cnv) {
								n_reads_supporting_deletions++;
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
					
					final double bounds_cov[] = IntStream.concat(
							IntStream.range(0, array_mid_start),
							IntStream.range(array_mid_end,coverage.length)
							).mapToDouble(IDX->coverage[IDX]).
						toArray();
					final double medianBound = Percentile.average().evaluate(bounds_cov);
					if(Double.isNaN(medianBound) || medianBound==0) {
						final Genotype gt2 = GenotypeBuilder.createMissing(gt.getSampleName(), 2);
						genotypes.add(gt2);
						continue;
						}
					
					// divide coverage per medianBound
					final double normalized_coverage[] = new double[coverage.length];
					for(int i=0;i< normalized_coverage.length;++i) {
						normalized_coverage[i] = coverage[i] / medianBound;
					}
					
					final double midDepth = Percentile.average().evaluate(normalized_coverage, array_mid_start, array_mid_end);
					final GenotypeBuilder gb;
					switch (getCNVTypeFromNormDepth(midDepth)) {
						case DEL2: gb = new GenotypeBuilder(gt.getSampleName(),Arrays.asList(delAllele,delAllele));
								alleles.add(delAllele);
								n_samples_with_cnv++;
								break;
						case DEL1: gb = new GenotypeBuilder(gt.getSampleName(),Arrays.asList(refAllele,delAllele));
								alleles.add(delAllele);
								n_samples_with_cnv++;
								break;
						case WILD: gb = new GenotypeBuilder(gt.getSampleName(),Arrays.asList(refAllele,refAllele));
								break;
						case DUP1: gb = new GenotypeBuilder(gt.getSampleName(),Arrays.asList(refAllele,dupAllele));
								alleles.add(dupAllele);
								n_samples_with_cnv++;
								break;
						case DUP2: gb = new GenotypeBuilder(gt.getSampleName(),Arrays.asList(dupAllele,dupAllele));
							n_samples_with_cnv++;
							alleles.add(dupAllele);break;
						default: throw new IllegalStateException();
						}
					
					gb.attribute(formatCN.getID(), midDepth);
					
					//final Coverage covL = calculateCoverage(coverage,0, array_mid_start);
					//final Coverage covM = calculateCoverage(coverage,array_mid_start, array_mid_end);
					//final Coverage covR = calculateCoverage(coverage,array_mid_end, coverage.length);
					//final Coverage coBounds = calculateCoverage(coverage,IDX->IDX< array_mid_start || IDX>array_mid_end);
										
					gb.attribute(nReadsSupportingDel.getID(), n_reads_supporting_deletions);
					
					/*
					final Function<Double, Integer> toInt = DBL->DBL==null?-1:DBL.intValue();
					gb.attribute(leftMedianDepth.getID(),toInt.apply(covL.median));
					gb.attribute(midMedianDepth.getID(),toInt.apply(covM.median));
					gb.attribute(rightMedianDepth.getID(),toInt.apply(covR.median));
					gb.attribute(leftVarianceDepth.getID(),toInt.apply(covL.variance));
					gb.attribute(midVarianceDepth.getID(),toInt.apply(covM.variance));
					gb.attribute(rightVarianceDepth.getID(),toInt.apply(covR.variance));
					*/
					/*
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
					if(covL.median!=null && covL.median.doubleValue() < this.min_depth) {
						gtFilters.add("LOWDPL");
						}
					//low dp right
					if(covR.median!=null && covR.median.doubleValue() < this.min_depth) {
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
					*/
					genotypes.add(gb.make());
					}
				
				
				if(alleles.size()<=1) continue;
				vcb.alleles(alleles);
				vcb.noID();
				vcb.genotypes(genotypes);
				vcb.attribute(infoSVSamples.getID(), n_samples_with_cnv);
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
			CloserUtil.close(out);
			sample2bam.values().forEach(F->CloserUtil.close(F));
			}
		}
	
	
	public static void main(final String[] args) {
		new ValidateCnv().instanceMainWithExit(args);
		}
	}
