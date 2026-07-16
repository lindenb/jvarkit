/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.expansionhunter;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.stream.Collectors;

import org.apache.commons.math3.stat.descriptive.AbstractUnivariateStatistic;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.rank.Median;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.math.stats.FisherCasesControls;
import com.github.lindenb.jvarkit.math.stats.FisherExactTest;
import com.github.lindenb.jvarkit.pedigree.CasesControls;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFReader;
import htsjdk.variant.vcf.VCFRecordCodec;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

/*
 
##FORMAT=<ID=ADFL,Number=1,Type=String,Description="Number of flanking reads consistent with the allele">
##FORMAT=<ID=ADIR,Number=1,Type=String,Description="Number of in-repeat reads consistent with the allele">
##FORMAT=<ID=ADSP,Number=1,Type=String,Description="Number of spanning reads consistent with the allele">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=LC,Number=1,Type=Float,Description="Locus coverage">
##FORMAT=<ID=REPCI,Number=1,Type=String,Description="Confidence interval for REPCN">
##FORMAT=<ID=REPCN,Number=1,Type=String,Description="Number of repeat units spanned by the allele">
##FORMAT=<ID=SO,Number=1,Type=String,Description="Type of reads that support the allele; can be SPANNING, FLANKING, or INREPEAT meaning that the reads span, flank, or are fully contained in the repeat">


##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">
##INFO=<ID=REF,Number=1,Type=Integer,Description="Reference copy number">
##INFO=<ID=REPID,Number=1,Type=String,Description="Repeat identifier as specified in the variant catalog">
##INFO=<ID=RL,Number=1,Type=Integer,Description="Reference length in bp">
##INFO=<ID=RU,Number=1,Type=String,Description="Repeat unit in the reference orientation">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=VARID,Number=1,Type=String,Description="Variant identifier as specified in the variant catalog">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">
##INFO=<ID=AC_Hom,Number=A,Type=Integer,Description="Allele counts in homozygous genotypes">
##INFO=<ID=AC_Het,Number=A,Type=Integer,Description="Allele counts in heterozygous genotypes">
##INFO=<ID=AC_Hemi,Number=A,Type=Integer,Description="Allele counts in hemizygous genotypes">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">


 +---------------------+---------+-------+------+---------+-------+---------+-------------+-------+-------------------+
 | Sample              | Type    | ADFL  | ADIR | ADSP    | GT    | LC      | REPCI       | REPCN | SO                |
 +---------------------+---------+-------+------+---------+-------+---------+-------------+-------+-------------------+
 | X10I9DM             | HOM_REF | 20/20 | 0/0  | 42/42   | 0/0   | 42.2212 | 8-8/8-8     | 8/8   | SPANNING/SPANNING |
 | X10I9DN             | HOM_REF | 15/15 | 0/0  | 28/28   | 0/0   | 41.8945 | 8-8/8-8     | 8/8   | SPANNING/SPANNING |
 | X10I9DO             | HOM_REF | 15/15 | 0/0  | 25/25   | 0/0   | 30.543  | 8-8/8-8     | 8/8   | SPANNING/SPANNING |
 | X10I9DP             | HOM_REF | 11/11 | 0/0  | 21/21   | 0/0   | 36.5863 | 8-8/8-8     | 8/8   | SPANNING/SPANNING |
 | X10I9DQ             | HOM_REF | 6/6   | 0/0  | 11/11   | 0/0   | 23.6014 | 8-8/8-8     | 8/8   | SPANNING/SPANNING |

 */

/**
BEGIN_DOC
 
# Input
 
Input is a list of indexed vcf files or one file with the '.list' suffix containing the path to the vcfs
 

END_DOC
*/
@Program(name="expansionhuntermerge",
description="Merge Vcf from ExpansionHunter.",
keywords= {"vcf","merge","ExpansionHunter"},
creationDate="20210210",
modificationDate="20260714",
jvarkit_amalgamion = true
)
public class ExpansionHunterMerge extends Launcher {
	private static final Logger LOG = Logger.of( ExpansionHunterMerge.class);
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-p","--percentile"},description="percentile to use 'average' or 'median'")
	private String percentile_name="median";
	@Parameter(names={"--skip-filtered"},description="Skip filtered variants")
	private boolean skip_filtered=false;
	@Parameter(names={"--factor"},description="multiple median/mean value of controls by 'factor'. if median value=100, then we count case having a size greater than 100*factor for the burden test")
	private double factor_value = 1.0;
	@Parameter(names={"--types"},description="types of call to consider using FORMAT/SO. There can be a bias of size depending of the nature of the call: SPANNING|FLANKING|INREPEAT")
	private String call_types_str= "SPANNING,FLANKING,INREPEAT";

	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection = new WritingSortingCollection();
	@ParametersDelegate
	private WritingVariantsDelegate writingVariants = new WritingVariantsDelegate();
	@ParametersDelegate
	private CasesControls casescontrols= new CasesControls();
	
	
	private boolean badVariance(final AbstractUnivariateStatistic percentile,final double[] array) {
		if(array.length<10) return false;
		final double mean = percentile.evaluate(array);
		int n_lt=0;
		int n_gt=0;
		for(double v:array) {
			if(v> mean*factor_value) n_gt++;
			if(v< mean/factor_value) n_lt++;
			}
		double n_samples= 0.1*array.length;//10% of sample with problems
		if(n_gt>=n_samples || n_lt>=n_samples) return true;
		return false;
		}
	
	private AbstractUnivariateStatistic getPercentile() {
		if(this.percentile_name.equalsIgnoreCase("median")) {
			return new Median();
			}
		else
			{
			return new Mean();
			}
		}
	
	@Override
		public int doWork(final List<String> args) {
		 	SortingCollection<VariantContext> sorter = null;
			try {
				final List<Path> inputs = IOUtils.unrollPaths(args);
				if(inputs.isEmpty()) {
					LOG.error("no input file.");
					return -1;
					}
				if(!this.percentile_name.matches("median|average|mean")) {
					LOG.error("Bad percentile name not (median|average|mean) "+this.percentile_name);
					return -1;
					}
			 	final AbstractUnivariateStatistic percentile= getPercentile();

				
				if(this.factor_value<1.0) {
					LOG.error("Factor should be >=1.0 "+this.factor_value);
					return -1;
					}
				
				final Set<String> call_types=new HashSet<String>(3);
				for(final String s: this.call_types_str.split("[, ;\\|]+")) {
					if(StringUtils.isBlank(s)) continue;
					if(!s.matches("(SPANNING|FLANKING|INREPEAT)")) {
						LOG.error("bad type "+s+" in "+this.call_types_str);
						return -1;
						}
					call_types.add(s);
					}
				if(call_types.isEmpty()) {
					LOG.error("no valid call-type");
					return -1;
					}
				
				this.casescontrols.load();
				
				final String REPID="REPID";
				final String RU="RU";
				final String REPCN = "REPCN";
				final String SO = "SO";
				final String fakeSample="___FAKE";
				SAMSequenceDictionary dict = null;
				final Set<String> samples = new TreeSet<>();
				final Set<VCFHeaderLine> metaData = new HashSet<>();
				/** loop over each vcf, register the samples, check required INFO exist */
				for(final Path path:inputs) {
					try(VCFReader r=VCFReaderFactory.makeDefault().open(path, false)) {
						final VCFHeader header = r.getHeader();
						if(header.getInfoHeaderLine(REPID)==null) {
							LOG.error("missing INFO/"+REPID);
							return -1;
							}
						// check only one sample
						if(header.getNGenotypeSamples()!=1) {
							LOG.error("expected one and only one genotyped sample in "+path);
							return -1;
							}
						final String sn = header.getSampleNamesInOrder().get(0);
						// check in case/control
						if(!casescontrols.isEmpty() && !casescontrols.contains(sn)) {
							LOG.warning("sample "+sn+" is not in case/controls list. Skipping "+path);
							continue;
							}
						
						final SAMSequenceDictionary d = header.getSequenceDictionary();
						if(d!=null) {
							if(dict==null) {
								dict = d;
								}
							else
								{
								SequenceUtil.assertSequenceDictionariesEqual(d, dict);
								}
							}
						
						if(sn.equals(fakeSample)) {
							LOG.error(sn+" cannot be named "+fakeSample+" in "+path);
							return -1;
							}
						if(samples.contains(sn)) {
							LOG.error("duplicate sample "+sn+" in "+path);
							return -1;
							}
						metaData.addAll(header.getMetaDataInInputOrder());
						samples.add(sn);
						}
					}
				
				if(!this.casescontrols.isEmpty()) {
					this.casescontrols.retain(samples);
					if(this.casescontrols.isEmpty()) {
						LOG.error("no common sample between vcf and case/controls");
						return -1;
						}
					if(!this.casescontrols.haveCases()) {
						LOG.error("no common cases between vcf and case/controls");
						return -1;
						}
					if(!this.casescontrols.haveControls()) {
						LOG.error("no common control between vcf and case/controls");
						return -1;
						}
					}
				
				final VCFHeader tmpHeader = new VCFHeader(metaData, Arrays.asList(fakeSample));
				final VCFInfoHeaderLine sampleInfo = new VCFInfoHeaderLine("SRCSAMPLE",1,VCFHeaderLineType.String,"SRC SAMPLE");
				final VCFInfoHeaderLine info_header_line_gt= new VCFInfoHeaderLine("FISHER_GT",1,VCFHeaderLineType.Float,"Fisher test. Cases have a higer number of repeats that the "+this.percentile_name+" of controls");
				final VCFInfoHeaderLine info_header_line_lt= new VCFInfoHeaderLine("FISHER_LT",1,VCFHeaderLineType.Float,"Fisher test. Cases have a lower number of repeats that the  "+this.percentile_name+" of controls");
				final VCFInfoHeaderLine info_header_line_min= new VCFInfoHeaderLine("FISHER",1,VCFHeaderLineType.Float,"min fisher found between INFO/FISHER_GT and INFO/FISHER_GT");
				final VCFInfoHeaderLine info_header_line_gt_desc = new VCFInfoHeaderLine("FISHER_GT_COUNT",1,VCFHeaderLineType.String,"Detail of Fisher test. Cases have a higer number of repeats that the  "+this.percentile_name+" of controls. controls case-ref|case-alt|ctrl-ref|ctrl-alt" );
				final VCFInfoHeaderLine info_header_line_lt_desc = new VCFInfoHeaderLine("FISHER_LT_COUNT",1,VCFHeaderLineType.String,"Detail of Fisher test. Cases have a lower number of repeats that the  "+this.percentile_name+" of controls. controls case-ref|case-alt|ctrl-ref|ctrl-alt");
				final VCFInfoHeaderLine info_header_median_ctrl_size = new VCFInfoHeaderLine(this.percentile_name.toUpperCase()+"_CTRL_SIZE",1,VCFHeaderLineType.Float, this.percentile_name+" size in controls");
				final VCFInfoHeaderLine info_header_median_case_size = new VCFInfoHeaderLine(this.percentile_name.toUpperCase()+"_CASE_SIZE",1,VCFHeaderLineType.Float, this.percentile_name+" size in cases");
				final VCFInfoHeaderLine info_header_case_histo = new VCFInfoHeaderLine("CASE_HISTO",1,VCFHeaderLineType.String,"cases histogram");
				final VCFInfoHeaderLine info_header_ctrl_histo = new VCFInfoHeaderLine("CTRL_HISTO",1,VCFHeaderLineType.String,"cpntrols histogram");
				final VCFFilterHeaderLine filter_header_bad_ctrls= new VCFFilterHeaderLine("VARIANCE_CTRLS","Too much variance of size in controls");
				//final VCFFilterHeaderLine filter_header_bad_cases= new VCFFilterHeaderLine("VARIANCE_CASES","Too much variance of size in cases");
				
				
				if(!this.casescontrols.isEmpty()) {
					metaData.add(info_header_line_gt);
					metaData.add(info_header_line_lt);
					metaData.add(info_header_line_min);
					metaData.add(info_header_line_gt_desc);
					metaData.add(info_header_line_lt_desc);
					metaData.add(info_header_median_ctrl_size);
					metaData.add(info_header_median_case_size);
					metaData.add(info_header_case_histo);
					metaData.add(info_header_ctrl_histo);
					metaData.add(filter_header_bad_ctrls);
					//metaData.add(filter_header_bad_cases);
					}
				
				VCFStandardHeaderLines.addStandardFormatLines(metaData, true,VCFConstants.GENOTYPE_FILTER_KEY);
				VCFStandardHeaderLines.addStandardInfoLines(metaData, true,VCFConstants.END_KEY);

				
				// contig sorter
				final Comparator<String> ctgComparator= dict==null?
						(A,B)->A.compareTo(B):
						new ContigDictComparator(dict)
						;
				// get required attribute
				final BiFunction<VariantContext, String, String> getAtt = (V,A)->{
					final String s1 = V.getAttributeAsString(A, "");
					if(StringUtils.isBlank(s1)) throw new IllegalStateException("INFO/"+A+" missing in "+V);
					return s1;
					};
				
				// variant sorter : sort on chromosome, pos, REPID
				final Comparator<VariantContext> comparator1= (V1,V2) ->{
					String s1 = V1.getContig();
					String s2 = V2.getContig();
					int i= ctgComparator.compare(s1,s2);
					if(i!=0) return i;
				
					i= Integer.compare(V1.getStart(), V2.getStart());
					if(i!=0) return i;
					
					s1 = getAtt.apply(V1,REPID);
					s2 = getAtt.apply(V2,REPID);
					return s1.compareTo(s2);
					};
				
				// second comparator, sort on sample identifier
				final Comparator<VariantContext> comparator2 = (V1,V2) ->{
					final int i= comparator1.compare(V1, V2);
					if(i!=0) return i;
					final String s1 = getAtt.apply(V1,sampleInfo.getID());
					final String s2 = getAtt.apply(V2,sampleInfo.getID());
					return s1.compareTo(s2);
					};
				
				/**
				 *  LOOP Over each VCF
				 * 
				 */
				tmpHeader.addMetaDataLine(sampleInfo);
				sorter = SortingCollection.newInstance(
                        VariantContext.class,
                        new VCFRecordCodec(tmpHeader, false),
                        comparator2,
                        this.writingSortingCollection.getMaxRecordsInRam(),
                        this.writingSortingCollection.getTmpPaths()
                        );
				sorter.setDestructiveIteration(true);
				for(final Path path:inputs) {
					LOG.info("adding "+path);
					try(VCFReader r=VCFReaderFactory.makeDefault().open(path, false)) {
						final VCFHeader header = r.getHeader();
						final String sn = header.getSampleNamesInOrder().get(0);
						if(!casescontrols.isEmpty() && !casescontrols.contains(sn)) {
							LOG.warning("sample "+sn+" is not in case/controls list. Skipping "+path);
							continue;
							}
						/* loop over each variant */
						try(CloseableIterator<VariantContext> iter = r.iterator()) {
							while(iter.hasNext()) {
								final VariantContext ctx = iter.next();
								/* rename the sample of the variant to fakeSample, same the original name in INFO */
								final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
								final Genotype gt0 = ctx.getGenotype(0);
								final Genotype gt = new GenotypeBuilder(gt0).name(fakeSample).make();
								vcb.genotypes(Collections.singletonList(gt));
								vcb.attribute(sampleInfo.getID(), gt0.getSampleName());
								sorter.add(vcb.make());
								}
							}
						}
					}
				sorter.doneAdding();
				
				
				
				final VCFHeader outputHeader = new VCFHeader(metaData, samples);
				if(dict!=null) outputHeader.setSequenceDictionary(dict);
				JVarkitVersion.getInstance().addMetaData(this, outputHeader);
				try(VariantContextWriter w= writingVariants.dictionary(dict).open(this.outputFile)) {
					w.writeHeader(outputHeader);
					try(CloseableIterator<VariantContext> iter = sorter.iterator()) {
						final EqualRangeIterator<VariantContext> eq = new EqualRangeIterator<>(iter,comparator1);
						while(eq.hasNext()) {
							final List<VariantContext> calls = eq.next();
							final VariantContext first = calls.get(0);
							final Map<String,VariantContext> sample2vc = new HashMap<>(calls.size());
							for(VariantContext ctx:calls) {
								final String sn=getAtt.apply(ctx,sampleInfo.getID());
								sample2vc.put(sn,ctx);
								}
							
							final Set<String> filters = new HashSet<String>();
							
							final Set<Allele> altAllelesSet = calls.stream().
									flatMap(V->V.getGenotypes().stream()).
									flatMap(GT->GT.getAlleles().stream()).
									filter(A->!(A.isReference() || A.isNoCall())).
									collect(Collectors.toSet());
								if(altAllelesSet.isEmpty()) continue;
								final List<Allele> altAllelesList = new ArrayList<>(altAllelesSet);
								final List<Allele> vcAlleles = new ArrayList<>(altAllelesList.size()+1);
								vcAlleles.add(first.getReference());
								vcAlleles.addAll(altAllelesList);

								final Map<String,Double> sample2mean_size= (this.casescontrols.isEmpty()?null: new HashMap<>(this.casescontrols.getTotalCount()));
								
								final List<Genotype> genotypes = new ArrayList<>(samples.size());
								for(final String sn : samples) {
									final VariantContext vcs = sample2vc.getOrDefault(sn, null);
									final Genotype gt;
									if(vcs==null) {
										gt=GenotypeBuilder.createMissing(sn,2);
										}
									else
										{
										final GenotypeBuilder gtb= new GenotypeBuilder(vcs.getGenotype(0)).
												name(sn);
										
										if(vcs.isFiltered()) {
											gtb.filter("LowQual");
											}
										
										gt = gtb.make();
										if(sample2mean_size!=null && 
												!gt.isNoCall() && 
												this.casescontrols.contains(sn) &&  
												gt.hasExtendedAttribute(REPCN) && 
												gt.hasExtendedAttribute(SO) && 
												(!this.skip_filtered ||!vcs.isFiltered())) {
											final String[] repcn=CharSplitter.SLASH.split(gt.getExtendedAttribute(REPCN, "").toString());
											final String[] so   =CharSplitter.SLASH.split(gt.getExtendedAttribute(SO   , "").toString());
											final List<Integer> values=new ArrayList<>(repcn.length);
											for(int i=0;i< repcn.length && i< so.length;++i) {
												if(repcn[i].equals(".")) continue;
												if(!call_types.contains(so[i])) continue;
												values.add(Integer.parseInt(repcn[i]));
												}
											
											if(!values.isEmpty()) {
													sample2mean_size.put(
														sn,
														values.stream()
															.mapToDouble(I->I.doubleValue())
															.average()
															.getAsDouble()
													);
												}
											}
										}
									genotypes.add(gt);
									}
								final VariantContextBuilder vcb = new VariantContextBuilder(null,
										first.getContig(),
										first.getStart(), first.getEnd(),
										vcAlleles
										);
								
								if(sample2mean_size!=null && !sample2mean_size.isEmpty()) {
									// create an histogram for for case and controls
									for(int side=0;side<2;++side) {
										final String histo = this.casescontrols.get(side).stream()
											.filter(S->sample2mean_size.containsKey(S))
											.map(S->sample2mean_size.get(S))
											.collect(Collectors.groupingBy(Function.identity(), Collectors.counting()))
											.entrySet()
											.stream()
											.map(KV->String.valueOf(KV.getKey())+":"+KV.getValue())
											.collect(Collectors.joining("|"))
											;
										if(!StringUtils.isBlank(histo)) {
											vcb.attribute((side==0?info_header_case_histo.getID():info_header_ctrl_histo.getID()), histo);
											}
										}
									
									
									
									
									
									/* median case size*/
										{
										final double[] case_sizes=	this.casescontrols.getCases().stream()
												.filter(S->sample2mean_size.containsKey(S))
												.mapToDouble(S->sample2mean_size.get(S))
												.sorted()
												.toArray();
										
										if(case_sizes.length>0) {
											final double median_case_size=  percentile.evaluate(case_sizes);
											vcb.attribute(info_header_median_case_size.getID(),median_case_size);
											
											/*
											if(badVariance(percentile,case_sizes)) {
												filters.add(filter_header_bad_cases.getID());
												}
											*/
											}
										}
										
									/* median ctrl size*/
									final double[] ctrl_mean_sizes = this.casescontrols.getControls().stream()
										.filter(S->sample2mean_size.containsKey(S))
										.mapToDouble(S->sample2mean_size.get(S))
										.sorted()
										.toArray();
									
									if(badVariance(percentile,ctrl_mean_sizes)) {
										filters.add(filter_header_bad_ctrls.getID());
										}
									
									/* check fisher for population lower or greater than reference population */
									if(ctrl_mean_sizes.length>0) {
										final double median_ctrl_size=  percentile.evaluate(ctrl_mean_sizes);
										vcb.attribute(info_header_median_ctrl_size.getID(),median_ctrl_size);
										
										
										// use only the case/control where we found a genotype/size
										final CasesControls casescontrols2 = this.casescontrols.clone().retain(sample2mean_size.keySet());
										if(!casescontrols2.isEmpty()) {
											final FisherCasesControls fisherFactory_gt = new FisherCasesControls(casescontrols2);
											final FisherCasesControls fisherFactory_lt = new FisherCasesControls(casescontrols2);
											
											for(String sn: sample2mean_size.keySet()) {
												double sample_size = sample2mean_size.get(sn);
												if(sample_size> median_ctrl_size*this.factor_value ) {
													fisherFactory_gt.accept(sn);
													}
												if(sample_size < median_ctrl_size/this.factor_value) {
													fisherFactory_lt.accept(sn);
													}
												}
											FisherExactTest fisher = fisherFactory_gt.getFisherExactTest();
											double min = fisher.getAsDouble();
											vcb.attribute(info_header_line_gt.getID(), fisher.getAsDouble());
											vcb.attribute(info_header_line_gt_desc.getID(), fisherFactory_gt.join("|"));
											
											
											fisher = fisherFactory_lt.getFisherExactTest();
											min = Math.min(min,fisher.getAsDouble());
											vcb.attribute(info_header_line_lt.getID(), fisher.getAsDouble());
											vcb.attribute(info_header_line_lt_desc.getID(), fisherFactory_lt.join("|"));
											
											vcb.attribute(info_header_line_min.getID(), min);
											}
										}
									}
								
								if(filters.isEmpty()) {
									vcb.passFilters();
									}
								else
									{
									vcb.filters(filters);
									}
								
								vcb.attribute(VCFConstants.END_KEY, first.getEnd());
								vcb.attribute(REPID, getAtt.apply(first, REPID));
								vcb.attribute(RU, getAtt.apply(first, RU));
								vcb.genotypes(genotypes);
								w.add(vcb.make());
							}
						eq.close();
						}
					}
				sorter.cleanup();
				return 0;
				}
			catch(final Throwable err) {
				LOG.error(err);
				return -1;
				}
			finally
				{
			
				}
			}
	
	public static void main(final String[] args)
		{
		new ExpansionHunterMerge().instanceMainWithExit(args);
		}

	}
