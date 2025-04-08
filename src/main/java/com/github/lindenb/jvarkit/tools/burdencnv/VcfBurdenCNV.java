/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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


History:
* 2014 creation
* 2015 moving to knime

*/
package com.github.lindenb.jvarkit.tools.burdencnv;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.OptionalDouble;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.Predicate;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.converter.FractionConverter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.stats.FisherExactTest;
import com.github.lindenb.jvarkit.pedigree.CasesControls;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.JexlVariantPredicate;
import com.github.lindenb.jvarkit.variant.vcf.VcfHeaderExtractor;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;

/**

BEGIN_DOC

Still in beta. Do not use.

## INPUT

input is the set of vcf.gz indexed files or
one file with the suffix '.list' containing the path to the vcfs

## Example

```
wget -O - "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.bed.gz" | gunzip -c |\
	awk -F '\t' '(NR>1 && $45="DEL" && $177>0.01 && $177!="NA") {L=int($3)-int($2);printf("%s\t%s\t%s\t%s\n",$1,$2,$3,(L>100?0.7:0.1));}' |\
	gzip > work/known.bed.gz


find MANTADIR -name "*.diploidSV.vcf.gz" > work/vcfs.list
java -jar jvarkit.jar vcfburdencnv \
 	--cases work/cases.txt \
 	--controls work/controls.txt \
 	--known work/known.bed.gz \
 	work/vcfs.list
```


END_DOC
*/
@Program(name="vcfburdencnv",
		description="Burden on CNV (experimental)",
		keywords={"vcf","burden","cnv"},
		creationDate = "20250404",
		modificationDate = "20250408",
		jvarkit_amalgamion = true
		)
public class VcfBurdenCNV
	extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfBurdenCNV.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@ParametersDelegate
	private CasesControls caseControls = new CasesControls();
	@Parameter(names={"--treshold"},description="do not display BED line if fisher > 'treshold'")
	private double p_value_treshold=1E-5;
	@Parameter(names={"--exclude-bed"},description="BED regions to exclude")
	private Path bedExcludePath=null;
	@Parameter(names={"--include-bed"},description="BED regions to include")
	private Path bedIncludePath=null;
	@Parameter(names={"--jexl"},description=JexlVariantPredicate.PARAMETER_DESCRIPTION)
	private String jexlExptr = null;
	@Parameter(names={"--known"},description="BED file containing the known frequent CNV. The 4th column must be a number between 0 and 1, which well be the mutual overlap between the known (e.g gnomad) variant and the user's variant. Default  value for column 4 uses --default-overlap")
	private Path bedKnownPath=null;
	@Parameter(names={"--default-overlap"}, description="default fraction overlap for option --known", converter = FractionConverter.class, splitter = NoSplitter.class)
	private  double default_fraction_overlap = 0.75;


	/** load bed for exclude.include bed files for given chrom */
	private final IntervalTree<Boolean> loadBed(final Path path,final SAMSequenceDictionary dict, final String contig) throws IOException {
		if(path==null) return null;
		final IntervalTree<Boolean> iTree=new IntervalTree<>();
		try(BedLineReader blr = new BedLineReader(path)) {
			blr.setContigNameConverter(ContigNameConverter.fromOneDictionary(dict));
			blr.stream().
				filter(B->B.getContig().equals(contig)).
				forEach(B->iTree.put(B.getStart(), B.getEnd(), Boolean.TRUE));
			}
		return iTree;
		}
	/** load bed for exclude.include bed files for given chrom */
	private final IntervalTree<Double> loadKnown(final SAMSequenceDictionary dict, final String contig) throws IOException {
		if(this.bedKnownPath==null) return null;
		final IntervalTree<Double> iTree=new IntervalTree<>();
		try(BedLineReader blr = new BedLineReader(this.bedKnownPath)) {
			blr.setContigNameConverter(ContigNameConverter.fromOneDictionary(dict));
			blr.stream().
				filter(B->B.getContig().equals(contig)).
				forEach(B->{
					final String column4 = B.getOrDefault(3, "");
					final double fract;
					if(StringUtils.isBlank(column4) || column4.equals(".")) {
						fract = this.default_fraction_overlap;
						}
					else
						{
						fract = Double.parseDouble(column4);
						}
					IntervalTree.Node<Double> n=iTree.find(B.getStart(), B.getEnd());
					if(n==null || n.getValue().compareTo(fract)>0) {
						iTree.put(B.getStart(), B.getEnd(),fract);
						}
				});
			}
		return iTree;
		}
	// TODO: give a chance to convert DELLY (no INFO/END) variants
	private Locatable variantToLocatable(final VariantContext ctx) {
		return ctx;
		}

	 /** return odd ratio https://en.wikipedia.org/wiki/Odds_ratio */
    private static OptionalDouble getOddRatio(	
    		final int case_ref,
    		final int case_alt, 
    		final int ctrl_ref, 
    		final int ctrl_alt 
    		) {
    	final int nCases = (case_ref+case_alt);
    	if(nCases==0) return OptionalDouble.empty();
    	final int nCtrls = (ctrl_ref+ctrl_alt);
    	if(nCtrls==0) return OptionalDouble.empty();
    	final double f1 = case_alt/(double)nCases;
    	final double f2 = ctrl_alt/(double)nCtrls;
    	if(f2==0) return OptionalDouble.empty();
    	return OptionalDouble.of(f1/f2);
    	}
	@Override
	public int doWork(final List<String> args) {
		try 
			{
			final List<Path> vcfPaths = new ArrayList<>();
			this.caseControls.load().checkHaveCasesControls();
			for(Path vcf:IOUtils.unrollPaths(args)) {
				IOUtil.assertFileIsReadable(vcf);
				vcfPaths.add(vcf);
				}
				
			SAMSequenceDictionary dict=null;
			// load every sample in each VCF
			final  Set<String> distinct_samples = new TreeSet<>();
			for(Path p : vcfPaths) {
				VCFHeader h= VcfHeaderExtractor.decode(p);
				if(h.getNGenotypeSamples()!=1) throw new IllegalArgumentException();
				String sn = h.getGenotypeSamples().get(0);
				if(!this.caseControls.contains(sn)) continue;
				SAMSequenceDictionary dict0 = SequenceDictionaryUtils.extractRequired(h);
				if(dict==null) {
					dict = dict0;
					}
				else
					{
					SequenceUtil.assertSequenceDictionariesEqual(dict, dict0);
					}
				distinct_samples.add(sn);
				}
			this.caseControls.retain(distinct_samples);
			this.caseControls.checkHaveCasesControls();
			final List<String> all_samples = new ArrayList<>(distinct_samples);
			final Map<String,Integer> sample2idx = new HashMap<>(all_samples.size());
			for(int i=0;i< all_samples.size();++i) {
				sample2idx.put(all_samples.get(i),i);
				}
			
			final BitSet is_case_bitset = new BitSet(all_samples.size());
			for(int i=0;i< all_samples.size();++i) {
				is_case_bitset.set(i,this.caseControls.isCase(all_samples.get(i)));
				}
			
			final Predicate<VariantContext> acceptVariant = JexlVariantPredicate.create(this.jexlExptr);
			
			try(PrintWriter out = super.openPathOrStdoutAsPrintWriter(outputFile)) {
				out.print("#CHROM");
				out.print("\t");
				out.print("start");
				out.print("\t");
				out.print("end");
				out.print("\t"); out.print("case_ref");
				out.print("\t"); out.print("case_alt");
				out.print("\t"); out.print("ctrl_ref");
				out.print("\t"); out.print("ctrl_alt");
				out.print("\t"); out.print("p_value");
				out.print("\t"); out.print("odd_ratio");
				out.println();
				
				
				for(SAMSequenceRecord ssr: dict.getSequences()) {
					final IntervalTree<BitSet> intervalTree=new IntervalTree<>();
					final IntervalTree<Double> knownIntervalTree = loadKnown(dict,ssr.getContig());

					for(Path p : vcfPaths) {
						try(VCFReader in = new VCFFileReader(p, true)) {
							final VCFHeader h= VcfHeaderExtractor.decode(p);
							final String sn = h.getGenotypeSamples().get(0);
							final Integer idx = sample2idx.get(sn);
							if(idx==null) continue;
							try(CloseableIterator<VariantContext> iter=in.query(ssr)) {
								while(iter.hasNext()) {
									final VariantContext ctx = iter.next();
									if(!acceptVariant.test(ctx)) continue;
									BitSet bitset = null;
									final Locatable ctx_loc = variantToLocatable(ctx);
									if(ctx_loc==null) continue;
									
									// overlap with known ?
									if(knownIntervalTree!=null) {
										boolean ok_with_known=true;
										for(Iterator<IntervalTree.Node<Double>> iter2=knownIntervalTree.overlappers(ctx_loc.getStart(), ctx_loc.getEnd());
											iter2.hasNext();) {
											final IntervalTree.Node<Double> node = iter2.next();
											final double shared_len = CoordMath.getLength(
													Math.max(ctx_loc.getStart(), node.getStart()),
													Math.min(ctx_loc.getEnd(), node.getEnd())
													);
											final double frac = node.getValue();
											if( CoordMath.getLength(ctx_loc.getStart(), ctx_loc.getEnd())/shared_len < frac ||
												ctx_loc.getLengthOnReference()/shared_len < frac)
												{
													ok_with_known=false;
													break;
												}
																							
											}
										if(!ok_with_known) continue;
										}
									
									final IntervalTree.Node<BitSet> samples = intervalTree.find(ctx_loc.getStart(), ctx_loc.getEnd());
									if(samples==null) {
										bitset = new BitSet(all_samples.size());
										intervalTree.put(ctx_loc.getStart(), ctx_loc.getEnd(),bitset);
										}
									else
										{
										bitset = samples.getValue();
										}
									bitset.set(idx, true);
									}
								}
							}
						}
					
					final IntervalTree<Boolean> excludeIntervalTree = loadBed(this.bedExcludePath,dict,ssr.getContig());
					final IntervalTree<Boolean> includeIntervalTree = loadBed(this.bedIncludePath,dict,ssr.getContig());
					for(int x=0;x< ssr.getLengthOnReference();x++) {
						if(excludeIntervalTree!=null && excludeIntervalTree.overlappers(x+1, x+1).hasNext()) continue;
						if(includeIntervalTree!=null && !includeIntervalTree.overlappers(x+1, x+1).hasNext()) continue;
						
						
						int case_ref = 0; 
						int case_alt = 0; 
						int ctrl_ref = 0; 
						int ctrl_alt = 0; 
						final  BitSet samples_at_x= new  BitSet(all_samples.size());
						final Iterator<IntervalTree.Node<BitSet>> iterator = intervalTree.overlappers(x+1,x+1);
			            while (iterator.hasNext()) {
			                final BitSet bitset = iterator.next().getValue();
			                samples_at_x.or(bitset);
			            	}
				            
					
						for(int i=0;i< all_samples.size();i++) {
		                	final boolean is_case = is_case_bitset.get(i);
		                	if(samples_at_x.get(i)) {
		                		if(is_case) case_alt++;
		                		else ctrl_alt++;
		                		}
		                	else
		                		{
		                		if(is_case) case_ref++;
		                		else ctrl_ref++;
		                		}
		                	}
				        final double p_value = FisherExactTest.compute(case_ref, case_alt, ctrl_ref, ctrl_alt).getAsDouble();
						if(p_value> p_value_treshold) continue;
						final OptionalDouble oddRatio = getOddRatio(case_ref, case_alt, ctrl_ref, ctrl_alt);
						out.print(ssr.getContig());
						out.print("\t");
						out.print(x);
						out.print("\t");
						out.print(x+1);
						out.print("\t"); out.print(case_ref);
						out.print("\t"); out.print(case_alt);
						out.print("\t"); out.print(ctrl_ref);
						out.print("\t"); out.print(ctrl_alt);
						out.print("\t"); out.print(p_value);
						out.print("\t"); out.print(oddRatio.isPresent()?String.valueOf(oddRatio.getAsDouble()):"NA");
						out.println();
						}
					
					} // end loop over ssr
				out.flush();
				}//end of write
			return 0;
			}			
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	
	public static void main(final String[] args)
		{
		new VcfBurdenCNV().instanceMainWithExit(args);
		}
	}
