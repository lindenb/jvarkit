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

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

/**
BEGIN_DOC



END_DOC

 */

@Program(name="mergesv",
generate_doc=false,
description="Merge SV results",
keywords= {"cnv","indel","sv"}
)
public class MergeStructuralVariants extends Launcher{
	private static final Logger LOG = Logger.build(MergeStructuralVariants.class).make();
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile=null;
	@Parameter(names={"-d","--distance"},description="Two intervals are the same if their bounds are distant by less than xxx bases. "+ DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class ,splitter=com.github.lindenb.jvarkit.util.jcommander.NoSplitter.class)
	private int max_distance = 100;
	@Parameter(names={"-m","--max-length"},description="ignore variant longer than 'x' bases. Ignore this parameter if 'x' <=0 ")
	private int max_variant_length = -1;
	
	private final Allele REF_ALLELE = Allele.create("N", true);

	private class CnvCall implements Locatable
		{
		private final VariantContext ctx;
		private final Interval _startInterval;
		private final Interval _endInterval;
		
		boolean echoed_flag = false;
		CnvCall(final VariantContext ctx) {
			this.ctx = ctx;
			this._startInterval = _getInterval(ctx.getStart(), "CIPOS");
			this._endInterval = _getInterval(ctx.getEnd(), "CIEND");
			}
		public String getSample() {
			return ctx.getGenotype(0).getSampleName();
			}
		@Override
		public String getContig() {
			return ctx.getContig();
			}
		
		private Interval _getInterval(final int pos,final String att)
			{
			int x0 = 0;
			int x1 = 0;
			// CIPOS95(lumpy) Description="Confidence interval (95%) around POS for imprecise variants
			// CIPOS Confidence interval around POS for imprecise variant
			for(final String suffix :new String[]{"95",""}) { 
				if(this.ctx.hasAttribute(att + suffix)) {
					try {
						final List<Integer> list = this.ctx.getAttributeAsIntList(att, 0);
						x0 = list.get(0);
						x1 = list.get(1);
						break;
						}
					catch(final Throwable err)
						{
						
						}
					}
				}
			return new Interval(this.getContig(),Math.max(1,pos+x0),pos+x1);
			}
		
		public Interval getStartInterval() {
			return this._startInterval;
			}
		
		public Interval getEndInterval() {
			return this._endInterval;
			}
		
		@Override
		public int getStart() {
			return Math.min(
				getStartInterval().getStart(),
				getEndInterval().getStart()
				);
			}
		@Override
		public int getEnd() {
			return Math.max(
				getStartInterval().getEnd(),
				getEndInterval().getEnd()
				);
			}
		
		StructuralVariantType getType() {
			return ctx.getStructuralVariantType();
		}
		
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + getContig().hashCode();
			result = prime * result + getStart();
			result = prime * result + getEnd();
			result = prime * result + getType().hashCode();
			// sample NON
			return result;
			}
	
		@Override
		public boolean equals(final Object obj) {
			if (this == obj)  return true;
			if (obj == null)  return false;
			if (!(obj instanceof CnvCall)) {
				return false;
				}
			final CnvCall other = (CnvCall) obj;
			if (getEnd() != other.getEnd()) {
				return false;
				}
			if (getStart() != other.getStart()) {
				return false;
				}
			if (!getContig().equals(other.getContig())) {
				return false;
				}
			if (!getType().equals(other.getType())) {
				return false;
				}
			// sample NON
			return true;
			}
		@Override
		public String toString() {
			return getContig()+":"+getStart()+":"+getEnd();
		}
		}
	
	private boolean testOverlapping(
				final StructuralVariantType svType,
				final CnvCall a,final CnvCall b ) {
		
		return  a.ctx.getStructuralVariantType().equals(b.ctx.getStructuralVariantType()) &&
				a.getContig().equals(b.getContig()) &&
				a.getStartInterval().withinDistanceOf(b.getStartInterval(), this.max_distance) &&
				a.getEndInterval().withinDistanceOf(b.getEndInterval(), this.max_distance)
				;
				
		/*
		return a.getContig().equals(b.getContig()) &&
			   !(a.getEnd()<b.getStart() || a.getStart()>b.getEnd()) && 
			   Math.abs(a.getStart()-b.getStart()) <= this.max_distance &&
			   Math.abs(a.getEnd()-b.getEnd()) <= this.max_distance
			   ; */
		}
	
	private int getExtendFor(final VariantContext ctx) {
		return 0;
	}
	
	@Override
	public int doWork(final List<String> args) {
		if(this.max_distance<0) {
			LOG.error("bad max_distance :" +this.max_distance);
			return -1;
			}
		VariantContextWriter out = null;
		try {
			final List<File> inputs = IOUtils.unrollFiles2018(args);
			if(inputs.isEmpty()) {
				LOG.error("input is empty");
				return -1;
				}
			SAMSequenceDictionary dict = null;
			
			final Set<VCFHeaderLine> metadata = new HashSet<>();
			final IntervalTreeMap<List<CnvCall>> all_calls = new IntervalTreeMap<>();
			final Set<String> all_samples = new TreeSet<>();
			
			for(final File input: inputs) {
				
				final VCFFileReader vcfFileReader = new VCFFileReader(input,false);
				final VCFHeader header = vcfFileReader.getFileHeader();
				metadata.addAll(header.getMetaDataInInputOrder());
				
				final SAMSequenceDictionary dict0 = header.getSequenceDictionary();
				if(dict0==null || dict0.isEmpty()) {
					vcfFileReader.close();
					throw new JvarkitException.VcfDictionaryMissing(input);
					}
				if(dict==null)
					{
					dict = dict0;
					}
				else if(!SequenceUtil.areSequenceDictionariesEqual(dict, dict0))
					{
					LOG.error(JvarkitException.DictionariesAreNotTheSame.getMessage(dict0, dict));
					vcfFileReader.close();
					return -1;
					}
				if(header.getNGenotypeSamples()!=1) {
					LOG.error("Expected one and only one genotyped sample in "+input);
					vcfFileReader.close();
					return -1;
					}
				final String sample = header.getGenotypeSamples().get(0);
				vcfFileReader.iterator().stream().
					filter(V->V.getStructuralVariantType()!=null).
					filter(V->this.max_variant_length<=0 || (V.getEnd()-V.getStart()+1)<=this.max_variant_length).
					map(V->new CnvCall(V)).
					forEach(V->{
						final int extend = getExtendFor(V.ctx);
						
						final Interval key = new Interval(V.getContig(), 
								Math.max(1,V.getStart()-extend),
								V.getEnd()+extend
								);
						List<CnvCall> callList = all_calls.get(key);
						if(callList==null) {
							callList = new ArrayList<>();
							all_calls.put(key, callList);
							}
						callList.add(V);
						});
				
				vcfFileReader.close();
				all_samples.add(sample);
			}
			
			final Comparator<String> contigComparator = new ContigDictComparator(dict);
			
			
			final Comparator<CnvCall> comparator = (A,B)->{
				int i = contigComparator.compare(A.getContig(), B.getContig());
				if(i!=0) return i;
				i = Integer.compare(A.getStart(), B.getStart());
				if(i!=0) return i;
				i = Integer.compare(A.getEnd(), B.getEnd());
				if(i!=0) return i;
				i = A.getType().compareTo(B.getType());
				return i;
				};
				
			final List<CnvCall> intervals_list = 
					all_calls.
					values().
					stream().
					flatMap(L->L.stream()).
					sorted(comparator).
					collect(Collectors.toList());
			
			
			metadata.add(new VCFInfoHeaderLine("SAMPLES",
					VCFHeaderLineCount.UNBOUNDED,
					VCFHeaderLineType.String,
					"Samples carrying the SV"
					));
			metadata.add(new VCFInfoHeaderLine("NSAMPLES",
					1,
					VCFHeaderLineType.Integer,
					"Number of Samples carrying the SV"
					));
			
			metadata.add(new VCFInfoHeaderLine("SVLEN",
					1,
					VCFHeaderLineType.Integer,
					"SV length"
					));
			metadata.add(new VCFInfoHeaderLine("CIPOS",
					2,
					VCFHeaderLineType.Integer,
					"Confidence interval around POS for imprecise variants"
					));
			metadata.add(new VCFInfoHeaderLine("CIEND",
					2,
					VCFHeaderLineType.Integer,
					"Confidence interval around END for imprecise variants"
					));
			
			metadata.add(new VCFInfoHeaderLine("IMPRECISE",
					0,
					VCFHeaderLineType.Flag,
					"Imprecise structural variation"
					));
			
			metadata.add(new VCFFormatHeaderLine(
					"OV",1,
					VCFHeaderLineType.Integer,
					"Number calls (with different sample) overlapping this genotype"
					));
			
			final VCFHeader header = new VCFHeader(
					metadata,
					all_samples
					);
			header.setSequenceDictionary(dict);
			JVarkitVersion.getInstance().addMetaData(this, header);
			
			out =  super.openVariantContextWriter(this.outputFile);
			out.writeHeader(header);
			for(final CnvCall interval:intervals_list)
				{
				final Set<Allele> alleles = new HashSet<>();
				alleles.add(REF_ALLELE);
				
				final List<CnvCall> overlappingList = all_calls.
					getOverlapping(interval).
					stream().
					flatMap(L->L.stream()).
					filter(C->C.getType().equals(interval.getType())).
					filter(C->!C.echoed_flag).
					collect(Collectors.toList());
				
				if(overlappingList.isEmpty()) continue;
				
				final StructuralVariantType svType =  overlappingList.get(0).getType();
				
				final CnvCall baseCall = overlappingList.stream().
						filter(C->C.equals(interval)).
						findFirst().
						orElse(null);
				if(baseCall==null)
					{
					continue;
					}
				final List<CnvCall> callsToPrint = overlappingList.stream().
						filter(C->testOverlapping(svType,C,baseCall)).
						collect(Collectors.toList())
						;
				
				if(callsToPrint.isEmpty()) {
					throw new IllegalStateException();
				}
				
				final VariantContextBuilder vcb= new VariantContextBuilder();
				vcb.chr(baseCall.getContig());
				vcb.start(baseCall.getStart());
				vcb.stop(baseCall.getEnd());
				vcb.attribute(VCFConstants.END_KEY, baseCall.getEnd());
				vcb.attribute(VCFConstants.SVTYPE, baseCall.getType().name());
				vcb.attribute("SVLEN", (1+baseCall.getEnd()-baseCall.getStart()));
				
				for(int side=0;side<2;side++)
					{
					final Function<CnvCall,Integer> extractor;
					if(side==0)
						{
						extractor = C->C.getStart();
						}
					else
						{
						extractor = C->C.getEnd();
						}
					final List<Integer> list = Arrays.asList(
						callsToPrint.stream().
							mapToInt(C->extractor.apply(C)-extractor.apply(baseCall)).
							min().
							orElse(0),
						callsToPrint.stream().
							mapToInt(C->extractor.apply(C)-extractor.apply(baseCall)).
							max().
							orElse(0)
						);
					vcb.attribute(
							side==0?"CIPOS":"CIEND", 
							list
							);
					}
				vcb.attribute("IMPRECISE", true);
				
				final Map<String,Genotype> sample2gt = new HashMap<>(callsToPrint.size());
				for(final CnvCall call: callsToPrint)
					{
					if(sample2gt.containsKey(call.getSample()))
						{
						LOG.warn("Sample "+call.getSample()+" exits twice at the same loc " + call+" could be two small SV overlapping a big one.");
						continue;
						}
					call.echoed_flag  = true;
					
					final GenotypeBuilder gb = new GenotypeBuilder(call.ctx.getGenotype(0));
					final List<Allele> gtAlleles=call.ctx.getGenotype(0).getAlleles().stream().map(A->A.isReference()?REF_ALLELE:A).collect(Collectors.toList());
					gb.alleles(gtAlleles);
					alleles.addAll(gtAlleles);
					
					
					gb.attribute("OV",
							all_calls.getOverlapping(call).
								stream().
								flatMap(col->col.stream()).
								filter(C->!C.getSample().equals(call.getSample())).
								count()
							);
					

					
					
					sample2gt.put(call.getSample(),gb.make());
					
					}
				vcb.attribute("SAMPLES",new ArrayList<>(sample2gt.keySet()));
				vcb.attribute("NSAMPLES",sample2gt.size());
				vcb.genotypes(sample2gt.values());
				alleles.remove(Allele.NO_CALL);
				vcb.alleles(alleles);
				out.add(vcb.make());
				}
			
			all_calls.values().stream().
				flatMap(C->C.stream()).
				filter(C->!C.echoed_flag).
				forEach(C->LOG.warn("Bug: Not printed "+C));
			
			out.close();
			out=null;
			return 0;
		} catch(final Throwable err) {
			LOG.error(err);
			return -1;
		} finally {
			CloserUtil.close(out);
		}
		}
	
	
public static void main(final String[] args) {
	new MergeStructuralVariants().instanceMainWithExit(args);
}
}
