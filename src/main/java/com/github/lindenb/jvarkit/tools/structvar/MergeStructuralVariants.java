/*
The MIT License (MIT)

Copyright (c) 2018 Pierre Lindenbaum

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
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
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
	@Parameter(names={"-d","--distance"},description="two intervals are the same if their bounds are distant by less than xxx bases")
	private int max_distance = 100;
	
	private final Allele REF_ALLELE = Allele.create("N", true);

	private class CnvCall implements Locatable
		{
		private final VariantContext ctx;
		boolean echoed_flag = false;
		CnvCall(final VariantContext ctx) {
			this.ctx = ctx;
			}
		public String getSample() {
			return ctx.getGenotype(0).getSampleName();
			}
		@Override
		public String getContig() {
			return ctx.getContig();
			}
		@Override
		public int getStart() {
			return ctx.getStart();
			}
		@Override
		public int getEnd() {
			return ctx.getEnd();
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
		}
	
	private boolean testOverlapping(final Locatable a,final Locatable b ) {		
		return a.getContig().equals(b.getContig()) &&
			   !(a.getEnd()<b.getStart() || a.getStart()>b.getEnd()) && 
			   Math.abs(a.getStart()-b.getStart()) <= this.max_distance &&
			   Math.abs(a.getEnd()-b.getEnd()) <= this.max_distance
			   ;
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
					map(V->new CnvCall(V)).
					forEach(V->{
						final Interval key = new Interval(V.getContig(), V.getStart(), V.getEnd());
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
				final CnvCall baseCall = overlappingList.stream().
						filter(C->C.equals(interval)).
						findFirst().
						orElse(null);
				if(baseCall==null)
					{
					continue;
					}
				final List<CnvCall> callsToPrint = overlappingList.stream().
						filter(C->testOverlapping(C,baseCall)).
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
					alleles.addAll(call.ctx.getGenotype(0).getAlleles());
					
					
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
