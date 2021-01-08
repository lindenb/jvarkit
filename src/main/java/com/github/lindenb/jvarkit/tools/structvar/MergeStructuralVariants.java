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
package com.github.lindenb.jvarkit.tools.structvar;

import java.io.File;
import java.nio.file.Path;
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
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFReader;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

/**
BEGIN_DOC



END_DOC

 */

@Program(name="mergesv",
generate_doc=false,
description="Merge SV results",
keywords= {"cnv","indel","sv"},
modificationDate="20190815"
)
public class MergeStructuralVariants extends Launcher{
	private static final Logger LOG = Logger.build(MergeStructuralVariants.class).make();
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile=null;
	@Parameter(names={"-d","--distance"},description="Two BND variants are the same if their bounds are distant by less than xxx bases. "+ DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class ,splitter=com.github.lindenb.jvarkit.util.jcommander.NoSplitter.class)
	private int max_distance = 100;
	@Parameter(names={"-f","--fraction"},description="Two CNV/DEL/.. variants are the same if they share 'x' fraction of their size.")
	private double max_fraction = 0.75;

	@Parameter(names={"-m","--max-length"},description="ignore variant longer than 'x' bases. Ignore this parameter if 'x' <=0 ")
	private int max_variant_length = -1;

	private class VcfInput {
		//final Path path;
		final SAMSequenceDictionary dict;
		final String sample;
		final IntervalTreeMap<CnvCall> callMap=new IntervalTreeMap<>();
		
		VcfInput(final Path path) {
			//this.path = path;
			final VCFReader vcfFileReader = VCFReaderFactory.makeDefault().open(path,false);
			final VCFHeader header = vcfFileReader.getHeader();
			
			this.dict = SequenceDictionaryUtils.extractRequired(header);
			
			if(!header.hasGenotypingData()) {
				this.sample = path.toString();
				}
			else if(header.getNGenotypeSamples()!=1) {
				CloserUtil.close(vcfFileReader);
				throw new IllegalArgumentException("Expected one and only one genotyped sample in "+path);
				}
			else
				{
				this.sample = header.getGenotypeSamples().get(0);
				}
			vcfFileReader.iterator().stream().
				filter(V->V.getStructuralVariantType()!=null).
				filter(V->max_variant_length<=0 || (V.getEnd()-V.getStart()+1)<= max_variant_length).
				map(V->new CnvCall(V,sample)).
				forEach(V->{
					final int extend = getExtendFor(V);
					
					final Interval key = new Interval(V.getContig(), 
							Math.max(1,V.getStart()-extend),
							V.getEnd()+extend
							);
					this.callMap.put(key,V);
					});
			CloserUtil.close(vcfFileReader);
			}
		}

	private class CnvCall implements Locatable
		{
		//private final VariantContext ctx;
		private final Interval _startInterval;
		private final Interval _endInterval;
		private final String sample;
		private final String svType;
		@SuppressWarnings("unused")
		private final GenotypeType genotypeType;
		boolean echoed_flag = false;
		CnvCall(final VariantContext ctx,final String sample) {
			//this.ctx = ctx;
			this.sample = sample;
			this.svType = ctx.getAttributeAsString(VCFConstants.SVTYPE, "");
			this._startInterval = _getInterval(ctx,ctx.getStart(), "CIPOS");
			this._endInterval = _getInterval(ctx,ctx.getEnd(), "CIEND");
			if(ctx.hasGenotypes()) {
				genotypeType = ctx.getGenotype(0).getType();
				}
			else
				{
				genotypeType = GenotypeType.HET;
				}
			}
		public String getSample() {
			return sample;
			}
		@Override
		public String getContig() {
			return this._startInterval.getContig();
			}
		
		private Interval _getInterval(final VariantContext ctx,final int pos,final String att)
			{
			int x0 = 0;
			int x1 = 0;
			// CIPOS95(lumpy) Description="Confidence interval (95%) around POS for imprecise variants
			// CIPOS Confidence interval around POS for imprecise variant
			for(final String suffix :new String[]{"95",""}) { 
				if(ctx.hasAttribute(att + suffix)) {
					try {
						final List<Integer> list = ctx.getAttributeAsIntList(att, 0);
						x0 = list.get(0);
						x1 = list.get(1);
						break;
						}
					catch(final Throwable err)
						{
						throw new IllegalArgumentException(err);
						}
					}
				}
			return new Interval(ctx.getContig(),Math.max(1,pos+x0),pos+x1);
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
		
		String getType() {
			return svType;
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
	
	private boolean testOverlapping(final CnvCall a,final CnvCall b ) {
		if(!a.svType.equals(b.svType)) return false;
		if(!a.getContig().equals(b.getContig())) return false;
		
		if(a.svType.equals("BND")) {
			return  a.getStartInterval().withinDistanceOf(b.getStartInterval(), this.max_distance) &&
					a.getEndInterval().withinDistanceOf(b.getEndInterval(), this.max_distance)
					;
			}
		else
			{
			final Interval interval1 = new Interval(a);
			final Interval interval2 = new Interval(b);
			if(!interval1.overlaps(interval2)) return false;
			int p1 = Math.max(interval1.getStart(),interval2.getStart());
			int p2 = Math.min(interval1.getEnd(),interval2.getEnd());
			double len = CoordMath.getLength(p1,p2);
			if(len/interval1.getLengthOnReference() < this.max_fraction ) return false; 
			if(len/interval2.getLengthOnReference() < this.max_fraction ) return false; 
			return true;
			}
		}
	
	private int getExtendFor(final CnvCall ctx) {
		return 0;
	}
	
	private final List<VcfInput> vcfFilesInput = new ArrayList<>();

	@SuppressWarnings("resource")
	@Override
	public int doWork(final List<String> args) {
		if(this.max_distance<0) {
			LOG.error("bad max_distance :" +this.max_distance);
			return -1;
			}
		VariantContextWriter out = null;
		try {
			final List<Path> inputPaths=(IOUtils.unrollPaths(args));
			if(inputPaths.isEmpty()) {
				LOG.error("input is empty");
				return -1;
				}
			SAMSequenceDictionary dict = null;
			
			final Set<VCFHeaderLine> metadata = new HashSet<>();
			
			for(int i=0;i< inputPaths.size();i++) {
				final Path input = inputPaths.get(i);
				LOG.info("reading ("+(i+1)+"/"+inputPaths.size()+") "+input);
				final VcfInput vcfInput = new VcfInput(input);
				this.vcfFilesInput.add(vcfInput);
				
				if(dict==null)
					{
					dict = vcfInput.dict;
					}
				else if(!SequenceUtil.areSequenceDictionariesEqual(dict, vcfInput.dict))
					{
					LOG.error(JvarkitException.DictionariesAreNotTheSame.getMessage(vcfInput.dict, dict));
					return -1;
					}
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
					this.vcfFilesInput.stream().
					flatMap(F->F.callMap.values().stream()).
					sorted(comparator).
					collect(Collectors.toList());
			
			VCFStandardHeaderLines.addStandardFormatLines(metadata, true, 
					VCFConstants.GENOTYPE_KEY
					);
			VCFStandardHeaderLines.addStandardInfoLines(metadata, true, 
					VCFConstants.END_KEY
					);
			
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
			
			/*metadata.add(new VCFFormatHeaderLine(
					"OV",1,
					VCFHeaderLineType.Integer,
					"Number calls (with different sample) overlapping this genotype"
					));*/
			
			metadata.add(new VCFInfoHeaderLine(
					VCFConstants.SVTYPE,1,
					VCFHeaderLineType.String,
					"SV type"
					));
			
			final VCFHeader header = new VCFHeader(
					metadata,
					(Set<String>)this.vcfFilesInput.stream().map(F->F.sample).collect(Collectors.toCollection(TreeSet::new))
					);
			
			header.setSequenceDictionary(dict);
			JVarkitVersion.getInstance().addMetaData(this, header);
			
			out =  super.openVariantContextWriter(this.outputFile);
			out.writeHeader(header);
			for(final CnvCall baseCall:intervals_list)
				{
				if(baseCall.echoed_flag) continue;
				
				final Allele ref = Allele.create("N",true); 
				final Allele alt = Allele.create("<"+baseCall.svType+">",false); 

				final List<Allele> alleles = Arrays.asList(ref,alt);
				
				
				final List<CnvCall> callsToPrint = this.vcfFilesInput.stream().
					flatMap(F->F.callMap.getOverlapping(baseCall).stream().
						filter(C->C.getType().equals(baseCall.getType())).
						filter(C->!C.echoed_flag).
						filter(C->testOverlapping(C, baseCall))
						).collect(Collectors.toList());
				
				if(callsToPrint.isEmpty()) {
					throw new IllegalStateException();
					}
				
				
				final VariantContextBuilder vcb= new VariantContextBuilder();
				vcb.chr(baseCall.getContig());
				vcb.start(baseCall.getStart());
				vcb.stop(baseCall.getEnd());
				vcb.attribute(VCFConstants.END_KEY, baseCall.getEnd());
				vcb.attribute(VCFConstants.SVTYPE, baseCall.getType());
				vcb.attribute("SVLEN", (1+baseCall.getEnd()-baseCall.getStart()));
				
				for(int side=0;side<2;side++)
					{
					final Function<CnvCall,Integer> coordExtractor;
					if(side==0)
						{
						coordExtractor = C->C.getStart();
						}
					else
						{
						coordExtractor = C->C.getEnd();
						}
					final List<Integer> list = Arrays.asList(
						callsToPrint.stream().
							mapToInt(C->coordExtractor.apply(C)-coordExtractor.apply(baseCall)).
							min().
							orElse(0),
						callsToPrint.stream().
							mapToInt(C->coordExtractor.apply(C)-coordExtractor.apply(baseCall)).
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
					
					//final List<Allele> alleles=Arrays.asList(ref,alt);
					
					final GenotypeBuilder gb = new GenotypeBuilder(call.sample,alleles);
					
					
										
					sample2gt.put(call.getSample(),gb.make());
					}
				vcb.attribute("SAMPLES",new ArrayList<>(sample2gt.keySet()));
				vcb.attribute("NSAMPLES",sample2gt.size());
				vcb.genotypes(sample2gt.values());
				vcb.alleles(alleles);
				out.add(vcb.make());				
				} //end of loop interval
			/*
			all_calls.values().stream().
				flatMap(C->C.stream()).
				filter(C->!C.echoed_flag).
				forEach(C->LOG.warn("Bug: Not printed "+C));
			 */			
			
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
