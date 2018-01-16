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
package com.github.lindenb.jvarkit.tools.lumpysv;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.math.stats.Percentile;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;


/**
BEGIN_DOC

## See also:



END_DOC
 */
@Program(name="lumpymerge",
description="Custom VCF merger for Lumpy.",
keywords={"lumpy","vcf","sort"},
generate_doc=false
)
public class LumpyMerge 
	 extends Launcher {
	
	private static final Logger LOG = Logger.build(LumpyMerge.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-f","--fraction"},description="Overlap fraction")
	private double fraction_overlap = 0.8;
	
	private class LumpyVar
		{
		private final VariantContext ctx;
		LumpyVar(final VariantContext ctx)
			{
			this.ctx = ctx;
			}
		LumpyVar next=null;
		
		private String getStrand() {
			String s = ctx.getAttributeAsString("STRANDS", ":");
			int colon = s.indexOf(":");
			return s.substring(0, colon);
		}
		
		boolean overlap(final LumpyVar other) {
			if(!this.ctx.getContig().equals(other.ctx.getContig())) return false;
			if(this.ctx.getStructuralVariantType()!=other.ctx.getStructuralVariantType()) return false;
			if(!this.getStrand().equals(other.getStrand())) return false;
			
			Interval L1 = LumpyConstants.getIntervalLeft(this.ctx);
			Interval L2 = LumpyConstants.getIntervalLeft(other.ctx);
			if(!LumpyMerge.this.overlap(L1,L2)) return false;
		
			if(this.ctx.getStructuralVariantType()==StructuralVariantType.BND) {
				L1 = LumpyConstants.getIntervalRight(this.ctx);
				L2 = LumpyConstants.getIntervalRight(other.ctx);
				if(!LumpyMerge.this.overlap(L1,L2)) return false;
				}
			
			
			
			if(this.next!=null)
				{
				return this.next.overlap(other);
				}
			return true;
			}
		void add(final LumpyVar L)
			{
			if(next==null) 
				{
				next=L;
				}
			else
				{
				next.add(L);
				}	
			}
		List<LumpyVar> asList() {
			final List<LumpyVar> L = new ArrayList<>();
			LumpyVar curr=this;
			while(curr!=null) { L.add(curr);curr=curr.next;}
			return L;
			}
		}
	
	private boolean overlap(final Interval i1,final Interval i2)
		{
		if(!i1.intersects(i2)) return false;
		int L1 = i1.length();
		int L2 = i2.length();
		int  L3 = i1.getIntersectionLength(i2);
		if(L3< (int)(this.fraction_overlap*L1)) return false;
		if(L3< (int)(this.fraction_overlap*L2)) return false;
		return true;
		}
	
	@Override
	public int doWork(final List<String> args) {
	final IntervalTreeMap<List<LumpyVar>> intervalTreeMap = new IntervalTreeMap<>();
	VariantContextWriter vcw = null;
	VCFFileReader vcfIn= null;
	final List<File> inputs = IOUtil.unrollFiles(
			args.stream().map(S->new File(S)).collect(Collectors.toList()),
			".vcf",".vcf.gz");
	if(inputs.isEmpty()) {
		LOG.error("empty vcf list");
		return -1;
		}
	try {
		final Set<VCFHeaderLine> metaData = new HashSet<>();
		//metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_KEY));

		final Set<String> sampleNames = new TreeSet<>();
		for(final File vcfFile : inputs)
			{
			int countVariant = 0;
			LOG.info("Read "+vcfFile);
			vcfIn  = new VCFFileReader(vcfFile,false);
			VCFHeader header = vcfIn.getFileHeader();
			if(!LumpyConstants.isLumpyHeader(header))
				{
				LOG.error("doesn't look like a Lumpy-SV vcf header "+vcfFile);
				return -1;
				}
			if(!header.hasGenotypingData()) {
				LOG.error("No sample in "+vcfFile);
				return -1;
				}
			for(final String sampleName : header.getSampleNamesInOrder())
				{
				if(sampleNames.contains(sampleName)) {
					LOG.error("Sample found twice "+sampleName+" in "+vcfFile);
					return -1;
					}
				sampleNames.add(sampleName);
				}
			metaData.addAll(
					header.getMetaDataInInputOrder().
					stream().filter(H->!H.getKey().equals("fileDate")).
						collect(Collectors.toSet())
					);

			final CloseableIterator<VariantContext> iter =vcfIn.iterator();
			while(iter.hasNext()) {
				final VariantContext ctx = new VariantContextBuilder(iter.next()).
					rmAttribute("PRPOS").
					rmAttribute("IMPRECISE").
					rmAttribute("SU").
					rmAttribute("PE").
					rmAttribute("SR").
					make();
				
				if(!LumpyConstants.isLumpyVariant(ctx))
					{
					LOG.error("doesn't look like a Lumpy-SV variant "+ctx);
					return -1;
					}
				++countVariant;
				final LumpyVar lCtx = new LumpyVar(ctx);
				final Interval rgn = LumpyConstants.getIntervalLeft(ctx);
				boolean found=false;
				for(final List<LumpyVar> L:intervalTreeMap.getOverlapping(rgn))
					{
					int x=0;
					while(!found && x<L.size())
						{
						final LumpyVar o = L.get(x);
						
						if(o.overlap(lCtx))
							{
							o.add(lCtx);
							found=true;
							break;
							}
						++x;
						}
					}
				if(!found)
					{
					List<LumpyVar> L = intervalTreeMap.get(rgn);
					if(L==null)
						{
						L=new ArrayList<>();
						intervalTreeMap.put(rgn,L);
						}
					L.add(lCtx);
					}
				}
			LOG.info(vcfFile+ " N:"+countVariant);
			iter.next();
			vcfIn.close();
			}
		final VCFHeader outHeader = new VCFHeader(metaData, sampleNames);
		
		final List<Allele> DIPLOID_NO_CALLS=Arrays.asList(Allele.NO_CALL,Allele.NO_CALL);
		vcw = super.openVariantContextWriter(this.outputFile);
		vcw.writeHeader(outHeader);
		for(final Interval rgn:intervalTreeMap.keySet() ) {
			final List<LumpyVar> chromvariants = intervalTreeMap.get(rgn);
			
			
			for(int i=0;i< chromvariants.size();i++)
				{
				List<LumpyVar> varl = chromvariants.get(i).asList();
				
				final LumpyVar first = varl.get(0);
				final Map<String,Genotype> sample2genotype = new HashMap<>(sampleNames.size());
				sampleNames.forEach(S->sample2genotype.put(S, new GenotypeBuilder(S,DIPLOID_NO_CALLS).
						attribute("SU",0).
						attribute("SR",0).
						attribute("PE",0).
						make()));
				final int variantStart = varl.stream().mapToInt(V->V.ctx.getStart()).min().getAsInt();
				//final int variantStart_x = varl.stream().mapToInt(V->V.ctx.getStart()).max().getAsInt();
				//final int variantEnd_x = varl.stream().mapToInt(V->V.ctx.getEnd()).min().getAsInt();
				final int variantEnd = varl.stream().mapToInt(V->V.ctx.getEnd()).max().getAsInt();
				final VariantContextBuilder vcb = new VariantContextBuilder(
						"lumpymerge",
						first.ctx.getContig(),
						variantStart,
						variantEnd,
						first.ctx.getAlleles()
						);
				vcb.attribute("END", variantEnd);
				vcb.attribute("SVTYPE", first.ctx.getAttribute("SVTYPE"));
				vcb.attribute("SVLEN", (int)Percentile.median().evaluate(varl.stream().mapToInt(V->V.ctx.getEnd()-V.ctx.getStart())));
				vcb.attribute("STRANDS",  first.getStrand());
				vcb.attribute("CIPOS",Arrays.asList(0,0));
				vcb.attribute("CIEND",Arrays.asList(0,0));

				
				
				
				varl.stream().flatMap(V->V.ctx.getGenotypes().stream()).forEach(G->{
					sample2genotype.put(G.getSampleName(), G);
				});
				
				vcb.genotypes(sample2genotype.values());
				vcw.add(vcb.make());
				}
		
			}
		vcw.close();vcw=null;
		
		return 0;
		}
	catch(final Exception err)
		{
		LOG.error(err);
		return -1;
		}
	finally
		{
		CloserUtil.close(vcfIn);
		CloserUtil.close(vcw);
		}
	}
	 
	public static void main(final String[] args) {
		new LumpyMerge().instanceMainWithExit(args);
	}
}
