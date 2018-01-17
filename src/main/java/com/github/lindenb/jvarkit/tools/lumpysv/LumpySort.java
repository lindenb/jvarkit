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

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.math.stats.Percentile;
import com.github.lindenb.jvarkit.util.iterator.ForwardPeekIterator;
import com.github.lindenb.jvarkit.util.iterator.ForwardPeekIteratorImpl;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.vcf.JexlVariantPredicate;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFEncoder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderVersion;


/**
BEGIN_DOC

## Example



END_DOC
 */
@Program(name="lumpysort",
description="sort and merge a set of Lumpy-SV VCF files.",
keywords={"lumpy","vcf","sort"},
generate_doc=false
)
public class LumpySort 
	 extends Launcher {
	
	private static final Logger LOG = Logger.build(LumpySort.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-f","--fraction"},description="Overlap fraction")
	private double fraction_overlap = 0.8;
	@Parameter(names={"-slop","--slop"},description="slop for 'x' bases in both directions.")
	private int slop_size = 0;
	@Parameter(names={"-secondary","--secondary"},description="keep SECONDARY record. \"Secondary breakend in a multi-line variants\".")
	private boolean keep_secondary = false;
	@Parameter(names={"-vf","--variant-filter"},description=JexlVariantPredicate.PARAMETER_DESCRIPTION,converter=JexlVariantPredicate.Converter.class)
	private Predicate<VariantContext> variantFilter = JexlVariantPredicate.create("");
	@Parameter(names={"-dm","--do-not-merge"},description="Do not merge genotypes, just sort")
	private boolean do_not_merge_ctx = false;
	@Parameter(names={"-gt","--genotype"},description="genotypes with SU>0 will be called with ALT (haploid)")
	private boolean do_genotype = false;
	@ParametersDelegate
	private Launcher.WritingSortingCollection sortingArgs = new Launcher.WritingSortingCollection();
	
	/** encoder for VariantCtx -> line */
	private VCFEncoder vcfEncoder = null;
	/** encoder for line -> VariantCtx */
	private final VCFCodec vcfCodec = new VCFCodec();
	
	/** VariantContext wrapper */
	private class LumpyVar
		{
		private final VariantContext ctx;
		boolean consummed = false;
		LumpyVar(final VariantContext ctx)
			{
			this.ctx = ctx;
			}
		private Set<String> getGenotypedSamples() {
			return ctx.getGenotypes().stream().
					filter(G->isAvailableGenotype(G)).
					map(G->G.getSampleName()).
					collect(Collectors.toSet());
			}	
		
		private Interval getInterval() {
			if(!ctx.hasAttribute("CIPOS")) throw new IllegalArgumentException("No CIPOS in "+ctx);
			final List<Integer> ciposL= ctx.getAttributeAsIntList("CIPOS",0);
			if(ciposL.size()!=2) throw new IllegalArgumentException("len(CIPOS)!=2 in "+ctx);
			if(!ctx.hasAttribute("CIEND")) throw new IllegalArgumentException("No CIEND in "+ctx);
			final List<Integer> ciendL= ctx.getAttributeAsIntList("CIEND",0);
			if(ciendL.size()!=2) throw new IllegalArgumentException("len(CIEND)!=2 in "+ctx);

			return new Interval(
				ctx.getContig(),
				Math.max(0,ctx.getStart() + ciposL.get(0) - LumpySort.this.slop_size),
				ctx.getEnd() + ciendL.get(1) + LumpySort.this.slop_size
				);
			}
		private Interval getBndInterval() {
			if(!ctx.hasAttribute("CIPOS")) throw new IllegalArgumentException("No CIPOS in "+ctx);
			final List<Integer> ciposL= ctx.getAttributeAsIntList("CIPOS",0);
			if(ciposL.size()!=2) throw new IllegalArgumentException("len(CIPOS)!=2 in "+ctx);
			if(!ctx.hasAttribute("CIEND")) throw new IllegalArgumentException("No CIEND in "+ctx);
			final List<Integer> ciendL= ctx.getAttributeAsIntList("CIEND",0);
			if(ciendL.size()!=2) throw new IllegalArgumentException("len(CIEND)!=2 in "+ctx);
			
			String cL;
			int pL;
			if(ctx.getStructuralVariantType()==StructuralVariantType.BND) {
				final  Map.Entry<String,Integer> entry = LumpyConstants.getBnDContigAndPos(ctx.getAlternateAllele(0).getDisplayString());
				cL = entry.getKey();
				pL = entry.getValue();
				} 
			else
				{
				cL = ctx.getContig();
				pL = ctx.getEnd();
				}
			
			return new Interval(
				cL,
				Math.max(0,pL+ciposL.get(0) - LumpySort.this.slop_size),
				pL+ciendL.get(1) + LumpySort.this.slop_size
				);	
			}
		
		boolean canMerge(final LumpyVar o)
			{
			if(this.consummed) return false;
			if(o.consummed) return false;
			// we cannot have common available variants between two ctx
			final Set<String> commonSamples= new HashSet<String>(this.getGenotypedSamples());
			commonSamples.retainAll(o.getGenotypedSamples());
			
			
			if(!commonSamples.isEmpty()) {
				return false;
			}
			if(this.compare1(o)!=0) return false;
			
			Interval L1 = this.getInterval();
			Interval L2 = o.getInterval();
			if(!LumpySort.this.overlap(L1,L2)) return false;
			if(this.ctx.getStructuralVariantType()==StructuralVariantType.BND) {
				L1 = this.getBndInterval();
				L2 = o.getBndInterval();
				if(!LumpySort.this.overlap(L1,L2)) return false;
				}
			
			return true;
			}
		
		private String getStrands()
			{
			String s= ctx.getAttributeAsString("STRANDS","");
			int colon = s.indexOf(':');
			return colon==-1?s:s.substring(0,2);
			}
			
		public int compare1(final LumpyVar o) {
			final StructuralVariantType st1= this.ctx.getStructuralVariantType();
			final StructuralVariantType st2= o.ctx.getStructuralVariantType();
			int i= st1.compareTo(st2);
			if(i!=0) return i;
			//
			String s1 = this.ctx.getContig();
			String s2 = o.ctx.getContig();
			i = s1.compareTo(s2);
			if(i!=0) return i;
			//
			if(st1.equals(StructuralVariantType.BND)) {
				s1 = LumpyConstants.getBnDContig(this.ctx.getAlternateAllele(0).getDisplayString());
				s2 = LumpyConstants.getBnDContig(o.ctx.getAlternateAllele(0).getDisplayString());
				i = s1.compareTo(s2);
				if(i!=0) return i;
				}
			//
			s1 = this.getStrands();
			s2 = o.getStrands();
			i = s1.compareTo(s2);
			if(i!=0) return i;
			return 0;
			}
		
		public int compare2(final LumpyVar o) {
			int i = compare1(o);
			if(i!=0) return i;
			i = this.ctx.getStart() - o.ctx.getStart();
			if(i!=0) return i;
			i = this.ctx.getEnd() - o.ctx.getEnd();
			if(i!=0) return i;
			
			String S1 = variantContextToLine(this.ctx);
			String S2 = variantContextToLine(o.ctx);
			return S1.compareTo(S2);
			}
		
		}
	
	/** codec for Sorting collection */
	private class LumpyVarCodec
		extends AbstractDataCodec<LumpyVar>
		{
		@Override
		public LumpyVar decode(final DataInputStream dis) throws IOException {
			String line;
			try
				{
				line = readString(dis);
				}
			catch(final Exception err)
				{
				return null;
				}
			return new LumpyVar(linetoVariantContext(line));
			}
		@Override
		public void encode(final DataOutputStream dos, final  LumpyVar v) throws IOException {
			writeString(dos,variantContextToLine(v.ctx));
			}
		@Override
		public LumpyVarCodec clone() {
			return new LumpyVarCodec();
			}
		}
	
	/** variant decoder */
	private VariantContext linetoVariantContext(final String line) {
		return this.vcfCodec.decode(line);
	}
	
	/** encoder variant decoder */
	private  String variantContextToLine(final VariantContext ctx) {
		return this.vcfEncoder.encode(ctx);
	}
	
	/** returns true two interval overlap with fraction_overlap  */
	private boolean overlap(final Interval i1,final Interval i2)
		{
		if(!i1.intersects(i2)) return false;
		final int L1 = i1.length();
		final int L2 = i2.length();
		final int  L3 = i1.getIntersectionLength(i2);
		if(L3< (int)(this.fraction_overlap*L1)) return false;
		if(L3< (int)(this.fraction_overlap*L2)) return false;
		return true;
		}
	
	/** returns true if there is a SU greater than 0 */
	private boolean isAvailableGenotype(final Genotype g)
		{
		if(!g.hasExtendedAttribute("SU")) {
			return false;
		}
		final Object su  = g.getExtendedAttribute("SU", 0);
		
		if(su==null ) return false;
		
		int suv=(su instanceof Integer ?
				Integer.class.cast(su).intValue():
				Integer.parseInt(su.toString())
				);
		if(suv<=0)
			{
			return false;
			}
		return true;
		}
	
	@Override
	public int doWork(final List<String> args) {
	VariantContextWriter vcw = null;
	LineIterator vcfIn= null;
	SortingCollection<LumpyVar> sorting = null;
	final List<File> inputs = IOUtil.unrollFiles(
			args.stream().map(S->new File(S)).collect(Collectors.toList()),
			".vcf",".vcf.gz");
	if(inputs.isEmpty()) {
		LOG.error("empty vcf list");
		return -1;
		}
	try {
		final Set<VCFHeaderLine> metaData = new HashSet<>();
		final Set<String> sampleNames = new TreeSet<>();

		
		for(final File vcfFile : inputs)
			{
			LOG.info("Reading Header of "+vcfFile);
			final VCFFileReader r  = new VCFFileReader(vcfFile,false);
			final VCFHeader header = r.getFileHeader();
			if(!LumpyConstants.isLumpyHeader(header))
				{
				LOG.error("doesn't look like a Lumpy-SV vcf header "+vcfFile);
				r.close();
				return -1;
				}
			
			if(!header.hasGenotypingData()) {
				LOG.error("No sample in "+vcfFile);
				r.close();
				return -1;
				}
			for(final String sampleName : header.getSampleNamesInOrder())
				{
				if(sampleNames.contains(sampleName)) {
					LOG.error("Sample found twice "+sampleName+" in "+vcfFile);
					r.close();
					return -1;
					}
				sampleNames.add(sampleName);
				}
			metaData.addAll(
					header.getMetaDataInInputOrder().
					stream().filter(H->!H.getKey().equals("fileDate")).
						collect(Collectors.toSet())
					);
			r.close();
			}
		final VCFHeader outHeader = new VCFHeader(
			metaData,
			sampleNames
			);
		final VCFHeaderVersion versions[]=VCFHeaderVersion.values();
		this.vcfEncoder = new VCFEncoder(outHeader, false, true);
		this.vcfCodec.setVCFHeader(
				outHeader,
				versions[versions.length-1]
				);
		sorting = SortingCollection.newInstance(
				LumpyVar.class,
				new LumpyVarCodec(),
				(V1,V2)->V1.compare2(V2),
				this.sortingArgs.getMaxRecordsInRam(),
				this.sortingArgs.getTmpPaths()
				);
		sorting.setDestructiveIteration(true);
		long total_variants = 0L;

		
		for(int idx=0;idx< inputs.size();++idx)
			{
			final File vcfFile = inputs.get(idx);
			LOG.info("Read "+(idx+1)+"/"+inputs.size()+" Variants of "+vcfFile);
			int nVariant = 0;
			final VCFFileReader r  = new VCFFileReader(vcfFile,false);
			
			final List<Genotype> missing =new ArrayList<>(sampleNames.size());
			for(final String sn:sampleNames)
				{
				if(r.getFileHeader().getSampleNamesInOrder().contains(sn)) continue;
				missing.add(GenotypeBuilder.createMissing(sn, 2));
				}
			
			final CloseableIterator<VariantContext> iter = r.iterator();
			while(iter.hasNext()) {
				VariantContext ctx = iter.next();
				if(!this.keep_secondary) {
					if(ctx.hasAttribute("SECONDARY")) continue;
					}
				if(!this.variantFilter.test(ctx)) continue;
				final List<Genotype> gtList  = new ArrayList<>(ctx.getGenotypes());
				if(this.do_genotype)
					{
					for(int gi=0;gi< gtList.size();gi++)
						{
						Genotype g= gtList.get(gi);
						if(!isAvailableGenotype(g)) continue;
						final GenotypeBuilder gb = new GenotypeBuilder(g.getSampleName(), ctx.getAlternateAlleles());
						gb.attributes(g.getExtendedAttributes());
						gtList.set(gi, gb.make());
						}
					}
				
				
				gtList.addAll(missing);
				
				ctx = new VariantContextBuilder(ctx).
						genotypes(gtList).
						rmAttribute("PRPOS").
						make();
				
				
				
				sorting.add(new LumpyVar(ctx));
				nVariant++;
				total_variants++;
				}
			iter.close();
			r.close();
			LOG.info("Read Variants of "+vcfFile+" N="+nVariant+" Total:"+total_variants);
			System.gc();
			}
			
		sorting.doneAdding();
		
		LOG.info("Writing output");
		final List<Allele> ALLELES_NO_CALLS=
				this.do_genotype
				? Collections.singletonList(Allele.NO_CALL)
				: Arrays.asList(Allele.NO_CALL,Allele.NO_CALL)
				;
		final CloseableIterator<LumpyVar> iter = sorting.iterator();
		final ForwardPeekIterator<LumpyVar> fwdIter =  new ForwardPeekIteratorImpl<LumpySort.LumpyVar>(iter);

		vcw = super.openVariantContextWriter(this.outputFile);
		vcw.writeHeader(outHeader);
		for(;;)
			{
			if(!fwdIter.hasNext()) break;
			final LumpyVar first = fwdIter.next();
			if(this.do_not_merge_ctx)
				{
				vcw.add(first.ctx);
				continue;
				}
			
			if(first.consummed) {
				continue;
			}
			final List<LumpyVar> buffer = new ArrayList<>();
			buffer.add(first);
			for(int x=0;;++x)
				{
				final LumpyVar lv = fwdIter.peek(x);
				if(lv==null){
					break;
				}
				if(first.compare1(lv)!=0) 
					{
					break;
					}
				if(lv.consummed) continue;
				if(lv.ctx.getStart()>first.ctx.getEnd()) {
					break;
					}
				if(first.canMerge(lv))
					{
					buffer.add(lv);
					}
				else
					{
					}
				}
			
			final int variantStart = buffer.stream().
					mapToInt(V->V.ctx.getStart()).
					min().getAsInt();
			final int variantEnd = buffer.stream().
					mapToInt(V->V.ctx.getEnd()).
					max().getAsInt();
			final VariantContextBuilder vcb = new VariantContextBuilder(
					"lumpymerge",
					first.ctx.getContig(),
					variantStart,
					variantEnd,
					first.ctx.getAlleles()
					);
			vcb.attribute("END", variantEnd);
			vcb.attribute("SVTYPE", first.ctx.getAttribute("SVTYPE"));
			vcb.attribute("SVLEN", (int)Percentile.median().evaluate(buffer.stream().mapToInt(V->V.ctx.getEnd()-V.ctx.getStart())));
			vcb.attribute("STRANDS",  first.getStrands());
			vcb.attribute("CIPOS",Arrays.asList(0,0));
			vcb.attribute("CIEND",Arrays.asList(0,0));

			final Map<String,Genotype> sample2genotype = new HashMap<>(sampleNames.size());
			
			
			buffer.stream().flatMap(V->V.ctx.getGenotypes().stream()).
				filter(G->isAvailableGenotype(G)).
				forEach(G->{
				sample2genotype.put(G.getSampleName(), G);
			});
			
			for(final String sn: sampleNames)
				{
				if(!sample2genotype.containsKey(sn))
					{
					sample2genotype.put(sn, new GenotypeBuilder(sn,ALLELES_NO_CALLS).
							attribute("SU",0).
							attribute("SR",0).
							attribute("PE",0).
							make());
					}	
				}
			
			vcb.genotypes(sample2genotype.values());
			vcw.add(vcb.make());

			buffer.stream().forEach(V->{V.consummed=true;});
			}
		vcw.close();vcw=null;
		fwdIter.close();
		iter.close();
		
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
		new LumpySort().instanceMainWithExit(args);
	}
}
