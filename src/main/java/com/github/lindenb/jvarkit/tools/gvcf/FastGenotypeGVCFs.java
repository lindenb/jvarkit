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
package com.github.lindenb.jvarkit.tools.gvcf;

import java.io.Closeable;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.SmartComparator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;
import com.github.lindenb.jvarkit.util.vcf.ContigPosRef;
import com.github.lindenb.jvarkit.util.vcf.VariantAttributesRecalculator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;


@Program(name="fastgenotypegvcfs",
	description="Fast Genotype Gvcfs",
	generate_doc=false
	)
public class FastGenotypeGVCFs extends Launcher {
	
	private static final Logger LOG = Logger.build(FastGenotypeGVCFs.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	
	
	private class GVCFVariantIterator
		implements Closeable
		{
		private final File gvcfFile;
		private final VCFFileReader vcfFileReader;
		private final CloseableIterator<VariantContext> iter;
		private final List<VariantContext> buffer = new ArrayList<>();
		private final List<String> samples;
		GVCFVariantIterator(final File vcf) {
			this.gvcfFile = vcf;
			this.vcfFileReader = new VCFFileReader(vcf,false);
			this.iter = this.vcfFileReader.iterator();
			this.samples = this.vcfFileReader.getFileHeader().getSampleNamesInOrder();
			}
		
		String getSource() {
			return this.gvcfFile.getPath();
		}
		
		@Override
		public void close() {
			CloserUtil.close(this.iter);
			CloserUtil.close(this.vcfFileReader);
			}
		
		private VariantContext cleanup(final VariantContext ctx) {
			return ctx;
			}
		
		private boolean isVariant(final VariantContext ctx) {
			if(ctx.getNAlleles()==2 && ctx.getAlleles().get(1).equals(Allele.NON_REF_ALLELE)) return false;
			return true;
			}
		
		ContigPosRef lookup() {
			for(int i=0;i< this.buffer.size();++i)
				{
				final VariantContext vc = this.buffer.get(i);
				if(!isVariant(vc)) continue;
				return new ContigPosRef(vc);
				}
			while(this.iter.hasNext())
				{
				final VariantContext vc =  this.iter.next();
				this.buffer.add(cleanup(vc));
				if(isVariant(vc)) return new ContigPosRef(vc);
				}
			return null;
			}
		
		private VariantContext makeNoCall(final ContigPosRef lookedUp) {
			final VariantContextBuilder vcb = new VariantContextBuilder(
					getSource(),
					lookedUp.getContig(), 
					lookedUp.getStart(),
					lookedUp.getEnd(), 
					Collections.singletonList(lookedUp.getReference())
					);
			vcb.genotypes(this.samples.stream().
				map(S->GenotypeBuilder.createMissing(S, 2)).
				collect(Collectors.toList()));
			return vcb.make();
			}
		
		
		VariantContext next(final ContigPosRef lookedUp) {
			int i=0;
			while(i< buffer.size())
				{
				final VariantContext vc = this.buffer.get(i);
				int diff = FastGenotypeGVCFs.this.contigComparator.compare(vc.getContig(), lookedUp.getContig());
				if(diff < 0) {
					//LOG.debug("remove "+vc+" for "+lookedUp);
					this.buffer.remove(i);
					continue;
					}
				else if(diff>0)
					{
					//LOG.debug("return nocall for "+lookedUp);
					return makeNoCall(lookedUp);
					}
				
				if(vc.getEnd()<lookedUp.getStart())
					{
					//LOG.debug("remove2 "+vc+" for "+lookedUp);
					this.buffer.remove(i);
					continue;
					}
				
				if(!isVariant(vc)) {
					if(vc.getStart()>lookedUp.getEnd())
						{
						//LOG.debug("return2 nocall for "+lookedUp);
						return makeNoCall(lookedUp);
						}
					if(vc.getStart()<=lookedUp.getPos() &&
						lookedUp.getPos() <= vc.getEnd())
						{
						//LOG.debug("return overlapping:\t "+vc+"\n\t(NOT_VARIANT) for "+lookedUp);
						final VariantContextBuilder vcb = new VariantContextBuilder(
								getSource(),
								lookedUp.getContig(), 
								lookedUp.getStart(),
								lookedUp.getEnd(), 
								Collections.singletonList(lookedUp.getReference())
								);
						final List<Genotype> genotypes= new ArrayList<>(this.samples.size());
						for(final Genotype gt:vc.getGenotypes())
							{
							
							final GenotypeBuilder gb = new GenotypeBuilder(gt.getSampleName(),
									gt.getAlleles().stream().map(A->A.isReference()?lookedUp.getReference():A).collect(Collectors.toList())
									);
							if(gt.hasAD()) gb.AD(gt.getAD());
							if(gt.hasDP()) gb.DP(gt.getDP());
							if(gt.hasGQ()) gb.GQ(gt.getGQ());
							if(gt.hasPL()) gb.PL(gt.getPL());
							genotypes.add(gb.make());
							}
						
						vcb.genotypes(genotypes);
						return vcb.make();
						}
					else
						{
						//LOG.warn("unknown case :"+lookedUp+" "+vc);
						return makeNoCall(lookedUp);
						}
					}
				else
					{					
					if(vc.getStart()>lookedUp.getEnd())
						{
						//LOG.debug("return3 nocall "+vc+" for "+lookedUp);
						return makeNoCall(lookedUp);
						}
					if(!vc.getReference().equals(lookedUp.getReference()))
						{
						//LOG.debug("skip "+vc+" for "+lookedUp);
						i++;
						continue;
						}
					if(vc.getStart()!=lookedUp.getPos())
						{	
						//LOG.warn("boom :"+lookedUp+" "+vc);
						return makeNoCall(lookedUp);
						}
					//https://gatkforums.broadinstitute.org/gatk/discussion/4216/non-ref-in-gvcf
					/* if(vc.getGenotypes().stream().anyMatch(G->G.getAlleles().contains(Allele.NON_REF_ALLELE))) {
						throw new RuntimeException("Boum "+vc+" "+lookedUp+" "+getSource());
						}*/
					
					final VariantContextBuilder vcb = new VariantContextBuilder(vc);
					vcb.alleles(vc.getAlleles().
							stream().
							filter(A->!A.equals(Allele.NON_REF_ALLELE)).
							collect(Collectors.toList())
							);
					vcb.genotypes(
							vc.getGenotypes().stream().
							map(G->G.getAlleles().contains(Allele.NON_REF_ALLELE)?
									GenotypeBuilder.createMissing(G.getSampleName(),2):
									G).
							collect(Collectors.toList())
							)
							;
					this.buffer.remove(i);
					return vcb.make();
					}
				}
			return makeNoCall(lookedUp);
			}
		}
	
	
	private SAMSequenceDictionary dictionary =null;
	
	private Comparator<String> contigComparator = null;
	
	private final Comparator<ContigPosRef> contigPosRefComparator = (S1,S2) ->{
		int i= 	contigComparator.compare(S1.getContig(), S2.getContig());
		if(i!=0) return i;
		i = Integer.compare(S1.getPos(), S2.getPos());
		if(i!=0) return i;
		return S1.getReference().compareTo(S2.getReference());
		};

		
	
	@Override
	public int doWork(final List<String> args) {
		VariantContextWriter w=null;
		
		try {
			final List<GVCFVariantIterator> gvcfSources =  
					IOUtil.unrollFiles(args.stream().map(F->new File(F)).collect(Collectors.toSet()),".g.vcf",".g.vcf.gz" ).
					stream().
					map(F->new GVCFVariantIterator(F)).
					collect(Collectors.toList())
					;
			if(args.isEmpty())
				{
				LOG.error("No gvcf file was given");
				return -1;
				}
			this.dictionary  = gvcfSources.get(0).vcfFileReader.getFileHeader().getSequenceDictionary();
			if(this.dictionary==null)
				{
				LOG.error("Dict missing in "+gvcfSources.get(0).gvcfFile);
				return -1;
				}
			this.contigComparator = new ContigDictComparator(this.dictionary);
			
			gvcfSources.stream().map(S->S.vcfFileReader.getFileHeader().getSequenceDictionary()).forEach(D->{
				if(D==null || !SequenceUtil.areSequenceDictionariesEqual(D, dictionary))
					{
					throw new JvarkitException.UserError("dict missing or dict are not the same");
					}
				});
			
			
			if(	gvcfSources.stream().
					flatMap(S->S.samples.stream()).
					collect(Collectors.groupingBy(Function.identity(),Collectors.counting())).
					entrySet().stream().anyMatch(P->P.getValue()!=1L))
				{
				LOG.error("Duplicate sample name. check input");
				return -1;
				}
			final Set<VCFHeaderLine> metaData=new HashSet<>();
			
			VCFStandardHeaderLines.addStandardFormatLines(metaData, true,
					VCFConstants.GENOTYPE_KEY,
					VCFConstants.GENOTYPE_ALLELE_DEPTHS,
					VCFConstants.DEPTH_KEY,
					VCFConstants.GENOTYPE_QUALITY_KEY,
					VCFConstants.GENOTYPE_PL_KEY
					);
			metaData.addAll(gvcfSources.stream().flatMap(S->S.vcfFileReader.getFileHeader().getFormatHeaderLines().stream()).collect(Collectors.toSet()));
			
			final VCFHeader header= new VCFHeader(
					metaData, 
					gvcfSources.stream().flatMap(S->S.samples.stream()).
					sorted(new SmartComparator()).
					collect(Collectors.toList())
					);
			final VariantAttributesRecalculator attCalc = new VariantAttributesRecalculator();
			attCalc.setHeader(header);
			
			w= super.openVariantContextWriter(outputFile);
			w.writeHeader(header);
			
			for(;;)
				{
				String id = null;
				ContigPosRef next = null;
				for(GVCFVariantIterator it:gvcfSources)
					{
					ContigPosRef cpr = it.lookup();
					if(cpr==null) continue;
					if(next==null || contigPosRefComparator.compare(cpr, next)<0)
						{
						next = cpr;
						}
					}
				if(next==null) break;
				final Set<Allele> alleles = new HashSet<>();
				final List<Genotype> genotypes = new ArrayList<>();
				alleles.add(next.getReference());
				for(final GVCFVariantIterator it:gvcfSources)
					{
					final VariantContext vc = it.next(next);
					if(vc.hasID()) id=vc.getID();
					Objects.requireNonNull(vc, "vc is null");
					alleles.addAll(
							vc.getGenotypes().
								stream().
								flatMap(G->G.getAlleles().stream()).
								filter(A->A.isCalled()).
								collect(Collectors.toSet())
								);
					genotypes.addAll(vc.getGenotypes());
					}
				if(alleles.size()<2) continue;
				
				final VariantContextBuilder vcb = new VariantContextBuilder(
						null,
						next.getContig(), 
						next.getStart(),
						next.getEnd(), 
						alleles
						);
				if(id!=null) vcb.id(id);
				vcb.genotypes(genotypes);
				w.add(attCalc.apply(vcb.make()));
				}
			
			
			
			for(final GVCFVariantIterator src:gvcfSources) src.close();
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}	
		finally
			{
			
			CloserUtil.close(w);
			}
		}
	public static void main(String[] args) {
		new FastGenotypeGVCFs().instanceMainWithExit(args);

	}

}
