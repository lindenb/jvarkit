package com.github.lindenb.jvarkit.tools.vcfcomposite;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParserFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParserFactory;

@Program(name="vcfcomposite",description="TODO",keywords={"vcf","disease","annotation","pedigree"})
public class VCFComposite extends Launcher {
	private static final Logger LOG= Logger.build(VCFComposite.class).make();
	@Parameter(names={"-ped","--pedigree"},description="Pedigree file",required=true)
	private File pedigreeFile=null;
	@Parameter(names={"-o","--out"},description="OUtput file")
	private File outputFile=null;
	@Parameter(names={"-acceptFiltered","--acceptFiltered"},description="Accept FILTERED variants")
	private boolean acceptFiltered=false;
	@Parameter(names={"-acceptID","--acceptID"},description="Accept variants with ID")
	private boolean acceptID=false;
	@Parameter(names={"-m","--model"},description="Model type",required=true)
	private Type modelType=null;
	@Parameter(names={"-models"},description="List the available models and exits",help=true)
	private boolean listModels=false;

	
	
	private Pedigree pedigree=null;
	
	
	private static class GeneAndVariant
		{
		String gene;
		String ctxLine;
		}	
	private static class GeneAndVariantCodec
		extends AbstractDataCodec<GeneAndVariant>
		{
		
		public GeneAndVariantCodec() {
			}
		@Override
		public GeneAndVariant decode(DataInputStream dis) throws IOException {
			final GeneAndVariant gav=new GeneAndVariant();
			try 
				{
				gav.gene=dis.readUTF();
				}
			catch(IOException err)
				{
				return null;
				}
			gav.ctxLine = super.readString(dis);
			return gav;
			}
		@Override
		public void encode(DataOutputStream dos, GeneAndVariant o) throws IOException {
			dos.writeUTF(o.gene);
			super.writeString(dos, o.ctxLine);
			}
		@Override
		public AbstractDataCodec<GeneAndVariant> clone() {
			return new GeneAndVariantCodec();
			}
		}	
	
	private final Predicate<Genotype> deNovoGenotype = new Predicate<Genotype>() {
		@Override
		public boolean test(Genotype t) {
			Pedigree.Person c = pedigree.getPersonById(t.getSampleName());
			if(c==null) return false;
			
			return false;
		}
	};
	
	
	public static enum Type {
		RecessiveHomVar,
		RecessiveComposite
	};
	
	private abstract class DiseaseModel
		{
		public List<Genotype> getGenotypesOfAffected(final VariantContext ctx)
			{
			return ctx.getGenotypes().stream().filter(G->{
				final Pedigree.Person p= pedigree.getPersonById(G.getSampleName());
				if(p==null) return false;
				if(!p.isAffected()) return false;
				return true;
				}).collect(Collectors.toList());
			}
		public List<Genotype> getGenotypesOfUnaffected(final VariantContext ctx)
			{
			return ctx.getGenotypes().stream().filter(G->{
				final Pedigree.Person p= pedigree.getPersonById(G.getSampleName());
				if(p==null) return false;
				if(!p.isUnaffected()) return false;
				return true;
				}).collect(Collectors.toList());
			}
		public abstract boolean accept(final VariantContext ctx);
		public abstract void scan(final String geneName,List<VariantContext> variants);
		}
	
	private class RecessiveHomVar extends DiseaseModel
		{
		@Override
		public boolean accept(final VariantContext ctx) {
			if(!getGenotypesOfAffected(ctx).stream().anyMatch(G->G.isHomVar())) return false;
			if(getGenotypesOfUnaffected(ctx).stream().anyMatch(G->G.isHomVar())) return false;
			return true;
			}
		public void scan(final String geneName,List<VariantContext> variants)
			{
			LOG.info(geneName);
			}
		}
	private class RecessiveComposite extends DiseaseModel
		{
		protected boolean isGenotypeForAffected(Genotype g)
			{
			if(g==null) return false;
			return g.isHet() || g.isHomVar();
			}
		
		@Override
		public boolean accept(final VariantContext ctx) {
			if(!getGenotypesOfAffected(ctx).stream().
					anyMatch(G->isGenotypeForAffected(G))) return false;			
			return true;
			}
		public void scan(final String geneName,List<VariantContext> variants)
			{

			if(variants.size()<2) return;
			for(Pedigree.Person c: pedigree.getAffected()) {
				for(int x=0;x+1< variants.size();++x)
					{
					final Genotype gcx = variants.get(x).getGenotype(c.getId());
					if(gcx==null || !isGenotypeForAffected(gcx)) continue;
					for(int y=x+1;y< variants.size();++y)
						{
						final Genotype gcy = variants.get(y).getGenotype(c.getId());
						if(gcy==null || !isGenotypeForAffected(gcy)) continue;
						boolean unaffected_are_ok=true;
						//check unaffected don't have same haplotype
						for(final Pedigree.Person unaffected: pedigree.getUnaffected()) {
							final Genotype gux = variants.get(x).getGenotype(unaffected.getId());
							final Genotype guy = variants.get(y).getGenotype(unaffected.getId());
							if(gux!=null && guy!=null &&
								gux.sameGenotype(gcx, true) &&
								guy.sameGenotype(gcy, true)
								)
								{
								LOG.debug("same in unaffected");
								unaffected_are_ok=false;
								break;
								}
							}
						if(unaffected_are_ok)
							{
							LOG.info(">>>>>>>>>>><<"+geneName);
							}
						}
					}
				}
			}
		}

	
	protected DiseaseModel createModel() {
		if(this.modelType==null) throw new NullPointerException();
		switch(this.modelType) {
			case RecessiveHomVar: return new RecessiveHomVar();
			case RecessiveComposite: return new RecessiveComposite();
			default: throw new IllegalStateException();
			}
		
		}
	@Override
	public int doWork(final List<String> args)
		{
		if(listModels)
			{
			for(final Type t:Type.values()) {
				stdout().print(t.name());
				stdout().println();
				}
			return 0;
			}

		try {
			this.pedigree = Pedigree.readPedigree(pedigreeFile);
			
			if(this.pedigree.getAffected().isEmpty()) {
				LOG.error("No Affected sample in "+this.pedigreeFile);
				return -1;
			}
			if(this.pedigree.getUnaffected().isEmpty()) {
				LOG.error("No Unaffected sample in "+this.pedigreeFile);
				return -1;
			}
			final DiseaseModel model = this.createModel();
			final String inputName= super.oneFileOrNull(args);
			final LineIterator r= (inputName==null?
					IOUtils.openStreamForLineIterator(stdin()):
					IOUtils.openURIForLineIterator(inputName)
					);
			final VCFCodec codec = new VCFCodec();
			final VCFHeader header = (VCFHeader)codec.readActualHeader(r);
			final AnnPredictionParser annParser=new AnnPredictionParserFactory(header).get();
			final VepPredictionParser vepParser=new VepPredictionParserFactory(header).get();
			
			
			//final VCFHeader h2=new VCFHeader(header.getMetaDataInInputOrder(),header.getSampleNamesInOrder());
			//h2.addMetaDataLine(new VCFInfoHeaderLine(this.TAG,1,VCFHeaderLineType.String,"Values from bigwig file: "+BIGWIG));
			
			
			SortingCollection<GeneAndVariant> sorting=null;
			String prevContig=null;	
			
			for(;;)
				{
				String line;
				final VariantContext ctx;
				if(r.hasNext())
					{
					line = r.next();
					ctx=codec.decode(line);
					}
				else
					{
					line=null;
					ctx=null;
					}
				
				if(ctx==null || !ctx.getContig().equals(prevContig))
					{
					if(sorting!=null)
						{
						LOG.debug("Dump contig "+prevContig);
						sorting.doneAdding();
						CloseableIterator<GeneAndVariant> iter2=sorting.iterator();
						EqualRangeIterator<GeneAndVariant> eqiter= new EqualRangeIterator<>(iter2,(A,B)->A.gene.compareTo(B.gene));
						while(eqiter.hasNext())
							{
							final List<GeneAndVariant> variants=eqiter.next();
							model.scan(
								variants.get(0).gene,
								variants.stream().
									map(L->codec.decode(L.ctxLine)).
									collect(Collectors.toList())
								);
							}
						eqiter.close();
						iter2.close();
						sorting.cleanup();
						}
					sorting=null;
					if(ctx==null) break;
					prevContig=ctx.getContig();
					
					}
				if(!ctx.isVariant()) continue;
				if(!acceptFiltered && ctx.isFiltered()) continue;
				if(!acceptID && ctx.hasID()) continue;
				if(!model.accept(ctx)) continue;
				final Set<String> geneKeys = new HashSet<>();
				for(final AnnPredictionParser.AnnPrediction pred: annParser.getPredictions(ctx)) {
					geneKeys.addAll(pred.getGeneKeys().stream().map(S->ctx.getContig()+"_"+S).collect(Collectors.toSet()));

					}
				for(final VepPredictionParser.VepPrediction pred: vepParser.getPredictions(ctx)) {
					geneKeys.addAll(pred.getGeneKeys().stream().map(S->ctx.getContig()+"_"+S).collect(Collectors.toSet()));
					}
				
				if(sorting==null) {
					sorting = SortingCollection.newInstance(GeneAndVariant.class,
							new GeneAndVariantCodec(),
							(A,B)->{int i=A.gene.compareTo(B.gene);if(i!=0) return i;return A.ctxLine.compareTo(B.ctxLine);},
							50000
							);
					sorting.setDestructiveIteration(true);
					}
				
				for(final String gk:geneKeys) 
					{
					final GeneAndVariant gav=new GeneAndVariant();
					gav.gene = gk;
					gav.ctxLine = line;
					sorting.add(gav);
					}
				
				}
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			
			}
		}
	public static void main(String[] args) {
		new VCFComposite().instanceMainWithExit(args);
	}
	
	}
