package com.github.lindenb.jvarkit.tools.gvcf;

import java.io.File;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFStandardHeaderLines;


@Program(name="fastgenotypegvcfs",description="Fast Genotype Gvcfs",generate_doc=false)
public class FastGenotypeGVCFs extends Launcher {
	
	private static final Logger LOG = Logger.build(FastGenotypeGVCFs.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	
	
	private final List<Source> gvcfSources = new ArrayList<>();
	private final static Allele NON_REF = Allele.create("<NON_REF>", false);
	private SAMSequenceDictionary dictionary =null;
	
	private final Function<String,Integer> contig2tid = S -> {
		final int tid1 = dictionary.getSequenceIndex(S);
		if(tid1==-1) throw new IllegalStateException();
		return tid1;
		};
	private final Comparator<String> contigComparator = (S1,S2) ->{
		return contig2tid.apply(S1) - contig2tid.apply(S2);
		};
	
	private final Comparator<VariantContext> variantContigPosRefComparator = (S1,S2) ->{
		int i= 	contigComparator.compare(S1.getContig(), S2.getContig());
		if(i!=0) return i;
		i = S1.getStart() - S2.getStart();
		if(i!=0) return i;
		return S1.getReference().compareTo(S2.getReference());
		};

		
	private abstract class GvcfThing implements Locatable
		{
		final VariantContext ctx;
		GvcfThing(final VariantContext ctx) {
			this.ctx = ctx;
			}
		public int getTid() {
			return contig2tid.apply(this.getContig());
			}
		@Override
		public String getContig() {
			return ctx.getContig();
			}
		@Override
		public int getStart() {
			return ctx.getStart();
			}
		public abstract boolean isVariant();
		public final boolean isBlock() { return !isVariant();}
		@Override
		public String toString() {
			return ctx.getContig()+":"+ctx.getStart()+"-"+ctx.getEnd()+":"+ctx.getAlleles().stream().map(A->A.toString()).collect(Collectors.joining("/"));
			}
		}
	private  class GvcfVariant extends GvcfThing
		{
		GvcfVariant(final VariantContext ctx) {
			super(ctx);
			}
		@Override
		public boolean isVariant() {
			return true;
			}
		@Override
		public int getEnd() {
			return ctx.getEnd();
			}
		}
	private  class GvcfBlock extends GvcfThing
		{
		GvcfBlock(final VariantContext ctx) {
			super(ctx);
			}
		
		@Override
		public int getEnd() {
			if(!ctx.hasAttribute("END")) throw new IllegalStateException("No \"END=\" in "+ctx);
			int end= ctx.getAttributeAsInt("END", -1);
			if(end<this.getStart())  throw new IllegalStateException("\"END=\""+end+"<=start in "+ctx);
			return end;
			}
		@Override
		public boolean isVariant() {
			return false;
			}
		}
	private  class Source
		{
		final File gvcfFile;
		private VCFFileReader vcfFileReader=null;
		private CloseableIterator<VariantContext> iter=null;
		private GvcfThing current=null;
		private String sampleName=null;
		Source(final File gvcfFile) {
			this.gvcfFile = gvcfFile;
			}
		void open()
			{	
			this.vcfFileReader = new VCFFileReader(this.gvcfFile, false);
			this.iter = vcfFileReader.iterator();
			final List<String> samples = vcfFileReader.getFileHeader().getSampleNamesInOrder();
			if(samples.size()!=1) {
				throw new IllegalArgumentException("Expected "+this.gvcfFile+" to contain only one sample but got "+samples);
				
				}
			this.sampleName=samples.get(0);
			}
		
		GvcfThing get()
			{
			if(this.current!=null) return this.current;
			if(!this.iter.hasNext()) return null;
			final VariantContext ctx= this.iter.next();
			final List<Allele> altAlleles = ctx.getAlternateAlleles();
			if(altAlleles.size()==1 && altAlleles.get(0).equals(NON_REF))
				{
				this.current = new GvcfBlock(ctx);
				}
			else
				{
				this.current = new GvcfVariant(ctx);
				if(this.current.ctx.getGenotype(0).getAlleles().contains(NON_REF)){
					LOG.warn("Ignoring variant "+ctx+" because a genotype is "+NON_REF+" (bug?)");
					this.current=null;
					return get();
					}
				}
			
			return this.current;
			}
		
		void close() {
			if(this.vcfFileReader!=null) {
				this.vcfFileReader.close();
				CloserUtil.close(this.iter);
				this.iter=null;
				this.vcfFileReader=null;
				}
			}
		
		}
	
	@Override
	public int doWork(final List<String> args) {
		VariantContextWriter w=null;
		try {
			this.gvcfSources.addAll(IOUtil.unrollFiles(args.stream().map(S->new File(S)).collect(Collectors.toSet()),
					".g.vcf",".g.vcf.gz"
					).stream().map(F->new Source(F)).collect(Collectors.toList()));
			if(args.isEmpty())
				{
				LOG.error("No gvcf file was given");
				return -1;
				}
			this.gvcfSources.stream().forEach(S->S.open());
			this.dictionary  = this.gvcfSources.get(0).vcfFileReader.getFileHeader().getSequenceDictionary();
			if(this.dictionary==null)
				{
				LOG.error("Dict missing in "+this.gvcfSources.get(0).gvcfFile);
				return -1;
				}
		
			
			this.gvcfSources.stream().map(S->S.vcfFileReader.getFileHeader().getSequenceDictionary()).forEach(D->{
				if(D==null || !SequenceUtil.areSequenceDictionariesEqual(D, dictionary))
					{
					throw new JvarkitException.UserError("dict missing or dict are not the same");
					}
				});
			
			
			if(	gvcfSources.stream().map(S->S.sampleName).collect(Collectors.toSet()).stream().count() != this.gvcfSources.size())
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
					gvcfSources.stream().map(S->S.sampleName).collect(Collectors.toList())
					);
			
			w= super.openVariantContextWriter(outputFile);
			w.writeHeader(header);
			
			
			int contigTid=0;
			while(contigTid< dictionary.size())
				{
				final SAMSequenceRecord ssr = this.dictionary.getSequence(contigTid);
				int pos=0;
				while(pos< ssr.getSequenceLength())
					{
					//LOG.debug(ssr.getSequenceName()+" "+pos+" "+this.gvcfSources.size());
					GvcfVariant variantAtThisPos = null;
					int minEnd=ssr.getSequenceLength();
					//cleanup
					for(final Source src:this.gvcfSources)
						{
						for(;;) {
							final GvcfThing gvcfthing = src.get();
							// LOG.debug(""+gvcfthing+" "+src.sampleName+" "+ssr.getSequenceName()+":"+pos);
							if(gvcfthing==null)
								{
								//no more variant avaialble
								break;
								}
							else if(contigTid> gvcfthing.getTid()) {
								//observed contig is after gvcfthing.contig
								src.current=null;
								continue;
								}
							else if(contigTid< gvcfthing.getTid()) {
								//observed contig is before gvcfthing.contig
								break;
								}
							else if( gvcfthing.getEnd() < pos) {
								// variant information is before observed pos
								src.current=null;
								continue;
								}
							else if( gvcfthing.getStart() > pos) {
								// variant information is after observed pos
								minEnd=Math.min(minEnd, gvcfthing.getStart()-1);
								break;
								}
							else if(gvcfthing.isVariant())
								{
								if(variantAtThisPos==null || 
									variantContigPosRefComparator.compare(GvcfVariant.class.cast(gvcfthing).ctx, variantAtThisPos.ctx)<0)
									{
									variantAtThisPos = GvcfVariant.class.cast(gvcfthing);
									}
								break;
								}
							else if(gvcfthing.isBlock())
								{
								minEnd=Math.min(minEnd, gvcfthing.getEnd());
								break;
								}
							else
								{
								LOG.debug("??");
								}
							}
						}
					
					if(variantAtThisPos==null)
						{
						pos = minEnd +1;
						}
					else
						{
						final VariantContext archetype = variantAtThisPos.ctx;
						final List<VariantContext> allVariants =
								this.gvcfSources.stream().map(S->S.get()).
								filter(G->G!=null && G.isVariant()).
								map(G->GvcfVariant.class.cast(G).ctx).
								filter(V->variantContigPosRefComparator.compare(V, archetype)==0).
								collect(Collectors.toList());
								;
						final Set<Allele> alleles = allVariants.stream().
								flatMap(V->V.getGenotypes().stream()).
								flatMap(G->G.getAlleles().stream()).
								filter(A->!(A.equals(NON_REF) || A.isNoCall())).
								collect(Collectors.toSet());
						
						alleles.add(archetype.getReference());
						
						
						final VariantContextBuilder vcb=new VariantContextBuilder(
								getClass().getName(),
								archetype.getContig(), 
								archetype.getStart(),
								archetype.getEnd(),
								alleles
								);
						
						if(archetype.hasID())
							{
							vcb.id(archetype.getID());
							}
						final List<Genotype> genotypes = new ArrayList<>(allVariants.size());
						for(VariantContext ctx: allVariants)
							{
							Genotype genotype = ctx.getGenotype(0);
							GenotypeBuilder gb = new GenotypeBuilder(genotype);
							genotypes.add(gb.make());
							}
						vcb.genotypes(genotypes);
						
						final VariantContext genotypedVariant ;
						
						try {
						genotypedVariant= vcb.make();
						
						} catch(Exception err2)
							{
							LOG.debug(allVariants);
							LOG.error(err2);
							return -1;
							}
						w.add(genotypedVariant);
						

						// reset source for the current variant
						for(final Source src:this.gvcfSources)
							{
							if(src.current!=null && 
								variantContigPosRefComparator.compare(src.current.ctx,archetype)==0)
								{
								src.current=null;
								}
							}
						}
					
					}
				++contigTid;
				}
			
			w.close();w=null;
			
			this.gvcfSources.stream().forEach(S->S.close());
			this.gvcfSources.clear();
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}	
		finally
			{
			for(final Source src:this.gvcfSources) src.close();
			CloserUtil.close(w);
			}
		}
	public static void main(String[] args) {
		new FastGenotypeGVCFs().instanceMainWithExit(args);

	}

}
