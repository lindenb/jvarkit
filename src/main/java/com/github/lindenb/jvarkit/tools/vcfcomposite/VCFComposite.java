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
package com.github.lindenb.jvarkit.tools.vcfcomposite;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.BiPredicate;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.AbstractVCFCodec;
import htsjdk.variant.vcf.VCFEncoder;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;
import com.github.lindenb.jvarkit.util.vcf.JexlGenotypePredicate;
import com.github.lindenb.jvarkit.util.vcf.JexlVariantPredicate;
import htsjdk.variant.vcf.VCFIterator;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParserFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParserFactory;
/*
BEGIN_DOC

## Input

input is a VCF file annotated with SNPEff or VEP.



END_DOC
*/
@Program(name="vcfcomposite",
	description="(in developpement) Finds Variants involved in a Het Composite Disease",
	keywords={"vcf","disease","annotation","pedigree"}
)
public class VCFComposite extends Launcher {
	private static final String INFO_TAG="COMPOSITE";
	private static final Logger LOG= Logger.build(VCFComposite.class).make();
	@Parameter(names={"-p","-ped","--pedigree"},description=Pedigree.OPT_DESCRIPTION)
	private File pedigreeFile=null;
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile=null;
	@Parameter(names={"-vf","--variant-filter"},description=JexlVariantPredicate.PARAMETER_DESCRIPTION,converter=JexlVariantPredicate.Converter.class)
	private Predicate<VariantContext> variantJexl = JexlVariantPredicate.create("");
	@Parameter(names={"-gf","--genotype-filter"},description=JexlGenotypePredicate.PARAMETER_DESCRIPTION,converter=JexlGenotypePredicate.Converter.class)
	private BiPredicate<VariantContext,Genotype> genotypeFilter = JexlGenotypePredicate.create("");
	@Parameter(names={"-max","--max","--max-variants"},description="[20180718] Max variants per gene. If different than -1, used to set a maximum number of variants per gene; The idea is to filter out the gene having a large number of variants.")
	private int max_number_of_variant_per_gene = -1;
	@Parameter(names={"--filter"},description="[20180718] set FILTER for the variants that are not part of a composite mutation.")
	private String filterTag = "NOT_COMPOSITE";

	
	
	//@Parameter(names={"-m","--model"},description="Model type",required=true)
	//private Type modelType=null;
	//@Parameter(names={"-models"},description="List the available models and exits",help=true)
	//private boolean listModels=false;
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection = new WritingSortingCollection();

	
	private AbstractVCFCodec vcfDecoder = null;
	private VCFEncoder vcfEncoder = null;
	private Pedigree pedigree = null;
	
	private class VariantLine
		{
		long id;
		VariantContext ctx;
		VariantLine(long id,VariantContext ctx) {
			this.id = id;
			this.ctx = ctx;
			}
		}
	private class VariantLineCodec extends AbstractDataCodec<VariantLine>
		{
		@Override
		public VariantLine decode(final DataInputStream dis) throws IOException {
			try {
				long n = dis.readLong();
				final VariantContext ctx = vcfDecoder.decode(IOUtils.readString(dis));
				return new VariantLine(n,ctx);
				} 
			catch(final EOFException err)
				{
				return null;
				}
			}

		@Override
		public void encode(DataOutputStream dos, VariantLine object) throws IOException {
			dos.writeLong(object.id);
			IOUtils.writeString(dos,vcfEncoder.encode(object.ctx));
			}

		@Override
		public VariantLineCodec clone() {
			return new VariantLineCodec();
			}

		}
	
	private class GeneIdentifier
		implements Comparable<GeneIdentifier>
		{
		String geneName;
		String contig="";//not null please
		String source;
		GeneIdentifier() {
			}
		GeneIdentifier(final String geneName,final String source) {
			this.geneName = geneName;
			this.source = source;
			}

		@Override
		public int hashCode() {
			return this.geneName.hashCode()*31 + this.source.hashCode();
			}
		@Override
		public boolean equals(final Object obj) {
			return this.compareTo(GeneIdentifier.class.cast(obj))==0;
			}
		@Override
		public int compareTo(final GeneIdentifier o) {
			int i = this.geneName.compareTo(o.geneName);
			if(i!=0) return i;
			i = this.contig.compareTo(o.contig);
			if(i!=0) return i;
			i = this.source.compareTo(o.source);
			return i;
			}
		@Override
		public String toString() {
			return geneName+"("+source+") on "+contig;
			}
		}

	private class GeneIdentifierCodec
		extends AbstractDataCodec<GeneIdentifier>
		{
		@Override
		public GeneIdentifier decode(DataInputStream in) throws IOException {
			final GeneIdentifier gi = new GeneIdentifier();
			try {
				gi.geneName = in.readUTF();
				gi.contig = in.readUTF();
				gi.source = in.readUTF();
				return gi;
			} catch(EOFException err) 
				{
				return null;
				}
			}
		
		@Override
		public void encode( final DataOutputStream dos, final GeneIdentifier object) throws IOException {
			dos.writeUTF(object.geneName);
			dos.writeUTF(object.contig);
			dos.writeUTF(object.source);
			}
		
		@Override
		public GeneIdentifierCodec clone() {
			return new GeneIdentifierCodec();
			}
		}
	
	private class GeneAndVariant
		{
		final GeneIdentifier gene;
		final VariantLine variant;
				
		public GeneAndVariant(final GeneIdentifier geneIdentifier, final VariantLine variant) {
			this.gene = geneIdentifier;
			this.variant = variant;
			}
		
		int compareGene(final GeneAndVariant other) {
			return this.gene.compareTo(other.gene);
			}
		
		int compareGeneThenIndex(final GeneAndVariant other) {
			final int i = compareGene(other);
			if(i!=0) return i;
			return Long.compare(this.variant.id, other.variant.id);
			}
		}
	
	private class GeneAndVariantCodec
		extends AbstractDataCodec<GeneAndVariant>
		{
		final GeneIdentifierCodec geneCodec = new GeneIdentifierCodec();
		final VariantLineCodec vcCodec = new VariantLineCodec();
		@Override
		public GeneAndVariant decode(DataInputStream dis) throws IOException {
			try {
				final GeneIdentifier gi  = geneCodec.decode(dis);
				final VariantLine vl = vcCodec.decode(dis);
				return new GeneAndVariant(gi, vl);
				}
			catch(final EOFException err)
				{
				return null;
				}
			}
		@Override
		public void encode(DataOutputStream dos, GeneAndVariant object) throws IOException {
			geneCodec.encode(dos,object.gene);
			vcCodec.encode(dos,object.variant);
			}
		@Override
		public GeneAndVariantCodec clone() {
			return new GeneAndVariantCodec();
			}
		}
	
	
	
	public static enum Type {
		//RecessiveHomVar(""),
		HetRecessiveComposite("child carry at least two HET variant in the same gene unit.")
		;
		
		private final String desc;
		Type(final String desc ) {
			this.desc = desc;
		}
		
		public String getDescription() {
			return desc;
		}
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
		
		
		public abstract void scan(
				final GeneIdentifier geneName,
				final List<VariantLine> variants
				);
		}
	
	
	/* affected individual contains at least two HET variants  */
	private class HetRecessiveComposite extends DiseaseModel
		{
		/** affected can be HET */
		private boolean isGenotypeForAffected(final Genotype g)
			{
			return g!=null && g.isHet() ;
			}
		
		@Override
		public boolean accept(final VariantContext ctx) {
			if(getGenotypesOfAffected(ctx).stream().
					noneMatch(G->isGenotypeForAffected(G))) return false;			
			if(getGenotypesOfUnaffected(ctx).stream().anyMatch(G->G.isHomVar())) return false;
			return true;
			}
		
		@Override
		public void scan(
				final GeneIdentifier geneKey,
				final List<VariantLine> variants)
			{
			if(variants.size()<2) {
				return;
				}
			// test if this gene is a big pool of variant, false positive. e.g: highly variables genes.
			if(max_number_of_variant_per_gene>=0) {
				if(variants.size()< max_number_of_variant_per_gene) {
					LOG.warn("Too many variants "+variants.size()+" for "+geneKey);
					return;
				}
			}
			
			/* loop over affected */
			for(final Pedigree.Person child: pedigree.getAffected()) {
				for(int x=0;x+1< variants.size();++x)
					{
					final VariantContext vcx = variants.get(x).ctx;
					final Genotype gcx = vcx.getGenotype(child.getId());
					// child variant  n. y  must be HOM_VAR or HET
					if(gcx==null || !isGenotypeForAffected(gcx)) continue;
					// filtered ?
					if(!genotypeFilter.test(vcx,gcx)) continue;
					// search for the second snp
					for(int y=x+1;y< variants.size();++y)
						{
						final VariantContext vcy = variants.get(y).ctx;
						final Genotype gcy = vcy.getGenotype(child.getId());
						// child variant n. y must be HOM_VAR or HET
						if(gcy==null || !isGenotypeForAffected(gcy)) continue;
						// filtered ?
						if(!genotypeFilter.test(vcy,gcy)) continue;
						
						boolean unaffected_are_ok=true;
						//check unaffected indididual don't have same haplotype
						for(final Pedigree.Person unaffected: pedigree.getUnaffected()) {
							final Genotype gux = variants.get(x).ctx.getGenotype(unaffected.getId());
							if(gux!=null && gux.isHomVar()) {
								unaffected_are_ok=false;
								break;
								}
							
							final Genotype guy = variants.get(y).ctx.getGenotype(unaffected.getId());
							
							if(guy!=null && guy.isHomVar()) {
								unaffected_are_ok=false;
								break;
								}
							
							if(gux!=null && guy!=null &&
								gux.sameGenotype(gcx, true) &&
								guy.sameGenotype(gcy, true)
								)
								{
								unaffected_are_ok=false;
								break;
								}
							}
						if(unaffected_are_ok)
							{
							for(int side=0;side<2;++side) 
								{
								final VariantContext vc=(side==0?vcx:vcy);
								final Set<String> set = getAnnotationsForVariant(vc);
								final VariantContextBuilder vcb = new VariantContextBuilder(vc);
								
								
								final StringBuilder sb=new StringBuilder();
								sb.append("gene|").append(geneKey.geneName);
								sb.append("|source|").append(geneKey.source);
								if(side==1) {
									sb.append("|pos|").append(vcx.getStart());
									sb.append("|ref|").append(vcx.getReference().getDisplayString());
									}
								else {
									sb.append("|pos|").append(vcy.getStart());
									sb.append("|ref|").append(vcy.getReference().getDisplayString());
									}
								sb.append("|sample|").append(child.getId());
								set.add(sb.toString());
								set.remove("");
								vcb.attribute(INFO_TAG, new ArrayList<>(set));
								
								if(side==0)
									{
									variants.get(x).ctx = vcb.make();
									}
								else
									{
									variants.get(y).ctx = vcb.make();
									}
								}
							}
						}
					}
				}
			}
		}

	
	protected DiseaseModel createModel() {
		return new HetRecessiveComposite();
		}
	
	private Set<String> getAnnotationsForVariant(VariantContext vc) {
		return new HashSet<>(vc.getAttributeAsStringList(INFO_TAG, ""));
		}
	
	@Override
	protected int doVcfToVcf(final String inputName,
		final VCFIterator iterin,
		final VariantContextWriter out) {
		final DiseaseModel model = this.createModel();

		final VCFHeader header = iterin.getHeader();
		if(this.pedigreeFile!=null)
			{
			try {
				this.pedigree = Pedigree.newParser().parse(this.pedigreeFile);
				}
			catch(final IOException err)
				{
				throw new RuntimeIOException(err);
				}
			}
		else
			{
			this.pedigree = Pedigree.newParser().parse(header);
			}
		
		if(this.pedigree==null) {
			LOG.error("pedigree missing");
			return -1;
			}
		
		if(this.pedigree.getAffected().isEmpty()) {
			LOG.error("No Affected sample in pedigre. " + this.pedigree);
			return -1;
			}
		
		if(this.pedigree.getUnaffected().isEmpty()) {
			LOG.error("No Unaffected sample in " + this.pedigree);
			return -1;
			}
		header.addMetaDataLine(new VCFInfoHeaderLine(INFO_TAG, VCFHeaderLineCount.UNBOUNDED,VCFHeaderLineType.String,"Variant of VCFComposite"));
		header.addMetaDataLine(new VCFFilterHeaderLine(this.filterTag, "Not a Variant fir VCFComposite"));
		final SAMSequenceDictionary dict = header.getSequenceDictionary();
		final Comparator<String> contigCmp;
		if(dict==null || dict.isEmpty())
			{
			contigCmp = (A,B)->A.compareTo(B); 
			}
		else
			{
			contigCmp = new ContigDictComparator(dict);
			}
		
		final Comparator<VariantContext> ctxComparator = (V1,V2)->{
			int i = contigCmp.compare(V1.getContig(), V2.getContig());
			if(i!=0) return i;
			i = Integer.compare(V1.getStart(), V2.getStart());
			if(i!=0) return i;
			return V1.getReference().compareTo(V2.getReference());
			};
		
		final Comparator<VariantLine> variantLineComparator = (V1,V2)->{
			final int i = ctxComparator.compare(V1.ctx,V2.ctx);
			if(i!=0) return i;
			return Long.compare(V1.id, V2.id);
			};


		long ID_GENERATOR = 0L;
		this.vcfDecoder = VCFUtils.createDefaultVCFCodec();//iterin.getCodec();
		this.vcfEncoder = new VCFEncoder(header, false, true);
		final AnnPredictionParser annParser=new AnnPredictionParserFactory(header).get();
		final VepPredictionParser vepParser=new VepPredictionParserFactory(header).get();
		SortingCollection<GeneAndVariant> sorting=null;
		SortingCollection<VariantLine> outputSorter = null;

		try
			{
			LOG.info("reading variants and genes");
			/* Gene and variant sorter */
			sorting = SortingCollection.newInstance(GeneAndVariant.class,
					new GeneAndVariantCodec(),
					GeneAndVariant::compareGeneThenIndex,
					this.writingSortingCollection.getMaxRecordsInRam(),
					this.writingSortingCollection.getTmpPaths()
					);
			sorting.setDestructiveIteration(true);
			/* Variant sorter */
			outputSorter = SortingCollection.newInstance(
					VariantLine.class,
					new VariantLineCodec(),
					variantLineComparator,
					this.writingSortingCollection.getMaxRecordsInRam(),
					this.writingSortingCollection.getTmpPaths()
					);
			outputSorter.setDestructiveIteration(true);
			
			/* read input */
			while(iterin.hasNext()) {
				final VariantContext ctx = iterin.next();
				final VariantLine variantLine = new VariantLine(++ID_GENERATOR,ctx);
				if(!this.variantJexl.test(ctx))
					{
					outputSorter.add(variantLine);
					continue;
					}
				
				if(!model.accept(ctx))
					{
					outputSorter.add(variantLine);
					continue;
					}
				
				final Set<GeneIdentifier> geneKeys = new HashSet<>();
				
				for(final AnnPredictionParser.AnnPrediction pred: annParser.getPredictions(ctx)) {
					if(pred.isIntergenicRegion())
						{
						continue;
						}
					
					if(!StringUtil.isBlank(pred.getGeneName())) 
						{
						geneKeys.add(new GeneIdentifier(pred.getGeneName(),"ANN_GeneName"));
						}
					if(!StringUtil.isBlank(pred.getGeneId())) 
						{
						geneKeys.add(new GeneIdentifier(pred.getGeneId(),"ANN_GeneId"));
						}
					if(!StringUtil.isBlank(pred.getFeatureId())) 
						{
						geneKeys.add(new GeneIdentifier(pred.getFeatureId(),"ANN_FeatureId"));
						}
					}
				
				for(final VepPredictionParser.VepPrediction pred: vepParser.getPredictions(ctx)) {
					if(!StringUtil.isBlank(pred.getGene())) 
						{
						geneKeys.add(new GeneIdentifier(pred.getGene(),"VEP_Gene"));
						}
					if(!StringUtil.isBlank(pred.getFeature())) 
						{
						geneKeys.add(new GeneIdentifier(pred.getFeature(),"VEP_Feature"));
						}
					if(!StringUtil.isBlank(pred.getSymbol())) 
						{
						geneKeys.add(new GeneIdentifier(pred.getSymbol(),"VEP_Symbol"));
						}
					if(!StringUtil.isBlank(pred.getHgncId())) 
						{
						geneKeys.add(new GeneIdentifier(pred.getHgncId(),"VEP_HgncId"));
						}
					}
				if(geneKeys.isEmpty()) {
					outputSorter.add(variantLine);
					continue;
					}
				
				for(final GeneIdentifier gk:geneKeys) 
					{
					final GeneAndVariant gav = new GeneAndVariant(gk,variantLine);
					gav.gene.contig = ctx.getContig();
					sorting.add(gav);
					}
				}
			sorting.doneAdding();
			
			LOG.info("compile per gene");
			//compile data
			CloseableIterator<GeneAndVariant> iter2=sorting.iterator();
			EqualRangeIterator<GeneAndVariant> eqiter= new EqualRangeIterator<>(iter2,(A,B)->A.gene.compareTo(B.gene));
			while(eqiter.hasNext())
				{
				final List<GeneAndVariant> variants=eqiter.next();
				model.scan(
					variants.get(0).gene,
					variants.stream().
						map(L->L.variant).
						collect(Collectors.toList())
					);
				for(final GeneAndVariant ga:variants) outputSorter.add(ga.variant);
				}
			eqiter.close();
			iter2.close();
			sorting.cleanup();
			//
			
			
			LOG.info("write variants");
			CloseableIterator<VariantLine> iter1 = outputSorter.iterator();
			EqualRangeIterator<VariantLine > eqiter1 = new EqualRangeIterator<>(iter1,variantLineComparator);
			out.writeHeader(header);
			while(eqiter1.hasNext())
				{
				final List<VariantLine> array = eqiter1.next();
				final VariantContext firstCtx = array.get(0).ctx;
				final Set<String> set= getAnnotationsForVariant(firstCtx);
				final VariantContext outCtx;
				final VariantContextBuilder vcb = new VariantContextBuilder(firstCtx);
				for(int y=1;y<array.size();++y) {
					set.addAll(getAnnotationsForVariant(array.get(y).ctx));
					}
				if(set.isEmpty())
					{
					vcb.filter(this.filterTag);
					}
				else
					{
					if(!firstCtx.isFiltered())
						{
						vcb.passFilters();
						}
					vcb.attribute(INFO_TAG, new ArrayList<>(set));
					}

				outCtx = vcb.make();
					
				out.add(outCtx);
				}
			outputSorter.cleanup();
			eqiter1.close();
			iter1.close();
			
			return 0;
			}
		catch(final Exception err) {
			LOG.error(err);
			return -1;
			}
		}
	
	
	@Override
	public int doWork(final List<String> args) {
		if(StringUtil.isBlank(this.filterTag)) 
			{
			LOG.error("bad name for FILTER");
			return -1;
			}
		try {
			/*
			if(this.listModels)
				{
				final PrintWriter out= super.openFileOrStdoutAsPrintWriter(this.outputFile);
				for(final Type t:Type.values()) {
					out.println(t.name());
					out.println("\t"+t.getDescription());
					}
				out.flush();
				out.close();
				return 0;
				}*/
			
			return doVcfToVcf(args, this.outputFile);
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		}
	public static void main(final String[] args) {
		new VCFComposite().instanceMainWithExit(args);
		}
	
	}
