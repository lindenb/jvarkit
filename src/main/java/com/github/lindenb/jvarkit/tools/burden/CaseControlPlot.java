package com.github.lindenb.jvarkit.tools.burden;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;
import java.util.function.Predicate;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.semontology.Term;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
/**
BEGIN_DOC

END_DOC
 */
@Program(
	name="casecontrolplot",
	description="Plot CASE/CTRL data from VCF files",
	keywords={"maf","burden","case","control","plot","chart","vcf"},
	terms=Term.ID_0000018
	)
public class CaseControlPlot extends Launcher
	{
	private static final Logger LOG=Logger.build(CaseControlPlot.class).make();
	@Parameter(names={"-o","--out"},description="Output Directory, or a filename ending with '.zip'",required=true)
	private File outputDirOrZip=null;
	@Parameter(names={"-tee","--tee"},description="Output the incoming VCF to stdout. Useful to capture intermediate results in pipelines.")
	private boolean teeToStdout=false;
	@Parameter(names={"-format","--format"},description="How to print doubles, printf-life precision format")
	private String precisionFormat="%.5f";
	@Parameter(names={"-prefix","--prefix"},description="Output files prefix")
	private String prefix="";
	@Parameter(names={"-ped","--pedigree"},description=Pedigree.OPT_DESCRIPTION+" If not defined, I will try to extract the pedigree from the VCFheader.")
	private File pedigreeFile;
	
	
	public static interface MafExtractor
		{
		public Double apply(final VariantContext ctx,Allele alt);
		}
	
	public static class AttributeMafExtractor implements MafExtractor
		{
		private String attribute = null;
		public AttributeMafExtractor() {this(null);}
		public AttributeMafExtractor(final String tag) {this.attribute = tag;}
		public void setAttribute(String tag)
			{
			this.attribute = tag;
			}
		public String getAttribute()
			{
			return attribute;
			}
		@Override
		public Double apply(final VariantContext ctx,final Allele alt) {
			final String att = this.getAttribute();
			if(att==null || att.isEmpty()) return null;
			int index=ctx.getAlleleIndex(alt);
			if( index<=0) return null;//can't be REF==0
			final List<Object> L = ctx.getAttributeAsList(att);
			if(index>=L.size()) return null;
			final Object o = L.get(index);
			if(o==null || ".".equals(o)) return null;
			try {
				double f = Double.parseDouble(String.valueOf(o));
				if(f<0.0 || f>1.0) return null;
				return f;
				}
			catch(NumberFormatException err) {
				return null;
				}
			}
		}
	
	public static class GenotypeMafExtractor implements MafExtractor
		{
		@Override
		public Double apply(final VariantContext ctx,final Allele alt) {
			return null;
			}
		}
	
	
	
	/** interface to extract X/Y data from a variant */
	public static class CaseControlExtractor
		implements Function<VariantContext, List<Point2D.Double>>
		{
		private String name;
		private Predicate<VariantContext> variantPredicate = C -> true;
		private MafExtractor caseExtractor = new GenotypeMafExtractor();
		private MafExtractor ctrlExtractor= new GenotypeMafExtractor();
		/** name of this extractor */
		public String getName() {
			return this.name;
			}
		/** Variant context filter: should I accept this variant ? */
		public Predicate<VariantContext> getFilter()
			{
			return this.variantPredicate;
			}
		/** case extractor */
		public MafExtractor getCaseMAFExtractor() {
			return this.caseExtractor;
			}
		/** control extractor */
		public MafExtractor getControlMAFExtractor() {
			return this.ctrlExtractor;
			}
		@Override
		public List<Point2D.Double> apply(final VariantContext vc)
			{
			if(vc==null || !vc.isVariant()) return  Collections.emptyList();
			if(getFilter()!=null)
				{
				if(!getFilter().test(vc)) return Collections.emptyList();
				}
			final List<Allele> alts = vc.getAlternateAlleles();
			final List<Point2D.Double> points = new ArrayList<>(alts.size());
			for(int altidx=0; altidx < alts.size();++altidx)
				{
				final Allele alt= alts.get(altidx);
				final Double casex = this.getCaseMAFExtractor().apply(vc, alt);
				if(casex==null || casex<0.0 || casex>1.0) continue;
				final Double ctrly = this.getControlMAFExtractor().apply(vc, alt);
				if(ctrly==null || ctrly<0.0 || ctrly>1.0) continue;
				points.add(new Point2D.Double(casex,ctrly));
				}
			return points;
			}
		}
	
	private class CaseControlHandler
		{
		final String baseName;
		CaseControlExtractor extractor;
		PrintWriter pw;
		
		CaseControlHandler(final String baseName) {
		this.baseName = baseName;
		}
		
		void visit(final VariantContext ctx) {
			final List<Point2D.Double> points = this.extractor.apply(ctx);
			for(final Point2D.Double pt : points ) {
				pw.printf(precisionFormat,pt.getX());
				pw.print('\t');
				pw.printf(precisionFormat,pt.getY());
				pw.println();
				}
			}
		public void close() {
			CloserUtil.close(pw);
			}
		}
	
	@Override
	public int doWork(final List<String> args)
		{
		ArchiveFactory archiveFactory=null;
		VcfIterator in = null;
		VariantContextWriter teeVariantWriter = null;
		final List<CaseControlHandler> handlers= new ArrayList<>();
		final List<CaseControlExtractor> excractors= new ArrayList<>();
		try
			{
			
			in = super.openVcfIterator(oneFileOrNull(args));
			final VCFHeader header= in.getHeader();
			final Pedigree pedigree;
			if( this.pedigreeFile!=null) {
				pedigree = Pedigree.newParser().parse(this.pedigreeFile);
				}
			else
				{
				pedigree = Pedigree.newParser().parse(header);
				}
			if(pedigree==null || pedigree.isEmpty()) {
				LOG.error("No pedigree defined , or it is empty");
				return -1;
				}
			if(!pedigree.getPersons().stream().
					filter(F->F.isAffected() && header.getSampleNamesInOrder().indexOf(F.getId())!=-1).
					findFirst().isPresent()){
					LOG.error("No Affected individuals in pedigree/header");
					return -1;
					}
			if(!pedigree.getPersons().stream().
					filter(F->F.isUnaffected() && header.getSampleNamesInOrder().indexOf(F.getId())!=-1).
					findFirst().isPresent()){
					LOG.error("No Unaffected individuals in pedigree/header");
					return -1;
					}
			if( this.teeToStdout)
				{
				teeVariantWriter = super.openVariantContextWriter(null);
				teeVariantWriter.writeHeader(header);
				}
			archiveFactory = ArchiveFactory.open(this.outputDirOrZip);
			for(final CaseControlExtractor extractor : excractors) {
				final CaseControlHandler handler = new CaseControlHandler(this.prefix + extractor.getName());
				handler.pw = archiveFactory.openWriter(handler.baseName);
				handlers.add(handler);
				}
			final SAMSequenceDictionaryProgress progress= new SAMSequenceDictionaryProgress(header);
			while(in.hasNext())
				{
				final VariantContext ctx = progress.watch(in.next());
				if(teeVariantWriter!=null) teeVariantWriter.add(ctx);
				for(final CaseControlHandler handler : handlers) {
					handler.visit(ctx);
					}
				}
			for(final CaseControlHandler handler : handlers) {
				handler.close();
				}
			progress.finish();
			if(teeVariantWriter!=null) {
				teeVariantWriter.close();
				teeVariantWriter=null;
				}
			in.close();in=null;
			archiveFactory.close();archiveFactory=null;
			return RETURN_OK;
			}
		catch (final Exception e)
			{
			LOG.error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(in);
			CloserUtil.close(archiveFactory);
			CloserUtil.close(teeVariantWriter);
			for(final CaseControlHandler handler : handlers) {
				handler.close();
				}
			}
		}
	
	public static void main(final String[] args)
		{
		new CaseControlPlot().instanceMainWithExit(args);
		}

	}
