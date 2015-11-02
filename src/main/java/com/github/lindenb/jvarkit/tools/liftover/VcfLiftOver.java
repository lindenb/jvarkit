package com.github.lindenb.jvarkit.tools.liftover;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryFactory;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class VcfLiftOver extends AbstractVcfLiftOver
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(VcfLiftOver.class);

	@Override
	public Command createCommand() {
		return new MyCommand();
		}

	static private class MyCommand extends AbstractVcfLiftOver.AbstractVcfLiftOverCommand
		{    
	private LiftOver liftOver=null;
	private SAMSequenceDictionary newDict=null;

	protected Allele revcomp(Allele a)
		{
		if(a.isNoCall()) return a;
		if(a.isSymbolic()) return a;
		String seq=a.getBaseString();
		StringBuilder sb=new StringBuilder(seq.length());
		for(int i=seq.length()-1;i>=0;--i)
			{
			sb.append(AcidNucleics.complement(seq.charAt(i)));
			}
		return Allele.create(sb.toString(), a.isReference());
		}

	@Override
		protected Collection<Throwable> doVcfToVcf(String inputName,
				VcfIterator in, VariantContextWriter out) throws IOException {
		final String TAG="LIFTOVER";
		VariantContextWriter failed=null;
		
		VCFHeader header=in.getHeader();
		

		if(this.failedFile!=null)
			{
			VCFHeader header2=new VCFHeader(header);
			addMetaData(header2);
			failed=VCFUtils.createVariantContextWriter(failedFile);
			failed.writeHeader(header2);
			}
		
		VCFHeader header3;
		
		if(newDict==null)
			{
			header3=new VCFHeader(header);
			LOG.warn("##contig files should be changed.");
			}
		else
			{
			Set<VCFHeaderLine> hh=new LinkedHashSet<VCFHeaderLine>();
			for(VCFHeaderLine h:header.getMetaDataInInputOrder())
				{
				if(!h.getKey().equals(VCFConstants.CONTIG_HEADER_KEY))
					{
					System.err.println(h);
					hh.add(h);
					}
				}
			hh.addAll(VCFUtils.samSequenceDictToVCFContigHeaderLine(newDict));
			
			header3=new VCFHeader(
					hh,header.getGenotypeSamples()
					);
			}
		addMetaData(header3);
		header3.addMetaDataLine(new VCFInfoHeaderLine(TAG,1,VCFHeaderLineType.String,"Chromosome|Position before liftOver."));
		out.writeHeader(header3);
		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
		while(in.hasNext())
			{
			VariantContext ctx=in.next();
			progress.watch(ctx.getContig(),ctx.getStart());
			
			Interval lifted=liftOver.liftOver(
					new Interval(ctx.getContig(),ctx.getStart(),ctx.getEnd(),
					false,//negative strand
					""));
			if(lifted==null )
				{
				if(failed!=null) failed.add(ctx);
				}
			else
				{
				VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				vcb.chr(lifted.getContig());
				vcb.start(lifted.getStart());
				vcb.stop(lifted.getEnd());
				vcb.attribute(TAG,ctx.getContig()+"|"+ctx.getStart() );
				
				if(lifted.isNegativeStrand())//strandess has changed
					{
					List<Genotype> genotypes=new ArrayList<Genotype>(header.getSampleNamesInOrder().size());
					Set<Allele> alleles=new HashSet<Allele>();
					alleles.add(revcomp(ctx.getReference()));
					for(String sample:header.getSampleNamesInOrder())
						{
						Genotype g=ctx.getGenotype(sample);
						GenotypeBuilder gb=new GenotypeBuilder(g);
						List<Allele> alleles2=new ArrayList<Allele>();
						for(Allele a0:g.getAlleles())
							{
							alleles2.add(revcomp(a0));
							}
						alleles.addAll(alleles2);
						gb.alleles(alleles2);
						genotypes.add(gb.make());
						}
					vcb.genotypes(genotypes);
					vcb.alleles(alleles);
					}
				out.add(vcb.make());
				}
			}
		CloserUtil.close(failed);	
		return RETURN_OK;
		}
	
	@Override
		protected Collection<Throwable> call(String inputName) throws Exception
			{
			try {
				if(liftOverFile==null)
				{
				return wrapException("LiftOver file is undefined.");
				}
			if(REF==null)
				{
				return wrapException("REF undefined.");
				}
			this.newDict=new SAMSequenceDictionaryFactory().load(REF);
			this.liftOver=new LiftOver(liftOverFile);
			this.liftOver.setLiftOverMinMatch(minMatch);
			return doVcfToVcf(inputName);
				} catch (Exception err) {
				return wrapException(err);
				}
			}
		}
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new VcfLiftOver().instanceMainWithExit(args);
		}

	}
