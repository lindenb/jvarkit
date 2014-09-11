package com.github.lindenb.jvarkit.tools.vcffixindels;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Pattern;

import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter2;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

import htsjdk.tribble.TribbleException;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;



public class VCFFixIndels extends AbstractVCFFilter2
	{
	private VCFFixIndels()
		{
		}
	
	@Override
	public String getProgramDescription() {
		return " Fix samtools indels (for @SolenaLS).";
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/VCFFixIndels";
		}
	
	
	@Override
	protected void doWork(VcfIterator r, VariantContextWriter w)
			throws IOException
		{
		long nChanged=0L;
		final String TAG="INDELFIXED";
		VCFHeader header=r.getHeader();
		
		VCFHeader h2=new VCFHeader(header.getMetaDataInInputOrder(),header.getSampleNamesInOrder());
		h2.addMetaDataLine(new VCFInfoHeaderLine(TAG,1,VCFHeaderLineType.String,"Fix Indels for @SolenaLS."));
		
		w.writeHeader(h2);
	
		final Pattern dna=Pattern.compile("[ATGCatgc]+");
		while(r.hasNext())
			{
			VariantContext ctx=r.next();
			VariantContextBuilder b=new VariantContextBuilder(ctx);
			List<Allele> alleles=ctx.getAlternateAlleles();
			if(alleles.size()!=1 || 
				!dna.matcher(ctx.getReference().getBaseString()).matches() ||
				!dna.matcher(alleles.get(0).getBaseString()).matches()
				)
				{
				w.add(ctx);
				continue;
				}
			StringBuffer ref=new StringBuffer(ctx.getReference().getBaseString().toUpperCase());
			StringBuffer alt=new StringBuffer(alleles.get(0).getBaseString().toUpperCase());
			int start=ctx.getStart();
			int end=ctx.getEnd();
			
			boolean changed=false;
			
			/**** we trim on the right side ****/
			//REF=TGCTGCGGGGGCCGCTGCGGGGG 	ALT=TGCTGCGGGGG
			while(	alt.length()>1 &&
					alt.length() < ref.length() &&
					ref.charAt(ref.length()-1)==alt.charAt(alt.length()-1)
					)
				{
				changed=true;
				ref.setLength(ref.length()-1);
				alt.deleteCharAt(alt.length()-1);
				end--;
				}
			
			//REF=TGCTGCGGGGG 	ALT= TGCTGCGGGGGCCGCTGCGGGGG
			while(	ref.length()>1 &&
					alt.length() > ref.length() &&
					ref.charAt(ref.length()-1)==alt.charAt(alt.length()-1)
					)
				{
				changed=true;
				ref.setLength(ref.length()-1);
				alt.deleteCharAt(alt.length()-1);
				end--;
				}
			
			
			
			/**** we trim on the left side ****/

			//REF=TGCTGCGGGGGCCGCTGCGGGGG 	ALT=TGCTGCGGGGG
			while(	alt.length()>1 &&
					alt.length() < ref.length() &&
					ref.charAt(0)==alt.charAt(0)
					)
				{
				changed=true;
				ref.deleteCharAt(0);
				alt.deleteCharAt(0);
				start++;
				}
			
			//REF=TGCTGCGGGGG 	ALT= TGCTGCGGGGGCCGCTGCGGGGG
			while(	ref.length()>1 &&
					alt.length() > ref.length() &&
					ref.charAt(0)==alt.charAt(0)
					)
				{
				changed=true;
				ref.deleteCharAt(0);
				alt.deleteCharAt(0);
				start++;
				}
			
			
			if(!changed)
				{
				w.add(ctx);
				continue;
				}
			
			
			/*
			LOG.info(line);
			LOG.info("ctx.getStart() "+ctx.getStart());
			LOG.info("ctx.getEnd() "+ ctx.getEnd());

			

			LOG.info("start " + start);
			LOG.info("end "+end);
			LOG.info("ref " + ref.toString());
			LOG.info("alt "+alt.toString());
			*/

			Allele newRef=Allele.create(ref.toString(),true);
			Allele newAlt=Allele.create(alt.toString(),false);
			
			
			Allele newalleles[]=new Allele[]{newRef,newAlt};
		
			b.attribute(TAG, ctx.getReference().getBaseString()+"|"+alleles.get(0).getBaseString()+"|"+ctx.getStart());
			b.start(start);
			b.stop(end);
			b.alleles(Arrays.asList(newalleles));
			
			nChanged++;
			
			
			
			VariantContext ctx2=b.make();
			try {
				w.add(ctx2);
				}
			catch(TribbleException err)
				{
				error(err,"Cannot convert new context:"+ctx2+" old context:"+ctx);
				w.add(ctx);
				}
			}
		
		info("indels changed:"+nChanged);
		}
	

	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+""))!=-1)
			{
			switch(c)
				{
				default:
					{
					switch(handleOtherOptions(c, opt,args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		return super.doWork(opt.getOptInd(), args);
		}
	
	public static void main(String[] args) throws IOException
		{
		new VCFFixIndels().instanceMainWithExit(args);
		}
}
