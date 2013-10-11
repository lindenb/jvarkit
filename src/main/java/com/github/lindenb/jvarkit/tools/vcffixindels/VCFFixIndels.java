package com.github.lindenb.jvarkit.tools.vcffixindels;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Pattern;

import net.sf.picard.cmdline.Usage;
import net.sf.picard.util.Log;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import org.broad.tribble.TribbleException;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter;


public class VCFFixIndels extends AbstractVCFFilter
	{
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+" Fix samtools indels (for @SolenaLS).";
	
	private static final Log LOG=Log.getInstance(VCFFixIndels.class);
	
	
	
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
				LOG.error(err,"Cannot convert new context:"+ctx2+" old context:"+ctx);
				w.add(ctx);
				}
			}
		
		LOG.info("indels changed:"+nChanged);
		}
	
	public static void main(String[] args) throws IOException
		{
		new VCFFixIndels().instanceMainWithExit(args);
		}
}
