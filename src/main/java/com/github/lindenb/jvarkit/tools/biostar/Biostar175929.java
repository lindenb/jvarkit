/*
The MIT License (MIT)

Copyright (c) 2016 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.biostar;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;


public class Biostar175929 extends AbstractBiostar175929
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(Biostar175929.class);
	private PrintWriter pw;
	
	private void recursive(
			final GenomicSequence chromosome,
			final List<VariantContext> variants,
			int index,
			StringBuilder title,
			StringBuilder sequence
			
			)
		{
		if(index==variants.size())
			{
			int lastPos0= (variants.get(index-1).getEnd()-1);
			for(int i=0;i< this.extendBases && lastPos0+i< chromosome.length();++i)
				{
				sequence.append(Character.toLowerCase(chromosome.charAt(lastPos0+i)));
				}
			pw.print(">");
			pw.print(title);
			for(int i=0;i< sequence.length();++i)
				{
				if(i%60==0) pw.println();
				pw.print(sequence.charAt(i));
				}
			pw.println();
			return;
			}
		if(index==0)
			{
			int firstPos0= (variants.get(0).getStart()-1);
			int chromStart=Math.max(0, firstPos0-extendBases);
			title.append(variants.get(0).getContig()+":"+chromStart);
			for(int i=Math.max(0, firstPos0-extendBases);i< firstPos0 ;++i)
				{
				sequence.append(Character.toLowerCase(chromosome.charAt(i)));
				}
			}
		else
			{
			int endPos0= (variants.get(index-1).getEnd());
			int begPos0= (variants.get(index).getStart()-1);
			while(endPos0< begPos0)
				{
				sequence.append(Character.toLowerCase(chromosome.charAt(endPos0)));
				endPos0++;
				}
			}
		final int  title_length= title.length();
		final int  sequence_length= sequence.length();
		final VariantContext ctx = variants.get(index);
		for(final Allele allele: ctx.getAlleles())
			{
			if(allele.isNoCall()) continue;
			if(allele.isSymbolic())  continue;
			title.setLength(title_length);
			sequence.setLength(sequence_length);
			title.append("|"+ctx.getContig()+":"+ctx.getStart()+"-"+ctx.getEnd()+"("+allele.getBaseString()+")");
			if(super.bracket);
			sequence.append(allele.getBaseString().toUpperCase());
			if(super.bracket);
			recursive(chromosome, variants, index+1, title, sequence);
			}
		
		
		}
	
	
	@Override
	protected Collection<Throwable> call(String inputName) throws Exception {
		if(super.faidx==null) 
			{
			return wrapException("Option -"+OPTION_FAIDX+" was not defined.");
			}
		IndexedFastaSequenceFile reference = null;
		VcfIterator iter=null;
		try {
			reference = new IndexedFastaSequenceFile(super.faidx);
			iter = super.openVcfIterator(inputName);
			this.pw = openFileOrStdoutAsPrintWriter();
			final List<VariantContext> variants = new ArrayList<>();
			for(;;)
				{
				VariantContext ctx = null;
				if(iter.hasNext()) {
					ctx = iter.next();
				}
				
				if( ctx == null || (!variants.isEmpty() && !ctx.getContig().equals(variants.get(0).getContig()))) {
					if(!variants.isEmpty()) 
						{
						LOG.info("chrom:" +variants.get(0).getContig()+ " N="+variants.size());
						final GenomicSequence genomic = new GenomicSequence(reference,variants.get(0).getContig());
						final StringBuilder title= new StringBuilder();
						final StringBuilder sequence= new StringBuilder();
						recursive(genomic,variants,0,title,sequence);
						variants.clear();
						}
					if( ctx == null) break;
					}
				variants.add(ctx);
				}
			iter.close();iter=null;
			this.pw.flush();
			this.pw.close();
			return RETURN_OK;
		} catch (Exception e) {
			return wrapException(e);
		}
			finally {
				CloserUtil.close(reference);
				CloserUtil.close(iter);
				CloserUtil.close(pw);
			}
		}
	
	public static void main(String[] args)throws Exception
		{
		new Biostar175929().instanceMain(args);
		}
	}
