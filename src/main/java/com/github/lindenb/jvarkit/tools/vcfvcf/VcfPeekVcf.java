/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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


History:
* 2015 creation

*/
package com.github.lindenb.jvarkit.tools.vcfvcf;

import java.io.IOException;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import com.github.lindenb.jvarkit.util.vcf.IndexedVcfFileReader;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

public class VcfPeekVcf extends AbstractVcfPeekVcf
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(VcfPeekVcf.class);
	private final Set<String> peek_info_tags=new HashSet<String>();
	private IndexedVcfFileReader indexedVcfFileReader=null;
	
	public VcfPeekVcf()
		{
		}
	
	
	@Override
	public Collection<Throwable> initializeKnime() {
		this.indexedVcfFileReader = null;
		try
			{
			this.indexedVcfFileReader = new IndexedVcfFileReader(TABIX);
			for(final String s: super.tagsAsString.split("[, \n]+")) {
				if(s.isEmpty()) continue;
				this.peek_info_tags.add(s);
			}
			} 
		catch(Exception err)
			{
			return wrapException(err);
			}
		return super.initializeKnime();
		}
	
	@Override
	public void disposeKnime() {
		CloserUtil.close(this.indexedVcfFileReader);
		this.indexedVcfFileReader=null;
		this.peek_info_tags.clear();
		super.disposeKnime();
		}
	
	/** public for knime */
	@Override
	public Collection<Throwable> doVcfToVcf(
			final String inputName, 
			final VcfIterator vcfIn,
			final VariantContextWriter out)
			throws IOException
		{
		final VCFHeader h = vcfIn.getHeader();
		final VCFHeader h2 = new VCFHeader(h);
		
		super.addMetaData(h2);
		
		for(final String key: this.peek_info_tags)
			{
			VCFInfoHeaderLine hinfo =this.indexedVcfFileReader.getHeader().getInfoHeaderLine(key);
			if(hinfo==null)
				{
				LOG.warn("INFO name="+key+" missing in "+this.TABIX);
				continue;
				}
			hinfo = VCFUtils.renameVCFInfoHeaderLine(hinfo, this.peekTagPrefix+key);
			if(h2.getInfoHeaderLine(hinfo.getID())!=null)
				{
				throw new IOException("key "+this.peekTagPrefix+key+" already defined in VCF header");
				}
			h2.addMetaDataLine(hinfo);;
			}
		
		out.writeHeader(h2);
		final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(h);
		while(vcfIn.hasNext())
			{
			final VariantContext ctx=progress.watch(vcfIn.next());
						
			final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
			CloseableIterator<VariantContext> iter= this.indexedVcfFileReader.iterator(
					ctx.getContig(),
					Math.max(0,ctx.getStart()-1),
					(ctx.getEnd()+1)
					);
			while(iter.hasNext())
				{
				VariantContext ctx2=iter.next();
				if(!ctx.getContig().equals(ctx2.getContig())) continue;
				if(ctx.getStart()!=ctx2.getStart()) continue;
				if(!ctx.getReference().equals(ctx2.getReference())) continue;
				
				if(this.altAlleleCheck)
					{
					boolean found_all_alt=true;
					for(Allele alt: ctx.getAlternateAlleles())
						{
						if(!ctx2.hasAlternateAllele(alt))
							{
							found_all_alt=false;
							break;
							}
						}
					if(!found_all_alt) continue;
					}
				if(this.peekId && ctx2.hasID())
					{
					vcb.id(ctx2.getID());
					}
				for(final String key: this.peek_info_tags)
					{
					if(!ctx2.hasAttribute(key)) continue;
					final Object o = ctx2.getAttribute(key);
					vcb.attribute(this.peekTagPrefix+key, o);
					}
				}
			iter.close();
			iter=null;
			
			out.add(vcb.make());
				
			if(out.checkError()) break;
			}
		progress.finish();
		return RETURN_OK;
		}
	
	@Override
	protected Collection<Throwable> call(String inputName) throws Exception {
		return doVcfToVcf(inputName);
		}
	
	public static void main(String[] args) throws IOException
		{
		new VcfPeekVcf().instanceMainWithExit(args);
		}
}
