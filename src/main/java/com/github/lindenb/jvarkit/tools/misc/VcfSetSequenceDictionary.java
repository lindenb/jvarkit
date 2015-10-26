/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;


import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.util.CloserUtil;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryFactory;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;


public class VcfSetSequenceDictionary extends AbstractVcfSetSequenceDictionary
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(VcfSetSequenceDictionary.class);

	@Override
	public  Command createCommand() {
			return new MyCommand();
		}
		 
	public  static class MyCommand extends AbstractVcfSetSequenceDictionary.AbstractVcfSetSequenceDictionaryCommand
		{		

	private SAMSequenceDictionary dict=null;
	private LinkedHashMap<String, Integer> newdict=null;
	

	@Override
		protected Collection<Throwable> doVcfToVcf(String inputName,
				VcfIterator in, VariantContextWriter out) throws IOException
		{
		VCFHeader header=in.getHeader();
		Set<VCFHeaderLine> meta2=new LinkedHashSet<VCFHeaderLine>();
		for(VCFHeaderLine L:header.getMetaDataInInputOrder())
			{
			if(!L.getKey().equals(VCFHeader.CONTIG_KEY))
				{
				meta2.add(L);
				}
			}
		addMetaData(meta2);
		
		if(dict!=null)
			{	
			meta2.addAll(VCFUtils.samSequenceDictToVCFContigHeaderLine(dict));
			}
		else
			{
			LOG.warn("No sequence dictionary was defined");
			}
		VCFHeader header2=new VCFHeader(meta2, header.getSampleNamesInOrder());

		out.writeHeader(header2);
		
		while(in.hasNext())
			{
			VariantContext ctx=in.next();
			if(dict!=null && dict.getSequenceIndex(ctx.getContig())==-1)
				{
				LOG.warn("Unknown chromosome "+ctx.getContig());
				}
			if(newdict!=null)
				{
				Integer length=this.newdict.get(ctx.getContig());
				if(length==null) length=0;
				if(ctx.getEnd()>length)
					{
					this.newdict.put(ctx.getContig(),ctx.getEnd());
					}
				}
			
			
			out.add(ctx);
			}
		return RETURN_OK;
		}
	
	@Override
	protected Collection<Throwable> call(String inputName) throws Exception
			{
			try {
				if( newDictOut !=null)
					{
					if(!newDictOut.getName().endsWith(".dict"))
						{
						return wrapException("dictionary should end with .dict :"+newDictOut);
						}
					this.newdict=new LinkedHashMap<String, Integer>();	
					}
				if(super.REF!=null)
				{
				return wrapException("undefined REFERENCE");
				}
			
			this.dict=new SAMSequenceDictionaryFactory().load(this.REF);
			
			Collection<Throwable> ret = super.doVcfToVcf(inputName);
			
			if(!ret.isEmpty()) return ret;
			
			if(newDictOut!=null)
				{
				LOG.info("Saving alt dictionary "+newDictOut);
				FileWriter out=null;
				try
					{
					List<SAMSequenceRecord> list=new ArrayList<SAMSequenceRecord>(newdict.size());
					for(String k:this.newdict.keySet())
						{
						list.add(new SAMSequenceRecord(k, this.newdict.get(k)));
						}
					SAMFileHeader sfh=new SAMFileHeader();
					sfh.setSequenceDictionary(
							new SAMSequenceDictionary(list)
							);
			        final SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
			        codec.setValidationStringency(htsjdk.samtools.ValidationStringency.SILENT);
					out=new FileWriter(newDictOut);
					codec.encode(out, sfh);
					out.flush();
					}
				catch(Exception err2)
					{
					return wrapException(err2);
					}
				finally
					{
					CloserUtil.close(out);
					}
				}
			return RETURN_OK;
				} catch (Exception e) {
				return wrapException(e);
				}
			}
	

		}
	
	public static void main(String[] args)
		{
		new VcfSetSequenceDictionary().instanceMainWithExit(args);
		}
	}
