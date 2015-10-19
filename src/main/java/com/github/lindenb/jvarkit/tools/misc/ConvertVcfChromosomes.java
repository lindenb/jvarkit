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

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class ConvertVcfChromosomes extends AbstractConvertVcfChromosomes
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(ConvertBedChromosomes.class);

	 @Override
	public  Command createCommand() {
			return new MyCommand();
		}
		 
	public  class MyCommand extends AbstractConvertVcfChromosomes.AbstractConvertVcfChromosomesCommand
	 	{		
		private Map<String,String> customMapping=new HashMap<String,String>();
		private Set<String> unmappedChromosomes=new HashSet<String>();
	
		
		private String convertName(String chrom)throws IOException
			{
			if(chrom==null) throw new NullPointerException();
			String newname=customMapping.get(chrom);
			if(newname==null)
				{
				if(!unmappedChromosomes.contains(chrom))
					{
					LOG.warn("unmapped chromosome "+chrom);
					unmappedChromosomes.add(chrom);
					}
				if(super.ignore_if_no_mapping) return null;
				
				if(super.use_original_chrom_name_if_no_mapping)
					{	
					return chrom;
					}
				throw new IOException("No mapping found to convert name of chromosome \""+chrom+"\"");
				}
			return newname;
			}
	
	
		private Collection<Throwable> doWork(VcfIterator in, VariantContextWriter out) throws IOException
			{
			VCFHeader header1=in.getHeader();
		
			Set<VCFHeaderLine> meta2=new LinkedHashSet<VCFHeaderLine>();
			for(VCFHeaderLine L:header1.getMetaDataInInputOrder())
				{
				if(!L.getKey().equals(VCFHeader.CONTIG_KEY))
					{
					meta2.add(L);
					}
				}
			VCFHeader header2=new VCFHeader(meta2,header1.getSampleNamesInOrder());
	
			if(header1.getSequenceDictionary()!=null)
				{
				List<SAMSequenceRecord> ssrs=new ArrayList<SAMSequenceRecord>();
	
				for(int i=0;i< header1.getSequenceDictionary().size();++i)
					{
					SAMSequenceRecord ssr=header1.getSequenceDictionary().getSequence(i);
					String newName=convertName(ssr.getSequenceName());
					if(newName==null)
						{
						//skip unknown chromosomes
						continue;
						}
					ssr=new SAMSequenceRecord(newName, ssr.getSequenceLength());
					ssrs.add(ssr);
					}
				header2.setSequenceDictionary(new SAMSequenceDictionary(ssrs));
				}
			addMetaData(header2);
			out.writeHeader(header2);
			
			while(in.hasNext())
				{
				VariantContext ctx=in.next();
				VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				String newName=convertName(ctx.getContig());
				if(newName==null)
					{
					//skip unknown chromosomes
					continue;
					}
				vcb.chr(newName);
				ctx=vcb.make();
	
				out.add(ctx);
				}
			
			if(!unmappedChromosomes.isEmpty())
				{
				LOG.warn("Unmapped chromosomes: "+unmappedChromosomes);
				}
			return RETURN_OK;
			}
				@Override
		protected Collection<Throwable> call(String inputName) throws Exception {
			if(super.mappingFile==null)
			{
			return wrapException("undefined mapping file");
			}
				
			BufferedReader in=null;
			try
				{
				this.customMapping =  super.loadCustomChromosomeMapping(getMappingFile());
				}
			catch(Exception err)
				{
				return wrapException(err);
				}
			finally
				{
				CloserUtil.close(in);
				}

			
			
			VcfIterator vcfin=null;
			VariantContextWriter w = null;
			try {
				
				
				
				vcfin= openVcfIterator(inputName);
				w= openVariantContextWriter();
				return doWork(vcfin, w);
			} catch (Exception e) {
				return wrapException(e);
				}
			finally
				{
				CloserUtil.close(vcfin);
				CloserUtil.close(w);
				}
			}
	 	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new ConvertVcfChromosomes().instanceMainWithExit(args);
		}
	}
