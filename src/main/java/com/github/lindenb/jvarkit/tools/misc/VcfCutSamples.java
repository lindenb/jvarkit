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
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;


/**
 * 
 * VcfCutSamples
 *
 */
public class VcfCutSamples
	extends AbstractVcfCutSamples
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(VcfBurden.class);

	@Override
	public  Command createCommand() {
			return new MyCommand();
		}
		 
	public  static class MyCommand extends AbstractVcfCutSamples.AbstractVcfCutSamplesCommand
		{		
		@Override
		protected Collection<Throwable> doVcfToVcf(
					String inputName,
					VcfIterator in,
					VariantContextWriter out
					) throws IOException
			{
			VCFHeader header=in.getHeader();
			final Set<String> samples1=new HashSet<String>(header.getSampleNamesInOrder());
			
			for(String my:this.getUserSamples())
				{
				if(!samples1.contains(my))
					{
					String msg="user sample "+my+" is not present in VCF Header : "+samples1;
					if(this.missing_sample_is_error)
						{
						throw new RuntimeException(msg);
						}
					else
						{
						LOG.warn(msg);
						}
					}
				}
			
			List<String> samples2=new ArrayList<String>();
	
			for(String sample: header.getSampleNamesInOrder())
				{
				if(this.getUserSamples().contains(sample))
					{
					if(!super.invert)
						{
						samples2.add(sample);
						}
					}
				else
					{
					if(super.invert)
						{
						samples2.add(sample);
						}
					}
				}
			
			
			VCFHeader header2=new VCFHeader(
					header.getMetaDataInInputOrder(),
					samples2
					);
			addMetaData(header2);
			
			out.writeHeader(header2);
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header);
			while(in.hasNext())
				{	
				VariantContext ctx=progress.watch(in.next());
				
				VariantContextBuilder vb=new VariantContextBuilder(ctx);
				List<Genotype> genotypes=new ArrayList<Genotype>();
				Set<Allele> alleles=new HashSet<Allele>();
				boolean only_no_call=true;
				for(String sample:samples2)
					{
					Genotype g=ctx.getGenotype(sample);
					if(g.isNoCall()) continue;
					alleles.addAll(g.getAlleles());
					genotypes.add(g);
					if(g.isCalled()) only_no_call=false;
					}
				
				if(removeCtxIfNoCall && only_no_call) continue;
				
				alleles.add(ctx.getReference());
				vb.alleles(alleles);
				vb.genotypes(genotypes);
				out.add(vb.make());
				if(out.checkError()) break;
				}
			progress.finish();
			return RETURN_OK;
			}
		
	@Override
	public Collection<Throwable> initializeKnime()
		{
		if(super.sampleFile!=null)
			{
			BufferedReader r=null;
			try
				{
				LOG.info("Reading "+super.sampleFile);
				r=IOUtils.openFileForBufferedReading(super.sampleFile);
				String line;
				while((line=r.readLine())!=null)
					{
					if(line.startsWith("#") || line.trim().isEmpty()) continue;
					this.userSamples.add(line);
					}
				}
			catch(Exception err)
				{
				return wrapException(err);
				}
			finally
				{
				CloserUtil.close(r);
				}
			}
		return super.initializeKnime();
		}
	@Override
	public void disposeKnime() {
		this.userSamples.clear();
		super.disposeKnime();
		}
	
	@Override
	protected Collection<Throwable> call(String inputName) throws Exception
		{
		return doVcfToVcf(inputName);
		}
	
	}
	
	public static void main(String[] args)
		{
		new VcfCutSamples().instanceMainWithExit(args);
		}
	}
