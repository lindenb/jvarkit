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
package com.github.lindenb.jvarkit.tools.onesamplevcf;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

/*
 * VcfMultiToOne
 */
public class VcfMultiToOne extends AbstractVcfMultiToOne
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(VcfMultiToOne.class);

	@Override
	public Command createCommand()
		{
		return new MyCommand();
		}
	
	static private class MyCommand extends AbstractVcfMultiToOne.AbstractVcfMultiToOneCommand
		{

	
	public static final String DEFAULT_VCF_SAMPLE_NAME="SAMPLE";
	public static final String DEFAULT_SAMPLE_TAGID="SAMPLENAME";
	public static final String DEFAULT_SAMPLE_FILETAGID="SAMPLESOURCE";
	public static final String SAMPLE_HEADER_DECLARATION="VcfMultiToOne.Sample";
	

	
	
	
	
	
	
	
	@Override
	protected Collection<Throwable> call(String inputName) throws Exception
		{
		List<String> arguments = getInputFiles();
		VariantContextWriter  out=null;
		Set<String> args= IOUtils.unrollFiles(arguments);
		List<VcfIterator> inputs=new ArrayList<>(args.size()+1);
		List<String> inputFiles=new ArrayList<>(args.size()+1);
		
		
		try
			{
			if(args.isEmpty() && arguments.isEmpty())
				{
				inputs.add(VCFUtils.createVcfIteratorFromStream(stdin()));
				inputFiles.add("stdin");
				}
			else if(args.isEmpty())
				{
				return wrapException(getMessageBundle("illegal.number.of.arguments"));
				}
			else
				{
				for(String vcfFile: args)
					{
					inputs.add(VCFUtils.createVcfIterator(vcfFile));
					inputFiles.add(VCFUtils.escapeInfoField(vcfFile));
					}
				}
			SAMSequenceDictionary dict=null;
			Set<String> sampleNames=new HashSet<String>();

			Set<VCFHeaderLine> metaData = new HashSet<VCFHeaderLine>();
			for(VcfIterator in:inputs)
				{
				VCFHeader header = in.getHeader();
				if(dict==null)
					{
					dict = header.getSequenceDictionary();
					}
				else if(header.getSequenceDictionary()==null)
					{
					return wrapException(getMessageBundle("no.dict.in.vcf"));
					}
				else if(!SequenceUtil.areSequenceDictionariesEqual(dict, header.getSequenceDictionary()))
					{
					return wrapException(getMessageBundle("not.the.same.sequence.dictionaries"));
					}
				metaData.addAll(in.getHeader().getMetaDataInInputOrder());
				sampleNames.addAll(in.getHeader().getSampleNamesInOrder());
				}
			
			
			final Comparator<VariantContext> comparator = (dict==null?
					VCFUtils.createChromPosRefComparator():
					VCFUtils.createTidPosRefComparator(dict)
					);
			
			addMetaData(metaData);
			metaData.add(new VCFInfoHeaderLine(
					DEFAULT_SAMPLE_TAGID,1,VCFHeaderLineType.String,
					"Sample Name from multi-sample vcf"
					));
			metaData.add(new VCFInfoHeaderLine(
					DEFAULT_SAMPLE_FILETAGID,1,VCFHeaderLineType.String,
					"Origin of sample"
					));
			
			for(String sample:sampleNames)
				{
				metaData.add(
					new VCFHeaderLine(
					SAMPLE_HEADER_DECLARATION,
					sample));
				}
			
			VCFHeader h2 = new VCFHeader(
					metaData,
					Collections.singleton(DEFAULT_VCF_SAMPLE_NAME)
					);
			
			out= openVariantContextWriter();
			out.writeHeader(h2);
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dict);
			for(;;)
				{
				if(out.checkError()) break;
				/* get 'smallest' variant */
				VariantContext smallest = null;
				int idx=0;
				int best_idx = -1;
				
				while(idx < inputs.size())
					{
					VcfIterator in= inputs.get(idx);
					if(!in.hasNext())
						{
						CloserUtil.close(in);
						inputs.remove(idx);
						inputFiles.remove(idx);
						}
					else
						{
						VariantContext ctx = in.peek();
						if( smallest==null ||
							comparator.compare(smallest,ctx)>0)
							{
							smallest = ctx;
							best_idx = idx;
							}
						++idx;
						}
					}
				
				if(smallest==null) break;
				
				VariantContext ctx = progress.watch(inputs.get(best_idx).next());
				
				
				if(ctx.getNSamples()==0)
					{
					if(keep_no_call)
						{
						VariantContextBuilder vcb = new VariantContextBuilder(ctx);
						vcb.attribute(DEFAULT_SAMPLE_FILETAGID,inputFiles.get(best_idx));
						vcb.genotypes(GenotypeBuilder.createMissing(DEFAULT_VCF_SAMPLE_NAME,2));
						out.add(vcb.make());
						}
					continue;
					}
				
				for(int i=0;i< ctx.getNSamples();++i)
					{
					Genotype g= ctx.getGenotype(i);
					String sample = g.getSampleName();
					
					if(!g.isCalled() && !keep_no_call) continue;
					if(!g.isAvailable() && !keep_non_available) continue;
					if(g.isHomRef() && !keep_hom_ref) continue;
					
					
					GenotypeBuilder gb=new GenotypeBuilder(g);
					gb.name(DEFAULT_VCF_SAMPLE_NAME);
					
					
					VariantContextBuilder vcb=new VariantContextBuilder(ctx);
					vcb.attribute(DEFAULT_SAMPLE_TAGID, sample);
					vcb.attribute(DEFAULT_SAMPLE_FILETAGID,inputFiles.get(best_idx));
					
					
					vcb.genotypes(gb.make());
					out.add(vcb.make());
					}
				}
			progress.finish();
			return RETURN_OK;
			}
		catch(Exception err)
			{
			return wrapException(err);
			}
		finally
			{
			CloserUtil.close(inputs);
			CloserUtil.close(out);
			}
		}
	
		
			
	
		}
	
	public static void main(String[] args)
		{
		new VcfMultiToOne().instanceMainWithExit(args);
		}
	}
