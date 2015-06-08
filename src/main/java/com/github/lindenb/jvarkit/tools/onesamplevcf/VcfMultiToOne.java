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

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

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
import com.github.lindenb.jvarkit.knime.AbstractKnimeApplication;
import com.github.lindenb.jvarkit.util.htsjdk.HtsjdkVersion;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

/*
 * VcfMultiToOne
 */
public class VcfMultiToOne extends AbstractKnimeApplication
	{
	private boolean keep_no_call=true;
	private boolean keep_hom_ref=true;
	private boolean keep_non_available=true;
	public static final String DEFAULT_VCF_SAMPLE_NAME="SAMPLE";
	public static final String DEFAULT_SAMPLE_TAGID="SAMPLENAME";
	public static final String DEFAULT_SAMPLE_FILETAGID="SAMPLESOURCE";
	public static final String SAMPLE_HEADER_DECLARATION="VcfMultiToOne.Sample";
	

	private int countFilteredVariants=0;
	
	/** return the number of variants in the output vcf */
	public int getVariantCount()
		{
		return this.countFilteredVariants;
		}

	
	public VcfMultiToOne()
		{
		}
	
	/** general utility for program using VCFMulti2One:
	 *  Extract SampleNames
	 */
	static Set<String> extractSampleNames(final VCFHeader header)
		{
		List<String> sample_list =header.getSampleNamesInOrder();
		if(sample_list.size()!=1 || !sample_list.get(0).equals(DEFAULT_VCF_SAMPLE_NAME))
			{
			throw new IllegalArgumentException("Not a VCF produced by VcfMultiToOne");
			}
		Set<String> samples = new TreeSet<String>();
		for(VCFHeaderLine h:header.getMetaDataInInputOrder())
			{
			if(h.getKey().equals(SAMPLE_HEADER_DECLARATION))
				{
				sample_list.add(h.getValue());
				}
			}
		return samples;
		}
	
	
	public void setKeepHomRef(boolean keep_hom_ref) {
		this.keep_hom_ref = keep_hom_ref;
	}
	public void setKeepNoCall(boolean keep_no_call) {
		this.keep_no_call = keep_no_call;
	}
	public void setKeepNonAvailable(boolean keep_non_available)
		{
		this.keep_non_available = keep_non_available;
		}
	
	
	@Override
	public String getProgramDescription() {
		return "Convert VCF with multiple samples to a VCF with one SAMPLE, duplicating variant and adding the sample name in the INFO column";
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return DEFAULT_WIKI_PREFIX+"VcfMultiToOne";
		}
	
	/** open VariantContextWriter */
	protected VariantContextWriter createVariantContextWriter()
		throws IOException
		{
		if(getOutputFile()==null)
			{
			return VCFUtils.createVariantContextWriterToStdout();
			}
		else
			{
			info("opening vcf writer to "+getOutputFile());
			return VCFUtils.createVariantContextWriter(getOutputFile());
			}

		}
	
	
	@Override
	public int executeKnime(List<String> arguments)
		{
		VariantContextWriter  out=null;
		Set<String> args= IOUtils.unrollFiles(arguments);
		List<VcfIterator> inputs=new ArrayList<>(args.size()+1);
		List<String> inputFiles=new ArrayList<>(args.size()+1);
		
		
		try
			{
			if(args.isEmpty() && arguments.isEmpty())
				{
				inputs.add(VCFUtils.createVcfIteratorStdin());
				inputFiles.add("stdin");
				}
			else if(args.isEmpty())
				{
				error(getMessageBundle("illegal.number.of.arguments"));
				return -1;
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
					error(getMessageBundle("no.dict.in.vcf"));
					return -1;
					}
				else if(!SequenceUtil.areSequenceDictionariesEqual(dict, header.getSequenceDictionary()))
					{
					error(getMessageBundle("not.the.same.sequence.dictionaries"));
					return -1;
					}
				metaData.addAll(in.getHeader().getMetaDataInInputOrder());
				sampleNames.addAll(in.getHeader().getSampleNamesInOrder());
				}
			
			
			final Comparator<VariantContext> comparator = (dict==null?
					VCFUtils.createChromPosRefComparator():
					VCFUtils.createTidPosRefComparator(dict)
					);
			metaData.add(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
			metaData.add(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
			metaData.add(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkVersion",HtsjdkVersion.getVersion()));
			metaData.add(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkHome",HtsjdkVersion.getHome()));
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
			
			out= createVariantContextWriter();
			out.writeHeader(h2);
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dict);
			for(;;)
				{
				if(checkOutputError()) break;
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
						++countFilteredVariants;
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
					++countFilteredVariants;
					}
				}
			progress.finish();
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(inputs);
			CloserUtil.close(out);
			}
		return 0;
		}
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -c discard if variant is no-call");
		out.println(" -r discard if variant is hom-ref");
		out.println(" -a discard if variant is non-available");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "crao:"))!=-1)
			{
			switch(c)
				{
				case 'o': this.setOutputFile(args[opt.getOptInd()]);
				case 'c': this.setKeepNoCall(false); break;
				case 'r': this.setKeepHomRef(false); break;
				case 'a': this.setKeepNonAvailable(false); break;
				default: 
					{
					switch(handleOtherOptions(c, opt, args))
						{
						case EXIT_FAILURE:return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		
		return mainWork(opt.getOptInd(), args);
		}

	public static void main(String[] args)
		{
		new VcfMultiToOne().instanceMainWithExit(args);
		}
	}
