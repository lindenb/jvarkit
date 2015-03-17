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
package com.github.lindenb.jvarkit.tools.vcfconcat;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.knime.AbstractKnimeApplication;
import com.github.lindenb.jvarkit.util.cli.GetOpt;
import com.github.lindenb.jvarkit.util.htsjdk.HtsjdkVersion;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class VcfConcat extends AbstractKnimeApplication
	{
	public static final String VARIANTSOURCE="VARIANTSOURCE";
	private Set<String> inputFiles=new HashSet<String>();
	private File outputfile=null;
	/** number of variants filtered */
	private int countFilteredVariants=0;
	
	public VcfConcat()
		{
		}
	
	@Override
	public void setOutputFile(File out)
		{
		this.outputfile=out;
		}
	public File getOutputFile()
		{
		return outputfile;
		}
	
	public Set<String> getInputFiles()
		{
		return inputFiles;
		}

	@Override
	protected String getOnlineDocUrl()
		{
		return DEFAULT_WIKI_PREFIX+"VcfConcat";
		}
	
	@Override
	public String getProgramDescription()
		{
		return "Concatenante sorted VCF with same sample, does NOT merge genotypes";
		}
	
	/** return the number of variants in the output vcf */
	public int getVariantCount()
		{
		return this.countFilteredVariants;
		}

	
	private int fromFiles(VariantContextWriter out) throws IOException
		{
		List<VcfIterator> inputs=new ArrayList<VcfIterator>(getInputFiles().size());
		List<String> inputFiles=new ArrayList<>(getInputFiles().size());
		List<String> samples=new ArrayList<>();
		SAMSequenceDictionary dict=null;
		this.countFilteredVariants=0;
		try
			{
			Set<VCFHeaderLine> metaData = new HashSet<VCFHeaderLine>();
			
			/* open each vcf file */
			for(String vcfFile:this.getInputFiles())
				{
				info("Opening "+vcfFile);
				VcfIterator r=VCFUtils.createVcfIterator(vcfFile);
				
				/* check VCF dict */
				VCFHeader header = r.getHeader();
				if(dict==null && inputs.isEmpty())
					{
					dict = header.getSequenceDictionary();
					}
				else if(!inputs.isEmpty() &&
					(
					(dict==null && header.getSequenceDictionary()!=null) ||
					(dict!=null && header.getSequenceDictionary()==null))
					)
					{
					error(getMessageBundle("not.the.same.sequence.dictionaries"));
					return -1;
					}
				else if(!inputs.isEmpty() && dict!=null && 
					!SequenceUtil.areSequenceDictionariesEqual(dict, header.getSequenceDictionary()))
					{
					error(getMessageBundle("not.the.same.sequence.dictionaries"));
					return -1;
					}
				/* check samples */
				if(inputs.isEmpty())
					{
					samples = header.getSampleNamesInOrder();
					}
				else if(!header.getSampleNamesInOrder().equals(samples))
					{
					error("No same samples");
					return -1;
					}
				
				metaData.addAll(header.getMetaDataInInputOrder());
				inputs.add(r);
				inputFiles.add(VCFUtils.escapeInfoField(vcfFile));
				}
			/* create comparator according to dict*/
			final Comparator<VariantContext> comparator = (dict==null?
					VCFUtils.createChromPosRefComparator():
					VCFUtils.createTidPosRefComparator(dict)
					);
			metaData.add(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
			metaData.add(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
			metaData.add(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkVersion",HtsjdkVersion.getVersion()));
			metaData.add(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkHome",HtsjdkVersion.getHome()));
			metaData.add(new VCFInfoHeaderLine(
					VARIANTSOURCE,1,VCFHeaderLineType.String,
					"Origin File of Varant"
					));
			VCFHeader h2 = new VCFHeader(
					metaData,
					samples
					);
			out.writeHeader(h2);
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dict);
			for(;;)
				{
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
				
				if(smallest==null || checkOutputError()) break;
				
				VariantContext ctx = progress.watch(inputs.get(best_idx).next());
				VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				vcb.attribute(VARIANTSOURCE, inputFiles.get(best_idx));
				out.add(vcb.make());
				++countFilteredVariants;
				}
			progress.finish();
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(inputs);
			}

		}
	
	@Override
	public int doWork(String[] args)
		{
		GetOpt getopt=new GetOpt();
		int c;
		while((c=getopt.getopt(args, super.getGetOptDefault()+""))!=-1)
			{
			switch(c)
				{
				default: switch(super.handleOtherOptions(c, getopt, args))
					{
					case EXIT_FAILURE: return -1;
					case EXIT_SUCCESS: return 0;
					default:break;
					}
				}
			}
				
		return mainWork(getopt.getOptInd(), args);
		}

	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new VcfConcat().instanceMainWithExit(args);
		}

	@Override
	public int executeKnime(List<String> args)
		{
		VariantContextWriter w=null;
		BufferedReader r=null;
		try
			{
			if(args.isEmpty())
				{
				info("Reading filenames from stdin");
				r= IOUtils.openStdinForBufferedReader();
				String line;
				while((line=r.readLine())!=null)
					{
					if(line.trim().isEmpty()|| line.startsWith("#"))
						continue;
					getInputFiles().add(line);
					}
				r.close();
				r=null;
				}
			else
				{
				for(String filename:args)
					{
					getInputFiles().addAll(IOUtils.unrollFiles(Arrays.asList(filename)));	
					}
				}
			if(getInputFiles().isEmpty())
				{
				error("No input");
				return -1;
				}
			if(getOutputFile()==null)
				{
				w= VCFUtils.createVariantContextWriterToStdout();
				}
			else
				{
				info("opening vcf writer to "+getOutputFile());
				w= VCFUtils.createVariantContextWriter(getOutputFile());
				}
			return fromFiles(w);
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(w);
			CloserUtil.close(r);
			}
		}

	}
