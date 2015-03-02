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
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.htsjdk.HtsjdkVersion;
import com.github.lindenb.jvarkit.util.picard.PicardException;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.knime.KnimeApplication;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;


/**
 * 
 * VcfCutSamples
 *
 */
public class VcfCutSamples
	extends AbstractCommandLineProgram
	implements KnimeApplication
	{
	/** output file */
	private File outputFile=null;
	/** selected sample */
	private Set<String> user_samples=new HashSet<String>();
	/** invert selection */
	private boolean invert=false;
	/** remove variant if no call at all */
	private boolean removeCtxIfNoCall=false;
	/** a missing sample is an error */
	private boolean missing_sample_is_error=true;
	/** count filtered variants */
	private int countFilteredVariants=0;
	
	
	public VcfCutSamples()
		{
		}
	
	public Set<String> getUserSamples() {
		return user_samples;
		}
	
	public void setInvert(boolean invert)
		{
		this.invert = invert;
		}
	
	public void setMissingSampleIsError(boolean missing_sample_is_error) {
		this.missing_sample_is_error = missing_sample_is_error;
	}

	@Override
	public String getProgramDescription() {
		return "Select/Exclude some samples from a VCF";
		}
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/VcfCutSamples";
		}
	
	
	@Override
	public void checkKnimeCancelled() {
		
		}
	
	@Override
	public void disposeKnime() {
		
		}
	
	@Override
	public int initializeKnime() {
		return 0;
		}
	
	public int getVariantCount()
		{
		return this.countFilteredVariants;
		}

	@Override
	public void setOutputFile(File out) {
		this.outputFile=out;
		}
	
	public File getOutputFile() {
		return outputFile;
		}
	
	public void setRemoveCtxIfNoCall(boolean removeCtxIfNoCall) {
		this.removeCtxIfNoCall = removeCtxIfNoCall;
		}

	
	private void filterVcfIterator(VcfIterator in, VariantContextWriter out)
			throws IOException
		{
		VCFHeader header=in.getHeader();
		this.countFilteredVariants=0;
		final Set<String> samples1=new HashSet<String>(header.getSampleNamesInOrder());
		
		for(String my:this.getUserSamples())
			{
			if(!samples1.contains(my))
				{
				String msg="user sample "+my+" is not present in VCF Header : "+samples1;
				if(this.missing_sample_is_error)
					{
					throw new PicardException(msg);
					}
				else
					{
					warning(msg);
					}
				}
			}
		
		List<String> samples2=new ArrayList<String>();

		for(String sample: header.getSampleNamesInOrder())
			{
			if(this.getUserSamples().contains(sample))
				{
				if(!invert)
					{
					samples2.add(sample);
					}
				}
			else
				{
				if(invert)
					{
					samples2.add(sample);
					}
				}
			}
		
		
		VCFHeader header2=new VCFHeader(
				header.getMetaDataInInputOrder(),
				samples2
				);
		header2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		header2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
		header2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkVersion",HtsjdkVersion.getVersion()));
		header2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkHome",HtsjdkVersion.getHome()));

		
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
			++countFilteredVariants;
			}
		progress.finish();
		}
	
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -v invert.");
		out.println(" -S (name) add sample.");
		out.println(" -f (file) read file containing sample names. Optional.");
		out.println(" -r remove variant if there is not any called genotype on the line. Optional.");
		out.println(" -E unknown user sample is not a fatal error . Optional.");
		out.println(" -o (file) File out. Default: stdout.");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "ErvS:f:o:"))!=-1)
			{
			switch(c)
				{
				case 'E': setMissingSampleIsError(false);break;
				case 'v': setInvert(true); break;
				case 'r': setRemoveCtxIfNoCall(true); break;
				case 'S': getUserSamples().add(opt.getOptArg()); break;
				case 'o': setOutputFile(new File(opt.getOptArg())); break;
				case 'f':
					BufferedReader r=null;
					try
						{
						r=IOUtils.openURIForBufferedReading(opt.getOptArg());
						String line;
						while((line=r.readLine())!=null)
							{
							if(line.startsWith("#") || line.trim().isEmpty()) continue;
							this.user_samples.add(line);
							}
						}
					catch(Exception err)
						{
						error(err);
						return -1;
						}
					finally
						{
						CloserUtil.close(r);
						}
					break;
				default: 
					{
					switch(handleOtherOptions(c, opt, null))
						{
						case EXIT_FAILURE:return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		
		this.initializeKnime();
		List<String> L=new ArrayList<String>();
		for(int i=opt.getOptInd();i<args.length;++i)
			{
			L.add(args[i]);
			}
		return this.executeKnime(L);
		}

	@Override
	public int executeKnime(List<String> args)
		{
		VariantContextWriter vcw=null;
		VcfIterator vcfIn=null;
		try
			{
			if(args.isEmpty())
				{
				vcfIn = VCFUtils.createVcfIteratorStdin();
				}
			else if(args.size()==1)
				{
				vcfIn= VCFUtils.createVcfIterator(args.get(0));
				}
			else
				{
				error(getMessageBundle("illegal.number.of.arguments"));
				return -1;
				}
			if(getOutputFile()==null)
				{
				vcw = VCFUtils.createVariantContextWriterToStdout();
				}
			else
				{
				vcw = VCFUtils.createVariantContextWriter(getOutputFile());
				}
			this.filterVcfIterator(vcfIn,vcw);
			vcw.close(); vcw=null;
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(vcfIn);
			CloserUtil.close(vcw);
			}
		
		}

	
	public static void main(String[] args)
		{
		new VcfCutSamples().instanceMainWithExit(args);
		}
	}
