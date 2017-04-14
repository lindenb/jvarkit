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
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.github.lindenb.jvarkit.util.htsjdk.HtsjdkVersion;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;


/**
 * 
 * VcfCutSamples
 *
 */
@Program(name="vcfcutsamples",
		description="Select/Exclude some samples from a VCF",
		keywords={"vcf"},
		deprecatedMsg="use bcftools or gatk SelectVariants"
		)
public class VcfCutSamples
	extends Launcher
	{
	private static final Logger LOG=Logger.build(VcfCutSamples.class).make();
	
	/** output file */
	@Parameter(names={"-o","--out"},description=" File out. Default: stdout")
	private File outputFile=null;
	/** selected sample */
	@Parameter(names={"-S","--samples"},description="Sample name")
	private Set<String> user_samples=new HashSet<String>();

	@Parameter(names={"-f"},description="read file containing sample names")
	private File sampleFile=null;
	
	@Parameter(names="--invert",description=" invert selection")
	private boolean invert=false;
	/** remove variant if no call at all */
	@Parameter(names="-r",description="remove variant if no call at all")
	private boolean removeCtxIfNoCall=false;
	
	@Parameter(names="-E",description=" a missing sample is an error")
	private boolean missing_sample_is_error=true;
	
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
	protected int doVcfToVcf(String inputName, VcfIterator in,
			VariantContextWriter out) {
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
					LOG.warning(msg);
					}
				}
			}
		
		final List<String> samples2=new ArrayList<String>();

		for(final String sample: header.getSampleNamesInOrder())
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
		
		
		final VCFHeader header2=new VCFHeader(
				header.getMetaDataInInputOrder(),
				samples2
				);
		header2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		header2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
		header2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkVersion",HtsjdkVersion.getVersion()));
		header2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkHome",HtsjdkVersion.getHome()));

		
		out.writeHeader(header2);
		final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header);
		while(in.hasNext())
			{	
			final VariantContext ctx=progress.watch(in.next());
			
			final VariantContextBuilder vb=new VariantContextBuilder(ctx);
			final List<Genotype> genotypes=new ArrayList<Genotype>();
			final Set<Allele> alleles=new HashSet<Allele>();
			boolean only_no_call=true;
			for(final String sample:samples2)
				{
				final Genotype g=ctx.getGenotype(sample);
				if(g.isNoCall()) continue;
				alleles.addAll(g.getAlleles());
				genotypes.add(g);
				if(g.isCalled()) only_no_call=false;
				}
			
			if(removeCtxIfNoCall && only_no_call) continue;
			
			alleles.add(ctx.getReference());
			vb.alleles(alleles);
			vb.genotypes(genotypes);
			out.add(VCFUtils.recalculateAttributes(vb.make()));
			}
		progress.finish();
		return 0;
		}
	
	@Override
	public int doWork(List<String> args) {
		if( this.sampleFile!=null)
			{
			BufferedReader r=null;
			try
				{
				r=IOUtils.openFileForBufferedReading(this.sampleFile);
				String line;
				while((line=r.readLine())!=null)
					{
					if(line.startsWith("#") || line.trim().isEmpty()) continue;
					this.user_samples.add(line);
					}
				}
			catch(Exception err)
				{
				LOG.error(err);
				return -1;
				}
			finally
				{
				CloserUtil.close(r);
				}
			}
		return this.doVcfToVcf(args,outputFile);
		}



	
	public static void main(String[] args)
		{
		new VcfCutSamples().instanceMainWithExit(args);
		}
	}
