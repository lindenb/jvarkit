package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;


import net.sf.picard.PicardException;
import net.sf.samtools.util.CloserUtil;

import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter2;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;


/**
 * 
 * VcfCutSamples
 *
 */
public class VcfCutSamples extends AbstractVCFFilter2
	{
	private Set<String> user_samples=new HashSet<String>();
	private boolean invert=false;
	private boolean removeCtxIfNoCall=false;
	private boolean missing_sample_is_error=true;
	private VcfCutSamples()
		{
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
	protected void doWork(VcfIterator in, VariantContextWriter out)
			throws IOException
		{
		VCFHeader header=in.getHeader();
		
		final Set<String> samples1=new HashSet<String>(header.getSampleNamesInOrder());
		
		for(String my:this.user_samples)
			{
			if(!samples1.contains(my))
				{
				String msg="user sample "+my+" is not present in VCF Header : "+samples1;
				if(missing_sample_is_error)
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
			if(this.user_samples.contains(sample))
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

		
		out.writeHeader(header2);
		while(in.hasNext())
			{	
			VariantContext ctx=in.next();
			
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
			}
			
		}
	
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -v invert.");
		out.println(" -S (name) add sample.");
		out.println(" -f (file) read file containing sample names. Optional.");
		out.println(" -r remove variant if there is not any called genotype on the line. Optional.");
		out.println(" -E unknown user sample is not a fatal error . Optional.");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "ErvS:f:"))!=-1)
			{
			switch(c)
				{
				case 'E': missing_sample_is_error=false;break;
				case 'v': invert=true; break;
				case 'r': removeCtxIfNoCall=true; break;
				case 'S': user_samples.add(opt.getOptArg()); break;
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
		
		return doWork(opt.getOptInd(), args);
		}

	public static void main(String[] args)
		{
		new VcfCutSamples().instanceMainWithExit(args);
		}
	}
