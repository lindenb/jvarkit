package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;


import net.sf.picard.PicardException;

import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.GenotypeBuilder;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter2;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;


public class VcfRenameSamples extends AbstractVCFFilter2
	{
	private Map<String,String> oldNameToNewName=new HashMap<String,String>();
	private boolean missing_user_name_is_error=true;
	
	private VcfRenameSamples()
		{
		}
	

	@Override
	public String getProgramDescription() {
		return "Rename the Samples in a VCF";
		}
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/VcfSampleRename";
		}
	
	

	@Override
	protected void doWork(VcfIterator in, VariantContextWriter out)
			throws IOException
		{
		VCFHeader header1=in.getHeader();
		final Set<String> samples1 = new LinkedHashSet<String>(header1.getSampleNamesInOrder());
		
		final List<String> newHeader=new ArrayList<String>(samples1);
		for(int i=0;i< newHeader.size();++i)
			{
			String destName=this.oldNameToNewName.get(newHeader.get(i));
			if(destName==null) continue;
			newHeader.set(i, destName);
			}
		if(newHeader.size()!= new HashSet<String>(newHeader).size())
			{
			throw new PicardException(
					"Error in input : there are some diplicates in the resulting new VCF header: "+newHeader);
			}
				
		for(String srcName:this.oldNameToNewName.keySet())
			{			
			if(!samples1.contains(srcName))
				{
				if(missing_user_name_is_error)
					{
					throw new IOException("Source Sample "+srcName+" missing in "+samples1+". Use option -E to ignore");
					}
				else
					{
					warning("Missing src-sample:"+srcName);
					}
				}
			}
			
		VCFHeader header2=new VCFHeader(header1.getMetaDataInInputOrder(), newHeader);
		header2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		header2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
		out.writeHeader(header2);
		
		
		while(in.hasNext() )
			{	
			VariantContext ctx=in.next();
			VariantContextBuilder b=new VariantContextBuilder(ctx);
			List<Genotype> genotypes=new ArrayList<Genotype>();
			for(String oldName:samples1)
				{
				Genotype g=ctx.getGenotype(oldName);
				
				String destName=this.oldNameToNewName.get(oldName);
				if(destName!=null)
					{
					GenotypeBuilder gb=new GenotypeBuilder(g);
					gb.name(destName);
					g=gb.make();
					}
				genotypes.add(g);
				}
			b.genotypes(genotypes);
			out.add(b.make());
			}
		
		}
	
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -f (filename). Tab delimited file containing old-name\\tnew-name");
		out.println(" -E ignore error like src sample missing in VCF.");
		super.printOptions(out);
		}
	private void parseNames(File f) throws IOException
		{
		Pattern tab=Pattern.compile("[\t]");
		BufferedReader in=IOUtils.openFileForBufferedReading(f);
		String line;
		while((line=in.readLine())!=null)
			{
			if(line.startsWith("#")) continue;
			if(line.trim().isEmpty()) continue;
			String tokens[]=tab.split(line);
			if(tokens.length<2) throw new IOException("Expected two columns in \""+line+"\" : "+f);
			tokens[0]=tokens[0].trim();
			tokens[1]=tokens[1].trim();
			
			if(tokens[0].isEmpty()) throw new IOException("Empty src-name in \""+line+"\" : "+f);
			if(tokens[1].isEmpty()) throw new IOException("Empty dest-name in \""+line+"\" : "+f);
			if(tokens[0].equals(tokens[1]))
				{
				warning("Same src/dest in in \""+line+"\" : "+f);
				continue;
				}
			if(this.oldNameToNewName.containsKey(tokens[0]))
				{
				String dest=this.oldNameToNewName.get(tokens[0]);
				if(dest.equals(tokens[1]))
					{
					warning(tokens[0]+" -> "+tokens[1]+" defined twice");
					continue;
					}
				else
					{
					throw new IOException("Empty src-name was first defined as "+
							dest+" and now as "+ tokens[1] + " in \""+line+"\" : "+f);
					}
				}
			this.oldNameToNewName.put(tokens[0], tokens[1]);
			}
		in.close();
		}
	@Override
	public int doWork(String[] args)
		{
		
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "f:E"))!=-1)
			{
			switch(c)
				{
				case 'E': missing_user_name_is_error=false;break;
				case 'f':
					{
					try
						{
						parseNames(new File(opt.getOptArg()));
						}
					catch(IOException err)
						{
						error(err);
						return -1;
						}
					break;
					}
				default: 
					{
					switch(handleOtherOptions(c, opt))
						{
						case EXIT_FAILURE:return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		if(this.oldNameToNewName.isEmpty())
			{
			warning("No replacement was defined");
			}
		if(new HashSet<String>(this.oldNameToNewName.values()).size()!=this.oldNameToNewName.size())
			{
			error("Some dest-name have been defined twice.");
			return -1;
			}
		return doWork(opt.getOptInd(), args);
		}

	public static void main(String[] args)
		{
		new VcfRenameSamples().instanceMainWithExit(args);
		}
	}
