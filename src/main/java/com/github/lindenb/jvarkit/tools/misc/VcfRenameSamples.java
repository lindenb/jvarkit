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
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

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
import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;


public class VcfRenameSamples extends AbstractVcfRenameSamples
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(VcfRenameSamples.class);

	@Override
	public  Command createCommand() {
			return new MyCommand();
		}
		 
	public  static class MyCommand extends AbstractVcfRenameSamples.AbstractVcfRenameSamplesCommand
		{		

	
	private Map<String,String> oldNameToNewName=new HashMap<String,String>();
	private boolean missing_user_name_is_error=true;
	
	@Override
		protected Collection<Throwable> doVcfToVcf(String inputName,
				VcfIterator in, VariantContextWriter out) throws IOException
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
			return wrapException( "Error in input : there are some diplicates in the resulting new VCF header: "+newHeader);
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
					return wrapException("Missing src-sample:"+srcName);
					}
				}
			}
			
		VCFHeader header2=new VCFHeader(header1.getMetaDataInInputOrder(), newHeader);
		addMetaData(header2);
		out.writeHeader(header2);
		
		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header1);
		while(in.hasNext() )
			{	
			VariantContext ctx=progress.watch(in.next());
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
			if(out.checkError()) break;
			}
		progress.finish();
		return RETURN_OK;
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
				LOG.warn("Same src/dest in in \""+line+"\" : "+f);
				continue;
				}
			if(this.oldNameToNewName.containsKey(tokens[0]))
				{
				String dest=this.oldNameToNewName.get(tokens[0]);
				if(dest.equals(tokens[1]))
					{
					LOG.warn(tokens[0]+" -> "+tokens[1]+" defined twice");
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
	protected Collection<Throwable> call(String inputName) throws Exception {
		return doVcfToVcf(inputName);
		}
	
	@Override
		public Collection<Throwable> initializeKnime()
		{
		if(assocFile==null)
			{
			return wrapException("undefined association file");
			}
		
		try {
			parseNames(this.assocFile);
			} 
		catch (Exception e) {
			return wrapException(e);
			}
		if(this.oldNameToNewName.isEmpty())
			{
			LOG.warn("No replacement was defined");
			}
		if(new HashSet<String>(this.oldNameToNewName.values()).size()!=this.oldNameToNewName.size())
			{
			return wrapException("Some dest-name have been defined twice.");
			}
		return super.initializeKnime();
		}
	
	}
	public static void main(String[] args)
		{
		new VcfRenameSamples().instanceMainWithExit(args);
		}
	}
