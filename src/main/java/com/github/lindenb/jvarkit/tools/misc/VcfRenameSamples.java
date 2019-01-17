/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.util.vcf.PostponedVariantContextWriter;
import htsjdk.variant.vcf.VCFIterator;


import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

/**
BEGIN_DOC

```bash
$ curl -s "https://raw.github.com/arq5x/gemini/master/test/test1.snpeff.vcf" |\
grep -A 2 CHROM  | cut -f 1-5,10-

#CHROM	POS	ID	REF	ALT	1094PC0005	1094PC0009	1094PC0012	1094PC0013
chr1	30860	.	G	C	0/0:7,0:7:15.04:0,15,177	0/0:2,0:2:3.01:0,3,39	0/0:6,0:6:12.02:0,12,143	0/0:4,0:4:9.03:0,9,119
chr1	69270	.	A	G	./.	./.	1/1:0,3:3:9.03:106,9,0	1/1:0,6:6:18.05:203,18,0

curl -s "https://raw.github.com/arq5x/gemini/master/test/test1.snpeff.vcf" |\
java -jar dist/vcfrenamesamples.jar -f <(echo -e "1094PC0005\tALPHA\n1094PC0012\tBETA\nSAMPLE\tGAMMA\n1094PC0009\tEPSILON") -E  |\
grep -A 2 "#CHROM" | cut -f 1-5,10-

#CHROM	POS	ID	REF	ALT	ALPHA	EPSILON	BETA	1094PC0013
chr1	30860	.	G	C	0/0:7,0:7:15:0,15,177	0/0:2,0:2:3:0,3,39	0/0:6,0:6:12:0,12,143	0/0:4,0:4:9:0,9,119
chr1	69270	.	A	G	./.	./.	1/1:0,3:3:9:106,9,0	1/1:0,6:6:18:203,18,0
```

END_DOC
*/

@Program(name="vcfrenamesamples",description="Rename the Samples in a VCF",
	keywords={"vcf","sample"}
	)
public class VcfRenameSamples extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfRenameSamples.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names="-f",description="Tab delimited file containing old-name\\tnew-name",required=true)
	private File mappingFile=null;
	@Parameter(names="-E",description= "error like src sample missing in VCF")
	private boolean missing_user_name_is_error=false;
	@ParametersDelegate
	private PostponedVariantContextWriter.WritingVcfConfig writingVcfArgs = new PostponedVariantContextWriter.WritingVcfConfig();

	private final Map<String,String> oldNameToNewName=new HashMap<>();
	
	public VcfRenameSamples()
		{
		}
	
	@Override
	protected int doVcfToVcf(
			final String inputName,
			final VCFIterator in,
			final VariantContextWriter out
			) {
		final VCFHeader header1=in.getHeader();
		final Set<String> samples1 = new LinkedHashSet<String>(header1.getSampleNamesInOrder());
		
		final List<String> newHeader=new ArrayList<String>(samples1);
		for(int i=0;i< newHeader.size();++i)
			{
			final String destName=this.oldNameToNewName.get(newHeader.get(i));
			if(destName==null) continue;
			newHeader.set(i, destName);
			}
		if(newHeader.size()!= new HashSet<String>(newHeader).size())
			{
			throw new RuntimeException(
					"Error in input : there are some diplicates in the resulting new VCF header: "+newHeader);
			}
				
		for(final String srcName:this.oldNameToNewName.keySet())
			{
			if(!samples1.contains(srcName))
				{
				if(missing_user_name_is_error)
					{
					LOG.error("Source Sample "+srcName+" missing in "+samples1+". Use option -E to ignore");
					return -1;
					}
				else
					{
					LOG.warning("Missing src-sample:"+srcName);
					}
				}
			}
			
		final VCFHeader header2=new VCFHeader(header1.getMetaDataInInputOrder(), newHeader);
		header2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		header2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
		JVarkitVersion.getInstance().addMetaData(getClass().getSimpleName(), header2);
		
		out.writeHeader(header2);
		
		final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header1);
		while(in.hasNext() )
			{	
			final VariantContext ctx=progress.watch(in.next());
			final VariantContextBuilder b=new VariantContextBuilder(ctx);
			final List<Genotype> genotypes=new ArrayList<Genotype>();
			for(final String oldName:samples1)
				{
				Genotype g=  ctx.getGenotype(oldName);
				
				final String destName=this.oldNameToNewName.get(oldName);
				if(destName!=null)
					{
					final GenotypeBuilder gb=new GenotypeBuilder(g);
					gb.name(destName);
					g=gb.make();
					}
				genotypes.add(g);
				}
			b.genotypes(genotypes);
			out.add(b.make());
			}
		progress.finish();
		return 0;
		}
	
	private void parseNames(final File f) throws IOException
		{
		final CharSplitter tab= CharSplitter.TAB;
		final BufferedReader in=IOUtils.openFileForBufferedReading(f);
		String line;
		while((line=in.readLine())!=null)
			{
			if(line.startsWith("#")) continue;
			if(StringUtil.isBlank(line)) continue;
			final String tokens[]=tab.split(line);
			if(tokens.length!=2) throw new IOException("Expected two columns in \""+line+"\" : "+f);
			tokens[0]=tokens[0].trim();
			tokens[1]=tokens[1].trim();
			
			if(StringUtil.isBlank(tokens[0])) throw new IOException("Empty src-name in \""+line+"\" : "+f);
			if(StringUtil.isBlank(tokens[1])) throw new IOException("Empty dest-name in \""+line+"\" : "+f);
			if(tokens[0].equals(tokens[1]))
				{
				LOG.warning("Same src/dest in in \""+line+"\" : "+f);
				continue;
				}
			if(this.oldNameToNewName.containsKey(tokens[0]))
				{
				final String dest=this.oldNameToNewName.get(tokens[0]);
				if(dest.equals(tokens[1]))
					{
					LOG.warning(tokens[0]+" -> "+tokens[1]+" defined twice");
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
	public int doWork(final List<String> args) {
		
		if(mappingFile==null) {
			LOG.error("undefined mapping file");
			return -1;
		}
		
		try
			{
			parseNames(mappingFile);
			}
		catch(IOException err)
			{
			LOG.error(err);
			return -1;
			}
		if(this.oldNameToNewName.isEmpty())
			{
			LOG.warning("No replacement was defined");
			}
		if(new HashSet<String>(this.oldNameToNewName.values()).size()!=this.oldNameToNewName.size())
			{
			LOG.error("Some dest-name have been defined twice.");
			return -1;
			}
		
		return doVcfToVcf(args, outputFile);
		}

	@Override
	protected VariantContextWriter openVariantContextWriter(final File outorNull) throws IOException {
		return new PostponedVariantContextWriter(this.writingVcfArgs,stdout(),this.outputFile);
		}
	
	public static void main(final String[] args)
		{
		new VcfRenameSamples().instanceMainWithExit(args);
		}
	}
