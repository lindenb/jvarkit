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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.vcfstripannot;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import com.github.lindenb.jvarkit.util.log.Logger;


import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import htsjdk.variant.vcf.VCFIterator;

/**
BEGIN_DOC

## Deprecated 

use `bcftools annotate` with option `-c`

## Example

```bash
$   java -jar dist/vcfreplacetag.jar -t INFO -k VDB/NEWNAME ~/jeter.vcf 
```

END_DOC
 */
@Program(
		name="vcfreplacetag",
		description="Replace the key for INFO/FORMAT/FILTER",
		keywords={"vcf"},
		deprecatedMsg="use `bcftools annotate` with option `-c`",
		modificationDate="20190321"
		)
public class VCFReplaceTag extends Launcher
	{
	private static final Logger LOG = Logger.build(VCFReplaceTag.class).make();	
	
	@Parameter(names={"-t","--type"},description="replace type: one of FORMAT,FILTER,INFO",required=true)
	private String replaceType = "INFO";

	@Parameter(names={"-k","--tag"},description="tag to replace . Format FROM/TO")
	private List<String> userTagList = new ArrayList<>();
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	
	private Map<String,String> transformMap=new HashMap<>();
	private int replaceTypeNo=-1;
	
	public VCFReplaceTag()
		{
		}
	
	@Override
	protected int doVcfToVcf(
			final String inputName,
			final VCFIterator r,
			final VariantContextWriter w)
			 {

			final VCFHeader header=r.getHeader();
			
			final HashSet<VCFHeaderLine> copyMeta= new HashSet<>(header.getMetaDataInInputOrder());
			
			for(final String key:this.transformMap.keySet())
				{
				switch(this.replaceTypeNo)
					{
					case 0://INFO
						{
						final VCFInfoHeaderLine info = header.getInfoHeaderLine(key);
						if(info!=null)
							{
							copyMeta.remove(info);
							copyMeta.add(VCFUtils.renameVCFInfoHeaderLine(info, this.transformMap.get(key)));
							}
						break;
						}
					case 1: //FORMAT
						{
							final VCFFormatHeaderLine fmt = header.getFormatHeaderLine(key);
						if(fmt!=null)
							{
							copyMeta.remove(fmt);
							copyMeta.add(VCFUtils.renameVCFFormatHeaderLine(fmt, this.transformMap.get(key)));
							}
						break;
						}
					case 2: //FILTER
						{
						final VCFFilterHeaderLine filter = header.getFilterHeaderLine(key);
						if(filter!=null)
							{
							copyMeta.remove(filter);
							copyMeta.add(VCFUtils.renameVCFFilterHeaderLine(filter, this.transformMap.get(key)));
							}
						break;
						}
					default: throw new IllegalStateException(""+this.replaceTypeNo);
					}
				}
			

			final VCFHeader h2=new VCFHeader(copyMeta,header.getSampleNamesInOrder());
			addMetaData(h2);
			
			
			final SAMSequenceDictionaryProgress progress= new SAMSequenceDictionaryProgress(h2);
			w.writeHeader(h2);
			while(r.hasNext())
				{
				VariantContext ctx=progress.watch(r.next());
				VariantContextBuilder b=new VariantContextBuilder(ctx);
				
					switch(this.replaceTypeNo)
						{
						case 0://INFO
							{
							for(String key:this.transformMap.keySet())
								{
								Object o = ctx.getAttribute(key);
								if(o!=null)
									{
									b.rmAttribute(key);
									b.attribute(this.transformMap.get(key), o);
									}
								}
							break;
							}
						case 1: //FORMAT
							{
							List<Genotype> newgenotypes=new ArrayList<>( ctx.getNSamples());
							for(int i=0;i< ctx.getNSamples();++i)
								{
								Genotype g= ctx.getGenotype(i);
								Map<String,Object> atts = g.getExtendedAttributes();
								GenotypeBuilder gb=new GenotypeBuilder(g);
								for(String key:this.transformMap.keySet())
									{
									Object o = atts.get(key);
									if(o!=null)
										{
										atts.remove(key);
										atts.put(this.transformMap.get(key),o);
										}
									}
								gb.attributes(atts);
								newgenotypes.add(gb.make());
								}
							b.genotypes(newgenotypes);
							break;
							}
						case 2: //FILTER
							{
							Set<String> filters=new HashSet<>(ctx.getFilters());
							for(String key:this.transformMap.keySet())
								{
								if(filters.contains(key))
									{
									filters.remove(key);
									filters.add(this.transformMap.get(key));
									}
								}
							b.filters(filters);
							break;
							}
						default: throw new IllegalStateException(""+this.replaceTypeNo);
						}
					
				w.add(b.make());
				if(w.checkError()) break;
				}	
			progress.finish();
			LOG.info("done");
			return 0;
			}
	
	@Override
	public int doWork(List<String> args) {
		if(this.replaceType==null)
			{
			throw new JvarkitException.CommandLineError("Undefined replaceType");
			}
		this.replaceType= this.replaceType.toUpperCase();
		if(this.replaceType.equals("FILTER")){replaceTypeNo=2;}
		else if(this.replaceType.equals("INFO")){replaceTypeNo=0;}
		else if(this.replaceType.equals("FORMAT")){replaceTypeNo=1;}
		else
			{	
			throw new JvarkitException.CommandLineError("Undefined replace type :"+this.replaceType);
			}
		for(final String s:this.userTagList)
			{
			int slash=s.indexOf('/');
			if(slash==-1)
				{
				throw new JvarkitException.CommandLineError("missing '/' in "+s);
				}
			if(slash==0 || slash+1==s.length())
				{
			 throw new JvarkitException.CommandLineError("bad '/' in "+s);
				}
			final String key = s.substring(0,slash);
			if(replaceTypeNo==1)
				{
				if(key.equals("DP") || key.equals("GT")|| key.equals("PL")||
					key.equals("AD")|| key.equals("GQ"))
					{
					throw new JvarkitException.CommandLineError("Cannot replace built-in FORMAT fields "+key);
					}
				}
			this.transformMap.put(key,s.substring(slash+1));
			}
		return doVcfToVcf(args,outputFile);
		}

	
	public static void main(String[] args) throws IOException
		{
		new VCFReplaceTag().instanceMainWithExit(args);
		}
	}
