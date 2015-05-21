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
package com.github.lindenb.jvarkit.tools.vcfstripannot;

import java.io.IOException;
import java.io.PrintStream;
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

import com.github.lindenb.jvarkit.util.htsjdk.HtsjdkVersion;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter3;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;


public class VCFReplaceTag extends AbstractVCFFilter3
	{
	private Set<String> userStrings=new HashSet<>();
	private Map<String,String> transformMap=new HashMap<>();
	private String replaceType=null;
	private int replaceTypeNo=-1;
	
	public VCFReplaceTag()
		{
		}
	
	@Override
	protected String getOnlineDocUrl()
		{
		return DEFAULT_WIKI_PREFIX+"VCFReplaceTag";
		}
	
	@Override
	protected String getProgramCommandLine()
		{
		return "Replace the key for INFO/FORMAT/FILTER";
		}
	
	public void setReplaceType(String replaceType) {
		this.replaceType = replaceType;
		}
	
	public Set<String> getUserReplace() {
		return userStrings;
		}
	
	
	@Override
	protected void doWork(String source,VcfIterator r, VariantContextWriter w)
			throws IOException
			{
			VCFHeader header=r.getHeader();
			
			HashSet<VCFHeaderLine> copyMeta= new HashSet<>(header.getMetaDataInInputOrder());
			
			for(String key:this.transformMap.keySet())
				{
				switch(this.replaceTypeNo)
					{
					case 0://INFO
						{
						VCFInfoHeaderLine info = header.getInfoHeaderLine(key);
						if(info!=null)
							{
							copyMeta.remove(info);
							copyMeta.add(VCFUtils.renameVCFInfoHeaderLine(info, this.transformMap.get(key)));
							}
						break;
						}
					case 1: //FORMAT
						{
						VCFFormatHeaderLine fmt = header.getFormatHeaderLine(key);
						if(fmt!=null)
							{
							copyMeta.remove(fmt);
							copyMeta.add(VCFUtils.renameVCFFormatHeaderLine(fmt, this.transformMap.get(key)));
							}
						break;
						}
					case 2: //FILTER
						{
						VCFFilterHeaderLine filter = header.getFilterHeaderLine(key);
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
			
			copyMeta.add(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
			copyMeta.add(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
			copyMeta.add(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkVersion",HtsjdkVersion.getVersion()));
			copyMeta.add(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkHome",HtsjdkVersion.getHome()));

			VCFHeader h2=new VCFHeader(copyMeta,header.getSampleNamesInOrder());
			
			
			
			SAMSequenceDictionaryProgress progress= new SAMSequenceDictionaryProgress(h2);

			
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
				this.incrVariantCount();
				if(this.checkOutputError()) break;
				}	
			progress.finish();
			}
	
	
	@Override
	public int initializeKnime() {
		if(this.replaceType==null)
			{
			error("Undefined replaceType");
			return -1;
			}
		this.replaceType= this.replaceType.toUpperCase();
		if(this.replaceType.equals("FILTER")){replaceTypeNo=2;}
		else if(this.replaceType.equals("INFO")){replaceTypeNo=0;}
		else if(this.replaceType.equals("FORMAT")){replaceTypeNo=1;}
		else
			{	
			error("Undefined replace type :"+this.replaceType);
			return -1;
			}
		for(String s:this.userStrings)
			{
			int slash=s.indexOf('/');
			if(slash==-1)
				{
				error("missing '/' in "+s);
				return -1;
				}
			if(slash==0 || slash+1==s.length())
				{
				error("bad '/' in "+s);
				return -1;
				}
			String key = s.substring(0,slash);
			if(replaceTypeNo==1)
				{
				if(key.equals("DP") || key.equals("GT")|| key.equals("PL")||
					key.equals("AD")|| key.equals("GQ"))
					{
					error("Cannot replace built-in FORMAT fields "+key);
					return -1;
					}
				}
			this.transformMap.put(key,s.substring(slash+1));
			}
		return super.initializeKnime();
		}
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -k (key-src/key-to) replace this tag. Multiple.");
		out.println(" -t (type) replace type: one of FORMAT,FILTER,INFO. Required");
		out.println(" -o (file) output file. Default stdout.");
		super.printOptions(out);
		}

	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "o:k:t:"))!=-1)
			{
			switch(c)
				{
				case 'k': this.getUserReplace().add(opt.getOptArg()); break;
				case 't': this.setReplaceType(opt.getOptArg()); break;
				case 'o': this.setOutputFile(opt.getOptArg());break;
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

	
	public static void main(String[] args) throws IOException
		{
		new VCFReplaceTag().instanceMainWithExit(args);
		}
	}
