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
* 2015 creation

*/
package com.github.lindenb.jvarkit.tools.vcfvcf;

import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.Set;

import com.github.lindenb.jvarkit.util.htsjdk.HtsjdkVersion;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter3;
import com.github.lindenb.jvarkit.util.vcf.IndexedVcfFileReader;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;


public class VcfPeekVcf extends AbstractVCFFilter3
	{
	private String TABIX=null;
	public Set<String> peek_info_tags=new HashSet<String>();
	private IndexedVcfFileReader indexedVcfFileReader=null;
	private String peekTagPrefix="";
	private boolean peekId=false;
	private boolean altAlleleCheck=false;
	

	
	public VcfPeekVcf()
		{
		}
	
	public Set<String> getPeekInfoTags() {
		return peek_info_tags;
		}
	
	@Override
	public int initializeKnime() {
		
		this.indexedVcfFileReader = null;
		try
			{
			this.indexedVcfFileReader = new IndexedVcfFileReader(TABIX);
			} 
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		return super.initializeKnime();
		}
	
	@Override
	public void disposeKnime() {
		CloserUtil.close(this.indexedVcfFileReader);
		this.indexedVcfFileReader=null;
		super.disposeKnime();
		}
	
	@Override
	protected void doWork(String inputSource, VcfIterator vcfIn,
			VariantContextWriter out) throws IOException
		{
		VCFHeader h = vcfIn.getHeader();
		VCFHeader h2 = new VCFHeader(h);
		
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkVersion",HtsjdkVersion.getVersion()));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkHome",HtsjdkVersion.getHome()));
		
		for(String key: this.peek_info_tags)
			{
			VCFInfoHeaderLine hinfo =this.indexedVcfFileReader.getHeader().getInfoHeaderLine(key);
			if(hinfo==null)
				{
				warning("INFO name="+key+" missing in "+this.TABIX);
				continue;
				}
			hinfo = VCFUtils.renameVCFInfoHeaderLine(hinfo, this.peekTagPrefix+key);
			if(h2.getInfoHeaderLine(hinfo.getID())!=null)
				{
				throw new IOException("key "+this.peekTagPrefix+key+" already defined in VCF header");
				}
			h2.addMetaDataLine(hinfo);;
			}
		
		out.writeHeader(h2);
		SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(h);
		while(vcfIn.hasNext())
			{
			VariantContext ctx=progress.watch(vcfIn.next());
						
			VariantContextBuilder vcb = new VariantContextBuilder(ctx);
			CloseableIterator<VariantContext> iter= this.indexedVcfFileReader.iterator(
					ctx.getContig(),
					Math.max(0,ctx.getStart()-1),
					(ctx.getEnd()+1)
					);
			while(iter.hasNext())
				{
				VariantContext ctx2=iter.next();
				if(!ctx.getContig().equals(ctx2.getContig())) continue;
				if(ctx.getStart()!=ctx2.getStart()) continue;
				if(!ctx.getReference().equals(ctx2.getReference())) continue;
				
				if(this.altAlleleCheck)
					{
					boolean found_all_alt=true;
					for(Allele alt: ctx.getAlternateAlleles())
						{
						if(!ctx2.hasAlternateAllele(alt))
							{
							found_all_alt=false;
							break;
							}
						}
					if(!found_all_alt) continue;
					}
				if(this.peekId && ctx2.hasID())
					{
					vcb.id(ctx2.getID());
					}
				for(String key: this.peek_info_tags)
					{
					if(!ctx2.hasAttribute(key)) continue;
					Object o = ctx2.getAttribute(key);
					vcb.attribute(this.peekTagPrefix+key, o);
					}
				}
			iter.close(); iter=null;
			
			out.add(vcb.make());
				
			incrVariantCount();
			if(checkOutputError()) break;
			}
		progress.finish();
		}
	
	public void setPeekTagPrefix(String peekTagPrefix) {
		this.peekTagPrefix = peekTagPrefix;
	}
	public void setPeekId(boolean peekId) {
		this.peekId = peekId;
	}
	
	public void setAltAlleleCheck(boolean altAlleleCheck) {
		this.altAlleleCheck = altAlleleCheck;
	}
	public void setIndexFile(String tABIX) {
		TABIX = tABIX;
		}
	@Override
	public String getProgramDescription() {
		return "Get the INFO from a VCF and use it for another VCF.";
		}
	
	@Override
    protected String getOnlineDocUrl() {
    	return DEFAULT_WIKI_PREFIX+"VcfPeekVcf";
    	}
	
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -f (file) The VCF file indexed with TABIX or tribble. Source of the annotations");
		out.println(" -t tag1,tag2,tag... the INFO keys to peek from the indexed file");
		out.println(" -i Replace the ID field if it exists.");
		out.println(" -a ALL alt allele must be found in indexed file.");
		out.println(" -p (prefix) prefix all database tags with this prefix to avoid collisions");
		out.println(" -o (out)  output file. default stdout");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"f:o:t:p:ia"))!=-1)
			{
			switch(c)
				{
				case 'f': setIndexFile(opt.getOptArg());break;
				case 'a': setAltAlleleCheck(true);break;
				case 'i': setPeekId(true); break;
				case 'p': setPeekTagPrefix(opt.getOptArg());break;
				case 't': for(String s: opt.getOptArg().split("[\\s;,]+"))
					{
					if(s.isEmpty()) continue;
					getPeekInfoTags().add(s);
					}
					break;
				case 'o': this.setOutputFile(opt.getOptArg());break;
				default:
					{
					switch(handleOtherOptions(c, opt,args))
						{
						case EXIT_FAILURE: return -1;
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
		new VcfPeekVcf().instanceMainWithExit(args);
		}
}
