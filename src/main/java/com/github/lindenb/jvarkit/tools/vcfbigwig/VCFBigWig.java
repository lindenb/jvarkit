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
* 2015 moved to htsjdk + knime

*/
package com.github.lindenb.jvarkit.tools.vcfbigwig;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;

import com.github.lindenb.jvarkit.util.htsjdk.HtsjdkVersion;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter3;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class VCFBigWig extends AbstractVCFFilter3
	{
	public boolean contained = true;
    private String biwWigFile;
	private BBFileReader bbFileReader=null;
	private String TAG=null;
	
	public VCFBigWig()
		{
		}
	
	public void setContained(boolean contained) {
		this.contained = contained;
		}
	
	private boolean isContained()
		{
		return contained;
		}
	
	@Override
	public int initializeKnime() {
		try
			{
			if(this.biwWigFile==null || this.biwWigFile.isEmpty())
				{
				error("Undefined BigWig file");
				return -1;
				}
			this.bbFileReader= new BBFileReader(this.biwWigFile);
			if(!this.bbFileReader.isBigWigFile())
				{
				this.bbFileReader=null;
				throw new IOException(this.biwWigFile+" is not a bigWIG file.");
				}
			
			if(this.TAG==null || this.TAG.isEmpty())
				{
				this.TAG=this.biwWigFile;
				int i=TAG.lastIndexOf(File.separator);
				if(i!=-1) TAG=TAG.substring(i+1);
				i=this.TAG.indexOf('.');
				this.TAG=this.TAG.substring(0,i);
				info("setting tag to "+this.TAG);
				}
			
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
		try
			{
			if(this.bbFileReader!=null)
				{
				CloserUtil.close(this.bbFileReader.getBBFis());
				}
			CloserUtil.close(this.bbFileReader);
			this.bbFileReader=null;
			}
		catch(Exception err)
			{
			error(err);
			}
		super.disposeKnime();
		}
	
	public void setBiwWigFile(String biwWigFile) {
		this.biwWigFile = biwWigFile;
		}
	
	public void setTag(String tAG) {
		this.TAG = tAG;
		}
	
	@Override
	protected void doWork(String inputSource, VcfIterator r,
			VariantContextWriter w) throws IOException
		{
		VCFHeader header=r.getHeader();
		
		VCFHeader h2=new VCFHeader(header.getMetaDataInInputOrder(),header.getSampleNamesInOrder());
		h2.addMetaDataLine(new VCFInfoHeaderLine(
				this.TAG,1,
				VCFHeaderLineType.Float,
				"Values from bigwig file: "+this.biwWigFile
				));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkVersion",HtsjdkVersion.getVersion()));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkHome",HtsjdkVersion.getHome()));
		w.writeHeader(h2);
		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header);
		
		List<Float> values=new ArrayList<Float>();
		while(r.hasNext())
			{
			
			VariantContext ctx = progress.watch(r.next());
			values.clear();
			
			BigWigIterator iter=this.bbFileReader.getBigWigIterator(
					ctx.getContig(),
					ctx.getStart()-1,
					ctx.getContig(),
					ctx.getEnd(),
					isContained()
					);
			while(iter!=null && iter.hasNext())
				{
				WigItem item=iter.next();
				float v=item.getWigValue();
				values.add(v);
				
				}
			if(values.isEmpty())
				{
				incrVariantCount();
				w.add(ctx);
				continue;
				}
			
			double total=0L;
			for(Float f:values) total+=f;
			VariantContextBuilder b=new VariantContextBuilder(ctx);
			b.attribute(this.TAG,(float)(total/values.size()));
			w.add(b.make());
			incrVariantCount();
			}
		progress.finish();
		}
	
	
	@Override
	public String getProgramDescription() {
		return "annotate a VCF with values from a bigwig file.";
		}
	
	@Override
    protected String getOnlineDocUrl() {
    	return DEFAULT_WIKI_PREFIX+"VCFBigWig";
    }
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -B (file) Path to the bigwig file.");
		out.println(" -T (name) of the INFO tag. default: name of the bigwig.");
		out.println(" -C Specifies wig values must be contained by region. if false: return any intersecting region values.");
		out.println(" -o (out)  output file. default stdout");
		super.printOptions(out);
		}

	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"B:Co:T:"))!=-1)
			{
			switch(c)
				{
				case 'B': this.setBiwWigFile(opt.getOptArg());break;
				case 'T': this.setTag(opt.getOptArg());break;
				case 'C': this.setContained(true); break;
				case 'o': this.setOutputFile(new File(opt.getOptArg()));break;
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
		new VCFBigWig().instanceMain(args);
		}
}
