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
package com.github.lindenb.jvarkit.tools.vcfbed;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.Set;

import com.github.lindenb.jvarkit.util.bio.bed.IndexedBedReader;
import com.github.lindenb.jvarkit.util.htsjdk.HtsjdkVersion;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter3;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;


/**
 * 
 * VCFBed
 *
 */
public class VCFBed extends AbstractVCFFilter3
	{
	private IndexedBedReader bedReader =null;
	private Chunk parsedFormat=null;
	private String FORMAT="${1}:${2}-${3}";
	private File TABIX;
	private String TAG="TAG";
	
	private static abstract class Chunk
		{
		public abstract String toString(IndexedBedReader.BedLine tokens);
		public Chunk next=null;
		}
	
	private static class PlainChunk extends Chunk
		{
		String s;
		PlainChunk(String s){this.s=s;}
		public String toString(IndexedBedReader.BedLine tokens)
			{
			return s+(next==null?"":next.toString(tokens));
			}
		}
	private static class ColChunk extends Chunk
		{
		int index;
		ColChunk(int index){ this.index=index;}
		public String toString(IndexedBedReader.BedLine tokens)
			{
			String s= tokens.get(index);
			if(s==null) s="";
			return s+(next==null?"":next.toString(tokens));
			}
		}

	
	private Chunk parseFormat(String s)
		{
		if(s==null || s.isEmpty()) return null;
		if(s.startsWith("${"))
			{
			int j=s.indexOf('}',2);
			if(j==-1) throw new IllegalArgumentException("bad format in \""+s+"\".");
			try
				{
				int col=Integer.parseInt(s.substring(2, j).trim());
				if(col<1) throw new IllegalArgumentException();
				ColChunk c=new ColChunk(col-1);
				c.next=parseFormat(s.substring(j+1));
				return c;
				}
			catch(Exception err)
				{
				 throw new IllegalArgumentException("bad format in \""+s+"\".",err);
				}
			}
		else if(s.startsWith("$"))
			{
			int j=1;
			while(j<s.length() && Character.isDigit(s.charAt(j)))
				{
				++j;
				}
			int col=Integer.parseInt(s.substring(1, j).trim());
			if(col<1) throw new IllegalArgumentException();
			ColChunk c=new ColChunk(col-1);
			c.next=parseFormat(s.substring(j));
			return c;
			}
		int i=0;
		StringBuilder sb=new StringBuilder();
		while(i< s.length() && s.charAt(i)!='$')
			{
			sb.append(s.charAt(i));
			i++;
			}
		PlainChunk c=new PlainChunk(sb.toString());
		c.next=parseFormat(s.substring(i));
		return c;
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return DEFAULT_WIKI_PREFIX+"VCFBed";
		}
	
	@Override
	public String getProgramDescription() {
		return " Cross information between a VCF and a BED .";
		}
	

	@Override
	protected void doWork(String inputSource,VcfIterator r, VariantContextWriter w)
			throws IOException
		{
		VCFHeader header=r.getHeader();

		VCFHeader h2=new VCFHeader(header.getMetaDataInInputOrder(),header.getSampleNamesInOrder());
		h2.addMetaDataLine(new VCFInfoHeaderLine(TAG, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "metadata added from "+TABIX+" . Format was "+FORMAT));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkVersion",HtsjdkVersion.getVersion()));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkHome",HtsjdkVersion.getHome()));

		
		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header);
		w.writeHeader(h2);
		while(r.hasNext())
			{
			VariantContext ctx= progress.watch(r.next());
			Set<String> annotations=new HashSet<String>();
			
			CloseableIterator<IndexedBedReader.BedLine> iter = this.bedReader.iterator(
					ctx.getContig(),
					ctx.getStart()-1,
					ctx.getEnd()+1
					);
			while(iter.hasNext())
				{
				IndexedBedReader.BedLine bedLine = iter.next();
				
				if(!ctx.getContig().equals(bedLine.getContig())) continue;
				if(ctx.getStart()-1 >= bedLine.getEnd() ) continue;
				if(ctx.getEnd()-1 < bedLine.getStart() ) continue;

				
				String newannot=this.parsedFormat.toString(bedLine);
				if(!newannot.isEmpty())
					annotations.add(VCFUtils.escapeInfoField(newannot));
				}
			CloserUtil.close(iter);
			
			if(annotations.isEmpty())
				{
				w.add(ctx);
				continue;
				}
			VariantContextBuilder vcb=new VariantContextBuilder(ctx);
			vcb.attribute(TAG, annotations.toArray());
			w.add(vcb.make());
			incrVariantCount();
			if(checkOutputError()) break;
			}
		progress.finish();
		}
	
	public void setFormat(String fORMAT) {
		FORMAT = fORMAT;
		}
	
	public void setTag(String tAG) {
		TAG = tAG;
		}
	
	public void setBedFile(File tABIX) {
		TABIX = tABIX;
		}
	
	public void printOptions(PrintStream out) {
		out.println(" -f (format). Field with ${number} will be replaced with the column of the BED. default:"+FORMAT);
		out.println(" -T (INFO). Key for the INFO field. default:"+TAG);
		out.println(" -B (path). BED file indexed with tribble/tabix");
		super.printOptions(out);
		}

	
	@Override
	public int initializeKnime() {
		try
			{
			this.info("parsing "+this.FORMAT);
			this.parsedFormat=parseFormat(this.FORMAT);
			if(this.parsedFormat==null) this.parsedFormat=new PlainChunk("");
			
			if(this.TABIX==null)
				{
				error("Undefined tabix file");
				return -1;
				}
			
			this.info("opening Bed "+this.TABIX);
			this.bedReader= new IndexedBedReader(this.TABIX);
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		return super.initializeKnime();
		}
	
	@Override
	public void disposeKnime()
		{
		CloserUtil.close(this.bedReader);
		this.bedReader = null;
		this.parsedFormat = null;
		super.disposeKnime();
		}
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"o:f:T:B:"))!=-1)
			{
			switch(c)
				{
				case 'f': this.setFormat(opt.getOptArg());break;
				case 'T': this.setTag(opt.getOptArg());break;
				case 'B': this.setBedFile(new File(opt.getOptArg()));break;
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

	
	
	public static void main(String[] args) throws Exception
		{
		new VCFBed().instanceMainWithExit(args);
		}
}
