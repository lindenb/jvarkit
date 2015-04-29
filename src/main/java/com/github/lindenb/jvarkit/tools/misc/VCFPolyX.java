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

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.htsjdk.HtsjdkVersion;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter3;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class VCFPolyX extends AbstractVCFFilter3
	{
    private File REF=null;
	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;
	
	
	public void setReference(File fasta)
		{
		this.REF=fasta;
		}
	
	@Override
	protected void doWork(String source,
			VcfIterator r, VariantContextWriter w)
			throws IOException
		{
		GenomicSequence genomicSequence=null;

		final String TAG="POLYX";
		VCFHeader header=r.getHeader();
		
		VCFHeader h2=new VCFHeader(header.getMetaDataInInputOrder(),header.getSampleNamesInOrder());
		h2.addMetaDataLine(new VCFInfoHeaderLine(TAG,1,
				VCFHeaderLineType.Integer,
				"number of repeated bases around REF")
				);
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkVersion",HtsjdkVersion.getVersion()));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkHome",HtsjdkVersion.getHome()));

		w.writeHeader(h2);
		SAMSequenceDictionaryProgress progress= new SAMSequenceDictionaryProgress(header);
		while(r.hasNext())
			{
			VariantContext ctx=progress.watch(r.next());
			VariantContextBuilder b=new VariantContextBuilder(ctx);
			if(genomicSequence==null || !ctx.getChr().equals(genomicSequence.getChrom()))
				{
				this.info("loading chromosome "+ctx.getChr());
				genomicSequence=new GenomicSequence(this.indexedFastaSequenceFile, ctx.getChr());
				}
			int pos0=ctx.getStart()-1;
			int count=1;
			char c0=Character.toUpperCase(genomicSequence.charAt(pos0));
			//go left
			pos0--;
			while(pos0>=0 && c0==Character.toUpperCase(genomicSequence.charAt(pos0)))
				{
				++count;
				pos0--;
				}
			//go right
			pos0=ctx.getEnd()-1;
			c0=Character.toUpperCase(genomicSequence.charAt(pos0));
			pos0++;
			while(pos0< genomicSequence.getSAMSequenceRecord().getSequenceLength()
				&& c0==Character.toUpperCase(genomicSequence.charAt(pos0)))
				{
				++count;
				++pos0;
				}
			b.attribute(TAG,count);
			w.add(b.make());
			
			if(this.checkOutputError()) break;
			}		
		progress.finish();
		genomicSequence=null;
		}

	@Override
	public String getProgramDescription() {
		return "Number of repeated REF bases around POS.";
		}
	
	@Override
    protected String getOnlineDocUrl() {
    	return DEFAULT_WIKI_PREFIX+"VCFPolyX";
    	}
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -R (fasta) fasta reference indexed with samtools faidx and picard.");
		out.println(" -o (out)  output file. default stdout");
		super.printOptions(out);
		}

	
	@Override
	public int initializeKnime() {
		try {
			if(this.REF==null)
				{
				error("Undefined Reference");
				return -1;
				}
			this.info("opening reference "+REF);
			this.indexedFastaSequenceFile=new IndexedFastaSequenceFile(REF);
			}
		catch (Exception e) {
			error(e);
			return -1;
			}
		return super.initializeKnime();
		}
	@Override
	public void disposeKnime() {
		CloserUtil.close(this.indexedFastaSequenceFile);
		this.indexedFastaSequenceFile=null;
		}
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"o:R:"))!=-1)
			{
			switch(c)
				{
				case 'R': this.setReference(new File(opt.getOptArg()));break;
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

	
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new VCFPolyX().instanceMainWithExit(args);
		}

}
