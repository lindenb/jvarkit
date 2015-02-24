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
package com.github.lindenb.jvarkit.tools.vcfbigwig;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;

import com.github.lindenb.jvarkit.util.picard.cmdline.Option;
import com.github.lindenb.jvarkit.util.picard.cmdline.Usage;
import htsjdk.samtools.util.Log;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class VCFBigWig extends AbstractVCFFilter
	{
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+" annotate a value from a bigwig file..";
	
	private static final Log LOG=Log.getInstance(VCFBigWig.class);
	
	@Option(shortName="BW",doc="Path to the bigwig file.",optional=false)
	public String BIGWIG;
	@Option(shortName="INFOTAG",doc="name of the INFO tag. default: name of the bigwig.",optional=true)
	public String TAG=null;

    @Option(shortName="CONT",doc=
    		"Specifies wig values must be contained by region. if false: return any intersecting region values.",optional=true)
	public boolean contained = true;
	
	private BBFileReader bbFileReader=null;
	
	private boolean isContained()
		{
		return contained;
		}
	
	
	@Override
	protected void doWork(VcfIterator r, VariantContextWriter w)
			throws IOException
		{
		LOG.info("opening bigwig: "+this.BIGWIG);
		this.bbFileReader=new BBFileReader(this.BIGWIG);
		if(!this.bbFileReader.isBigWigFile())
			{
			throw new IOException(BIGWIG+" is not a bigWIG file.");
			}
		
		if(TAG==null)
			{
			TAG=this.BIGWIG;
			int i=TAG.lastIndexOf(File.separator);
			if(i!=-1) TAG=TAG.substring(i+1);
			i=TAG.indexOf('.');
			TAG=TAG.substring(0,i);
			LOG.info("setting tag to "+this.TAG);
			}
		
		
		VCFHeader header=r.getHeader();
		
		VCFHeader h2=new VCFHeader(header.getMetaDataInInputOrder(),header.getSampleNamesInOrder());
		h2.addMetaDataLine(new VCFInfoHeaderLine(this.TAG,1,VCFHeaderLineType.Float,"Values from bigwig file: "+getCommandLine()));
		
		w.writeHeader(h2);
		
		
		List<Float> values=new ArrayList<Float>();
		while(r.hasNext())
			{
			
			VariantContext ctx=r.next();
			

			values.clear();
			
			BigWigIterator iter=this.bbFileReader.getBigWigIterator(
					ctx.getChr(),
					ctx.getStart()-1,
					ctx.getChr(),
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
				w.add(ctx);
				continue;
				}
			
			double total=0L;
			for(Float f:values) total+=f;

			
			VariantContextBuilder b=new VariantContextBuilder(ctx);
			
			b.attribute(this.TAG,(float)(total/values.size()));
			w.add(b.make());
			}
		}
	
	public static void main(String[] args) throws IOException
		{
		new VCFBigWig().instanceMain(args);
		}
}
