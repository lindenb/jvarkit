package com.github.lindenb.jvarkit.tools.vcfbigwig;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.util.Log;

import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;
import org.broad.tribble.readers.LineReader;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFCodec;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter;

public class VCFBigWig extends AbstractVCFFilter
	{
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+" annotate a value from a bigwig file..";
	
	private static final Log LOG=Log.getInstance(VCFBigWig.class);
	
	@Option(shortName="BW",doc="Path to the bigwig file.",optional=false)
	public String BIGWIG;
	@Option(shortName="TAG",doc="name of the INFO tag .",optional=true)
	public String TAG=null;

	
	private BBFileReader bbFileReader=null;
	private boolean contained=true;
	
	@Override
	protected void doWork(LineReader in, VariantContextWriter w)
			throws IOException
		{
		LOG.info("opening "+this.BIGWIG);
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
		
		
		VCFCodec codeIn=new VCFCodec();		
		VCFHeader header=(VCFHeader)codeIn.readHeader(in);
		
		VCFHeader h2=new VCFHeader(header.getMetaDataInInputOrder(),header.getSampleNamesInOrder());
		h2.addMetaDataLine(new VCFInfoHeaderLine(this.TAG,1,VCFHeaderLineType.String,"Values from bigwig file: "+BIGWIG));
		
		w.writeHeader(h2);
		
		
		String line;
		List<Float> values=new ArrayList<Float>();
		while((line=in.readLine())!=null)
			{
			
			VariantContext ctx=codeIn.decode(line);
			

			values.clear();
			
			BigWigIterator iter=this.bbFileReader.getBigWigIterator(
					ctx.getChr(),
					ctx.getStart()-1,
					ctx.getChr(),
					ctx.getEnd(),
					this.contained
					);
			while(iter.hasNext())
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
