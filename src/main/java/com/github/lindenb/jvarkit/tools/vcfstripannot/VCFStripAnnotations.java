package com.github.lindenb.jvarkit.tools.vcfstripannot;

import java.io.File;
import java.io.IOException;
import java.util.EnumSet;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import org.broad.tribble.readers.AsciiLineReader;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.variantcontext.writer.Options;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriterFactory;
import org.broadinstitute.variant.vcf.VCFCodec;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.util.Log;

import com.github.lindenb.jvarkit.util.picard.IOUtils;


public class VCFStripAnnotations extends CommandLineProgram
	{
	private static final Log LOG=Log.getInstance(VCFStripAnnotations.class);
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+" Removes one or more field from the INFO column of a VCF.";
    
	@Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="VCF file to process. Default stdin. ",optional=true)
	public File IN=null;
	@Option(shortName= StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="VCF file to generate. Default stdout. ",optional=true)
	public File OUT=null;
	@Option(shortName="K", doc="remove this INFO key",minElements=0)	
	public Set<String> KEY=new HashSet<String>();

	@Option(shortName="NF", doc="Reset the FILTER column",optional=true)	
	public boolean RESET_FILTER=false;

	
	@Override
	protected int doWork()
		{
		VariantContextWriter w=null;
		AsciiLineReader r=null;
		try {
			VCFCodec codeIn=new VCFCodec();
			
		
			if(IN==null)
				{	
				LOG.info("reading from stdin");
				r=new  AsciiLineReader(System.in);
				}
			else
				{
				LOG.info("reading from "+IN);
				r=new  AsciiLineReader(IOUtils.openFileForReading(IN));
				}
			
			VCFHeader header=(VCFHeader)codeIn.readHeader(r);
			
			EnumSet<Options> options = EnumSet.noneOf(Options.class);
			if(OUT==null) 
				{
				LOG.info("writing to stdout");
				w=VariantContextWriterFactory.create(System.out,
						null, options);
				}
			else
				{
				LOG.info("writing to "+OUT);
				w=VariantContextWriterFactory.create(OUT,
						null, options);
				}
			VCFHeader h2=new VCFHeader(header);
			
			for(Iterator<VCFInfoHeaderLine> h=h2.getInfoHeaderLines().iterator();
					h.hasNext();)
				{
				VCFInfoHeaderLine vih=h.next();
				if(this.KEY.contains(vih.getID()))
					h.remove();
				}
			h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName(),"REMOVED: "+this.KEY.toString()));
			w.writeHeader(h2);
			String line;
			
			while((line=r.readLine())!=null)
				{
				VariantContext ctx=codeIn.decode(line);
				VariantContextBuilder b=new VariantContextBuilder(ctx);
				for(String key:KEY) b.rmAttribute(key);
				if(RESET_FILTER) b.unfiltered();
				w.add(b.make());
				}
			return 0;
		} catch (Exception e) {
			LOG.error(e);
			return -1;
		}
		finally
		{
			if(r!=null) r.close();
			if(w!=null) w.close();
		}
	}
	
	
	public static void main(String[] args) throws IOException
		{
		new VCFStripAnnotations().instanceMainWithExit(args);
		}
}
