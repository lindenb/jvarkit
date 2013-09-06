package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.EnumSet;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import net.sf.picard.vcf.VcfIterator;
import net.sf.samtools.util.BlockCompressedOutputStream;

import org.broad.tribble.readers.LineIterator;
import org.broad.tribble.readers.LineIteratorImpl;
import org.broad.tribble.readers.LineReaderUtil;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.writer.Options;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriterFactory;
import org.broadinstitute.variant.vcf.VCFCodec;
import org.broadinstitute.variant.vcf.VCFConstants;
import org.broadinstitute.variant.vcf.VCFFilterHeaderLine;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;
import org.broadinstitute.variant.vcf.VCFHeaderLineCount;
import org.broadinstitute.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.picard.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.IOUtils;

public class FixVCF extends AbstractCommandLineProgram
	{
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+"  Fix a VCF if INFO or FILTER are missing";
	@Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="VCF file/URL to process. Default stdin. ",optional=true)
	public String IN=null;
	@Option(shortName= StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="VCF file to generate. Default stdout. ",optional=true)
	public File OUT=null;
	
	
	private static final Log LOG=Log.getInstance(FixVCF.class);

	
	
	@Override
	public String getVersion() {
		return "1.0";
		}
	
	private void doWork(
			InputStream vcfStream, VariantContextWriter w)
			throws IOException
		{
		
		VCFCodec vcfCodec = new VCFCodec();
		
		LineIterator r= new LineIteratorImpl(LineReaderUtil.fromBufferedStream(vcfStream));
		VCFHeader header=(VCFHeader) vcfCodec.readActualHeader(r);
		
		VCFHeader h2=new VCFHeader(header.getMetaDataInInputOrder(),header.getSampleNamesInOrder());
		
		File tmp=IoUtil.newTempFile("tmp", ".vcf.gz", super.TMP_DIR.toArray(new File[super.TMP_DIR.size()]));
		tmp.deleteOnExit();

		PrintWriter pw=new PrintWriter(new GZIPOutputStream(new FileOutputStream(tmp)));
		while(r.hasNext())
			{
			String line=r.next();
			pw.println(line);
			VariantContext ctx=vcfCodec.decode(line);
			for(String f:ctx.getFilters())
				{
				if(h2.getFilterHeaderLine(f)!=null) continue;
				//if(f.equals(VCFConstants.PASSES_FILTERS_v4)) continue; hum...
				if(f.isEmpty() || f.equals(VCFConstants.UNFILTERED)) continue;
 				LOG.info("Fixing missing Filter:"+f);
				h2.addMetaDataLine(new VCFFilterHeaderLine(f));
				}
			for(String tag:ctx.getAttributes().keySet())
				{
				if(h2.getInfoHeaderLine(tag)!=null) continue;
				LOG.info("Fixing missing INFO:"+tag);
				h2.addMetaDataLine(new VCFInfoHeaderLine(tag, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "undefined. Saved by "+getClass()));
				}
			}
		pw.flush();
		pw.close();
		pw=null;
		
		LOG.info("re-reading VCF frm tmpFile:" +tmp);
		//reopen tmp file
		VcfIterator in=new VcfIterator(new GZIPInputStream(new FileInputStream(tmp)));
		
		w.writeHeader(h2);
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName(),"Saved VCF FILTER AND INFO from IN="+this.IN));
		while(in.hasNext())
			{
			w.add(in.next());
			}
		in.close();
		tmp.delete();
		}
	
	
	protected VariantContextWriter createVariantContextWriter() throws IOException
	{
	if(OUT==null)
		{
		LOG.info("writing to stdout");
		return VariantContextWriterFactory.create(System.out,null,EnumSet.noneOf(Options.class));
		}
	else if(OUT.getName().endsWith(".gz"))
		{
		LOG.info("writing to "+OUT+" as bgz file.");
		BlockCompressedOutputStream bcos=new BlockCompressedOutputStream(OUT);
		return VariantContextWriterFactory.create(bcos,null,EnumSet.noneOf(Options.class));
		}
	else
		{
		LOG.info("writing to "+OUT);
		return  VariantContextWriterFactory.create(OUT,null,EnumSet.noneOf(Options.class));
		}
	}

	
	@Override
	protected int doWork()
		{
		InputStream r=null;
		VariantContextWriter w=null;
		try
			{
			if(IN==null)
				{
				LOG.info("reading from stdin");
				r=System.in;
				}
			else
				{
				LOG.info("reading from "+IN);
				r=IOUtils.openURIForReading(IN);
				}
			w=this.createVariantContextWriter();
			doWork(r,w);
			}
		catch (Exception e)
			{
			LOG.error(e);
			testRemoteGit();
			return -1;
			}
		finally
			{
			if(w!=null) w.close();
			if(r!=null) try{r.close();}catch(IOException err){}
			}	
		return 0;
		}

	
	
	
	public static void main(String[] args)
		{
		new FixVCF().instanceMainWithExit(args);
		}
	}
