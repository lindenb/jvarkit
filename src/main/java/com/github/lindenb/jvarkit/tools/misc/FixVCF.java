package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.EnumSet;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import net.sf.picard.io.IoUtil;
import net.sf.samtools.util.BlockCompressedOutputStream;
import net.sf.samtools.util.CloserUtil;

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

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class FixVCF
	extends com.github.lindenb.jvarkit.util.AbstractCommandLineProgram
	{
	private File tmpDir=null;

	@Override
	public String getProgramDescription() {
		return "Fix a VCF if INFO or FILTER are missing";
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/FixVCF";
		}
	
	@Override
	public void printOptions(PrintStream out) {
		out.println(" -o (filenameout) optional. default stdout.");
		out.println(" -T (dir) tmp directory. Optional");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		File fileout=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"o:T:"))!=-1)
			{
			switch(c)
				{
				case 'o':fileout=new File(args[opt.getOptInd()]);break;
				case 'T':tmpDir=new File(args[opt.getOptInd()]);break;
				default:
					{
					switch(handleOtherOptions(c, opt))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		
		if(tmpDir==null)
			{
			tmpDir=new File(System.getProperty("java.io.tmpdir"));
			}
		
		VariantContextWriter w=null;
		try
			{
			if(fileout==null)
				{
				this.info("writing to stdout");
				w= VariantContextWriterFactory.create(System.out,null,EnumSet.noneOf(Options.class));
				}
			else if(fileout.getName().endsWith(".gz"))
				{
				this.info("writing to "+fileout+" as bgz file.");
				BlockCompressedOutputStream bcos=new BlockCompressedOutputStream(fileout);
				w= VariantContextWriterFactory.create(bcos,null,EnumSet.noneOf(Options.class));
				}
			else
				{
				this.info("writing to "+fileout);
				w=  VariantContextWriterFactory.create(fileout,null,EnumSet.noneOf(Options.class));
				}
			

			
			if(opt.getOptInd()==args.length)
				{
				info("Reading from stdin");
				doWork("stdin",System.in,w);
				}
			else if(opt.getOptInd()+1==args.length)
				{
				String filename=args[opt.getOptInd()];
				info("Reading from "+filename);
				InputStream in=IOUtils.openURIForReading(filename);
				doWork(filename,System.in,w);
				CloserUtil.close(in);
				}
			else
				{
				error("Illegal number of arguments.");
				return -1;
				}
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(w);
			}
		}
	
	private int doWork(
			String filenameIn,
			InputStream vcfStream, VariantContextWriter w)
			throws IOException
		{
		
		VCFCodec vcfCodec = new VCFCodec();
		
		LineIterator r= new LineIteratorImpl(LineReaderUtil.fromBufferedStream(vcfStream));
		VCFHeader header=(VCFHeader) vcfCodec.readActualHeader(r);
		
		VCFHeader h2=new VCFHeader(header.getMetaDataInInputOrder(),header.getSampleNamesInOrder());
		
		File tmp=IoUtil.newTempFile("tmp", ".vcf.gz",new File[]{tmpDir});
		tmp.deleteOnExit();

		PrintWriter pw=new PrintWriter(new GZIPOutputStream(new FileOutputStream(tmp)));
		while(r.hasNext())
			{
			String line=r.next();
			pw.println(line);
			VariantContext ctx=null;
			
			try
				{
				ctx=vcfCodec.decode(line);
				}
			catch(Exception err)
				{
				pw.close();
				error(line);
				error(err);
				return -1;
				}
			for(String f:ctx.getFilters())
				{
				if(h2.getFilterHeaderLine(f)!=null) continue;
				//if(f.equals(VCFConstants.PASSES_FILTERS_v4)) continue; hum...
				if(f.isEmpty() || f.equals(VCFConstants.UNFILTERED)) continue;
 				info("Fixing missing Filter:"+f);
				h2.addMetaDataLine(new VCFFilterHeaderLine(f));
				}
			for(String tag:ctx.getAttributes().keySet())
				{
				if(h2.getInfoHeaderLine(tag)!=null) continue;
				info("Fixing missing INFO:"+tag);
				h2.addMetaDataLine(new VCFInfoHeaderLine(tag, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "undefined. Saved by "+getClass()));
				}
			}
		pw.flush();
		pw.close();
		pw=null;
		
		info("re-reading VCF frm tmpFile:" +tmp);
		//reopen tmp file
		
		VcfIterator in=new VcfIterator(new GZIPInputStream(new FileInputStream(tmp)));
		
		w.writeHeader(h2);
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName(),
				"Saved VCF FILTER AND INFO from "+filenameIn
				));
		while(in.hasNext())
			{
			w.add(in.next());
			}
		in.close();
		tmp.delete();
		return 0;
		}
	

	
	public static void main(String[] args)
		{
		new FixVCF().instanceMainWithExit(args);
		}
	}
