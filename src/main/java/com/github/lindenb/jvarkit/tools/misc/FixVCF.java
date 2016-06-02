package com.github.lindenb.jvarkit.tools.misc;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.SequenceInputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.SynchronousLineReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.AbstractVCFCodec;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.jvarkit.util.vcf.VcfIteratorImpl;

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
					switch(handleOtherOptions(c, opt, args))
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
				w= VCFUtils.createVariantContextWriterToStdout();
				}
			else
				{
				this.info("writing to "+fileout);
				w= VCFUtils.createVariantContextWriter(fileout);
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
		
		final AbstractVCFCodec vcfCodec = VCFUtils.createDefaultVCFCodec();
		
		LineIterator r= new LineIteratorImpl(new SynchronousLineReader(vcfStream));
		final VCFHeader header=(VCFHeader) vcfCodec.readActualHeader(r);
		
		//samples names have been changed by picard api and reordered !!!
		//re-create the original order
		List<String> sampleNamesInSameOrder=new ArrayList<String>(header.getSampleNamesInOrder().size());
		for(int col=0;col< header.getSampleNamesInOrder().size();++col )
			{
			for(String sample: header.getSampleNameToOffset().keySet())
				{
				if(header.getSampleNameToOffset().get(sample)==col)
					{
					sampleNamesInSameOrder.add(sample);
					break;
					}
				}
			}
		if(sampleNamesInSameOrder.size()!=header.getSampleNamesInOrder().size())
			{
			throw new IllegalStateException();
			}
		
		VCFHeader h2=new VCFHeader(
				header.getMetaDataInInputOrder(),
				sampleNamesInSameOrder
				);
		
		File tmp=IOUtil.newTempFile("tmp", ".vcf.gz",new File[]{tmpDir});
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
		
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName(),
				"Saved VCF FILTER AND INFO from "+filenameIn
				));

		
		//save header in memory
		ByteArrayOutputStream baos=new ByteArrayOutputStream();
		VariantContextWriter w2= VCFUtils.createVariantContextWriterToOutputStream(baos);
		w2.writeHeader(h2);
		w2.close();
		baos.close();
		 
		//reopen tmp file

		@SuppressWarnings("resource")
		VcfIterator in=new VcfIteratorImpl(new SequenceInputStream(
				new ByteArrayInputStream(baos.toByteArray()),
				new GZIPInputStream(new FileInputStream(tmp)))
				);
		
		w.writeHeader(h2);

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
