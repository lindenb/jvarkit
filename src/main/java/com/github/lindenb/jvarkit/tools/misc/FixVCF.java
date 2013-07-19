package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;

import org.broad.tribble.readers.LineReader;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFCodec;
import org.broadinstitute.variant.vcf.VCFConstants;
import org.broadinstitute.variant.vcf.VCFFilterHeaderLine;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;
import org.broadinstitute.variant.vcf.VCFHeaderLineCount;
import org.broadinstitute.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.AbstractVCFFilter;

public class FixVCF extends AbstractVCFFilter
	{
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+"  Fix a VCF if INFO or FILTER are missing";
	
	
	private static final Log LOG=Log.getInstance(FixVCF.class);

	@Override
	public String getVersion() {
		return "1.0";
		}
	
	@Override
	protected void doWork(LineReader r, VariantContextWriter w)
			throws IOException
		{
		VCFCodec codeIn=new VCFCodec();		
		VCFHeader header=(VCFHeader)codeIn.readHeader(r);
		String line;
		VCFHeader h2=new VCFHeader(header.getMetaDataInInputOrder(),header.getSampleNamesInOrder());
		
		File tmp=IoUtil.newTempFile("tmp", ".vcf.gz", super.TMP_DIR.toArray(new File[super.TMP_DIR.size()]));
		tmp.deleteOnExit();

		PrintWriter pw=new PrintWriter(new GZIPOutputStream(new FileOutputStream(tmp)));
		while((line=r.readLine())!=null)
			{
			pw.println(line);
			VariantContext ctx=codeIn.decode(line);
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
		BufferedReader in=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(tmp))));
		
		w.writeHeader(h2);
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName(),"Saved VCF FILTER AND INFO from IN="+super.IN));
		while((line=in.readLine())!=null)
			{
			VariantContext ctx=codeIn.decode(line);
			w.add(ctx);
			}
		in.close();
		tmp.delete();
		}
	
	
	
	
	public static void main(String[] args)
		{
		new FixVCF().instanceMainWithExit(args);
		}
	}
