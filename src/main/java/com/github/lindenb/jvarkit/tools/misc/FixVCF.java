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

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.io.SequenceInputStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.LineReaderUtil;
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
import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class FixVCF
	extends AbstractFixVCF
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(FixVCF.class);

	
	
	@Override
	public  Command createCommand() {
			return new MyCommand();
		}
		 
	public  static class MyCommand extends AbstractFixVCF.AbstractFixVCFCommand
	 	{		

	@Override
		protected Collection<Throwable> call(String inputName) throws Exception {
		
		if(super.tmpDir==null)
			{
			super.tmpDir=new File(System.getProperty("java.io.tmpdir"));
			}
		
		VariantContextWriter w=null;
		InputStream in = null;
		try
			{
			w = openVariantContextWriter();
			
			if(inputName==null)
				{
				in= stdin();
				doWork("stdin",in,w);
				}
			else
				{
				in=IOUtils.openURIForReading(inputName);
				doWork(inputName,in,w);
				}
			w.close();
			in.close();
			return RETURN_OK;
			}
		catch(Exception err)
			{
			return wrapException(err);
			}
		finally
			{
			CloserUtil.close(in);
			CloserUtil.close(w);
			}
		}
	
	private void doWork(
			String filenameIn,
			InputStream vcfStream,
			VariantContextWriter w)
			throws IOException
			{
			
			AbstractVCFCodec vcfCodec = VCFUtils.createDefaultVCFCodec();
			
			LineIterator r= new LineIteratorImpl(LineReaderUtil.fromBufferedStream(vcfStream));
			VCFHeader header=(VCFHeader) vcfCodec.readActualHeader(r);
			
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
			
			File tmp=IOUtil.newTempFile("tmp", ".vcf.gz",new File[]{super.tmpDir});
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
					LOG.error(line);
					LOG.error(err);
					throw err;
					}
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
			VcfIterator in=new VcfIterator(new SequenceInputStream(
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
			}
	
	 	}
	
	public static void main(String[] args)
		{
		new FixVCF().instanceMainWithExit(args);
		}
	}
