/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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
* 2015 creation

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

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;


/**
BEGIN_DOC

## Example

```bash
$ java -jar dist/fixvcf.jar < bad.vcf > ok.vcf
```


END_DOC

*/
@Program(name="fixvcf",description="Fix a VCF if INFO or FILTER are missing")
public class FixVCF
	extends Launcher
	{
	private static final Logger LOG = Logger.build(FixVCF.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-T","--tmpDir"},description="tmp directory")
	private File tmpDir = IOUtils.getDefaultTmpDir();
	

	@Override
	public int doWork(List<String> args) {
		VariantContextWriter w=null;
		InputStream in = null;
		try
			{
			String fname=super.oneFileOrNull(args);
			
			in = super.openInputStream(fname);
			w = super.openVariantContextWriter(this.outputFile);
			
			doWork((fname==null?"stdin":fname),in,w);
				
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(w);
			CloserUtil.close(in);
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
				LOG.error(line);
				LOG.error(err);
				return -1;
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
		VCFIterator in= new VCFIteratorBuilder().open(new SequenceInputStream(
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
