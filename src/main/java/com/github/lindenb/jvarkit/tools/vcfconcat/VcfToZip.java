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
package com.github.lindenb.jvarkit.tools.vcfconcat;

import java.io.FileOutputStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class VcfToZip extends AbstractVcfToZip
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(VcfToZip.class);

	public VcfToZip()
		{
		}
	
	@Override
	public Collection<Throwable> call() throws Exception {
		
		if(getOutputFile()!=null && !getOutputFile().getName().endsWith(".zip")) {
			return wrapException("Filename must end with '.zip' "+getOutputFile());
		}
		LineIterator lr=null;
		VcfIterator in = null;
		ZipOutputStream zout =null;
		FileOutputStream fout=null;
		int num_vcfs= 0;
		VariantContextWriter vcw = null;
		List<String> args = new ArrayList<>(IOUtils.unrollFiles(super.getInputFiles()));
		try {
			int optind=0;
			
			if(getOutputFile()!=null) {
				fout = new FileOutputStream(getOutputFile());
				zout = new ZipOutputStream(fout);
			} else
				{
				zout = new ZipOutputStream(stdout());
				}

			
			do {
				
				if(args.isEmpty()) {
					LOG.info("reading concatenated vcf from stdin");
					lr = IOUtils.openStreamForLineIterator(stdin());
					}
				else
					{
					LOG.info("reading concatenated vcf from "+args.get(optind));
					lr = IOUtils.openURIForLineIterator(args.get(optind));
					}
				
				
				while(lr.hasNext()) {
					++num_vcfs;
					in = VCFUtils.createVcfIteratorFromLineIterator(lr,true);
					final VCFHeader header= in.getHeader();
					String filename=null;
					if(super.titleHeaderStr!=null && !super.titleHeaderStr.isEmpty()) {
						final VCFHeaderLine h = header.getOtherHeaderLine(super.titleHeaderStr);
						if(h!=null && !h.getValue().trim().isEmpty()) filename=h.getValue().trim();
					}	
					if(filename==null || filename.trim().isEmpty()) {
						//create title
						filename= String.format("vcf2zip.%05d.vcf", num_vcfs);
						//set name in header
						if(super.titleHeaderStr!=null && !super.titleHeaderStr.isEmpty()) {
							header.addMetaDataLine(new VCFHeaderLine(super.titleHeaderStr.trim(),filename));
						}
					}
					if(!filename.endsWith(".vcf")) {
						filename+=".vcf";
					}
					if(!super.zipPrefix.isEmpty())
						{
						filename= super.zipPrefix+
								(super.zipPrefix.endsWith("/")?"":"/")+
								filename;
						}
					LOG.info(filename);
					final ZipEntry entry = new ZipEntry(filename);
					entry.setComment("Created with "+getName());
					zout.putNextEntry(entry);
					vcw = VCFUtils.createVariantContextWriterToOutputStream(
							IOUtils.uncloseableOutputStream(zout)
							);
					vcw.writeHeader(header);
					final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(header);
					while(in.hasNext()) {
						vcw.add(progress.watch(in.next()));
					}
					vcw.close();
					progress.finish();
					zout.closeEntry();
					in.close();
				}
				
				
				++optind;
			} while(optind< args.size());
			
			zout.finish();
			zout.flush();
			zout.close();zout=null;
			CloserUtil.close(fout);

			
			LOG.info("done. Number of VCFs:"+num_vcfs);
			return RETURN_OK;
		} catch (final Exception e) {
			return wrapException(e);
		} finally {
			CloserUtil.close(in);
			CloserUtil.close(lr);
			CloserUtil.close(zout);
			CloserUtil.close(fout);
		}
		}
	
	public static void main(String[] args) {
		new VcfToZip().instanceMainWithExit(args);
	}
	
	}
