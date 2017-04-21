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

import java.io.File;
import java.io.FileOutputStream;

import java.util.ArrayList;
import java.util.List;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

@Program(name="vcf2zip",description="Reads a stream of concatenated VCFs and insert them into a Zip file")
public class VcfToZip extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfToZip.class).make();


	@Parameter(names={"-o","--output"},description="Output zip file.")
	private File outputFile = null;


	@Parameter(names={"-p","--prefix"},description="Prefix all zip entries with this prefix")
	private String zipPrefix = "VCF";

	@Parameter(names={"-t","--title"},description="Try to find ##(TITLE)=abcdefghijk in the VCF header and use it as the name of the inserted VCF file")
	private String titleHeaderStr = "";

	public VcfToZip()
		{
		}
	
	@Override
	public int doWork(List<String> args) {
		
		if(this.outputFile!=null && this.outputFile.getName().endsWith(".zip")) {
			LOG.error("Filename must end with '.zip' "+outputFile);
			return -1;
		}
		LineIterator lr=null;
		VcfIterator in = null;
		ZipOutputStream zout =null;
		FileOutputStream fout=null;
		int num_vcfs= 0;
		VariantContextWriter vcw = null;
		args = new ArrayList<>(IOUtils.unrollFiles(args));
		try {
			int optind=0;
			
			if(this.outputFile!=null) {
				fout = new FileOutputStream(this.outputFile);
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
					if(this.titleHeaderStr!=null && !this.titleHeaderStr.isEmpty()) {
						final VCFHeaderLine h = header.getOtherHeaderLine(this.titleHeaderStr);
						if(h!=null && !h.getValue().trim().isEmpty()) filename=h.getValue().trim();
					}	
					if(filename==null || filename.trim().isEmpty()) {
						//create title
						filename= String.format("vcf2zip.%05d.vcf", num_vcfs);
						//set name in header
						if(this.titleHeaderStr!=null && !this.titleHeaderStr.isEmpty()) {
							header.addMetaDataLine(new VCFHeaderLine(this.titleHeaderStr.trim(),filename));
						}
					}
					if(!filename.endsWith(".vcf")) {
						filename+=".vcf";
					}
					if(!this.zipPrefix.isEmpty())
						{
						filename= this.zipPrefix+
								(this.zipPrefix.endsWith("/")?"":"/")+
								filename;
						}
					LOG.info(filename);
					final ZipEntry entry = new ZipEntry(filename);
					entry.setComment("Created with "+getProgramName());
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
			LOG.error(e);
			return -1;
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
