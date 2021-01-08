/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import htsjdk.variant.vcf.VCFIterator;

/**

BEGIN_DOC




### Motivation

This tool was used to create a zip from the output of VCFburdensplitter which is a stream of VCFs.



### Example



```

$ cat ~/input.vcf ~/input.vcf ~/input.vcf | java -jar dist/vcf2zip.jar -o jeter.zip
[main] INFO jvarkit - Command Line args : -o jeter.zip
[main] INFO jvarkit - Executing as lindenb@kaamelot-master01 on Linux 2.6.32-431.17.1.el6.x86_64 amd64; Java HotSpot(TM) 64-Bit Server VM 1.8.0_60-b27
[main] INFO jvarkit - reading concatenated vcf from stdin
[main] INFO jvarkit - VCF/vcf2zip.00001.vcf
[main] INFO jvarkit - Count: 499 Elapsed: 10 seconds(0.05%) Remains: 6 hours(99.95%) Last: 1:1431105
[main] INFO jvarkit - done: N=870
[main] INFO jvarkit - VCF/vcf2zip.00002.vcf
[main] INFO jvarkit - Count: 530 Elapsed: 10 seconds(0.05%) Remains: 5 hours(99.95%) Last: 1:1510577
[main] INFO jvarkit - done: N=870
[main] INFO jvarkit - VCF/vcf2zip.00003.vcf
[main] INFO jvarkit - Count: 530 Elapsed: 10 seconds(0.05%) Remains: 5 hours(99.95%) Last: 1:1510577
[main] INFO jvarkit - done: N=870
[main] INFO jvarkit - done. Number of VCFs:3
[main] INFO jvarkit - End JOB  [Mon May 02 12:30:24 CEST 2016] VcfToZip done. Elapsed time: 0.85 minutes.
$ unzip -t jeter.zip 
Archive:  jeter.zip
    testing: VCF/vcf2zip.00001.vcf    OK
    testing: VCF/vcf2zip.00002.vcf    OK
    testing: VCF/vcf2zip.00003.vcf    OK
No errors detected in compressed data of jeter.zip.

```





END_DOC
*/


@Program(name="vcf2zip",
description="Reads a stream of concatenated VCFs and insert them into a Zip file",
keywords= {"vcf","zip"},
modificationDate="20191121"
)
public class VcfToZip extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfToZip.class).make();


	@Parameter(names={"-o","--output"},description=ArchiveFactory.OPT_DESC,required=true)
	private Path outputFile = null;


	@Parameter(names={"-p","--prefix"},description="Prefix all zip entries with this prefix")
	private String zipPrefix = "VCF";

	@Parameter(names={"-t","--title"},description="Try to find ##(TITLE)=abcdefghijk in the VCF header and use it as the name of the inserted VCF file")
	private String titleHeaderStr = "";

	public VcfToZip()
		{
		}
	
	@Override
	public int doWork(List<String> args) {
		
	
		LineIterator lr=null;
		VCFIterator in = null;
		ArchiveFactory archive =null;
		FileOutputStream fout=null;
		int num_vcfs= 0;
		VariantContextWriter vcw = null;
		args = new ArrayList<>(IOUtils.unrollFiles(args));
		try {
			int optind=0;
			final Path tmpVcf = Files.createTempFile("tmp.",".vcf");
			
			do {
				
				if(args.isEmpty()) {
					lr = IOUtils.openStreamForLineIterator(stdin());
					}
				else
					{
					lr = IOUtils.openURIForLineIterator(args.get(optind));
					}
				archive  = ArchiveFactory.open(this.outputFile);
				
				while(lr.hasNext()) {
					++num_vcfs;
					in = VCFUtils.createVCFIteratorFromLineIterator(lr,true);
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
					vcw = VCFUtils.createVariantContextWriterToPath(tmpVcf);
					vcw.writeHeader(header);
					final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(header);
					while(in.hasNext()) {
						vcw.add(progress.watch(in.next()));
					}
					vcw.close();
					progress.finish();
					archive.copyTo(tmpVcf, filename);
					Files.delete(tmpVcf);
				}
				
				
				++optind;
			} while(optind< args.size());
			
			archive.close();archive=null;
			CloserUtil.close(fout);

			
			LOG.info("done. Number of VCFs:"+num_vcfs);
			return RETURN_OK;
		} catch (final Exception e) {
			LOG.error(e);
			return -1;
		} finally {
			CloserUtil.close(in);
			CloserUtil.close(lr);
			CloserUtil.close(archive);
			CloserUtil.close(fout);
		}
		}
	
	public static void main(String[] args) {
		new VcfToZip().instanceMainWithExit(args);
	}
	
	}
