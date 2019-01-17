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
* 2014 creation
* 2015 moving to knime

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFBuffer;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import htsjdk.variant.vcf.VCFIterator;

/**

BEGIN_DOC




### Output

Output filename (option -o) MUST contain the word __GROUPID__.



### Example



```
$ 

```






```


```

cat input.vcf | java -jar dist/forkvcf.jar -n 3 -o "_tmp.__GROUPID__.vcf"
[main] INFO jvarkit - opening VCF file "_tmp.00001.vcf" for writing
[main] INFO jvarkit - opening VCF file "_tmp.00002.vcf" for writing
[main] INFO jvarkit - opening VCF file "_tmp.00003.vcf" for writing

$ wc _tmp.0000*
   226   6819 143947 _tmp.00001.vcf
   226   6819 140792 _tmp.00002.vcf
   225   6161 125219 _tmp.00003.vcf
   
   
   


### See also


 *  https://github.com/lindenb/jvarkit/wiki/SplitVcf





END_DOC
*/


@Program(name="forkvcf",description="Fork a VCF.")
public class ForkVcf
	extends Launcher
	{
	private static final Logger LOG = Logger.build(ForkVcf.class).make();


	@Parameter(names={"-o","--output"},description="Output file Must contains "+REPLACE_GROUPID,required=true)
	private File outputFile = null;


	@Parameter(names={"-n","--count"},description="number of vcf files to generate")
	private int number_of_files = 2 ;

	@Parameter(names={"-c","--splitbychunk"},description="When this option is used, the variant are first saved in a temporary file, the number of variant is dividided by 'count' and the output files are lineray produced. The default is to dispatch the variants as they are coming in the stream.")
	private boolean split_by_chunk = false;

	@Parameter(names={"-m","--manifest"},description="optional save produced vcf filenames in this file.")
	private File manifestFile = null;


    @Parameter(names={"-T","--tmpDir"},description="mp directory")
    private File tmpDir = IOUtils.getDefaultTmpDir();

    @Parameter(names={"-maxRecordsInRam","--maxRecordsInRam"},description="Max records in RAM")
    private int maxRecordsInRam =50000;

	private final static String REPLACE_GROUPID="__GROUPID__";

	
	/** spit group */
	private class SplitGroup implements Closeable
		{
		final int groupId;
		VCFHeader header=null;
		VariantContextWriter _writer;
		
		SplitGroup(final int groupId)
			{
			this.groupId=groupId;
			}
		

		@Override
		public void close() {
			if(_writer!=null) CloserUtil.close(_writer);
			_writer=null;
			
		}
		

		
		public File getFile()
			{
			return new File(
					ForkVcf.this.outputFile.getParentFile(),
					ForkVcf.this.outputFile.getName().replaceAll(
							ForkVcf.REPLACE_GROUPID,
							String.valueOf(this.groupId)
							));
			}
		
		public void open(final VCFHeader src)
			{	
	
			final File fileout=getFile();
			LOG.info("opening VCF file \""+fileout+"\" for writing");
			final File parent=fileout.getParentFile();
			if(parent!=null) {
				parent.mkdirs();
			}
	
			this.header= new VCFHeader(src);
			this.header.addMetaDataLine(new VCFHeaderLine("ForkVcf.GroupId", String.valueOf(this.groupId)));
			try {
				this._writer = VCFUtils.createVariantContextWriter(fileout);
			} catch (IOException e) {
				throw new RuntimeIOException(e);
			}
			this._writer.writeHeader(this.header);
			}
		
		@Override
		public String toString() {
			return getFile().getPath();
			}
		}
	
	
		
	public ForkVcf()
		{
		}
	
	
	@Override
	public int doWork(List<String> args) {
	if (this.outputFile==null || !this.outputFile.getName().contains(REPLACE_GROUPID)) {
			LOG.error("Output file pattern undefined or doesn't contain " + REPLACE_GROUPID + " : "
					+ this.outputFile);
			return -1;
		}
		if (!(this.outputFile.getName().endsWith(".vcf") || this.outputFile.getName().endsWith(".vcf.gz"))) {
			LOG.error("output file must end with '.vcf' or '.vcf.gz'");
			return -1;
		}
		
		if(this.number_of_files<=0) {
			LOG.error("Bad value for number of files:"+this.number_of_files);
			return -1;
		}
		
		BufferedReader r=null;
		VCFIterator in =null;
		PrintWriter manifestWriter=null;
		final List<SplitGroup> groups = new ArrayList<>();
		VCFBuffer vcfBuffer=null;
		try 
			{
			in = openVCFIterator(oneFileOrNull(args));
			manifestWriter = (this.manifestFile==null?
					new PrintWriter(new NullOuputStream()):
					IOUtils.openFileForPrintWriter(this.manifestFile)
					);
			final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(in.getHeader());
			
			
			if (!this.split_by_chunk)
				{
				while(groups.size()<this.number_of_files)
					{
					final SplitGroup sg = new SplitGroup(groups.size()+1);
					sg.open(in.getHeader());
					manifestWriter.println(sg.getFile().getPath());
					groups.add(sg);
					}
				int idx=0;
				while(in.hasNext()) {
					final VariantContext ctx = progress.watch(in.next());
					groups.get(idx%this.number_of_files)._writer.add(ctx);
					++idx;
					}
				in.close();
				} 
			else {
				long count_variants=0;
				vcfBuffer=new VCFBuffer(this.maxRecordsInRam,this.tmpDir);
				vcfBuffer.writeHeader(in.getHeader());
				while(in.hasNext()) {
					final VariantContext ctx = progress.watch(in.next());
					vcfBuffer.add(ctx);
					++count_variants;
					}
				in.close();
				final long variant_per_file= Math.max(1L,(long)Math.ceil(count_variants/(double)this.number_of_files));

				LOG.info("done buffering. n="+count_variants+" now forking "+variant_per_file+" variants for "+this.number_of_files+" files.");
				VCFIterator iter2=vcfBuffer.iterator();
				long count_ctx=0L;
				while(iter2.hasNext()) {
					if(groups.isEmpty() || count_ctx>=variant_per_file)
						{
						if(!groups.isEmpty()) groups.get(groups.size()-1).close();
						final SplitGroup last = new SplitGroup(groups.size()+1);
						last.open(in.getHeader());
						manifestWriter.println(last.getFile().getPath());
						groups.add(last);
						count_ctx=0;
						}
					
					final VariantContext ctx = iter2.next();
					groups.get(groups.size()-1)._writer.add(ctx);
					count_ctx++;
					}
				iter2.close();
				vcfBuffer.close();
				vcfBuffer.dispose();
				vcfBuffer=null;
				
				//save remaining empty VCFs
				while(groups.size()<this.number_of_files)
					{
					LOG.info("creating empty vcf");
					final SplitGroup sg = new SplitGroup(groups.size()+1);
					sg.open(in.getHeader());
					manifestWriter.println(sg.getFile().getPath());
					sg.close();
					groups.add(sg);
					}
				}
			
			progress.finish();
			
	
			for (final SplitGroup g : groups) {
				g.close();
			}
			
			manifestWriter.flush();
			manifestWriter.close();
			manifestWriter=null;
			return RETURN_OK;
			}
		catch(final Exception err) {
			LOG.error(err);
			for (final SplitGroup g : groups) {
				CloserUtil.close(g);
				if(in!=null) g.getFile().delete();
			}
			return -1;
		} finally {
			if(vcfBuffer!=null) vcfBuffer.dispose();
			CloserUtil.close(r);
			CloserUtil.close(in);
			IOUtils.flush(manifestWriter);
			CloserUtil.close(manifestWriter);
		}
	}
	 	
	
	public static void main(String[] args)
		{
		new ForkVcf().instanceMainWithExit(args);
		}

	}
