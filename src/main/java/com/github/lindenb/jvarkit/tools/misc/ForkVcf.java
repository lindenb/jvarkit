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
* 2015 moving to knime

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFBuffer;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;


public class ForkVcf
	extends AbstractForkVcf
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(ForkVcf.class);

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
					ForkVcf.this.getOutputFile().getParentFile(),
					ForkVcf.this.getOutputFile().getName().replaceAll(
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
	protected Collection<Throwable> call(final String inputName) throws Exception {
		if (getOutputFile()==null || !getOutputFile().getName().contains(REPLACE_GROUPID)) {
			return wrapException("Output file pattern undefined or doesn't contain " + REPLACE_GROUPID + " : "
					+ this.getOutputFile());
		}
		if (!(getOutputFile().getName().endsWith(".vcf") || getOutputFile().getName().endsWith(".vcf.gz"))) {
			return wrapException("output file must end with '.vcf' or '.vcf.gz'");
		}
		
		if(super.number_of_files<=0) {
			return wrapException("Bad value for -"+OPTION_NUMBER_OF_FILES+":"+this.number_of_files);
		}
		
		BufferedReader r=null;
		VcfIterator in =null;
		PrintWriter manifestWriter=null;
		final List<SplitGroup> groups = new ArrayList<>();
		VCFBuffer vcfBuffer=null;
		try 
			{
			in = openVcfIterator(inputName);
			manifestWriter = (super.manifestFile==null?
					new PrintWriter(new NullOuputStream()):
					IOUtils.openFileForPrintWriter(super.manifestFile)
					);
			final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(in.getHeader());
			
			
			if (!super.split_by_chunk)
				{
				while(groups.size()<super.number_of_files)
					{
					final SplitGroup sg = new SplitGroup(groups.size()+1);
					sg.open(in.getHeader());
					manifestWriter.println(sg.getFile().getPath());
					groups.add(sg);
					}
				int idx=0;
				while(in.hasNext()) {
					final VariantContext ctx = progress.watch(in.next());
					groups.get(idx%super.number_of_files)._writer.add(ctx);
					++idx;
					}
				in.close();
				} 
			else {
				long count_variants=0;
				vcfBuffer=new VCFBuffer(super.getMaxRecordsInRam(), super.getTmpdir());
				vcfBuffer.writeHeader(in.getHeader());
				while(in.hasNext()) {
					final VariantContext ctx = progress.watch(in.next());
					vcfBuffer.add(ctx);
					++count_variants;
					}
				in.close();
				final long variant_per_file= Math.max(1L,(long)Math.ceil(count_variants/(double)super.number_of_files));

				LOG.info("done buffering. n="+count_variants+" now forking "+variant_per_file+" variants for "+super.number_of_files+" files.");
				VcfIterator iter2=vcfBuffer.iterator();
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
				while(groups.size()<super.number_of_files)
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
		catch(Exception err) {
			for (final SplitGroup g : groups) {
				CloserUtil.close(g);
				if(in!=null) g.getFile().delete();
			}
			return wrapException(err);
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
