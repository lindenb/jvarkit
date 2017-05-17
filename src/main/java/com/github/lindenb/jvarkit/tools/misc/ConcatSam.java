/*
The MIT License (MIT)

Copyright (c) 2016 Pierre Lindenbaum

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
* 2016 creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamFileHeaderMerger;
import htsjdk.samtools.util.CloserUtil;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;


/**

BEGIN_DOC


END_DOC

 */
@Program(name="concatsam",description="concat sam files",keywords={"sam","bam"})
public class ConcatSam extends Launcher
	{
	private static final Logger LOG = Logger.build(ConcatSam.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@ParametersDelegate
	private WritingBamArgs writingBamArgs = new WritingBamArgs();

	@Override
	public int doWork(List<String> args) {
	
		SAMFileWriter out=null;
		try
			{
			final Set<String> set = IOUtils.unrollFiles(args);
			if(set.isEmpty()) {
				LOG.info("Reading from stdin");
				final SamReader r = openSamReader(null);
				out = this.writingBamArgs.openSAMFileWriter(outputFile,r.getFileHeader(), true);
				final SAMRecordIterator iter = r.iterator();
				while(iter.hasNext()) out.addAlignment(iter.next());
				iter.close();
				out.close();
				r.close();
				}
			else
				{
				final List<SAMFileHeader> headers=new ArrayList<>();
				for(final String fname:set)
					{
					LOG.info("Reading header of " + fname);
					final SamReader r = openSamReader(fname);
					headers.add(r.getFileHeader());
					r.close();
					}
				
				final SAMFileHeader header= new SamFileHeaderMerger(
						SAMFileHeader.SortOrder.unsorted,
						headers,
						false
						).getMergedHeader();
				out = this.writingBamArgs.openSAMFileWriter(outputFile,header, true);
				for(final String fname:set)
					{
					LOG.info("Reading bam " + fname);
					final SamReader r = openSamReader(fname);
					final SAMRecordIterator iter = r.iterator();
					while(iter.hasNext()) out.addAlignment(iter.next());
					iter.close();
					r.close();
					}
				out.close();
				}
			
			
			return RETURN_OK;
			}
		catch(Exception err) 
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(out);
			}
		}
	

	public static void main(String[] args) {
		new ConcatSam().instanceMainWithExit(args);
		}
	}
