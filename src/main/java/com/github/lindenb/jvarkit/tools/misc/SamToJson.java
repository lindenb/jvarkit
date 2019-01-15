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

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.PrintWriter;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.SamJsonWriterFactory;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.util.CloserUtil;
/**
BEGIN_DOC


END_DOC
 */
@Program(name="sam2json",keywords={"sam","bam","json"},description="Convert a SAM input to JSON")
public class SamToJson extends Launcher {
	private static final Logger LOG = Logger.build(SamToJson.class).make();
	
	@Parameter(names={"-H","--header"},description="don't print SAM HEADER")
	private boolean print_header = false;

	@Parameter(names={"-name","--name"},description="do not print read name")
	private boolean disable_readName = false;

	@Parameter(names={"-atts","--atts"},description="do not print attributes")
	private boolean disable_atts = false;

	@Parameter(names={"-flag","--flag"},description="expand SAm Flags")
	private boolean expflag = false;

	@Parameter(names={"-cigar","--cigar"},description="expand cigar")
	private boolean excigar = false;

	@Parameter(names={"-out","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File output = null;

	
	private int call(final String inputName) throws Exception {
		PrintWriter out=null;
		SamReader sfr=null;
		final SamJsonWriterFactory factory =SamJsonWriterFactory.newInstance().
				printHeader(this.print_header).
				printReadName(!this.disable_readName).
				printAttributes(!this.disable_atts).
				expandFlag(this.expflag).
				expandCigar(this.excigar)
				;
		SAMFileWriter swf=null;
		try
			{
			sfr = super.openSamReader(inputName);
			out = super.openFileOrStdoutAsPrintWriter(this.output);
			swf = factory.open(sfr.getFileHeader(), out);
			final SAMRecordIterator iter=sfr.iterator();
			while(iter.hasNext() && !out.checkError())
				{
				swf.addAlignment(iter.next());
				}
			iter.close();
			out.flush();
			swf.close();swf=null;
			LOG.info("done");
			return 0;
			}
		catch(final Exception err) {
			LOG.error(err);
			return -1;
			}
		finally 	{
			CloserUtil.close(sfr);
			CloserUtil.close(swf);
			}
		}
	@Override
	public int doWork(List<String> args) {
		try {
			return call(oneFileOrNull(args));
		} catch (Exception e) {
			LOG.error(e);
			return -1;
		}
		}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new SamToJson().instanceMainWithExit(args);

	}

}
