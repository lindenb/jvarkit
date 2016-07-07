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

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.PrintWriter;
import java.util.Collection;

import com.github.lindenb.jvarkit.util.samtools.SamJsonWriterFactory;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.util.CloserUtil;


public class SamToJson extends AbstractSamToJson {
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(SamToJson.class);

	@Override
	protected Collection<Throwable> call(final String inputName) throws Exception {
		PrintWriter out=null;
		SamReader sfr=null;
		final SamJsonWriterFactory factory =SamJsonWriterFactory.newInstance().
				printHeader(super.print_header).
				printReadName(!super.disable_readName).
				printAttributes(!super.disable_atts).
				expandFlag(super.expflag).
				expandCigar(super.excigar)
				;
		SAMFileWriter swf=null;
		try
			{
			sfr = super.openSamReader(inputName);
			out = super.openFileOrStdoutAsPrintWriter();
			swf = factory.open(sfr.getFileHeader(), out);
			final SAMRecordIterator iter=sfr.iterator();
			while(iter.hasNext() && !out.checkError())
				{
				swf.addAlignment(iter.next());
				}
			iter.close();
			out.flush();
			swf.close();
			LOG.info("done");
			return RETURN_OK;
			}
		catch(final Exception err) {
			return wrapException(err);
			}
		finally 	{
			CloserUtil.close(sfr);
			}
		}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new SamToJson().instanceMainWithExit(args);

	}

}
