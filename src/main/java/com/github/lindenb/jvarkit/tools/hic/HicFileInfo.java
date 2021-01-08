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


*/
package com.github.lindenb.jvarkit.tools.hic;

import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.hic.HicReader;
import com.github.lindenb.jvarkit.hic.HicReaderFactory;
import com.github.lindenb.jvarkit.io.CustomSeekableStreamFactory;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.seekablestream.ISeekableStreamFactory;
import htsjdk.samtools.util.CloserUtil;
/**
BEGIN_DOC

### Example

```
$ java -jar dist/hicfileinfo.jar   "https://s3.amazonaws.com/hicfiles/external/goodell/tcell.hic" 
SAMSequenceDictionary:( sequences:26 length:3098789674  md5:d0494b28329e1e20d05d74611b913fa6)
# https://s3.amazonaws.com/hicfiles/external/goodell/tcell.hic
 HiC Version           : 8
 Build                 : hg19
 Base Pair Resolutions : 5000 10000 25000 50000 100000 250000 500000 1000000 2500000
 Fragments Resolutions : 1 2 5 20 50 100 200 500


```

END_DOC
 */
@Program(name="hicfileinfo",
	description="Prints information about a HI-C File/URL",
	keywords={"hic"},
	creationDate="20190606",
	modificationDate="20190606"
	)
public class HicFileInfo extends Launcher {
	private static final Logger LOG = Logger.build(HicFileInfo.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-A","--attributes"},description="print attributes")
	private boolean print_attributes = false;
	@Parameter(names={"-D","--dict"},description="print dictionary")
	private boolean print_dict = false;

	@Override
	public int doWork(final List<String> args) {
		PrintWriter out=null;
		try {
			final ISeekableStreamFactory seekableStreamFactory = new CustomSeekableStreamFactory().
					setUsingHttpHead(false).
					setNormalizeURI(false)
					;
			
			final HicReaderFactory hrf = new HicReaderFactory();
			hrf.setSeekableStreamFactory(seekableStreamFactory);
			
			out = super.openPathOrStdoutAsPrintWriter(this.outputFile);
			for(final String path: args) {
				try(final HicReader r = hrf.open(path)) {
					out.println("# "+ path);
					out.println(" HiC Version           : "+r.getVersion());
					out.println(" Build                 : "+r.getBuild());
					out.println(" Base Pair Resolutions : " + r.getBasePairResolutions().stream().sorted().map(I->String.valueOf(I)).collect(Collectors.joining(" ")));
					out.println(" Fragments Resolutions : " + r.getFragmentResolutions().stream().sorted().map(I->String.valueOf(I)).collect(Collectors.joining(" ")));
					
					if ( this.print_attributes ) {
						out.println(" Attributes : ");
						for(final Map.Entry<String, String> e : r.getAttributes().entrySet()) {
							out.println("   "+ e.getKey() +" : " + e.getValue());
							}
						}
					
					if( this.print_dict) {
						out.println(" Dictionary           : ");
						for(final SAMSequenceRecord ssr : r.getDictionary().getSequences()) {
							out.printf("   %2d %5s %7d\n",ssr.getSequenceIndex(),ssr.getSequenceName(),ssr.getSequenceLength());
							}
						}
					out.println();
					}
				}
			
			out.flush();
			out.close();
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally {
			CloserUtil.close(out);
			}
		}
	
	public static void main(String[] args) {
		new HicFileInfo().instanceMainWithExit(args);
	}

}
