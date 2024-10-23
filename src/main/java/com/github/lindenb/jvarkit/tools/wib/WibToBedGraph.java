/*
The MIT License (MIT)

Copyright (c) 2024 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.wib;

import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.samtools.util.IntervalParser;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.wib.WibReader;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;

/**
BEGIN_DOC

Example:

```
$ java -jar dist/jvarkit.jar wib2bedgraph  --tabix ~/phastCons100way.txt.gz --wib ~/phastCons100way.wib -r "chr2:1-125168903" | head

chr2	11391	11392	0.14969291
chr2	11392	11393	0.14667717
chr2	11393	11394	0.14064567
chr2	11394	11395	0.13461417
chr2	11395	11396	0.12858267
chr2	11396	11397	0.13461417
chr2	11397	11398	0.14366141
chr2	11398	11399	0.14969291
chr2	11399	11400	0.1557244
chr2	11400	11401	0.15874015

$ java -jar  dist/jvarkit.jar wib2bedgraph  --tabix ~/phastCons100way.txt.gz --wib ~/phastCons100way.wib -r "chr2:1-125168903" --wig | head

variableStep chrom=chr2 span=1
11392	0.14969291
11393	0.14667717
11394	0.14064567
11395	0.13461417
11396	0.12858267
11397	0.13461417
11398	0.14366141
11399	0.14969291
11400	0.1557244
```
 
END_DOC
*/
@Program(name="wib2bedgraph",
description="Extract Wib files to bedgraph or wig",
keywords={"wib","wig","bed"},
creationDate="20230819",
modificationDate="20230819",
jvarkit_amalgamion = true
)
public class WibToBedGraph extends Launcher {
	private static final Logger LOG = Logger.build(WibToBedGraph.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputPath = null;
	@Parameter(names={"-r","--regions"},description=IntervalParser.OPT_DESC,required=true)
	private String region="";
	@Parameter(names={"--tabix"},description=WibReader.TABIX_DESC,required=true)
	private String tabixURI="";
	@Parameter(names={"--wib"},description=WibReader.WIB_DESC,required=true)
	private String wibURI="";
	@Parameter(names={"--wig"},description="output WIG format instead of BED graph")
	private boolean output_wig = false;

	@Override
	public int doWork(List<String> args) {
		if(!args.isEmpty()) {
			LOG.error("too many parameters");
			return -1;
			}
	
		try(WibReader reader=new  WibReader(this.tabixURI, this.wibURI)) {
			final Locatable loc =new  IntervalParser().apply(this.region).orElseThrow(IntervalParser.exception(this.region));
			try(CloseableIterator<WibReader.WigItem> iter = reader.query(loc)) {
				try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(outputPath)) {
					if(output_wig) {
						pw.println("variableStep chrom=" + loc.getContig()+" span=1");
						while(iter.hasNext()) {
							final WibReader.WigItem item = iter.next();
							for(int x=item.getStart(); x<=item.getEnd();++x) {
								pw.print(x);
								pw.print("\t");
								pw.println((float)item.getValue());
								}
							}
						}
					else {
						while(iter.hasNext()) {
							final WibReader.WigItem item = iter.next();
							pw.print(item.getContig());
							pw.print('\t');
							pw.print(item.getStart()-1);
							pw.print('\t');
							pw.print(item.getEnd());
							pw.print('\t');
							pw.println((float)item.getValue());
							}
						}
					pw.flush();
					}
				}
			return 0;
			}
		catch(Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	
public static void main(String[] args) {
	new WibToBedGraph().instanceMainWithExit(args);
	}
}
