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
* 2016

*/
package com.github.lindenb.jvarkit.tools.biostar;

import java.io.OutputStream;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloserUtil;


/**

BEGIN_DOC

### Example

```
java -jar dist/biostar178713.jar -d 100000 -o out.zip in1.bed in2.bed 
```

END_DOC
*/


@Program(name="biostar178713",
	description="split bed file into several bed files where each region is separated of any other by N bases",
	biostars=178713,
	keywords="bed",
	creationDate="20160226",
	modificationDate="20200818"
	)
public class Biostar178713 extends Launcher
	{
	private static final Logger LOG = Logger.build(Biostar178713.class).make();


	@Parameter(names={"-o","--output"},description="Output file.zip .",required=true)
	private Path outputFile = null;

	@Parameter(names={"-d","--distance"},description="Distance between bed features." + DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=com.github.lindenb.jvarkit.util.jcommander.NoSplitter.class)
	private int distancebed = 100 ;

	
	@Override
	public int doWork(final List<String> args) {
		if(this.outputFile==null || !outputFile.getFileName().toString().endsWith(".zip")) {
			LOG.error("output file option  must be declared and must en with .zip");
			return -1;
		}
		
		final Set<String> inputs = IOUtils.unrollFiles(args);
		final List<BedLine> bedLines=new ArrayList<>();
		OutputStream fos = null;
		ZipOutputStream zout=null;

		try {
			if(inputs.isEmpty()) 
				{
				try(BedLineReader r = new BedLineReader(stdin(),"stdin")) {
					bedLines.addAll(r.stream().collect(Collectors.toList()));
					}
				}
			else for(final String input:inputs)
				{
				try(BedLineReader r = new BedLineReader(IOUtils.openURIForBufferedReading(input),input)) {
					bedLines.addAll(r.stream().collect(Collectors.toList()));
					}
				}
			LOG.info("sorting "+bedLines.size());
			Collections.sort(bedLines,new Comparator<BedLine>() {
				@Override
				public int compare(BedLine o1, BedLine o2) {
					int i= o1.getContig().compareTo(o2.getContig());
					if(i!=0) return i;
					i = o1.getStart() - o2.getStart();
					if(i!=0) return i;
					i = o1.getEnd() - o2.getEnd();
					return i;
					}
			});
			
			if(bedLines.isEmpty())
				{
				LOG.error("no bed line found");
				return -1;
				}
			
			LOG.info("creating zip "+this.outputFile);
			fos = Files.newOutputStream(this.outputFile);
			zout = new ZipOutputStream(fos);
			int chunk=0;
			while(!bedLines.isEmpty())
				{
				++chunk;
				
				final ZipEntry entry = new ZipEntry("bed"+String.format("%03d", chunk)+".bed");
				LOG.info("creating "+entry.getName());
				zout.putNextEntry(entry);
				PrintWriter pw = new PrintWriter(zout);
				BedLine prev= bedLines.get(0);
				bedLines.remove(0);
				pw.println(prev.join());
				int n_in_entry=1;
				int i=0;
				while(i<bedLines.size())
					{
					final BedLine curr = bedLines.get(i);
					final double distance  = curr.getStart()- prev.getEnd();
					if( !prev.getContig().equals(curr.getContig()) ||
						distance > this.distancebed) {
						pw.println(curr.join());
						prev=curr;
						bedLines.remove(i);
						++n_in_entry;
						}
					else
						{
						++i;
						}
					}
				pw.flush();
				zout.closeEntry();
				LOG.info("closing "+entry.getName()+" N="+n_in_entry);
				}
			
			zout.finish();
			zout.close();
			return 0;
		} catch (Throwable e) {
			LOG.error(e);
			return -1;
		} finally {
			CloserUtil.close(fos);
		}
	}
	public static void main(String[] args) {
		new Biostar178713().instanceMainWithExit(args);
	}

}
