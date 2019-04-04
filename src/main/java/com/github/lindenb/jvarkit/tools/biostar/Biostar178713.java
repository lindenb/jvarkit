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
* 2016

*/
package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.tribble.readers.LineIterator;


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
	keywords="bed"
	)
public class Biostar178713 extends Launcher
	{
	private static final Logger LOG = Logger.build(Biostar178713.class).make();


	@Parameter(names={"-o","--output"},description="Output file.zip .",required=true)
	private File outputFile = null;

	@Parameter(names={"-d","--distance"},description="Distance between bed features." + DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=com.github.lindenb.jvarkit.util.jcommander.NoSplitter.class)
	private int distancebed = 100 ;

	
	
	private final void readBed( Collection<BedLine> bed,final LineIterator in) {
		final BedLineCodec codec = new BedLineCodec();
		codec.readActualHeader(in);
		while(in.hasNext()){
		final BedLine line = codec.decode(in);	
		if(line==null) continue;
		bed.add(line);
		}
		CloserUtil.close(in);
		}
	
	@Override
	public int doWork(final List<String> args) {
		if(this.outputFile==null || !outputFile.getName().endsWith(".zip")) {
			LOG.error("output file option  must be declared and must en with .zip");
			return -1;
		}
		
		final Set<String> inputs = IOUtils.unrollFiles(args);
		List<BedLine> bedLines=new ArrayList<>();
		FileOutputStream fos = null;
		ZipOutputStream zout=null;

		try {
			if(inputs.isEmpty()) 
				{
				LOG.info("reading bed from stdin");
				LineIterator r = IOUtils.openStreamForLineIterator(stdin());
				this.readBed(bedLines, r);
				CloserUtil.close(r);
				}
			else for(final String input:inputs)
				{
				LOG.info("reading bed from "+input);
				LineIterator r = IOUtils.openURIForLineIterator(input);
				this.readBed(bedLines, r);
				CloserUtil.close(r);
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
			fos = new FileOutputStream(this.outputFile);
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
			return RETURN_OK;
		} catch (Exception e) {
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
