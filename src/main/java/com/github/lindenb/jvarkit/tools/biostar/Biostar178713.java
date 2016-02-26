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
* 2016

*/
package com.github.lindenb.jvarkit.tools.biostar;

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

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.tribble.readers.LineIterator;



public class Biostar178713 extends AbstractBiostar178713
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(Biostar178713.class);

	@Override
	public Collection<Throwable> initializeKnime() {
		if(super.outputFile==null || !outputFile.getName().endsWith(".zip")) {
			return wrapException("output file option -"+OPTION_OUTPUTFILE+" must be declared and must en with .zip");
		}
		return super.initializeKnime();
		}
	
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
	public Collection<Throwable> call() throws Exception {
		final Set<String> inputs = IOUtils.unrollFiles(this.getInputFiles());
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
				return wrapException("no bed line found");
				}
			
			LOG.info("creating zip "+super.outputFile);
			fos = new FileOutputStream(super.outputFile);
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
						distance > super.distancebed) {
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
		return wrapException(e);
		} finally {
			CloserUtil.close(fos);
		}
	}
	public static void main(String[] args) {
		new Biostar178713().instanceMainWithExit(args);
	}

}
