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
package com.github.lindenb.jvarkit.tools.biostar;

import java.nio.file.Path;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.iterator.EqualIterator;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.tools.pcr.ReadClipper;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

/**
 * 
 BEGIN_DOC

## Description

input must be sorted on read name using `samtools sort -n` or ` samtools collate`

## Example

```
samtools collate -O input.bam| java -jar dist/biostar480685.jar
```

 END_DOC
 *
 */

@Program(name="biostar480685",
keywords={"sam","bam","clip"},
description="paired-end bam clip bases outside insert range",
biostars= {480685},
creationDate="20201223",
modificationDate="20200220"
)
public class Biostar480685 extends Launcher {			
private static final Logger LOG = Logger.build(Biostar480685.class).make();
@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
private Path outputFile = null;
@Parameter(names= {"-R","--reference"},description=CRAM_INDEXED_REFENCE)
private Path faidx;
@ParametersDelegate
private WritingBamArgs writingBamArgs=new WritingBamArgs();


@Override
public int doWork(final List<String> args) {
	SamReader in=null;
	SAMFileWriter out=null;
	try {
		final SamReaderFactory srf = super.createSamReaderFactory();
		if(this.faidx!=null) {
			srf.referenceSequence(this.faidx);
			writingBamArgs.setReferencePath(this.faidx);
			}
		final String input = oneFileOrNull(args);
		if(input==null) {
			in = srf.open(SamInputResource.of(stdin()));
			}
		else
			{
			in = srf.open(SamInputResource.of(input));
			}
		final SAMFileHeader header = in.getFileHeader();
		if(!(header.getSortOrder().equals(SAMFileHeader.SortOrder.unsorted) || header.getSortOrder().equals(SAMFileHeader.SortOrder.queryname))) {
			LOG.error("input should be sorted with 'samtools sort -n' or 'samtools collate' but got " + header.getSortOrder());
			return -1;
			}
		final ReadClipper clipper = new ReadClipper();
		header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
		final SAMProgramRecord prg =  header.createProgramRecord();
		prg.setCommandLine(this.getProgramCommandLine());
		prg.setProgramName(this.getProgramName());
		prg.setProgramVersion(this.getGitHash());		
		JVarkitVersion.getInstance().addMetaData(this, header);
		out = this.writingBamArgs.openSamWriter(this.outputFile, header, true);
		try(CloseableIterator<List<SAMRecord>> iter = new EqualIterator<>(
				in.iterator(),
				(A,B)->A.getReadName().compareTo(B.getReadName()))
				) {
			while(iter.hasNext()) {
				final List<SAMRecord> buffer = iter.next();
				int read1_idx = -1;
				int read2_idx = -1;
				for(int i=0;i< buffer.size();i++) {
					final SAMRecord rec = buffer.get(i);
					if(!rec.getReadPairedFlag()) continue;
					if(rec.getReadUnmappedFlag()) continue;
					if(rec.getMateUnmappedFlag()) continue;
					if(rec.isSecondaryOrSupplementary()) continue;
					if(rec.getFirstOfPairFlag()) {
						read1_idx = i;
					}
					else if(rec.getSecondOfPairFlag()) {
						read2_idx = i;
					}
				}
				
				if(read1_idx==-1 || read2_idx==-1 || read1_idx==read2_idx)  continue;
				final SAMRecord rec1a = buffer.get(read1_idx);
				final SAMRecord rec2a = buffer.get(read2_idx);
				if(!rec1a.overlaps(rec2a)) continue;
				final int chromStart = Math.max(rec1a.getStart(),rec2a.getStart());
				final int chromEnd  = Math.min(rec1a.getEnd(),rec2a.getEnd());
				if(chromStart > chromEnd)  continue;
				
				final SimpleInterval rgn = new SimpleInterval(rec1a.getContig(), chromStart, chromEnd);
				final SAMRecord rec1b = clipper.clip(rec1a, rgn);
				if(rec1b==null || rec1b.getReadUnmappedFlag()) continue;
				final SAMRecord rec2b = clipper.clip(rec2a, rgn);
				if(rec2b==null || rec2b.getReadUnmappedFlag()) continue;
				
				rec1b.setAttribute("PG", prg.getId());
				rec2b.setAttribute("PG", prg.getId());
				
				
				rec1b.setAlignmentStart(chromStart);
				rec1b.setMateAlignmentStart(rec2b.getAlignmentStart());
				rec2b.setAlignmentStart(chromStart);
				rec2b.setMateAlignmentStart(rec1b.getAlignmentStart());
				
				rec1b.setAttribute("MC", rec2b.getCigarString());
				rec2b.setAttribute("MC", rec1b.getCigarString());
				rec1b.setAttribute("NM", null);
				rec2b.setAttribute("NM", null);
	
				
				buffer.set(read1_idx, rec1b);
				buffer.set(read2_idx, rec2b);
				for(SAMRecord rec: buffer) {
					out.addAlignment(rec);
					}
				}
		}
		in.close();
		in=null;
		out.close();
		out=null;
		
		return 0;
	} catch(final Throwable err) {
		LOG.error(err);
		return -1;
	} finally {
		CloserUtil.close(in);
		CloserUtil.close(out);
	}
}

public static void main(final String[] args) {
	new Biostar480685().instanceMainWithExit(args);
}
}
