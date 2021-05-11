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

import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;

import com.github.lindenb.jvarkit.jcommander.OnePassBamLauncher;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;

/** 

BEGIN_DOC

##Example

```bash
 java -jar dist/biostar9469733.jar --regions "RF01:1-1000" src/test/resources/S*.bam
```

END_DOC

*/

@Program(name="biostar9469733",
description="Extract reads mapped within chosen intronic region from BAM file",
keywords= {"sam","bam","rnaseq","bed"},
creationDate="20210511",
modificationDate="20210511",
biostars=9469733
)
public class Biostar9469733 extends OnePassBamLauncher
	{
	private static final Logger LOG = Logger.build(Biostar9469733.class).make();
	private IntervalTreeMap<Locatable> intervalTreeMap = null;
	@Override
	protected int beforeSam()
		{
		if(super.regionFiles==null) {
			LOG.error("option --regions must be defined");
			return -1;
			}
		this.intervalTreeMap = super.regionFiles.toIntervalTreeMap();
		if(this.intervalTreeMap.isEmpty()) {
			LOG.warn("no region defined.");
			}
		return super.beforeSam();
		}
	
	/** we just check any getAlignmentBlocks overlaps the intervalTreeMap */
	private boolean testCigar(final SAMRecord rec) {
		if(rec.getReadUnmappedFlag()) return false;
		return rec.getAlignmentBlocks().
				stream().
				map(B->new SimpleInterval(rec.getContig(),B.getReferenceStart(),B.getReferenceStart()+B.getLength()-1)).
				anyMatch(B->this.intervalTreeMap.containsOverlapping(B));
		}
	
	@Override
	protected Function<SAMRecord, List<SAMRecord>> createSAMRecordFunction()
		{
		return R->testCigar(R)?
				Collections.singletonList(R):
				Collections.emptyList();
		}
		
	public static void main(final String[] args) throws IOException
		{
		new Biostar9469733().instanceMainWithExit(args);
		}
		

	}
