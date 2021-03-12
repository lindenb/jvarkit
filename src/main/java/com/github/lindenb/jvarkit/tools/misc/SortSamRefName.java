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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.util.Comparator;

import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.jcommander.OnePassBamLauncher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordQueryNameComparator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;

/**
BEGIN_DOC

## Example

```bash
$  java -jar dist/sortsamrefname.jar /commun/data/packages/samtools/1.2/samtools/examples/toy.sam  2> /dev/null 
@HD	VN:1.4	SO:unsorted
@SQ	SN:ref	LN:45
@SQ	SN:ref2	LN:40
@CO	SortSamRefName 1c7bc5e674136947586779a2aac53e576db4a67f /commun/data/packages/samtools/1.2/samtools/examples/toy.sam
r001	83	ref	37	30	9M	=	7	-39	CAGCGCCAT	*
r001	163	ref	7	30	8M4I4M1D3M	=	37	39	TTAGATAAAGAGGATACTG	*	XX:B:S,12561,2,20,112
r002	0	ref	9	30	1S2I6M1P1I1P1I4M2I	*	0	0	AAAAGATAAGGGATAAA	*
r003	0	ref	9	30	5H6M	*	0	0	AGCTAA	*
r003	16	ref	29	30	6H5M	*	0	0	TAGGC	*
r004	0	ref	16	30	6M14N1I5M	*	0	0	ATAGCTCTCAGC	*
x1	0	ref2	1	30	20M	*	0	0	AGGTTTTATAAAACAAATAA	????????????????????
x2	0	ref2	2	30	21M	*	0	0	GGTTTTATAAAACAAATAATT	?????????????????????
x3	0	ref2	6	30	9M4I13M	*	0	0	TTATAAAACAAATAATTAAGTCTACA	??????????????????????????
x4	0	ref2	10	30	25M	*	0	0	CAAATAATTAAGTCTACAGAGCAAC	?????????????????????????
x5	0	ref2	12	30	24M	*	0	0	AATAATTAAGTCTACAGAGCAACT	????????????????????????
x6	0	ref2	14	30	23M	*	0	0	TAATTAAGTCTACAGAGCAACTA	???????????????????????
```

END_DOC
 */

@Program(
		name="sortsamrefname",
		description="Sort a BAM on chromosome/contig and then on read/querty name",
		keywords={"sam","sort"},
		biostars= {154220,483658},
		creationDate="20150812",
		modificationDate="20210312"
		)
public class SortSamRefName extends OnePassBamLauncher {
	private static final Logger LOG = Logger.build(SortSamRefName.class).make();

	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection=new WritingSortingCollection();
	
	private static class RefNameComparator implements
		Comparator<SAMRecord>
		{
		private final SAMRecordQueryNameComparator nameCmp=new SAMRecordQueryNameComparator();
		
		@Override
		public int compare(final SAMRecord o1, final SAMRecord o2)
			{
	        final int refIndex1 = o1.getReferenceIndex();
	        final int refIndex2 = o2.getReferenceIndex();
	        final int cmp = refIndex1 - refIndex2;
	        
	        if (cmp != 0)
	        	{
	        	if (refIndex1 == -1 ) return 1;
	        	if (refIndex2 == -1 ) return -1;
	            return cmp;
	        	}
			return nameCmp.compare(o1, o2);
			}
		
		}
	
	@Override
	protected Logger getLogger()
		{
		return LOG;
		}

	@Override
	protected SAMFileHeader createOutputHeader(SAMFileHeader headerIn) {
		final SAMFileHeader hOut =  super.createOutputHeader(headerIn);
		hOut.setSortOrder(SAMFileHeader.SortOrder.unsorted);
		return hOut;
		}
	
	@Override
	protected void scanIterator(
			final SAMFileHeader headerIn,
			final CloseableIterator<SAMRecord> iter,
			final SAMFileWriter out)
		{
		
		SortingCollection<SAMRecord> sorter=null;
		try
			{
			final BAMRecordCodec bamRecordCodec=new BAMRecordCodec(headerIn);
			final RefNameComparator refNameComparator=new RefNameComparator();
			sorter =SortingCollection.newInstance(
					SAMRecord.class,
					bamRecordCodec,
					refNameComparator,
					this.writingSortingCollection.getMaxRecordsInRam(),
					this.writingSortingCollection.getTmpPaths()
					);
			sorter.setDestructiveIteration(true);
			
			while(iter.hasNext())
				{
				sorter.add(iter.next());
				}
			sorter.doneAdding();
			
			
			try(CloseableIterator<SAMRecord> iter2 = sorter.iterator()) {
				while(iter2.hasNext())
					{
					out.addAlignment(iter2.next());
					}
				}
			sorter.cleanup();
			sorter=null;
			}
		catch(final Throwable err)
			{
			throw new SAMException(err);
			}
		finally
			{
			CloserUtil.close(iter);
			}
		}
	
	public static void main(final String[] args) throws IOException
		{
		new SortSamRefName().instanceMainWithExit(args);
		}
	}
