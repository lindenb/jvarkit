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

import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMRecordQueryNameComparator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;

import java.io.File;
import java.io.IOException;
import java.util.Comparator;
import java.util.List;

import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

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
		description="Sort a BAM of contig and then on name",
		keywords={"sam","sort"},
		biostars=154220
		)
public class SortSamRefName extends Launcher
	{
	private static final Logger LOG = Logger.build(SortSamRefName.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection=new WritingSortingCollection();
	@ParametersDelegate
	private WritingBamArgs writingBamArgs=new WritingBamArgs();


	
	
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
	
	public SortSamRefName()
		{
		}

	@Override
	public int doWork(final List<String> args) {
		SamReader in=null;
		SAMFileWriter out=null;
		SAMRecordIterator iter=null;
		CloseableIterator<SAMRecord> iter2=null;
		SortingCollection<SAMRecord> sorter=null;
		try
			{
			in  = openSamReader(oneFileOrNull(args));
			final SAMFileHeader header= in.getFileHeader();
			
			final BAMRecordCodec bamRecordCodec=new BAMRecordCodec(header);
			final RefNameComparator refNameComparator=new RefNameComparator();
			sorter =SortingCollection.newInstance(
					SAMRecord.class,
					bamRecordCodec,
					refNameComparator,
					this.writingSortingCollection.getMaxRecordsInRam(),
					this.writingSortingCollection.getTmpPaths()
					);
			sorter.setDestructiveIteration(true);
			
			final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header).logger(LOG);
			iter = in.iterator();
			while(iter.hasNext())
				{
				sorter.add( progress.watch(iter.next()));
				}
			in.close();in=null;
			sorter.doneAdding();
			
			final SAMFileHeader header2=header.clone();
			header2.addComment(getProgramName()+" "+getVersion()+" "+getProgramCommandLine());
			header2.setSortOrder(SortOrder.unsorted);
			out = this.writingBamArgs.openSAMFileWriter(outputFile,header2, true);
			iter2 = sorter.iterator();
			while(iter2.hasNext())
				{
				out.addAlignment(iter2.next());
				}
			out.close();out=null;
			sorter.cleanup();
			progress.finish();
			return RETURN_OK;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter2);
			CloserUtil.close(iter);
			CloserUtil.close(out);
			CloserUtil.close(in);
			}
		}
	
	public static void main(final String[] args) throws IOException
		{
		new SortSamRefName().instanceMainWithExit(args);
		}
	}
