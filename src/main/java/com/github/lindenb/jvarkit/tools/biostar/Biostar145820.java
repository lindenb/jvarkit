/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
* 2015 creation

*/
package com.github.lindenb.jvarkit.tools.biostar;

import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMRecordQueryNameComparator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.BinaryCodec;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;

import java.io.InputStream;
import java.io.OutputStream;
import java.nio.file.Path;
import java.util.Comparator;
import java.util.List;
import java.util.Random;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.samtools.SamRecordJEXLFilter;

/*
BEGIN_DOC

## Example

```bash
$ java -jar jvarkit.jar biostar145820 -n 10  -o out.bam  in.bam
```

## CIted In:

 * "MED25 connects enhancer-promoter looping and MYC2-dependent activation of jasmonate signalling. Wang et al. Nature Plants 5, 616-625 (2019)  https://doi.org/10.1038/s41477-019-0441-9 
 * Makamure, C.E., Justinen, S., Martínez, D.E. et al. Cool temperature inhibits binary fission and results in phenotypic and transcriptomic changes that suggest inducible aging in Diadumene lineata. BMC Res Notes 18, 293 (2025). https://doi.org/10.1186/s13104-025-07378-x 
```
END_DOC

 */
@Program(name="biostar145820",
	description="subsample/shuffle BAM to fixed number of alignments.",
	biostars=145820,
	keywords= {"sam","bam","shuffle"},
	creationDate="20150615",
	modificationDate="20211005",
	references={"MED25 connects enhancer-promoter looping and MYC2-dependent activation of jasmonate signalling. Wang et al. Nature Plants 5, 616-625 (2019)  https://doi.org/10.1038/s41477-019-0441-9 "},
	jvarkit_amalgamion =  true,
	menu="BAM Visualization"
	)
public class Biostar145820 extends Launcher
	{
	private static final Logger LOG = Logger.of(Biostar145820.class);

	@Parameter(names={"-f","--filter","--jexl"},description = SamRecordJEXLFilter.FILTER_DESCRIPTION,converter=SamRecordJEXLFilter.StringConverter.class)
	private SamRecordFilter filter  = SamRecordJEXLFilter.buildAcceptAll();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-n"},description=" number of reads. negative: all reads, shuffle output.")
	private long count=-1L;
	@Parameter(names={"--reference","-R"},description=CRAM_INDEXED_REFENCE)
	private Path refPath = null;
	@Parameter(names={"--seed"},description="Random seed. -1 == current time")
	private long seed = -1L;

	
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection = new WritingSortingCollection();
	@ParametersDelegate
	private WritingBamArgs writingBamArgs = new WritingBamArgs();

	
	private static class RandSamRecord
		{
		long l_rand_index;
		SAMRecord samRecord;
		}
	private static class RandSamRecordComparator
		implements Comparator<RandSamRecord>
		{
		final SAMRecordQueryNameComparator secondCompare =new SAMRecordQueryNameComparator();
		@Override
		public int compare(final RandSamRecord o1,final RandSamRecord o2)
			{
			final int i = Long.compare(o1.l_rand_index , o2.l_rand_index);
			if(i!=0) return i;
			return secondCompare.compare(o1.samRecord, o2.samRecord);
			}
		}
	
	private static class RandSamRecordCodec
		implements SortingCollection.Codec<RandSamRecord>
		{
	    private final BinaryCodec binaryCodec = new BinaryCodec();
	    private final BAMRecordCodec bamRecordCodec;
	    private final SAMFileHeader header;
		RandSamRecordCodec(final SAMFileHeader header)
			{
			this.header=header;
			this.bamRecordCodec=new BAMRecordCodec(this.header);
			}
		
		@Override
		public RandSamRecord decode()
			{
			long r;
			  try {
		          r = this.binaryCodec.readLong();
		        }
		     catch (final Exception e) {
		            return null;
		        }
			  final RandSamRecord o =new RandSamRecord();
			  o.l_rand_index=r;
			  o.samRecord = this.bamRecordCodec.decode();
			  return o;
			}
		
		@Override
		public void encode(RandSamRecord val) {
			this.binaryCodec.writeLong(val.l_rand_index);
			this.bamRecordCodec.encode(val.samRecord);
			}
		@Override
		public void setInputStream(InputStream is) {
			this.binaryCodec.setInputStream(is);
			this.bamRecordCodec.setInputStream(is);
			}
		@Override
		public void setOutputStream(OutputStream os) {
			this.binaryCodec.setOutputStream(os);
			this.bamRecordCodec.setOutputStream(os);
			}
		
		@Override
		public RandSamRecordCodec clone()  {
			return new RandSamRecordCodec(this.header);
			}
		}
	
	@Override
	public int doWork(final List<String> args) {
		SamReader samReader=null;
		SAMRecordIterator iter=null;
		SAMFileWriter samWriter=null;
		final Random random=new Random(this.seed==-1L?System.currentTimeMillis():this.seed);
		CloseableIterator<RandSamRecord> iter2=null;
		try
			{
			final String input = oneFileOrNull(args);
			final SamReaderFactory srf = super.createSamReaderFactory().validationStringency(ValidationStringency.LENIENT);
			if(this.refPath!=null) srf.referenceSequence(this.refPath);
			samReader = srf.open(input==null?SamInputResource.of(stdin()):SamInputResource.of(input));
			
			final SAMFileHeader header=samReader.getFileHeader().clone();
						
			header.setSortOrder(SortOrder.unsorted);
			header.addComment("Processed with "+getProgramName()+" : "+getProgramCommandLine());
			
			
			final ProgressFactory.Watcher<SAMRecord> progress=ProgressFactory.newInstance().dictionary(samReader.getFileHeader()).logger(LOG).build();
			iter=samReader.iterator();

			
			final SortingCollection<RandSamRecord> sorter=SortingCollection.newInstance(
					RandSamRecord.class,
					new RandSamRecordCodec(header),
					new RandSamRecordComparator(), 
					this.writingSortingCollection.getMaxRecordsInRam(),
					this.writingSortingCollection.getTmpPaths()
					);
			sorter.setDestructiveIteration(true);
			while(iter.hasNext())
				{
				final SAMRecord rec = progress.apply(iter.next());
				if(this.filter.filterOut(rec)) {
					continue;
				}
				final RandSamRecord r=new RandSamRecord();
				r.l_rand_index  = random.nextLong();
				r.samRecord =  rec;

				sorter.add(r);
				}
			iter.close();iter=null;
			sorter.doneAdding();
			iter2=sorter.iterator();
			
			samWriter = writingBamArgs.openSamWriter(this.outputFile, header, true);

			final SAMFileWriter finalw = samWriter;
			
			iter2.stream().
				limit(this.count<0L?Long.MAX_VALUE-1:this.count).
				map(R->R.samRecord).
				forEach(R->finalw.addAlignment(R));
				
			iter2.close();iter2=null;
			sorter.cleanup();
			progress.close();
			}
		catch(final Throwable e)
			{
			LOG.error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(iter2);
			CloserUtil.close(samReader);
			CloserUtil.close(samWriter);
			}
		return 0;
		}

	
	public static void main(final String[] args)
		{
		new Biostar145820().instanceMainWithExit(args);
		}
		

	}
