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


*/
package com.github.lindenb.jvarkit.tools.samrmdupnames;


import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

/**

BEGIN_DOC


### Examples


#### Example 1


```

java -jar  dist/samgrep.jar -R r001  -- samtools-0.1.18/examples/toy.sam 

@HD     VN:1.4  SO:unsorted
@SQ     SN:ref  LN:45
@SQ     SN:ref2 LN:40
@PG     ID:0    PN:com.github.lindenb.jvarkit.tools.samgrep.SamGrep     VN:dac03b80e9fd88a15648b22550e57d10c9bed725     CL:-R r001 samtools-0.1.18/examples/toy.sam
r001    163     ref     7       30      8M4I4M1D3M      =       37      39      TTAGATAAAGAGGATACTG     *       XX:B:S,12561,2,20,112
r001    83      ref     37      30      9M      =       7       -39     CAGCGCCAT       *

```


#### Example 4


```

java -jar  dist/samgrep.jar -R r001 -- -n 1 samtools-0.1.18/examples/toy.sam 

@HD     VN:1.4  SO:unsorted
@SQ     SN:ref  LN:45
@SQ     SN:ref2 LN:40
@PG     ID:0    PN:com.github.lindenb.jvarkit.tools.samgrep.SamGrep     VN:dac03b80e9fd88a15648b22550e57d10c9bed725     CL:-R r001 -n 1 samtools-0.1.18/examples/toy.sam
r001    163     ref     7       30      8M4I4M1D3M      =       37      39      TTAGATAAAGAGGATACTG     *       XX:B:S,12561,2,20,112

```






END_DOC
*/


import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.StringUtil;


/**

BEGIN_DOC

### Motivation

I got a BAM file with the same read duplicated. They have the same position, the same flags


END_DOC
*/


@Program(name="samrmdupnames",
	description="remove duplicated names in sorted BAM",
		keywords={"sam","bam"},
		modificationDate="20221207",
		creationDate="20240405",
		jvarkit_amalgamion = true,
		biostars = {9558284,9610986},
		menu="BAM Manipulation"
		)
public class SamRemoveDuplicatedNames extends Launcher {
	private static final Logger LOG = Logger.build(SamRemoveDuplicatedNames.class).make();
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	protected Path outputFile=null;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION)
	protected Path faidxPath =null;
	@Parameter(names={"--validation-stringency"},description="SAM Reader Validation Stringency")
	protected ValidationStringency validationStringency = ValidationStringency.LENIENT;

	@ParametersDelegate
	protected WritingBamArgs writingBamArgs = new WritingBamArgs();
	@Parameter(names={"--report"},description="Only report duplicate names in the output bam file")
	private boolean report_duplicate_names = false;
	
	/*
	private boolean sameUnmappedRead(final SAMRecord rec1, final SAMRecord rec2) {
		if(!(rec1.getReadUnmappedFlag() && rec2.getReadUnmappedFlag())) throw new IllegalStateException(rec1.getReadName()+" "+rec2.getReadName());
		int i =  rec1.getReadName().compareTo(rec2.getReadName());
		if(i!=0) return false;
		if(rec1.getReadPairedFlag()!=rec2.getReadPairedFlag())  throw new IllegalStateException(rec1.getReadName());
		if(rec1.getReadPairedFlag()) {
			if(rec1.getFirstOfPairFlag() && rec2.getFirstOfPairFlag()) return true;
			if(rec2.getSecondOfPairFlag() && rec2.getSecondOfPairFlag()) return true;
			return false;
			}
		return true;
		}*/
	
	private boolean sameRead(final SAMRecord rec1, final SAMRecord rec2) {
		if(rec1.isSecondaryOrSupplementary() && rec2.isSecondaryOrSupplementary()) return false;
		int i =  Integer.compare(rec1.getFlags(), rec2.getFlags());
		if(i!=0) return false;
		i = rec1.getReadName().compareTo(rec2.getReadName());
		return i==0;
	}
	 
	
	private boolean samePos(final SAMRecord rec1, final SAMRecord rec2) {
		final int tid1 = rec1.getReferenceIndex();
		final int tid2 = rec2.getReferenceIndex();
		int i = Integer.compare(tid1, tid2);
		if(i!=0) return false;
		final int p1 = rec1.getAlignmentStart();
		final int p2 = rec2.getAlignmentStart();
		i = Integer.compare(p1, p2);
		return i==0;
		}
	
	private long dump(final SAMFileWriter w,final List<SAMRecord> buffer) {
		long n_removed = 0L;
		if(report_duplicate_names) {
			for(int x=0; x+1 < buffer.size();++x) {
				final SAMRecord rec1 = buffer.get(x);
				boolean saved1=false;
				for(int y=x+1; y < buffer.size();++y) {
					final SAMRecord rec2 = buffer.get(y);
					if(sameRead(rec1,rec2)) {
						if(!saved1) {
							w.addAlignment(rec1);
							n_removed++;
							saved1 = true;
							}
						w.addAlignment(rec2);
						n_removed++;
						}
					
					}
				}
			}
			
		else
			{
			int x=0;
			while(x +1 < buffer.size()) {
				SAMRecord best_mapq = buffer.get(x);
				int y = x+1;
				while(y < buffer.size()) {
					final SAMRecord rec2 = buffer.get(y);
					final boolean same_name = sameRead(best_mapq,rec2);
					if(same_name) {
						//update 20250405, read at[y] as better	MAPQ than read[x], just keep that read[y]
						if(best_mapq.getMappingQuality() < rec2.getMappingQuality()) {
							best_mapq  = rec2;
							buffer.set(x, best_mapq);
							}
						buffer.remove(y);
						n_removed++;
						}
					else
						{
						y++;
						}
					}
				x++;
				}
			
			
			for(SAMRecord rec2: buffer) {
				w.addAlignment(rec2);
				}
			}
		buffer.clear();
		return n_removed;
		}
	
	@Override
	public int doWork(final List<String> args) {
		try {
			long n_removed = 0L;
			final String input = oneFileOrNull(args);
			final SamReaderFactory swf  = super.createSamReaderFactory().
					referenceSequence(this.faidxPath).
					validationStringency(this.validationStringency);
			
			try(SamReader sr=StringUtil.isBlank(input)?swf.open(SamInputResource.of(stdin())):swf.open(Paths.get(input))) {
				final SAMFileHeader header= sr.getFileHeader();
				if(!header.getSortOrder().equals(SAMFileHeader.SortOrder.coordinate)) {
					throw new JvarkitException.BamBadSortOrder(SAMFileHeader.SortOrder.coordinate, header.getSortOrder());
					}
				JVarkitVersion.getInstance().addMetaData(this, header);
				
				try(CloseableIterator<SAMRecord> iter =  sr.iterator() ) {
					try(SAMFileWriter w = this.writingBamArgs.
								setReferencePath(this.faidxPath).
								setProgressLogger(LOG).
								openSamWriter(this.outputFile, header, true)) {
						final List<SAMRecord> buffer = new ArrayList<>();
						for(;;) {
							SAMRecord rec= iter.hasNext()?iter.next():null;
							
							
							/* both reads unmaped, we are in the non-mappend section at the end of the SAM file */
							if(rec!=null && rec.getReadUnmappedFlag() &&
								(!rec.getReadPairedFlag() || 
								rec.getMateUnmappedFlag())) {
								n_removed += dump(w,buffer);
								if(!this.report_duplicate_names) {
									w.addAlignment(rec);
									}
								continue;
								} /* end of section for unmapped reads */
							
							
							if(rec==null || (!buffer.isEmpty() && !samePos(buffer.get(0),rec))) {
								n_removed += dump(w,buffer);
								if(rec==null) break;
								}
							
							buffer.add(rec);
							}
						}
					}
				}
			LOG.info("N duplicate found: "+n_removed);
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
    
    public static void main(final String[] argv)
		{
	    new SamRemoveDuplicatedNames().instanceMainWithExit(argv);
		}	

	}
