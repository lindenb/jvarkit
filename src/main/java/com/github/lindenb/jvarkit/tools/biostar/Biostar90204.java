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
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.biostar;

import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.SamRecordJEXLFilter;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.util.CloseableIterator;


/**
BEGIN_DOC

##Example

```bash
$ java -jar dist/jvarkit.jar biostar90204 -m bam.manifest -n 3 -a 5 samtools-0.1.18/examples/toy.sam

$ cat bam.manifest
_splitbam.00001.bam	1	3
_splitbam.00002.bam	4	6
_splitbam.00003.bam	7	9
_splitbam.00004.bam	10	12

$ samtools-0.1.18/samtools view -h _splitbam.00003.bam 
@HD	VN:1.4	SO:unsorted
@SQ	SN:ref	LN:45
@SQ	SN:ref2	LN:40
@PG	ID:0	PN:com.github.lindenb.jvarkit.tools.biostar.Biostar90204	VN:7e17f8bd273cf081d4415bc4f579cd34e2c681d1	CL:-m bam.manifest -n 3 -a 5 samtools-0.1.
18/examples/toy.sam
@CO	SPLIT:3
@CO	SPLIT:Starting from Read7
x1	0	ref2	1	30	20M	*	0	0	AGGTTTTATAAAACAAATAA	????????????????????
x2	0	ref2	2	30	21M	*	0	0	GGTTTTATAAAACAAATAATT	?????????????????????
x3	0	ref2	6	30	9M4I13M	*	0	0	TTATAAAACAAATAATTAAGTCTACA	??????????????????????????
```

END_DOC
*/
@Program(name="biostar90204",
	keywords={"sam","bam","split","util"},
	description="Bam version of linux split.",
	modificationDate = "20250408",
	biostars=90204,
	jvarkit_amalgamion =  true,
	menu="Biostars"
	)
public class Biostar90204 extends Launcher
	{
	private static final Logger LOG = Logger.of(Biostar90204.class);
	@Parameter(names= {"-R"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private Path referencePath =null;
	@Parameter(names= {"-p","--prefix"},description="(prefix) output file prefix.")
	private String prefix="_splitbam";
	@Parameter(names= {"-a","--padding"},description="'0' padding length")
	private int suffix_length=2;
	@Parameter(names= {"-M","--manifest","-o"},description="Manifest file. Optional")
	private Path manifestFile=null;
	@Parameter(names={"--filter"},description=SamRecordJEXLFilter.FILTER_DESCRIPTION,converter=SamRecordJEXLFilter.StringConverter.class)
	private SamRecordFilter samRecordFilter = SamRecordJEXLFilter.buildAcceptAll();
	@Parameter(names= {"-n","--count"},description="Number of records per file.")
	private long record_per_file=-1L;

	
	@ParametersDelegate
	private WritingBamArgs writingBamArgs=new WritingBamArgs();
	
	@Override
	public int doWork(final List<String> args) {
		
		if(this.suffix_length<0)
			{
			LOG.error("Bad value of suffix_length:"+ this.suffix_length);
			return -1;
			}
		if(this.record_per_file<1L)
			{
			LOG.error("Bad value of record_per_file:"+this.record_per_file);
			return -1;
			}
		
		
		SAMFileWriter sfw=null;
		
		try
			{
			final SamReaderFactory srf = super.createSamReaderFactory();
			srf.referenceSequence(this.referencePath);
			final String input = oneFileOrNull(args);
			try(SamReader samFileReader = srf.open(input==null?SamInputResource.of(stdin()):SamInputResource.of(Paths.get(input)))) {
				final SAMFileHeader header=samFileReader.getFileHeader();
				int split_file_number=0;
				long nReads=0L;
				try(PrintWriter manifest=(this.manifestFile==null? new PrintWriter(new NullOuputStream()):IOUtils.openPathForPrintWriter(manifestFile))) {
					try(CloseableIterator<SAMRecord> iter=samFileReader.iterator()) {
						while(iter.hasNext())
							{
							final SAMRecord rec=iter.next();
							if(this.samRecordFilter.filterOut(rec)) continue;
							++nReads;
							if(sfw==null)
								{
								split_file_number++;
								final String pathname=(this.prefix.isEmpty()?"":this.prefix+".")+String.format("%0"+suffix_length+"d", split_file_number)+".bam";
								final Path out= Paths.get(pathname);
								manifest.write(pathname);
								manifest.write("\t"+(nReads)+"\t");
								
								final SAMFileHeader header2=header.clone();
								header2.addComment("SPLIT:"+split_file_number);
								header2.addComment("SPLIT:Starting from Read"+nReads);
								
								sfw = this.writingBamArgs.openSamWriter(out,header2, true);
								}
							sfw.addAlignment(rec);
							
							if(nReads%record_per_file==0)
								{
								sfw.close();
								manifest.write((nReads)+"\n");
								sfw=null;
								}
							}
						}
					if(sfw!=null)
						{
						sfw.close();
						manifest.write((nReads)+"\n");
						}
					manifest.flush();
					}
				}
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			if(sfw!=null) try {sfw.close();}catch(Throwable err) {}
			}
		return 0;
		}

	public static void main(final String[] args) {
		new Biostar90204().instanceMainWithExit(args);
	}

}
