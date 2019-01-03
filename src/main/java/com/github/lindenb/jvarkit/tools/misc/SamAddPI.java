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

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.math.stats.Percentile;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.samtools.SamRecordJEXLFilter;


/**

BEGIN_DOC


### Example


```
$ samtools view -h ~/src/gatk-ui/testdata/S1.bam | head -n 10
@HD	VN:1.5	GO:none	SO:coordinate
@SQ	SN:rotavirus	LN:1074
@RG	ID:S1	SM:S1
@PG	ID:bwa	PN:bwa	VN:0.7.12-r1044	CL:../bwa/bwa mem -R @RG\tID:S1\tSM:S1 ref.fa S1_01_R1.fq.gz S1_01_R2.fq.gz
@PG	ID:bwa.1	PN:bwa	VN:0.7.12-r1044	CL:../bwa/bwa mem -R @RG\tID:S1\tSM:S1 ref.fa S1_02_R1.fq.gz S1_02_R2.fq.gz
@PG	ID:bwa.2	PN:bwa	VN:0.7.12-r1044	CL:../bwa/bwa mem -R @RG\tID:S1\tSM:S1 ref.fa S1_03_R1.fq.gz S1_03_R2.fq.gz
rotavirus_1_317_5:0:0_7:0:0_2de	99	rotavirus	1	60	70M	=	248	317	GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAATATGGCGTCAACTCAGCAGATGGTCAGCTCTAATATT	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	MD:Z:33G4A3T14A7T4	RG:Z:S1	NM:i:5	AS:i:45	XS:i:0
rotavirus_1_535_4:0:0_4:0:0_1a6	163	rotavirus	1	60	70M	=	466	535	GGCTTTTACTGCTTTTCAGTGGTTGCTTCTCAAGATGGAGTGTACTCATCAGATGGTAAGCTCTATTATT	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	MD:Z:8A18G13C6G21	RG:Z:S1	NM:i:4	AS:i:50	XS:i:0
rotavirus_1_543_5:0:0_11:0:0_390	163	rotavirus	1	60	70M	=	487	530	GGCTTTTAATGCTTTTCATTTGATGCTGCTCAAGATGGAGTCTACACAGCAGATGGTCAGCTCTATTATT	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	MD:Z:18G1G1T22T11A12	RG:Z:S1	NM:i:5	AS:i:45XS:i:0
rotavirus_1_578_3:0:0_7:0:0_7c	99	rotavirus	1	60	70M	=	509	578	GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAAGATGGAGTCTCCTGAGCAGCTGGTAAGCTCTATTATT	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	MD:Z:43A2C5A17	RG:Z:S1	NM:i:3	AS:i:55	XS:i:0



$ java -jar dist/samaddpi.jar S1.bam | head -n 10
@HD	VN:1.5	GO:none	SO:coordinate
@SQ	SN:rotavirus	LN:1074
@RG	ID:S1	SM:S1	PI:498
@PG	ID:bwa	PN:bwa	VN:0.7.12-r1044	CL:../bwa/bwa mem -R @RG\tID:S1\tSM:S1 ref.fa S1_01_R1.fq.gz S1_01_R2.fq.gz
@PG	ID:bwa.1	PN:bwa	VN:0.7.12-r1044	CL:../bwa/bwa mem -R @RG\tID:S1\tSM:S1 ref.fa S1_02_R1.fq.gz S1_02_R2.fq.gz
@PG	ID:bwa.2	PN:bwa	VN:0.7.12-r1044	CL:../bwa/bwa mem -R @RG\tID:S1\tSM:S1 ref.fa S1_03_R1.fq.gz S1_03_R2.fq.gz
@CO	Processed with SamAddPI /home/lindenb/src/gatk-ui/testdata/S1.bam
rotavirus_1_317_5:0:0_7:0:0_2de	99	rotavirus	1	60	70M	=	248	317	GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAATATGGCGTCAACTCAGCAGATGGTCAGCTCTAATATT	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	MD:Z:33G4A3T14A7T4	RG:Z:S1	NM:i:5	AS:i:45	XS:i:0
rotavirus_1_535_4:0:0_4:0:0_1a6	163	rotavirus	1	60	70M	=	466	535	GGCTTTTACTGCTTTTCAGTGGTTGCTTCTCAAGATGGAGTGTACTCATCAGATGGTAAGCTCTATTATT	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	MD:Z:8A18G13C6G21	RG:Z:S1	NM:i:4	AS:i:50	XS:i:0
rotavirus_1_543_5:0:0_11:0:0_390	163	rotavirus	1	60	70M	=	487	530	GGCTTTTAATGCTTTTCATTTGATGCTGCTCAAGATGGAGTCTACACAGCAGATGGTCAGCTCTATTATT	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	MD:Z:18G1G1T22T11A12	RG:Z:S1	NM:i:5	AS:i:45XS:i:0

```
END_DOC
*/
@Program(
		name="samaddpi",
		description="Add predicted median insert size 'PI' to SAM Read groups (RG).",
		keywords= {"sam","bam"}
		)
public class SamAddPI extends Launcher
	{
	private static final Logger LOG = Logger.build(SamAddPI.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-w","--overwrite"},description="Overwrite median insert size if it already exists")
	private boolean overwrite_existing=false;
	@Parameter(names={"-N","--num-reads"},description="Number of reads to test. Negative=all = memory consuming.")
	private int num_read_to_test = 1_000_000;
	@Parameter(names={"--filter"},description=SamRecordJEXLFilter.FILTER_DESCRIPTION,converter=SamRecordJEXLFilter.StringConverter.class)
	private SamRecordFilter samRecordFilter = SamRecordJEXLFilter.buildDefault();
	
	@ParametersDelegate
	private WritingBamArgs writingBamArgs=new WritingBamArgs();
	
	
	@Override
	public int doWork(final List<String> args) {
		
		final Map<String,List<Integer>> rg2insertsize= new HashMap<>(); 
		SamReader sfr=null;
		SamReader sfrTmp =null;
		SAMFileWriter sfw=null;
		File tmpBam = null; 
		SAMFileWriter tmpBamWriter = null;
		SAMFileWriter outWriter = null;
		CloseableIterator<SAMRecord> iter = null;
		CloseableIterator<SAMRecord> iterTmp = null;
		try
			{
			sfr = openSamReader(oneFileOrNull(args));
			SAMFileHeader header=sfr.getFileHeader();
			for(final SAMReadGroupRecord rg:header.getReadGroups())
				{
				if(!overwrite_existing &&
				   rg.getPredictedMedianInsertSize()!=null)
					{
					continue;
					}
				rg2insertsize.put(rg.getId(), new ArrayList<>(
						num_read_to_test<1L?10000:num_read_to_test
						));
				}
			
			tmpBam = File.createTempFile("__addpi", ".bam");
			
			tmpBamWriter = this.writingBamArgs.openSAMFileWriter(tmpBam, header,true);
			iter = sfr.iterator();
			int n_processed=0;
			final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(header).logger(LOG);
			while(iter.hasNext() && (this.num_read_to_test<0 || n_processed< this.num_read_to_test))
				{
				final SAMRecord rec=progress.watch(iter.next());
				tmpBamWriter.addAlignment(rec);
				final SAMReadGroupRecord rg = rec.getReadGroup();
				if(rg==null || rg.getId()==null) continue;
				final List<Integer> insertlist = rg2insertsize.get(rg.getId());
				if(insertlist==null) continue;
				
				if(rec.getReadUnmappedFlag()) continue;
				if(!rec.getReadPairedFlag()) continue;
				if(!rec.getFirstOfPairFlag()) continue;
				if(rec.getMateUnmappedFlag()) continue;
				if(!rec.getReferenceName().equals(rec.getMateReferenceName())) continue;
				if(this.samRecordFilter.filterOut(rec)) continue;
				
				final int len = rec.getInferredInsertSize();
				if(len==0) continue;
				insertlist.add(Math.abs(len));
				++n_processed;
				}
			tmpBamWriter.close();tmpBamWriter=null;
			//reopen tmp file
			sfrTmp = super.createSamReaderFactory().open(tmpBam);
			iterTmp = sfrTmp.iterator();
			//update dMedianInsertSize
			for(final SAMReadGroupRecord rg:header.getReadGroups())
				{
				final List<Integer> insertlist = rg2insertsize.get(rg.getId());
				if(insertlist==null || insertlist.isEmpty()) continue;
				rg.setPredictedMedianInsertSize((int)Percentile.median().
						evaluate(insertlist.stream().mapToDouble(I->I.doubleValue())));
				}
			header.addComment("Processed with "+getClass().getSimpleName()+" "+getProgramCommandLine());
			outWriter =  this.writingBamArgs.openSAMFileWriter(this.outputFile, header,true);
			while(iterTmp.hasNext())
				{
				outWriter.addAlignment(iterTmp.next());
				}
			iterTmp.close();iterTmp=null;
			sfrTmp.close();sfrTmp=null;
			tmpBam.delete();
			//finish writing original input
			while(iter.hasNext())
				{
				outWriter.addAlignment(progress.watch(iter.next()));
				}
			progress.finish();
			iter.close();iter=null;
			sfr.close();sfr=null;
			
			outWriter.close();
			return RETURN_OK;
			}
		catch(final Exception err)
			{			
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(tmpBamWriter);
			if(tmpBam!=null) tmpBam.delete();
			CloserUtil.close(outWriter);
			CloserUtil.close(sfr);
			CloserUtil.close(sfw);
			}
		}
	
	public static void main(final String[] args) {
		new SamAddPI().instanceMainWithExit(args);

	}

}
