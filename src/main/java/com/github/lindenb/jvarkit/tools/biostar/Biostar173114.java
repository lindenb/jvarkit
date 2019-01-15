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
package com.github.lindenb.jvarkit.tools.biostar;
import java.io.File;
import java.util.List;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;

/**

BEGIN_DOC


### Example

```
 $ java -jar dist/biostar173114.jar --keepSequence    my.bam  

@HD	VN:1.5	GO:none	SO:coordinate
@SQ	SN:rotavirus	LN:1074
R0	0	rotavirus	1	60	70M	*	0	0	GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAATATGGCGTCAACTCAGCAGATGGTCAGCTCTAATATT	*
R1	0	rotavirus	1	60	70M	*	0	0	GGCTTTTACTGCTTTTCAGTGGTTGCTTCTCAAGATGGAGTGTACTCATCAGATGGTAAGCTCTATTATT	*
R2	0	rotavirus	1	60	70M	*	0	0	GGCTTTTAATGCTTTTCATTTGATGCTGCTCAAGATGGAGTCTACACAGCAGATGGTCAGCTCTATTATT	*
R3	0	rotavirus	1	60	70M	*	0	0	GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAAGATGGAGTCTCCTGAGCAGCTGGTAAGCTCTATTATT	*
R4	0	rotavirus	1	60	70M	*	0	0	GGCATTTAATGCTTAACAGTGGTTGCTGCTCAAGATGGAGTCTACTCAGCAGATGGTAAGCTCTCTTATT	*
R5	0	rotavirus	2	60	70M	*	0	0	GCTTTTAAAGCTTTTCAGTTGTTGCTGCTCAAGATGGAGTCTACTCAGCAGATGTTACGCTCTATTATTA	*
R6	0	rotavirus	2	60	51M19S	*	0	0	GCTTTTAATGCTTTTCAGTTGTTGCTGCACAAGATGGAGTCTACACAGCAGCTGTTCATCTCTCTTCATC	*
R7	0	rotavirus	2	60	70M	*	0	0	GCTTTTAATGCTTTTCAGTGGTTTCTTCTCACGATGGAGTCTACTCAGCAGAAGGTAAGCACTATTATTA	*
R8	0	rotavirus	2	60	70M	*	0	0	GCTTTTAAAGCATTACAGTTGTTGCAGCTCAAGAAGGAGACTACTCAGCAGATGGTAAGCTCTATAATTA	*
R9	0	rotavirus	2	60	70M	*	0	0	GCTTTTAATTCTATTCAGTGGTTGCTGCTCCAGAAGGAGTCTACTCAGGAGATGGTACGCTCTCTTATTA	*
Ra	0	rotavirus	2	60	70M	*	0	0	GCTTTTAATGCTTTTCAGTGGTTGCTGCTCAAGATGGAGTCTACTCAGCAGATGGTAAGCTCAATTATTA	*
Rb	0	rotavirus	2	60	70M	*	0	0	GCTTTTAATGCTTTTCAGTTGTAGCTGCTCAAGATGGAGTCTACTCATCAGATGGTAAGCTCTCTTCTTA	*
Rc	0	rotavirus	2	60	63M7S	*	0	0	TCTTTAAATGCTTTTCAGTGTTTGCTGCTCAAGATGGAGTCTACTCAGCAGAAGGTAAGCTCTCTTAAAC	*
Rd	0	rotavirus	2	60	66M4S	*	0	0	GCATTTAATGCTTTTCAGTGGTTGCTGCACAAGATGGAGTCTACTCAGCAGATTGTAAGCTCTATTCTAA	*
Re	0	rotavirus	3	60	70M	*	0	0	CTTTTAATGCTTTTCAGTGGTTGCTGCTCAAGAAGGCGTCTCCTGATGAGATGGTAAGCTCTATTATTAA	*
Rf	0	rotavirus	3	60	70M	*	0	0	CTTTTAATGGTTATGAGTGGTTGGTGCACAAGATGGAGTCTACTCAGCAGATGGTACTCTCTATAATTAA	*
R10	0	rotavirus	3	60	70M	*	0	0	CTTTTAAAGCTTTTCAGTGGTTGCTGCTCAAGATGGAGTCTACTCAGCAGATGGTACTCTCTATTCTTAA	*
R11	0	rotavirus	3	60	70M	*	0	0	CTTTTAAAGCTTTTCAGAGGTTGCTGCTCAAGATGTAGTCTACTCAGGAGATTGTAAGCTCTATTATTAA	
```

## See also

  * https://bioinformatics.stackexchange.com/a/3866/71


END_DOC
*/


@Program(name="biostar173114",
	description="make a bam file smaller by removing unwanted information see also https://www.biostars.org/p/173114/",
	keywords= {"sam","bam"})
public class Biostar173114 extends Launcher
	{
	private static final Logger LOG = Logger.build(Biostar173114.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@Parameter(names={"-keepSeq","--keepSequence"},description="keep read sequence")
	private boolean keepSequence = false;
	
	@Parameter(names={"-keepName","--keepName","--name"},description="keep Read Name, do not try to create a shorter name")
	private boolean keepName = false;
	@Parameter(names={"-keepAtt","--keepAttributes"},description="keep Attributes")
	private boolean keepAttributes = false;
	@Parameter(names={"-keepRG","--keepReadGroup"},description="if attributes are removed, keep the RG")
	private boolean keepReadGroup = false;
	@Parameter(names={"-mate","--mate"},description="keep Mate/Paired Information")
	private boolean keepMate = false;
	@Parameter(names={"-keepQuals","--keepQualities"},description="keep base qualities")
	private boolean keepQualities = false;
	@Parameter(names={"-keepCigar","--keepCigar"},description="keep cigar : don't remove hard clip")
	private boolean keepCigar = false;
	
	
	@ParametersDelegate
	private WritingBamArgs writingBamArgs =new WritingBamArgs().setBestCompression();
	
	@Override
	public int doWork(final List<String> args) {
		if(keepQualities) keepSequence=true;
		SamReader sfr=null;
		SAMFileWriter sfw=null;
		SAMRecordIterator iter=null;
		try
			{
			sfr = super.openSamReader(oneFileOrNull(args));
			sfw = this.writingBamArgs.openSAMFileWriter(this.outputFile,sfr.getFileHeader(), true);
			
			iter = sfr.iterator();
			final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(sfr.getFileHeader()).logger(LOG);
			long nReads = 0;
			while (iter.hasNext()) {
				final SAMRecord record = progress.watch(iter.next());

				if(!this.keepAttributes)
					{
					final SAMReadGroupRecord g=record.getReadGroup();
					record.clearAttributes();
					if(g!=null && this.keepReadGroup)
						{
						record.setAttribute("RG",g.getId());
						}
					}
					
					
				record.setReadName(this.keepName?
						record.getReadName():
						"R" + Long.toHexString(nReads++));
				
				if(!this.keepMate && record.getReadPairedFlag()) {
					record.setReadPairedFlag(false);
					record.setMateReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
					record.setMateAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
					record.setMateUnmappedFlag(false); 
					record.setMateNegativeStrandFlag(false);
					record.setInferredInsertSize(0);
					record.setProperPairFlag(false);
					}
				
				
				if (!this.keepCigar && !record.getReadUnmappedFlag() && record.getCigar() != null) {
					record.setCigar(new Cigar(record.getCigar().
							getCigarElements().stream().
							filter(C->!C.getOperator().equals(CigarOperator.H)).
							collect(Collectors.toList())));
					}
				
				if(!this.keepSequence) {
					record.setReadBases(SAMRecord.NULL_SEQUENCE);
					}
				if(!this.keepQualities) {
					record.setBaseQualities(SAMRecord.NULL_QUALS);
					}
				sfw.addAlignment(record);
				}
			progress.finish();
			sfw.close();
			sfw = null;
			LOG.info("done");
			return RETURN_OK;
		}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(sfr);
			CloserUtil.close(sfw);
			}
		}
	
	public static void main(final String[] args)throws Exception
		{
		new Biostar173114().instanceMainWithExit(args);
		}
	}
