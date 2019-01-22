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
import java.util.ArrayList;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloserUtil;


/**

BEGIN_DOC

## Example

Before fixing:
```bash
samtools view src.bam | grep "HWI-1KL149:6:C0KTCACXX:5:2107:2283:35906"
HWI-1KL149:6:C0KTCACXX:5:2107:2283:35906        177     3       1264832 37      101M    =       1264940 109     AGGTGGTGAAGCATGAGATGTAGGGAGAGCTGCTTTAAAACCCAGCACAAGGCTGGTT
GTACTGGCTCACACCTGTAATCCCAGGTCTTTGGGAGGCTGAG    """#""""#"""#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""#""""""#""##"""""""#"#"   X0:i:1  X1:i:0  BD:
Z:ABBACCABBABAAABABAABAACBBABAA@AAABAABBAABB@AAAAAAABA@ABAA@BAA@@BAAAAAAAA@@BAAAABA@ABAAABAACBBACBAABAA       MD:Z:0C0C0C1C0C0G1G0T1T0T0G0C0C0T0T0C0T0G0T0A0C0A2T2C0A0T0G0
C1T0G0T0G2G0T0T0T0G1G0T1T1C0C0A0A0G0T0G0C0G0A1T0G1G0C0T1T0T2C0G0T0G0T0G1C1C0A1G1G0T1C0C1T0T2C0A0     RG:Z:idp63088   XG:i:0  BI:Z:ABBAEDCCCBCBBABABAAAA@CBBAB@A@@A@AAAAAAA
BA@@AAA@A@BA@@BAA@A@@@@BAA@A@@@A@@AAAAABA@@BAAABBACBBACBBABAA       AM:i:37 NM:i:75 SM:i:37 XM:i:1  XO:i:0  MQ:i:37 XT:A:U
HWI-1KL149:6:C0KTCACXX:5:2107:2283:35906        113     3       1264940 37      101M    =       1264832 -109    GGTATCTCCATGCTCGAAGCCCTGACCTACTGTATTGCCCCGAAAGTCTTCCCTGCTG
TGGCTGCATCTTTTCCACGTGGATAATCTTGGTTCATCTCTAG    """##"""""""""""""""""""""""""#"""""""""""""""""""""""""""#"""""""""""""""#""""""""""""#""#""##"#"##"   X0:i:1  X1:i:0  BD:
Z:BBAABBBBAAABBBCBAABCBA@BAAAAAAABAAAAACCCBABAAAAAAACBAAAAABABA@AA@AAABBAAAAACB@BBAAAAAAAABBBAABBBBAAAA       MD:Z:0T0T1G0C0A1G0T1C0C0A1G0T0G0C0A0T0G0T0G0T0G2A0T4G0C0A0A0
T0G0T0G0C1G0G0T0G1C0A0G0T0T0G0C0A4C1A0T0G0C0G0T0G2G0G1C0G0T0G0A1C0G0T0G1G0C2T2T0C0G0T0G0T0A0T1       RG:Z:idp63088   XG:i:0  BI:Z:BABADDCCBBBCBBCBAABCBA@AABAAA@AAAAA@BBBB
BAAA@AA@AABA@@A@@A@BA@@A@AA@AAAAAAABB@BAAAAAAAA@CBAAABBBBAAAA       AM:i:37 NM:i:74 SM:i:37 XM:i:0  XO:i:0  MQ:i:37 XT:A:U
```
Fixing:
```bash
 java -jar dist/biostar76892.jar ONLYSAVEFIXED=true \
 	IN=src.bam \
 	OUT=fix.bam  \
 	VALIDATION_STRINGENCY=LENIENT
```
result

```bash
samtools view fix.bam | grep "HWI-1KL149:6:C0KTCACXX:5:2107:2283:35906"
HWI-1KL149:6:C0KTCACXX:5:2107:2283:35906        163     3       1264832 37      101M    =       1264940 109     AGGTGGTGAAGCATGAGATGTAGGGAGAGCTGCTTTAAAACCCAGCACAAGGCTGGTT
GTACTGGCTCACACCTGTAATCCCAGGTCTTTGGGAGGCTGAG    """#""""#"""#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""#""""""#""##"""""""#"#"   X0:i:1  X1:i:0  BD:
Z:ABBACCABBABAAABABAABAACBBABAA@AAABAABBAABB@AAAAAAABA@ABAA@BAA@@BAAAAAAAA@@BAAAABA@ABAAABAACBBACBAABAA       MD:Z:0C0C0C1C0C0G1G0T1T0T0G0C0C0T0T0C0T0G0T0A0C0A2T2C0A0T0G0
C1T0G0T0G2G0T0T0T0G1G0T1T1C0C0A0A0G0T0G0C0G0A1T0G1G0C0T1T0T2C0G0T0G0T0G1C1C0A1G1G0T1C0C1T0T2C0A0     RG:Z:idp63088   XG:i:0  BI:Z:ABBAEDCCCBCBBABABAAAA@CBBAB@A@@A@AAAAAAA
BA@@AAA@A@BA@@BAA@A@@@@BAA@A@@@A@@AAAAABA@@BAAABBACBBACBBABAA       AM:i:37 NM:i:75 SM:i:37 XM:i:1  XO:i:0  MQ:i:37 XT:A:U  rv:i:1
HWI-1KL149:6:C0KTCACXX:5:2107:2283:35906        83      3       1264940 37      101M    =       1264832 -109    GGTATCTCCATGCTCGAAGCCCTGACCTACTGTATTGCCCCGAAAGTCTTCCCTGCTG
TGGCTGCATCTTTTCCACGTGGATAATCTTGGTTCATCTCTAG    """##"""""""""""""""""""""""""#"""""""""""""""""""""""""""#"""""""""""""""#""""""""""""#""#""##"#"##"   X0:i:1  X1:i:0  BD:
Z:BBAABBBBAAABBBCBAABCBA@BAAAAAAABAAAAACCCBABAAAAAAACBAAAAABABA@AA@AAABBAAAAACB@BBAAAAAAAABBBAABBBBAAAA       MD:Z:0T0T1G0C0A1G0T1C0C0A1G0T0G0C0A0T0G0T0G0T0G2A0T4G0C0A0A0
T0G0T0G0C1G0G0T0G1C0A0G0T0T0G0C0A4C1A0T0G0C0G0T0G2G0G1C0G0T0G0A1C0G0T0G1G0C2T2T0C0G0T0G0T0A0T1       RG:Z:idp63088   XG:i:0  BI:Z:BABADDCCBBBCBBCBAABCBA@AABAAA@AAAAA@BBBB
BAAA@AA@AABA@@A@@A@BA@@A@AA@AAAAAAABB@BAAAAAAAA@CBAAABBBBAAAA       AM:i:37 NM:i:74 SM:i:37 XM:i:0  XO:i:0  MQ:i:37 XT:A:U  rv:i:1
```

END_DOC

*/

@Program(name="biostar76892",
description="fix strand of two paired reads close but on the same strand. ",
biostars=76892,
keywords={"sam","bam"})
public class Biostar76892 extends Launcher
	{
	private static final Logger LOG = Logger.build(Biostar76892.class).make();
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	
	@Parameter(names={"-osf","--osf"},description="only save pairs of reads which have been corrected by this program")
	private boolean onlySaveFixed = false;
	@Parameter(names={"-d","--maxc"},description="distance beween two reads."+DistanceParser.OPT_DESCRIPTION,splitter=com.github.lindenb.jvarkit.util.jcommander.NoSplitter.class,converter=DistanceParser.StringConverter.class)
	private int distance = 30 ;

	@ParametersDelegate
	private WritingBamArgs writingBamArgs = new WritingBamArgs();
    
	@Override
	public int doWork(final List<String> args) {
		SamReader sfr=null;
		SAMFileWriter sfw=null;
		try
			{
			sfr= super.openSamReader(oneFileOrNull(args));
			sfw = writingBamArgs.openSAMFileWriter(this.outputFile,sfr.getFileHeader(), true);
			long nRecords=0;
			final List<SAMRecord> buffer=new ArrayList<SAMRecord>();
			SAMRecordIterator iter=sfr.iterator();
			for(;;)
				{
				SAMRecord rec=null;
				//get next record
				if(iter.hasNext())
					{
					rec=iter.next();
					++nRecords;
					if(nRecords%1000000==0) LOG.info("records: "+nRecords);
					if(!rec.getReadPairedFlag() ||
						rec.getReadUnmappedFlag() ||
						rec.getMateUnmappedFlag() ||
						rec.getProperPairFlag() ||
						rec.getReferenceIndex()!=rec.getMateReferenceIndex() ||
						rec.getReadNegativeStrandFlag()==!rec.getMateNegativeStrandFlag()
						)
						{
						if(!onlySaveFixed) sfw.addAlignment(rec);
						continue;
						}
					}

				if(rec!=null)
					{
					int i=0;
					//cleanup buffer
					int mate_index=-1;
					while(i<buffer.size())
						{
						SAMRecord prev=buffer.get(i);
						if(prev.getReferenceIndex()!=rec.getReferenceIndex() ||
							prev.getAlignmentEnd() + distance < rec.getAlignmentStart())
							{
							if(!onlySaveFixed) sfw.addAlignment(prev);
							buffer.remove(i);
							}
						else if(prev.getReadName().equals(rec.getReadName()) && 
								(	(prev.getFirstOfPairFlag() && rec.getSecondOfPairFlag()) ||
									(rec.getFirstOfPairFlag() && prev.getSecondOfPairFlag()))
								)
							{
							mate_index=i;
							++i;
							}
						else
							{
							++i;
							}
						}
					
					if(mate_index==-1)
						{
						buffer.add(rec);
						}
					else
						{
						final SAMRecord mate=buffer.get(mate_index);
						buffer.remove(mate_index);
						LOG.info("changing "+rec+" "+mate);
						
						if(mate.getReadNegativeStrandFlag())
							{
							mate.setReadNegativeStrandFlag(false);
							rec.setMateNegativeStrandFlag(mate.getReadNegativeStrandFlag());
							}
						else
							{
							rec.setReadNegativeStrandFlag(false);
							mate.setMateNegativeStrandFlag(rec.getReadNegativeStrandFlag());
							}
						if(!mate.getReadFailsVendorQualityCheckFlag() && !rec.getReadFailsVendorQualityCheckFlag())
							{
							mate.setProperPairFlag(true);
							rec.setProperPairFlag(true);
							}
						mate.setAttribute("rv",1);
						rec.setAttribute("rv", 1);
						sfw.addAlignment(mate);
						sfw.addAlignment(rec);
						}
					
					}
				else
					{
					for(final SAMRecord r:buffer)
						{
						if(!onlySaveFixed)  sfw.addAlignment(r);
						}
					break;
					}
				}
			LOG.info("done");
			sfw.close();
			return RETURN_OK;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(sfw);
			CloserUtil.close(sfr);
			}
		}
	
	public static void main(final String[] args)throws Exception
		{
		new Biostar76892().instanceMainWithExit(args);
		}
	}
