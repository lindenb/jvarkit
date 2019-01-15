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
package com.github.lindenb.jvarkit.tools.mem;

import java.io.File;
import java.util.ArrayList;
import java.util.List;


import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.DefaultSAMRecordFactory;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecord.SAMTagAndValue;
import htsjdk.samtools.SAMRecordFactory;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SAMValidationError;
import htsjdk.samtools.util.CloserUtil;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

/**
BEGIN_DOC
## Example

The original SAM output:

```bash
$ samtools view file.bam | grep -F "HWI-1KL149:61:D2C11ACXX:4:1205:18030:70811"
HWI-1KL149:61:D2C11ACXX:4:1205:18030:70811	163	8	833200	60	62M39S	=	833209	62	ATCCCAGCAAGTATTGTACTCAGGCGTCTAGAAGCTCCAATCGCAGAATCCACAGAAAGAGCGCCTGTAGCAGAAGGACTTGATTGATGTTGAATGCAACA	D@CCFFFFFGHGFHJGIHGGGIIIJGIGGHGHGDHIJJJJIJIIIJJJIGIGCDEHIJEHEDDEFFEECACDCDDDDD@CCDDDDDDDDDDDBDDDDDD>A	NM:i:0	MD:Z:62AS:i:62	XS:i:19	SA:Z:8,836189,+,62S34M1D5M,49,1;
HWI-1KL149:61:D2C11ACXX:4:1205:18030:70811	83	8	833209	60	53M48S	=	833200	-62	AGTATTGTACTCAGGCGTCTAGAAGCTCCAATCGCAGAATCCACAGAAAGAGCGCCTGTAGCAGAAGGACTTGATTGATGTTGAATGCAACAGGGCGCCCC	DBDBBB?BCCDDDCC=5CCCFFDFB?ACHECGGGIHCGGCCHGBGGJHGHGGGJJIFBJJIIHFGGIIJIGIGIIIGHGCIIGIHGIIHGHHHFFFFFCCC	NM:i:0	MD:Z:53AS:i:53	XS:i:19	SA:Z:8,836189,-,53S34M1D14M,60,1;
HWI-1KL149:61:D2C11ACXX:4:1205:18030:70811	419	8	836189	49	62H34M1D5M	=	833209	-2929	GCCTGTAGCAGAAGGACTTGATTGATGTTGAATGCAACA	DEFFEECACDCDDDDD@CCDDDDDDDDDDDBDDDDDD>A	NM:i:1	MD:Z:34^T5	AS:i:35	XS:i:26	SA:Z:8,833200,+,62M39S,60,0;
HWI-1KL149:61:D2C11ACXX:4:1205:18030:70811	339	8	836189	60	53H34M1D14M	=	833200	-3038	GCCTGTAGCAGAAGGACTTGATTGATGTTGAATGCAACAGGGCGCCCC	JJIFBJJIIHFGGIIJIGIGIIIGHGCIIGIHGIIHGHHHFFFFFCCC	NM:i:1	MD:Z:34^A14	AS:i:41	XS:i:26	SA:Z:8,833209,-,53M48S,60,0;
```

Using BWAMemNOp:

```
$ java -jar dist/bwamemnop.jar -S  file.bam | grep -F "HWI-1KL149:61:D2C11ACXX:4:1205:18030:70811"
HWI-1KL149:61:D2C11ACXX:4:1205:18030:70811	163	8	833200	60	62M2927N34M1D5M	=	833209	62	ATCCCAGCAAGTATTGTACTCAGGCGTCTAGAAGCTCCAATCGCAGAATCCACAGAAAGAGCGCCTGTAGCAGAAGGACTTGATTGATGTTGAATGCAACA	D@CCFFFFFGHGFHJGIHGGGIIIJGIGGHGHGDHIJJJJIJIIIJJJIGIGCDEHIJEHEDDEFFEECACDCDDDDD@CCDDDDDDDDDDDBDDDDDD>A	MD:Z:62NM:i:0	AS:i:62	XS:i:19
HWI-1KL149:61:D2C11ACXX:4:1205:18030:70811	83	8	833209	60	53M2927N34M1D14M	=	833200	-62	AGTATTGTACTCAGGCGTCTAGAAGCTCCAATCGCAGAATCCACAGAAAGAGCGCCTGTAGCAGAAGGACTTGATTGATGTTGAATGCAACAGGGCGCCCC	DBDBBB?BCCDDDCC=5CCCFFDFB?ACHECGGGIHCGGCCHGBGGJHGHGGGJJIFBJJIIHFGGIIJIGIGIIIGHGCIIGIHGIIHGHHHFFFFFCCC	MD:Z:53	NM:i:0	AS:i:53	XS:i:19

```

converted to BED12 with bedtools/**bamToBed** and visualized in the UCSC genome browser... :


![http://i.imgur.com/NccxWdT.jpg](http://i.imgur.com/NccxWdT.jpg)


END_DOC
 */
@Program(
		name="bwamemnop",
		description="Merge the other BWA-MEM alignements with its SA:Z:* attributes to an alignment containing a cigar string with 'N' (  Skipped region from the reference.)",
		keywords={"bwa","sam","bam"}
		)
public class BWAMemNOp extends Launcher
	{
	private static final Logger LOG=Logger.build(BWAMemNOp.class).make();
	
	@Parameter(names={"-o"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile=null;
	@Parameter(names={"-m"},description="min soft clip length")
	private int min_soft_clip_length=10;
	@Parameter(names={"-S"},description=" only print converted reads")
	private boolean print_only_spit_read=false;

	@ParametersDelegate
	private WritingBamArgs writingBamArgs=new WritingBamArgs();
	
	static class CigarEvt
		{
		CigarOperator op;
		int read0;
		int ref1;
		}
	private static List<CigarEvt> cigarEvents(int read0,int ref1,Cigar cigar)
		{
		List<CigarEvt> L=new ArrayList<CigarEvt>();
		for(CigarElement ce:cigar.getCigarElements())
			{
			for(int i=0;i< ce.getLength();++i)
				{
				CigarEvt evt=new CigarEvt();
				evt.op=ce.getOperator();
				evt.read0=read0;
				evt.ref1=ref1;
				L.add(evt);
				switch(ce.getOperator())
					{
					case P:break;
					case H:break;
					case I:
					case S:
						{
						read0++;
						break;
						}
					case D:
					case N:
						{
						ref1++;
						break;
						}
					case X:case EQ:case M:
						{
						read0++;
						ref1++;
						break;
						}
					default:throw new IllegalStateException();
					}
				}
			}
		return L;
		}
	
	
	@Override
	public int doWork(List<String> args) {
		SamReader r=null;
		SAMFileWriter w=null;
		try
			{
			r = super.openSamReader(oneFileOrNull(args));
			SAMFileHeader header=r.getFileHeader();
			w = writingBamArgs.openSAMFileWriter(outputFile, header, true);
			
			SAMRecordFactory samRecordFactory=new DefaultSAMRecordFactory();
			SAMRecordIterator iter=r.iterator();
			while(iter.hasNext())
				{
				SAMRecord rec=iter.next();
				if(rec.getSupplementaryAlignmentFlag())
					{
					continue;
					}
				if(rec.getReadUnmappedFlag())
					{
					if(!print_only_spit_read) w.addAlignment(rec);
					continue;
					}
				Cigar cigar1=rec.getCigar();
				if(cigar1==null ||
					cigar1.isEmpty() ||
					!(cigar1.getCigarElement(cigar1.numCigarElements()-1).getOperator().equals(CigarOperator.S) ||
					  cigar1.getCigarElement(0).getOperator().equals(CigarOperator.S)
					 )
					) //last or first is soft clipping
					{
					if(!print_only_spit_read) w.addAlignment(rec);
					continue;
					}
				rec.getAlignmentStart();
				final List<SAMRecord> xps=SAMUtils.getOtherCanonicalAlignments(rec);
				if(xps.isEmpty())
					{
					if(!print_only_spit_read) w.addAlignment(rec);
					continue;
					}
				
				boolean found_one=false;
				for(final SAMRecord  xp:xps)
					{
					if(!rec.getReferenceName().equals(xp.getReferenceName())) continue;
					if(xp.getReadNegativeStrandFlag()!=rec.getReadNegativeStrandFlag() ) continue;
					
					Cigar cigar2=xp.getCigar();
					if(cigar2==null || cigar2.isEmpty() )
						{
						continue;
						}
					
					SAMRecord newrec=null;
					List<CigarEvt> L1=null;
					List<CigarEvt> L2=null;
					if( cigar1.getCigarElement(cigar1.numCigarElements()-1).getOperator().equals(CigarOperator.S) &&
						cigar1.getCigarElement(cigar1.numCigarElements()-1).getLength() >= this.min_soft_clip_length &&
						cigar2.getCigarElement(0).getOperator().equals(CigarOperator.S)  &&
						cigar2.getCigarElement(0).getLength() >= this.min_soft_clip_length &&
					   rec.getAlignmentEnd()< xp.getAlignmentStart()
					   )
						{
						newrec=samRecordFactory.createSAMRecord(header);
						int ref1=rec.getAlignmentStart();
						newrec.setAlignmentStart(ref1);
						L1=cigarEvents(0, ref1, cigar1);
						L2=cigarEvents(0, xp.getAlignmentStart(), cigar2);
						}
					else  if(
							cigar2.getCigarElement(cigar2.numCigarElements()-1).getOperator().equals(CigarOperator.S) &&
							cigar2.getCigarElement(cigar2.numCigarElements()-1).getLength() >= this.min_soft_clip_length &&
							cigar1.getCigarElement(0).getOperator().equals(CigarOperator.S) &&
							cigar1.getCigarElement(0).getLength() >= this.min_soft_clip_length &&
							xp.getAlignmentEnd()< rec.getAlignmentStart()
							)
						{
						newrec=samRecordFactory.createSAMRecord(header);
						int ref1=xp.getAlignmentStart();
						newrec.setAlignmentStart(ref1);
						L1=cigarEvents(0, ref1, cigar2);
						L2=cigarEvents(0, rec.getAlignmentStart(), cigar1);
						}
					
					if(newrec==null) continue;
					
					newrec.setFlags(rec.getFlags());
					newrec.setReadName(rec.getReadName());
					newrec.setReadBases(rec.getReadBases());
					newrec.setMappingQuality(rec.getMappingQuality());
					newrec.setReferenceIndex(rec.getReferenceIndex());
					newrec.setBaseQualities(rec.getBaseQualities());
					if(found_one)
						{
						newrec.setNotPrimaryAlignmentFlag(true);
						}
					
					found_one=true;
					
					for(final SAMTagAndValue tav: rec.getAttributes())
						{
						if(tav.tag.equals("SA")) continue;
						if(tav.tag.equals("NM")) continue;
						newrec.setAttribute(tav.tag, tav.value);
						}
					if(rec.getReadPairedFlag() && !rec.getMateUnmappedFlag())
						{
						newrec.setMateAlignmentStart(rec.getMateAlignmentStart());
						newrec.setMateReferenceIndex(rec.getMateReferenceIndex());
						newrec.setInferredInsertSize(rec.getInferredInsertSize());
						}
					while(!L1.isEmpty() &&
						(L1.get(L1.size()-1).op.equals(CigarOperator.S) ||
						 L1.get(L1.size()-1).op.equals(CigarOperator.D) ||
						 L1.get(L1.size()-1).op.equals(CigarOperator.H)))
						{
						L1.remove(L1.size()-1);
						}
					while(!L2.isEmpty() && L2.get(0).read0<=L1.get(L1.size()-1).read0)
						{
						L2.remove(0);
						}
					
					List<CigarElement> cigarElements=new ArrayList<CigarElement>();
					int i=0;
					while(i< L1.size())
						{
						int j=i+1;
						while(j< L1.size() && L1.get(i).op.equals(L1.get(j).op))
							{
							j++;
							}
						cigarElements.add(new CigarElement(j-i, L1.get(i).op));
						i=j;
						}
					
					//add 'N'
					cigarElements.add(
							new CigarElement((L2.get(0).ref1-L1.get(L1.size()-1).ref1)-1,
							CigarOperator.N)
							);
					
					
					i=0;
					while(i< L2.size())
						{
						int j=i+1;
						while(j< L2.size() && L2.get(i).op.equals(L2.get(j).op))
							{
							j++;
							}
						cigarElements.add(new CigarElement(j-i, L2.get(i).op));
						i=j;
						}
					
					//cleanup : case where  'S' is close to 'N'
					i=0;
					while(i+1< cigarElements.size())
						{
						CigarElement ce1=cigarElements.get(i);
						CigarElement ce2=cigarElements.get(i+1);
						
						if( i>0 &&
							ce1.getOperator().equals(CigarOperator.S) && 
							ce2.getOperator().equals(CigarOperator.N) )
							{
							cigarElements.set(i, new CigarElement(
								ce1.getLength(),
								CigarOperator.X));
							}
						else if(i+2 < cigarElements.size() &&
								ce1.getOperator().equals(CigarOperator.N) && 
								ce2.getOperator().equals(CigarOperator.S) )
							{
							cigarElements.set(i+1, new CigarElement(
								ce2.getLength(),
								CigarOperator.X));
							}
						i++;
							
						}

					
					
					newrec.setCigar(new Cigar(cigarElements));
					List<SAMValidationError> validations=newrec.isValid();
					if(validations!=null)
						{
						for(SAMValidationError err:validations)
							{
							LOG.warning(err.getType()+":"+err.getMessage());
							}
						}
					w.addAlignment(newrec);	
					}
				if(!found_one)
					{
					if(!print_only_spit_read)  w.addAlignment(rec);	
					}
				}
			iter.close();
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(r);
			CloserUtil.close(w);
			}
		}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new BWAMemNOp().instanceMainWithExit(args);
	}

}
