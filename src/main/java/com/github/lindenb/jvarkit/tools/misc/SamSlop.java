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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMValidationError;
import htsjdk.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;


/**

BEGIN_DOC




### Example



```

$ java -jar dist/samslop.jar -m 20 -M 10 -r ~/src/gatk-ui/testdata/ref.fa  ~/src/gatk-ui/testdata/S1.bam
 
@HD     VN:1.5  GO:none SO:unsorted
@SQ     SN:rotavirus    LN:1074
@RG     ID:S1   SM:S1
@PG     ID:bwa  PN:bwa  VN:0.7.12-r1044 CL:../bwa/bwa mem -R @RG\tID:S1\tSM:S1 ref.fa S1_01_R1.fq.gz S1_01_R2.fq.gz
@PG     ID:bwa.1        PN:bwa  VN:0.7.12-r1044 CL:../bwa/bwa mem -R @RG\tID:S1\tSM:S1 ref.fa S1_02_R1.fq.gz S1_02_R2.fq.gz
@PG     ID:bwa.2        PN:bwa  VN:0.7.12-r1044 CL:../bwa/bwa mem -R @RG\tID:S1\tSM:S1 ref.fa S1_03_R1.fq.gz S1_03_R2.fq.gz
rotavirus_1_317_5:0:0_7:0:0_2de 99      rotavirus       1       60      140M    =       248     317     GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAATATGGCGTCAACTCAGCAGATGGTCAGCTCTAATATTAATACTTCTT     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##########        MD:Z:33G4A3T14A7T4      RG:Z:S1 NM:i:5  AS:i:45 XS:i:0
rotavirus_1_535_4:0:0_4:0:0_1a6 163     rotavirus       1       60      140M    =       466     535     GGCTTTTACTGCTTTTCAGTGGTTGCTTCTCAAGATGGAGTGTACTCATCAGATGGTAAGCTCTATTATTAATACTTCTT     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##########        MD:Z:8A18G13C6G21       RG:Z:S1 NM:i:4  AS:i:50 XS:i:0
rotavirus_1_543_5:0:0_11:0:0_390        163     rotavirus       1       60      140M    =       487     530     GGCTTTTAATGCTTTTCATTTGATGCTGCTCAAGATGGAGTCTACACAGCAGATGGTCAGCTCTATTATTAATACTTCTT     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##########        MD:Z:18G1G1T22T11A12    RG:Z:S1 NM:i:5  AS:i:45      XS:i:0
rotavirus_1_578_3:0:0_7:0:0_7c  99      rotavirus       1       60      140M    =       509     578     GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAAGATGGAGTCTCCTGAGCAGCTGGTAAGCTCTATTATTAATACTTCTT     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##########        MD:Z:43A2C5A17  RG:Z:S1 NM:i:3  AS:i:55 XS:i:0
rotavirus_1_497_4:0:0_5:0:0_2f6 163     rotavirus       1       60      140M    =       432     497     GGCATTTAATGCTTAACAGTGGTTGCTGCTCAAGATGGAGTCTACTCAGCAGATGGTAAGCTCTCTTATTAATACTTCTT     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##########        MD:Z:3T10T0T48A5        RG:Z:S1 NM:i:4  AS:i:51 XS:i:0
rotavirus_2_509_4:0:0_8:0:0_68  99      rotavirus       1       60      4M      =       440     478     GGCTTTTAAAGCTTTTCAGTTGTTGCTGCTCAAGATGGAGTCTACTCAGCAGATGTTACGCTCTATTATTAATACTTCTTT    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##########       MD:Z:8T10G34G2A12       RG:Z:S1 NM:i:4  AS:i:50 XS:i:0
rotavirus_2_441_11:0:0_4:0:0_298        99      rotavirus       1       60      4M      =       372     440     GGCTTTTAATGCTTTTCAGTTGTTGCTGCACAAGATGGAGTCTACACAGCAGATGGTAAGCT       #+++++++++++++++++++++++++++++++++++++++++++++++++++##########  MD:Z:19G8T15T6  RG:Z:S1 NM:i:3  AS:i:36 XS:i:0
rotavirus_2_538_5:0:0_6:0:0_6b  99      rotavirus       1       60      4M      =       469     537     GGCTTTTAATGCTTTTCAGTGGTTTCTTCTCACGATGGAGTCTACTCAGCAGAAGGTAAGCACTATTATTAATACTTCTTT    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##########       MD:Z:23G2G4A20T7T9      RG:Z:S1 NM:i:5  AS:i:45 XS:i:0
rotavirus_2_459_8:0:0_3:0:0_388 99      rotavirus       1       60      4M      =       390     458     GGCTTTTAAAGCATTACAGTTGTTGCAGCTCAAGAAGGAGACTACTCAGCAGATGGTAAGCTCTATAATTAATACTTCTTT    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##########       MD:Z:8T2T2T4G5T8T4T25T4 RG:Z:S1 NM:i:8  AS:i:30 XS:i:0
rotavirus_2_577_7:0:0_6:0:0_64  163     rotavirus       1       60      4M      =       513     576     GGCTTTTAATTCTATTCAGTGGTTGCTGCTCCAGAAGGAGTCTACTCAGGAGATGGTACGCTCTCTTATTAATACTTCTTT    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##########       MD:Z:9G2T17A3T13C8A5A6  RG:Z:S1 NM:i:7  AS:i:35 XS:i:0
rotavirus_2_500_1:0:0_14:0:0_88 163     rotavirus       1       60      4M      =       453     492     GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAAGATGGAGTCTACTCAGCAGATGGTAAGCTCAATTATTAATACTTCTTT    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##########       MD:Z:62T7       RG:Z:S1 NM:i:1  AS:i:65 XS:i:0
rotavirus_2_445_5:0:0_10:0:0_c9 163     rotavirus       1       60      4M      =       377     434     GGCTTTTAATGCTTTTCAGTTGTAGCTGCTCAAGATGGAGTCTACTCATCAGATGGTAAGCTCTCTTCTTAATACTTCTTT    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##########       MD:Z:19G2T24G15A2A3     RG:Z:S1 NM:i:5  AS:i:48 XS:i:0
rotavirus_2_488_8:0:0_9:0:0_343 99      rotavirus       1       60      4M      =       419     487     GTCTTTAAATGCTTTTCAGTGTTTGCTGCTCAAGATGGAGTCTACTCAGCAGAAGGTAAGCTCTATTATTAATA   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##########      MD:Z:0G4T14G31T10       RG:Z:S1 NM:i:4  AS:i:47 XS:i:0
rotavirus_2_478_5:0:0_5:0:0_35c 163     rotavirus       1       60      4M      =       409     477     GGCATTTAATGCTTTTCAGTGGTTGCTGCACAAGATGGAGTCTACTCAGCAGATTGTAAGCTCTATTATTAATACTT        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##########   MD:Z:2T25T24G12 RG:Z:S1 NM:i:3  AS:i:53 XS:i:0
(...)
```



END_DOC
*/
@Program(name="SamSlop",description="extends sam by 'x' bases")
public class SamSlop extends Launcher
	{
	private static final Logger LOG = Logger.build(SamSlop.class).make();


	@Parameter(names={"-o","--output"},description="Output file. Optional . Default: stdout")
	private File outputFile = null;


	@Parameter(names={"-m","--extend5"},description="num bases to extends on 5'")
	private int extend5 = 0 ;

	@Parameter(names={"-M","--extend3"},description="num bases to extends on 3'")
	private int extend3 = 0 ;

	@Parameter(names={"-q","--defaultQual"},description="default base quality (one letter)")
	private String defaultQual = "#";

	@Parameter(names={"-c","--rmClip"},description="remove clipped bases")
	private boolean removeClip = false;

	@Parameter(names={"-r","--reference"},description="Indexed fasta Reference")
	private File faidx = null;
	
	
	@ParametersDelegate
	private WritingBamArgs writingBamArgs=new WritingBamArgs();
	
	
	@Override
	public int doWork(List<String> args) {
		if(this.faidx==null)
			{
			LOG.error("Reference was not specified.");
			return -1;
			}
		if(this.defaultQual.length()!=1)
			{
			LOG.error("default quality should have length==1 "+this.defaultQual);
			}
		GenomicSequence genomicSequence=null;
		SamReader sfr=null;
		SAMFileWriter sfw=null;
		IndexedFastaSequenceFile indexedFastaSequenceFile=null;
		final char defaultQUAL=this.defaultQual.charAt(0);
		try
			{
			final String inputName= oneFileOrNull(args);
			LOG.info("Loading reference");
			indexedFastaSequenceFile=new IndexedFastaSequenceFile(faidx);
			sfr = openSamReader(inputName);
			final SAMFileHeader header=sfr.getFileHeader();
			header.setSortOrder(SortOrder.unsorted);
			sfw = writingBamArgs.openSAMFileWriter(outputFile,header, true);
			final SAMSequenceDictionaryProgress progress= new SAMSequenceDictionaryProgress(header);
			final SAMRecordIterator iter=sfr.iterator();
			while(iter.hasNext())
				{
				final SAMRecord rec=progress.watch(iter.next());
				final Cigar cigar=rec.getCigar();
				if( rec.getReadUnmappedFlag() ||
					cigar==null ||
					cigar.isEmpty() ||
					rec.getReadBases()==SAMRecord.NULL_SEQUENCE ||
					(this.extend5<=0 && this.extend3<=0)
					)
					{
					sfw.addAlignment(rec);
					continue;
					}
				
				final StringBuilder sbs = new StringBuilder(rec.getReadString());
				final StringBuilder sbq = new StringBuilder(rec.getBaseQualityString());

				
				if(genomicSequence==null ||
					genomicSequence.getSAMSequenceRecord().getSequenceIndex()!=rec.getReferenceIndex())
					{
					genomicSequence=new GenomicSequence(indexedFastaSequenceFile, rec.getReferenceName());
					}
				
				final int refPos1=(this.removeClip?rec.getAlignmentStart():rec.getUnclippedStart());
				final int endAlignmend1= (this.removeClip?rec.getAlignmentEnd():rec.getUnclippedEnd());
				final List<CigarElement> cl = new ArrayList<>(cigar.getCigarElements());
				
				if(!this.removeClip)
					{
					//replace clip S/H by match M
					for(int i=0;i< cl.size();++i)
						{
						final CigarElement ce = cl.get(i);
						if(ce.getOperator()==CigarOperator.S || ce.getOperator()==CigarOperator.H)
							{
							cl.set(i,new CigarElement(ce.getLength(), CigarOperator.M));
							}
						}
					}
				
				if(this.extend5>0 && refPos1>1)
					{
					
					if(this.removeClip)
						{
						///remove hard + soft clip 5'
						while(!cl.isEmpty())
							{
							//first
							final CigarElement ce = cl.get(0);
							//not a clip
							if(!(ce.getOperator()==CigarOperator.S || ce.getOperator()==CigarOperator.H))
								{
								break;
								}
							 if( ce.getOperator()==CigarOperator.S)
							 	{
								sbs.replace(0,ce.getLength(),"");
								sbq.replace(0,ce.getLength(),"");
							 	}
							cl.remove(0);//remove first
							}
						}
				
					final StringBuilder prefix= new StringBuilder(this.extend5);
					///append + soft clip 5'
					for(int i= 0;i< this.extend5;++i)
						{
						int x1 = (refPos1-1)-i;
						if(x1<1) break;//break if out of genome
						prefix.insert(0,genomicSequence.charAt(x1-1));
						}
					sbs.insert(0, prefix.toString());
					for(int i= 0;i<prefix.length();++i) sbq.insert(0, defaultQUAL);//preprend quality
					cl.add(0, new CigarElement(prefix.length(), CigarOperator.M));//prepend cigar
					rec.setAlignmentStart(refPos1-prefix.length());//update start pos
					}
				
				if(this.extend3>0 && rec.getAlignmentEnd()< genomicSequence.length())
					{
					if(this.removeClip)
						{
						///remove hard + soft clip 3'
						while(!cl.isEmpty())
							{
							//last
							final CigarElement ce = cl.get(cl.size()-1);
							//not a clip
							if(!(ce.getOperator()==CigarOperator.S || ce.getOperator()==CigarOperator.H))
								{
								break;
								}
							 if( ce.getOperator()==CigarOperator.S)
							 	{
								sbs.setLength(sbs.length()-ce.getLength());
								sbq.setLength(sbq.length()-ce.getLength());
							 	}
							 //remove last
							cl.remove(cl.size()-1);
							}
						}
					int extend=0;
					for(int pos1= endAlignmend1+1;
							pos1<= (endAlignmend1+this.extend3) && pos1<= genomicSequence.length() ;
							++pos1)
						{
						sbs.append(genomicSequence.charAt(pos1-1));
						sbq.append(defaultQUAL);
						++extend;
						}
					cl.add(new CigarElement(extend, CigarOperator.M));//append cigar
					}
				//simplify cigar
				int idx=0;
				while(idx+1<cl.size())
					{
					final CigarElement ce1 = cl.get(idx);
					final CigarElement ce2 = cl.get(idx+1);
					if(ce1.getOperator()==ce2.getOperator())
						{
						cl.set(idx,new CigarElement(ce1.getLength()+ce2.getLength(), ce1.getOperator()));
						cl.remove(idx+1);
						}
					else
						{
						idx++;
						}
					}
				
				
				rec.setCigar(new Cigar(cl));
				rec.setReadString(sbs.toString());
				rec.setBaseQualityString(sbq.toString());
				final List<SAMValidationError> errors = rec.isValid();
				if(errors!=null && !errors.isEmpty()) {
					for(SAMValidationError err:errors)
						{
						LOG.error(err.getMessage());
						}
				}
				
				//info("changed "+rec.getCigarString()+" to "+newCigarStr+" "+rec.getReadName()+" "+rec.getReadString());
				
				
				sfw.addAlignment(rec);
				}
			progress.finish();
			return RETURN_OK;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(indexedFastaSequenceFile);
			CloserUtil.close(sfr);
			CloserUtil.close(sfw);
			}
		}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new SamSlop().instanceMainWithExit(args);

	}

}
