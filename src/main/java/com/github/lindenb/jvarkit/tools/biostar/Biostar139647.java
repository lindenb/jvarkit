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
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.DefaultSAMRecordFactory;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.tribble.readers.LineIterator;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

/**
BEGIN_DOC 
 
## Example

```bash
$ curl -sL "https://raw.githubusercontent.com/suryasaha/Pred_cutoff/60a6f980c9940dfb6e381c5394918f27cb14564f/data/Xylella-RpoH.aln" |\
  java -jardist/biostar139647.jar

@HD	VN:1.4	SO:unsorted
@SQ	SN:chrUn	LN:42
@PG	ID:0	VN:3a0c4ccb05e7492382e00328ac60951f215d9400	CL:(empty)	PN:Biostar139647
1	0	chrUn	1	60	42M	*	0	0	CATACTTGGTCATCGGTCGTGTCCTTGAAAGTGACTTGTTAA	*
2	0	chrUn	1	60	42M	*	0	0	TCTCTGAACCCCCTTGAAACCCCTACACTCAGCCATATATGC	*
3	0	chrUn	1	60	42M	*	0	0	TACCTTCGGGTCCTTGAAAATAGCGTCGCCGTGCTTATCTGT	*
4	0	chrUn	1	60	5M2D35M	*	0	0	TTGACAGCCGCTTGAGCAGGCGTCGGTCATCCCCACATTC	*
5	0	chrUn	1	60	18M1D9M1D13M	*	0	0	ATGCCTGGGTGGCTTGAAAGCTGGCGGCTTGCCCACATAC	*
6	0	chrUn	1	60	20M1D21M	*	0	0	TCAGTTTTATCGCTTGATATTCACTGAGACTGGCCACACAT	*

```

```
$ curl -sL "https://raw.github.com/biopython/biopython/master/Tests/Clustalw/opuntia.aln" 2> /dev/null | java -jar dist-1.128/biostar139647.jar 
@HD	VN:1.4	SO:unsorted
@SQ	SN:chrUn	LN:156
@PG	ID:0	VN:3a0c4ccb05e7492382e00328ac60951f215d9400	CL:(empty)	PN:Biostar139647
gi|6273285|gb|AF191659.1|AF191	0	chrUn	1	60	56M10D90M	*	0	0	TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATAATATATTTCA
AATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCCATTGATTTAGTGTACCAGA*
gi|6273284|gb|AF191658.1|AF191	0	chrUn	1	60	58M8D90M	*	0	0	TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATAATATATTT
CAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA	*
gi|6273289|gb|AF191663.1|AF191	0	chrUn	1	60	60M6D90M	*	0	0	TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATAATATAT
TTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTATACCAGA	*
gi|6273291|gb|AF191665.1|AF191	0	chrUn	1	60	156M	*	0	0	TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATATATATAATATATTT
CAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA	*
gi|6273287|gb|AF191661.1|AF191	0	chrUn	1	60	56M10D90M	*	0	0	TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATAATATATTTCA
AATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA*
gi|6273286|gb|AF191660.1|AF191	0	chrUn	1	60	56M10D90M	*	0	0	TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATAATATATTTAT
AATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA*
gi|6273290|gb|AF191664.1|AF191	0	chrUn	1	60	60M6D90M	*	0	0	TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATAATATAT
TTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA	*
```

END_DOC
 */

@Program(name="biostar139647",
	description="Convert alignment in Fasta/Clustal format to SAM/BAM file",
	biostars= 139647,
	keywords={"msa","sam","bam","clustal"}
	)
public class Biostar139647 extends Launcher
	{
	private static final Logger LOG = Logger.build(Biostar139647.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-R","--refname"},description="reference name. Optional")
	private String REF = "chrUn";
	
	@ParametersDelegate
	private WritingBamArgs writingBamArgs=new WritingBamArgs();

	private static final char CLIPPING=' ';
	private int align_length=0;
	private Map<String, AlignSequence> sample2sequence=new HashMap<String, AlignSequence>();
	//private AbstractSequence consensus=null;
	private enum Format{None,Clustal,Fasta};
	
	private abstract class AbstractSequence
		{
		abstract char at(int index);
		@Override
		public String toString()
			{
			final StringBuilder b=new StringBuilder(align_length);
			for(int i=0;i< align_length;++i) b.append(at(i));
			return b.toString();
			}
		}
	
	
	private abstract class Sequence extends AbstractSequence
		{
		final StringBuilder seq=new StringBuilder();
		char at(int index)
			{
			return(index< 0 || index >=seq.length()?CLIPPING:Character.toUpperCase(seq.charAt(index)));
			}
		}
	
	private class AlignSequence extends Sequence
		{
		String name;
		}
	
	
	@Override
	public int doWork(final List<String> args) {
		SAMFileWriter w=null;
		LineIterator r=null;
		try
			{
			final String filename = super.oneFileOrNull(args);
			if(filename==null)
				{
				LOG.info("Reading from stdin");
				r=IOUtils.openStreamForLineIterator(stdin());
				}
			else
				{
				LOG.info("Reading from "+filename);
				r=IOUtils.openURIForLineIterator(filename);
				}
			
			Format format=Format.None;

			
			while(r.hasNext() && format==Format.None)
				{
				String line=r.peek();
				if( line.trim().isEmpty()) { r.next(); continue;}
				if(line.startsWith("CLUSTAL"))
					{
					format=Format.Clustal;
					r.next();//consume
					break;
					}
				else if(line.startsWith(">"))
					{
					format=Format.Fasta;
					break;
					}
				else
					{
					LOG.error("MSA format not recognized in "+line);
					return -1;
					}
				}
			LOG.info("format : "+format);
			if(Format.Fasta.equals(format))
				{
				//this.consensus=new FastaConsensus();
				AlignSequence curr=null;
				while(r.hasNext())
					{
					String line=r.next();
					if(line.startsWith(">"))
						{
						curr=new AlignSequence();
						curr.name=line.substring(1).trim();
						if(sample2sequence.containsKey(curr.name))
							{
							LOG.error("Sequence ID "+curr.name +" defined twice");
							return -1;
							}
						sample2sequence.put(curr.name, curr);
						}
					else if(curr!=null)
						{
						curr.seq.append(line.trim());
						this.align_length=Math.max(this.align_length, curr.seq.length());
						}
					}
				}
			else if(Format.Clustal.equals(format))
				{
				AlignSequence curr=null;
				int columnStart=-1;
				while(r.hasNext())
					{
					String line=r.next();
					
					if( line.trim().isEmpty() || line.startsWith("CLUSTAL W"))
						{
						columnStart=-1;
						continue;
						}
					if(line.charAt(0)==' ')
						{
						if(columnStart==-1)
							{
							LOG.error("illegal consensus line for "+line);
							return -1;
							}	
						}
					else
						{
						 if(columnStart==-1)
							 {
							columnStart=line.indexOf(' ');
							if(columnStart==-1)
								{
								LOG.error("no whithespace in "+line);
								return -1;
								}
							while(columnStart< line.length() && line.charAt(columnStart)==' ')
								{
								columnStart++;
								}
							}
						String seqname=line.substring(0, columnStart).trim();
						curr=this.sample2sequence.get(seqname);
						if(curr==null)
							{
							curr=new AlignSequence();
							curr.name=seqname;
							this.sample2sequence.put(curr.name, curr);
							}
						curr.seq.append(line.substring(columnStart));
						this.align_length=Math.max(align_length, curr.seq.length());
						}
					}
				}
			else
				{
				LOG.error("Undefined input format");
				return -1;
				}
			final SAMFileHeader header=new SAMFileHeader();
			final SAMSequenceDictionary dict=new SAMSequenceDictionary();
			dict.addSequence(new SAMSequenceRecord(REF,this.align_length));
			header.setSortOrder(SortOrder.unsorted);
			header.setSequenceDictionary(dict);
			final SAMProgramRecord pgr=header.createProgramRecord();
			pgr.setProgramName(getProgramName());
			pgr.setProgramVersion(getVersion());
			if(getProgramCommandLine().trim().isEmpty())
				{
				pgr.setCommandLine("(empty)");
				}
			else
				{
				pgr.setCommandLine(getProgramCommandLine());
				}
			w = this.writingBamArgs.openSAMFileWriter(this.outputFile, header, false);
			
			final DefaultSAMRecordFactory samRecordFactory = new DefaultSAMRecordFactory();
			for(final  String seqName: this.sample2sequence.keySet())
				{
				final  AlignSequence shortRead = this.sample2sequence.get(seqName);
				final  SAMRecord rec = samRecordFactory.createSAMRecord(header);
				
				
				rec.setReadName(seqName);
				rec.setReadString(shortRead.seq.toString().replaceAll("[\\*\\- ]+", ""));
				int start=0;
				while(start< shortRead.seq.length() && !Character.isLetter(shortRead.at(start)))
					{
					start++;
					}
				rec.setAlignmentStart(start+1);
				
				int end=shortRead.seq.length()-1;
				while(end>0 && !Character.isLetter(shortRead.at(end)))
					{
					--end;
					}
				//rec.setAlignmentEnd(end+1); not supported
				
				
				int n=start;
				StringBuilder dna=new StringBuilder();
				final List<CigarElement> cigars =new ArrayList<>();
				while(n<=end )
					{
					if( !Character.isLetter(shortRead.at(n)))
						{
						cigars.add(new CigarElement(1,CigarOperator.D));
						}
					else
						{
						cigars.add(new CigarElement(1,CigarOperator.M));
						dna.append(shortRead.at(n));
						}
					n++;
					}
				
				//simplify cigar string
				n=0;
				while(n+1< cigars.size())
					{
					final CigarElement c1= cigars.get(n);
					final CigarElement c2= cigars.get(n+1);
					
					if(c1.getOperator().equals(c2.getOperator()))
						{
						cigars.set(n, new CigarElement(c1.getLength()+c2.getLength(), c1.getOperator()));
						cigars.remove(n+1);
						}
					else
						{
						++n;
						}
					}
				
				rec.setReadPairedFlag(false);
				rec.setMappingQuality(60);
				rec.setReferenceName(REF);
				rec.setReadString(dna.toString());
				rec.setCigar(new Cigar(cigars));
				rec.setAttribute("PG", pgr.getId());
				
				w.addAlignment(rec);
				}
			
			return RETURN_OK;
			}
		catch(final Exception err)
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
	
	
	public static void main(final String[] args) {
		new Biostar139647().instanceMainWithExit(args);
	}
}
