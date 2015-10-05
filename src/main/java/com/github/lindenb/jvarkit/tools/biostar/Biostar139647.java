/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

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

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
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
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;

import htsjdk.tribble.readers.LineIterator;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.command.Command;

public class Biostar139647 extends AbstractBiostar139647
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(Biostar139647.class);

	
	private static final char CLIPPING=' ';
	@Override
	public Command createCommand()
		{
		return new MyCommand();
		}
	
	static private class MyCommand extends AbstractBiostar139647.AbstractBiostar139647Command
				{
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
				StringBuilder b=new StringBuilder(align_length);
				for(int i=0;i< align_length;++i) b.append(at(i));
				return b.toString();
				}
			}
		
		
		
		private abstract class Sequence extends AbstractSequence
			{
			StringBuilder seq=new StringBuilder();
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
		protected Collection<Throwable> call(final String filename) throws Exception
			{
			SAMFileWriter w=null;
			LineIterator r=null;
			try
				{
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
						return wrapException("MSA format not recognized in "+line);
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
								return wrapException("Sequence ID "+curr.name +" defined twice");
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
								return wrapException("illegal consensus line for "+line);
								}	
							}
						else
							{
							 if(columnStart==-1)
								 {
								columnStart=line.indexOf(' ');
								if(columnStart==-1)
									{
									return wrapException("no whithespace in "+line);
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
					return wrapException("Undefined input format");
					}
				SAMFileHeader header=new SAMFileHeader();
				SAMSequenceDictionary dict=new SAMSequenceDictionary();
				dict.addSequence(new SAMSequenceRecord(super.REF,this.align_length));
				header.setSortOrder(SortOrder.unsorted);
				header.setSequenceDictionary(dict);
				SAMProgramRecord pgr=header.createProgramRecord();
				pgr.setProgramName(getName());
				pgr.setProgramVersion(getVersion());
				if(getProgramCommandLine().trim().isEmpty())
					{
					pgr.setCommandLine("(empty)");
					}
				else
					{
					pgr.setCommandLine(getProgramCommandLine());
					}
				
				w = openSAMFileWriter(header, false);
				DefaultSAMRecordFactory samRecordFactory = new DefaultSAMRecordFactory();
				for(String seqName: this.sample2sequence.keySet())
					{
					AlignSequence shortRead = this.sample2sequence.get(seqName);
					SAMRecord rec = samRecordFactory.createSAMRecord(header);
					
					
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
					List<CigarElement> cigars =new ArrayList<>();
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
						CigarElement c1= cigars.get(n);
						CigarElement c2= cigars.get(n+1);
						
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
					
					w.addAlignment(rec);
					}
				
				LOG.info("Done");
				return Collections.emptyList();
				}
			catch(Exception err)
				{
				return wrapException(err);
				}
			finally
				{
				CloserUtil.close(r);
				CloserUtil.close(w);
				}	
			}
		}
	public static void main(String[] args) {
		new Biostar139647().instanceMainWithExit(args);
	}
}
