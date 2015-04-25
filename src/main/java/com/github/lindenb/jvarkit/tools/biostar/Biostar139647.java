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
import com.github.lindenb.jvarkit.knime.AbstractKnimeApplication;

public class Biostar139647 extends AbstractKnimeApplication
	{
	private static final char CLIPPING=' ';
	private String REF="chrUn";
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
	protected String getOnlineDocUrl() {
		return DEFAULT_WIKI_PREFIX+ "Biostar139647";
		}
	@Override
	public String getProgramDescription() {
		return "Convert alignment in Fasta/Clustal format to SAM/BAM file see https://www.biostars.org/p/139647/";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println("-R (name) reference name. Optional");
		out.println("-o (file) output file. default: stdout");
		super.printOptions(out);
		}
	
	
	@Override
	public int executeKnime(List<String> args)
		{
		SAMFileWriter w=null;
		LineIterator r=null;
		SAMFileWriterFactory sfwf=null;
		try
			{
			if(args.isEmpty())
				{
				info("Reading from stdin");
				r=IOUtils.openStdinForLineIterator();
				}
			else if(args.size()==0)
				{
				String filename=args.get(0);
				info("Reading from "+filename);
				r=IOUtils.openURIForLineIterator(filename);
				}
			else
				{
				error(getMessageBundle("illegal.number.of.arguments"));
				return -1;
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
					error("MSA format not recognized in "+line);
					return -1;
					}
				}
			info("format : "+format);
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
							error("Sequence ID "+curr.name +" defined twice");
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
							error("illegal consensus line for "+line);
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
								error("no whithespace in "+line);
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
				error("Undefined input format");
				return -1;
				}
			SAMFileHeader header=new SAMFileHeader();
			SAMSequenceDictionary dict=new SAMSequenceDictionary();
			dict.addSequence(new SAMSequenceRecord(REF,this.align_length));
			header.setSortOrder(SortOrder.unsorted);
			header.setSequenceDictionary(dict);
			SAMProgramRecord pgr=header.createProgramRecord();
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
			
			sfwf = new SAMFileWriterFactory();
			if(getOutputFile()==null)
				{
				w = sfwf.makeSAMWriter(header, false, System.out);
				}
			else
				{
				w = sfwf.makeSAMOrBAMWriter(header, false, getOutputFile());
				}
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
			
			info("Done");
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(r);
			CloserUtil.close(w);
			}	
		}
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"R:o:"))!=-1)
			{
			switch(c)
				{
				case 'R': REF=opt.getOptArg();break;
				case 'o': setOutputFile(opt.getOptArg());break;
				default:
					{
					switch(handleOtherOptions(c, opt,args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		return mainWork(opt.getOptInd(), args);
		}
	
	public static void main(String[] args) {
		new Biostar139647().instanceMainWithExit(args);
	}
}
