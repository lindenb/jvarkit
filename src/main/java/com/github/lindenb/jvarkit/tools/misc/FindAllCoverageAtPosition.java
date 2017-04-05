/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collectors;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

/**
 BEGIN_DOC
 
## Example

```
$ find ./  -type f -name "*.bam" |\
   java -jar dist/findallcoverageatposition.jar -p "chr2:1234" 


#File   CHROM   POS SAMPLE  DEPTH   M   I   D   N   S   H   P   EQ  X
/path/to/Sample1.bam    2   1234    SAMPLE1 10  10  0   1   0   0   0   0   0   0
/path/to/Sample2.bam    2   1234    SAMPLE2 10  0   0   0   1   0   0   0   5   5
/path/to/Sample3.bam    2   1234    SAMPLE3 10  10  0   0   0   0   0   0   0   0
```

## See also

 * https://twitter.com/pjacock/status/538300664334798848
 * https://twitter.com/yokofakun/status/538300434109456385
 * https://twitter.com/pjacock/status/538299549455233024
 * FindAVariation

## History

 * 2017: moved to jcommander

END_DOC
 */
@Program(name="findallcoverageatposition",keywords={"bam","coverage","search","depth"},description="Find depth at specific position in a list of BAM files. My colleague Estelle asked: in all the BAM we sequenced, can you give me the depth at a given position ?")
public class FindAllCoverageAtPosition extends Launcher
	{
	private static final Logger LOG = Logger.build(FindAllCoverageAtPosition.class).make();
	private static final char INSERTION_CHAR='^';
	private static final char DELETION_CHAR='-';
	private static final char BASES_To_PRINT[]=new char[]{'A','C','G','T','N',INSERTION_CHAR,DELETION_CHAR};

	@Parameter(names={"-p","--position"},description="-p chrom:pos . Multiple separated by space. Add this chrom/position. Required")
	private String positionStr = "";

	@Parameter(names={"-f","--posfile"},description="File containing positions. if file suffix is '.bed': all positions in the range will be scanned.")
	private File positionFile = null;

	@Parameter(names={"-o","--out"},description="output file. Default: stdout")
	private File outputFile = null;

	
	private static class Mutation implements Comparable<Mutation>
		{
		final String chrom;
		final int pos;
		Mutation(final String s) {
			int colon=s.indexOf(':');
			if(colon==-1 || colon+1==s.length())
				{
				throw new IllegalArgumentException("Bad chrom:pos "+s);
				}
			
			this.chrom=s.substring(0,colon).trim();
			if(chrom.isEmpty())
				{
				throw new IllegalArgumentException("Bad chrom:pos "+s);
				}
			this.pos=Integer.parseInt(s.substring(colon+1));
			}
		
		Mutation(String chrom,int pos)
			{
			this.chrom=chrom;
			this.pos=pos;
			}
		
		@Override
		public int compareTo(Mutation o) {
			int i=this.chrom.compareTo(o.chrom);
			if(i!=0) return i;
			i=this.pos-o.pos;
			return i;
			}
		
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + chrom.hashCode();
			result = prime * result + pos;
			return result;
			}
		@Override
		public boolean equals(Object obj) {
			if (this == obj)return true;
			Mutation other = (Mutation) obj;
			if (pos != other.pos) return false;
			 if (!chrom.equals(other.chrom))
				return false;
			
			return true;
		}
		
		@Override
		public String toString() {
			return chrom+":"+pos;
			}
		
		}
	
	private static class CigarAndBases
		{
		Counter<CigarOperator> operators = new Counter<>();
		Counter<Character> bases = new Counter<>();
		}
	private PrintWriter out=null;
	private SamReaderFactory samReaderFactory;
    public FindAllCoverageAtPosition()
    	{
    	samReaderFactory=  SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
    	}		
    
    private Mutation convertFromSamHeader(File f,SAMFileHeader h,final Mutation src)
    	{
    	SAMSequenceDictionary dict=h.getSequenceDictionary();
    	if(dict==null)
    		{
    		LOG.warn("No dictionary in "+h);
    		return null;
    		}
    	SAMSequenceRecord rec=dict.getSequence(src.chrom);
    	if(rec!=null) return src;
    	String chromName=src.chrom;
		if(chromName.startsWith("chr"))
			{
			rec=dict.getSequence(chromName.substring(3));
			if(rec!=null) return new Mutation(rec.getSequenceName(),src.pos);
			}
		else
			{
			rec=dict.getSequence("chr"+chromName);
			if(rec!=null) return new Mutation(rec.getSequenceName(),src.pos);
			}
	
		if(chromName.equals("MT") &&  dict.getSequence("chrM")!=null)
			{
			return new Mutation("chrM",src.pos);
			}
		if(chromName.equals("chrM") &&  dict.getSequence("MT")!=null)
			{
			return new Mutation("MT",src.pos);
			}
    	return null;
    	}

    private void scan(final BufferedReader in,final Set<Mutation> mutations) throws Exception
    	{

    	final String DEFAULT_SAMPLE_NAME="(undefined)";
    	String line;
    	while((line=in.readLine())!=null)
			{
    		if(out.checkError()) break;
			if(line.isEmpty() || line.startsWith("#")) continue;
			File f=new File(line);
			if(!f.exists()) continue;
			if(!f.isFile()) continue;
			if(!f.canRead()) continue;
			String filename=f.getName();
			if(filename.endsWith(".cram"))
				{
				LOG.warn("Sorry CRAM is not supported "+filename);
				continue;
				}
			if(!filename.endsWith(".bam")) continue;
			
    			
			SamReader samReader=null;
			SAMRecordIterator iter=null;
			try
				{
				samReader = this.samReaderFactory.open(f);
				if(!samReader.hasIndex())
					{
					LOG.warn("no index for "+f);
					continue;
					}
				final SAMFileHeader header=samReader.getFileHeader();
				for(final Mutation src:mutations)
					{
					final Map<String, CigarAndBases> sample2count=new TreeMap<String,CigarAndBases>();
					for(SAMReadGroupRecord rg:header.getReadGroups())
						{
						if(rg!=null)
							{
							String sn= rg.getSample();
							if(sn!=null && !sn.trim().isEmpty())
								{
								sample2count.put(sn, new CigarAndBases());
								}
							}
						}
					
					if(sample2count.isEmpty())
						{
						sample2count.put(DEFAULT_SAMPLE_NAME, new CigarAndBases());
						}

					
					
					final Mutation m = convertFromSamHeader(f,header,src);
					if(m==null) continue;
					iter=samReader.query(m.chrom, m.pos-1, m.pos+1,	false);
					while(iter.hasNext())
						{
						final SAMRecord rec=iter.next();
						if(rec.getReadUnmappedFlag()) continue;
						if(rec.getReadFailsVendorQualityCheckFlag()) continue;
						if(rec.getNotPrimaryAlignmentFlag()) continue;
						if(rec.isSecondaryOrSupplementary()) continue;
						if(rec.getDuplicateReadFlag()) continue;
						final Cigar cigar=rec.getCigar();
						if(cigar==null) continue;
						final String readString = rec.getReadString().toUpperCase();
						String sampleName=DEFAULT_SAMPLE_NAME;
						final SAMReadGroupRecord rg=rec.getReadGroup();
						if(rg!=null)
							{
							String sn= rg.getSample();
							if(sn!=null && !sn.trim().isEmpty())
								{
								sampleName=sn;
								}
							}
						CigarAndBases counter= sample2count.get(sampleName);
						if(counter==null)
							{
							counter=new CigarAndBases();
							sample2count.put(sampleName, counter);
							}	
						
						
						int ref= rec.getUnclippedStart();
						int readPos = 0;
						for(int k=0;k<cigar.numCigarElements() && ref< m.pos+1;++k)
							{
							final CigarElement ce=cigar.getCigarElement(k);
							final CigarOperator op=ce.getOperator();
							switch(op)
								{
								case P: break;
								case I: 
									{
									if(ref==m.pos)
										{
										counter.operators.incr(op);
										counter.bases.incr(INSERTION_CHAR);
										}
									readPos += ce.getLength();
									break;
									}
								case D:case N:
								case M: case X: case EQ: 
								case H:
								case S:
									{
									for(int i=0;i< ce.getLength();++i )
										{
										if(ref==m.pos)
											{
											counter.operators.incr(op);
											switch(op)
												{
												case M:case X:case EQ:
													counter.bases.incr(readString.charAt(readPos));
													break;
												case D:case N:
													counter.bases.incr(DELETION_CHAR);
													break;
												default:break;
												}
											break;
											}	
										if(op.consumesReadBases()) ++readPos;
										ref++;
										}
									break;
									}
								default: throw new RuntimeException("unknown operator:"+op);
								}
							}
						}
					iter.close();
					iter=null;
						
					
					for(final String sample:sample2count.keySet())
						{
						final CigarAndBases counter= sample2count.get(sample);
						
						out.print(f);
						out.print('\t');
						out.print(m.chrom);
						out.print('\t');
						out.print(m.pos);
						out.print('\t');
						out.print(sample);
						out.print('\t');
						out.print(
								counter.operators.count(CigarOperator.M)+
								counter.operators.count(CigarOperator.EQ)+
								counter.operators.count(CigarOperator.X)
								);
						for(CigarOperator op:CigarOperator.values())
							{
							out.print('\t');
							out.print(counter.operators.count(op));
							}
						for(char c:BASES_To_PRINT)
							{
							out.print('\t');
							out.print(counter.bases.count(c));
							}
						
						out.println();
						}
					}//end of loop over mutations
				}
			catch(Exception err)
				{
				LOG.error(err);
				throw err;
				}
			finally
				{
				CloserUtil.close(iter);
				CloserUtil.close(samReader);
				}    				
			}		
	    
    	}
    
    @Override
    public int doWork(final List<String> args) {
		final Set<Mutation> mutations=new TreeSet<>();

		
		BufferedReader r = null;
		try
			{
			mutations.addAll(
					Arrays.asList(this.positionStr.split("[  ]")).
						stream().filter(S->!S.trim().isEmpty())
						.map(S->new Mutation(S)).
						collect(Collectors.toSet())
					);
			
			for(final String s:positionStr.split("[  ]"))
				{
				if(s.trim().isEmpty()) continue;
				mutations.add(new Mutation(s));
				}
			if(this.positionFile!=null) {
				String line;
				r= IOUtils.openFileForBufferedReading(positionFile);
				if(positionFile.getName().endsWith(".bed")) {
					final BedLineCodec codec = new BedLineCodec();
					while((line=r.readLine())!=null) {
						final BedLine bedLine =codec.decode(line);
						if(bedLine==null) continue;
						for(int x=bedLine.getStart();x<=bedLine.getEnd();++x)
							{
							mutations.add(new Mutation(bedLine.getContig(),x));
							}
						}
					}
				else
					{
					while((line=r.readLine())!=null) {
						if(line.trim().isEmpty() || line.startsWith("#")) continue;
						mutations.add(new Mutation(line));
					}
					}
				
				r.close();
				r=null;
			}
		
		
			if(mutations.isEmpty())
				{
				LOG.fatal( "undefined position \'str\'");
				return -1;
				}
		
			LOG.info("number of mutations "+mutations.size());
			
			
			this.out=this.openFileOrStdoutAsPrintWriter(this.outputFile);
			
			out.print("#File");
			out.print('\t');
			out.print("CHROM");
			out.print('\t');
			out.print("POS");
			out.print('\t');
			out.print("SAMPLE");
			out.print('\t');
			out.print("DEPTH");
			for(final CigarOperator op:CigarOperator.values())
				{
				out.print('\t');
				out.print(op.name());
				}
			for(char c:BASES_To_PRINT)
				{
				out.print('\t');
				out.print("Base("+c+")");
				}
			
			out.println();

			
			if(args.isEmpty())
				{
				LOG.info("Reading from stdin");
				r = new BufferedReader(new InputStreamReader(stdin()));
				scan(r,mutations);
				r.close();
				r=null;
				}
			else
				{				
				for(final String filename: args)
					{
					LOG.info("Reading from "+filename);
					r=IOUtils.openURIForBufferedReading(filename);
					scan(r,mutations);
					r.close();
					r=null;
					}
				}
			this.out.flush();
			return 0;
			}
		catch(Exception err)
			{
			LOG.severe(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.out);
			CloserUtil.close(r);
			}
		}
	
	/**
	 * main
	 */
	public static void main(String[] args) {
		new FindAllCoverageAtPosition().instanceMainWithExit(args);

	}

}
