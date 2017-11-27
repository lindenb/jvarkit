package com.github.lindenb.jvarkit.tools.splitread;

import java.io.File;
import java.io.PrintStream;
import java.util.Iterator;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMUtils;

@Program(name="splitread",
	description="TODO",keywords={"sam","bam"},generate_doc=false)
public class SplitRead extends Launcher{
	private static final Logger LOG=Logger.build(SplitRead.class).make();
    @Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
    private File outputFile = null;

    private PrintStream out=null;
	private float maxFractionCommon=0.1f;
	
	private class Fragment
		implements Comparable<Fragment>
		{
		String chrom;
		int pos;
		char strand;
		String cigar;
		
		@Override public int compareTo(final Fragment other)
			{
			int i=chrom.compareTo(other.chrom);
			if(i!=0) return i;
			return pos-other.pos;
			}
		
		void print()
			{
			SplitRead.this.out.print(chrom+"\t"+pos+"\t"+strand+"\t"+cigar);
			}
		}
	
	
	private void scanRecord(final SAMRecord record) throws Exception
		{
		if(record.getReadUnmappedFlag()) return;
		final List<SAMRecord> xpAlnList = SAMUtils.getOtherCanonicalAlignments(record);
		if(xpAlnList.isEmpty()) return;
		Cigar cigar1=record.getCigar();
		int readPos=0;
		// positions in the reads that are mapped by the cigar string
		final int readMap1[]=new int[record.getReadLength()];
		for(final CigarElement ce:cigar1.getCigarElements())
			{
			switch(ce.getOperator())
				{
				case I: case S:
					{
					readPos+=ce.getLength();
					break;
					}
				case M:case X:case EQ:
					{
					for(int i=0;i< ce.getLength();++i)
						{
						readMap1[readPos]+=1;
						readPos++;
						}
					break;
					}
				case P: case H: case D: case N: break;
				default: throw new RuntimeException("cigar operator not handled:"+ce.getOperator());
				}
			}
		// lopp over alternate reads
		for(final SAMRecord xpAln:xpAlnList)
			{
			readPos=0;
			float common=0f;
			for(final CigarElement ce:xpAln.getCigar())
				{
				switch(ce.getOperator())
					{
					case I: case S:
						{
						readPos+=ce.getLength();
						break;
						}
					case M:case X:case EQ:
						{
						//get the common region between the main alignemt and the
						// this canonical align
						for(int i=0;i< ce.getLength();++i)
							{
							if(readMap1[readPos]==1)//was in main align
								{
								common++;
								}
							readPos++;
							}
						break;
						}
					case P: case H: case D: case N: break;
					default: throw new RuntimeException("cigar operator not handled:"+ce.getOperator());
					}
				}
			float fraction=common/readMap1.length;
			if(  fraction > this.maxFractionCommon)
				{
				continue;
				}

			
			final Fragment f1=new Fragment();
			f1.chrom=record.getReferenceName();
			f1.pos=record.getAlignmentStart();
			f1.strand=(record.getReadNegativeStrandFlag()?'-':'+');
			f1.cigar=record.getCigarString();
			
			final Fragment f2=new Fragment();
			f2.chrom=xpAln.getReferenceName();
			f2.pos=xpAln.getAlignmentStart();
			f2.strand=xpAln.getReadNegativeStrandFlag()?'-':'+';
			f2.cigar=xpAln.getCigarString();

			this.out.print(
				record.getReadName()+"\t"+
				(record.getFirstOfPairFlag()?'F':'R')+"\t"
				);
			if(f1.compareTo(f2)<0)
				{
				f1.print();
				this.out.print("\t");
				f2.print();
				}
			else
				{
				f2.print();
				this.out.print("\t");
				f1.print();
				}	
			this.out.println("\t"+fraction);
			}
		}
	

	private void scan(final SamReader reader) throws Exception
		{
		long nrecords=0L;
		for(Iterator<SAMRecord> iter=reader.iterator();
				iter.hasNext(); )
			{
			SAMRecord record=iter.next();
			++nrecords;
			if(nrecords%1E6==0)
				{
				LOG.info("nRecord:"+nrecords);
				System.out.flush();
				}
			scanRecord(record);
			}
		}
	
	@Override
	public int doWork(final List<String> args) {		
		
		try
			{
			this.out = openFileOrStdoutAsPrintStream(this.outputFile);
			scan(super.openSamReader(oneFileOrNull(args)));
			this.out.flush();
			this.out.close();
			this.out=null;
			return 0;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(out);
			}
		}
		
	public static void main(String[] args) throws Exception
		{
		new SplitRead().instanceMainWithExit(args);
		}
	
	}
