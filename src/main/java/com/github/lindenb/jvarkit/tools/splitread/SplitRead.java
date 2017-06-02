package com.github.lindenb.jvarkit.tools.splitread;

import java.util.Iterator;
import java.util.List;

import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.OtherCanonicalAlign;
import com.github.lindenb.jvarkit.util.picard.OtherCanonicalAlignFactory;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMRecord;

public class SplitRead extends Launcher{
	private static final Logger LOG=Logger.build(SplitRead.class).make();
	
	
	private float maxFractionCommon=0.1f;
	
	private class Fragment
		{
		String chrom;
		int pos;
		char strand;
		String cigar;
		
		public int compareTo(Fragment other)
			{
			int i=chrom.compareTo(other.chrom);
			if(i!=0) return i;
			return pos-other.pos;
			}
		
		void print()
			{
			System.out.print(chrom+"\t"+pos+"\t"+strand+"\t"+cigar);
			}
		}
	
	
	private void scanRecord(final SAMRecord record,final OtherCanonicalAlignFactory xPalignFactory) throws Exception
		{
		if(record.getReadUnmappedFlag()) return;
		String xp=record.getStringAttribute("XP");
		if(xp==null) return;
		LOG.info(xp);
		Cigar cigar1=record.getCigar();
		int readPos=0;
		int readMap1[]=new int[record.getReadLength()];
		for(CigarElement ce:cigar1.getCigarElements())
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
		for(OtherCanonicalAlign xpAln:xPalignFactory.getXPAligns(record))
			{
			
			readPos=0;
			float common=0f;
			for(CigarElement ce:xpAln.getCigarElements())
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
							if(readMap1[readPos]==1)
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

			
			Fragment f1=new Fragment();
			f1.chrom=record.getReferenceName();
			f1.pos=record.getAlignmentStart();
			f1.strand=(record.getReadNegativeStrandFlag()?'-':'+');
			f1.cigar=record.getCigarString();
			
			Fragment f2=new Fragment();
			f2.chrom=xpAln.getReferenceName();
			f2.pos=xpAln.getAlignmentStart();
			f2.strand=xpAln.getReadNegativeStrandFlag()?'-':'+';
			f2.cigar=xpAln.getCigarString();

			System.out.print(
				record.getReadName()+"\t"+
				(record.getFirstOfPairFlag()?'F':'R')+"\t"
				);
			if(f1.compareTo(f2)<0)
				{
				f1.print();
				System.out.print("\t");
				f2.print();
				}
			else
				{
				f2.print();
				System.out.print("\t");
				f1.print();
				}	
			System.out.println("\t"+fraction);
			}
		}
	

	private void scan(SamReader reader) throws Exception
		{
		OtherCanonicalAlignFactory xpalignFactory=new OtherCanonicalAlignFactory(reader.getFileHeader());
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
			scanRecord(record,xpalignFactory);
			}
		}
	
	@Override
	public int doWork(List<String> args) {		
		int optind=0;
	
		try
			{
			scan(super.openSamReader(oneFileOrNull(args)));
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		}
		
	public static void main(String[] args) throws Exception
		{
		new SplitRead().instanceMainWithExit(args);
		}
	
	}
