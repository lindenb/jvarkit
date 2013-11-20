package com.github.lindenb.jvarkit.tools.fastq;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import net.sf.picard.PicardException;
import net.sf.picard.fastq.BasicFastqWriter;
import net.sf.picard.fastq.FastqRecord;
import net.sf.picard.fastq.FastqWriter;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.CloserUtil;
import net.sf.samtools.util.SortingCollection;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;

public class BamToFastq
	extends AbstractCommandLineProgram
	{
	
	@Override
	public String getProgramDescription() {
		return "Implementation of https://twitter.com/DNAntonie/status/402909852277932032 " +
				"Shrink your FASTQ.bz2 files by 40+% using this one weird tip -> order them by alignment to reference before compression";
		}
	
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/BamToFastq";
		}
	
	private static class MappedFastq
		{
		int tid;
		int pos;
		byte side;
		String name;
		String seq;
		String qual;
		@Override
		public String toString() {
			return "("+name+" tid:"+tid+" pos:"+pos+" side:"+(int)side+")";
			}
		}
	
	private static class MappedFastqComparator
		implements Comparator<MappedFastq>
		{
		@Override
		public int compare(MappedFastq o1, MappedFastq o2)
			{
			int i= o1.tid - o2.tid;
			if(i!=0) return i;
			i= o1.pos - o2.pos;
			if(i!=0) return i;
			return o1.name.compareTo(o2.name);
			}
		
		}
	
	private static class MappedFastqCodec extends AbstractDataCodec<MappedFastq>
		{
		@Override
		public void encode(DataOutputStream dos, MappedFastq o)
				throws IOException
			{
			dos.writeInt(o.tid);
			dos.writeInt(o.pos);
			dos.writeByte(o.side);
			dos.writeUTF(o.name);
			dos.writeUTF(o.seq);
			dos.writeUTF(o.qual);
			}
		@Override
		public MappedFastq decode(DataInputStream dis) throws IOException
			{
			MappedFastq m=new MappedFastq();
			try {
				m.tid=dis.readInt();
			} catch (IOException e) {
				return null;
				}
			m.pos=dis.readInt();
			m.side=dis.readByte();
			m.name=dis.readUTF();
			m.seq=dis.readUTF();
			m.qual=dis.readUTF();
			
			return m;
			}
		
		@Override
		public AbstractDataCodec<MappedFastq> clone()
			{
			return new MappedFastqCodec();
			}
		}
	
	
	
	private void echo(FastqWriter fqw,MappedFastq rec)
		{
		fqw.write(new FastqRecord(
				rec.name+" "+((int)rec.side)+":N:0:1",
				rec.seq,
				new String(""),
				rec.qual
				));
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -h get help (this screen)");
		out.println(" -v print version and exit.");
		out.println(" -L (level) log level. One of java.util.logging.Level . currently:"+getLogger().getLevel());
		out.println(" -t (dir) set temporary directory . Optional.");
		out.println(" -F (fastq) Save fastq_R1 to file (default: stdout) . Optional.");
		out.println(" -R (fastq) Save fastq_R2 to file (default: interlaced with forward) . Optional.");
		out.println(" -r  repair: insert missing read");
		out.println(" -N (int) max records in memory. Optional.");
		}

	
	
	@Override
	public int doWork(String[] args)
		{
		boolean repair_missing_read=false;
		int maxRecordsInRAM=500000;
		List<File> tmpDirs=new ArrayList<File>();
		SortingCollection<MappedFastq> fastqCollection=null;
		File forwardFile=null;
		File reverseFile=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args, "hvL:F:R:N:r"))!=-1)
			{
			switch(c)
				{
				case 'h': printUsage();return 0;
				case 'v': System.out.println(getVersion());return 0;
				case 'L': getLogger().setLevel(java.util.logging.Level.parse(opt.getOptArg()));break;
				case 'F': forwardFile=new File(opt.getOptArg());break;
				case 'R': reverseFile=new File(opt.getOptArg());break;
				case 't': tmpDirs.add(new File(opt.getOptArg()));break;
				case 'N': maxRecordsInRAM=Math.max(Integer.parseInt(opt.getOptArg()),100);break;
				case 'r': repair_missing_read=true;break;
				case ':': System.err.println("Missing argument for option -"+opt.getOptOpt());return -1;
				default: System.err.println("Unknown option -"+opt.getOptOpt());return -1;
				}
			}
		if(tmpDirs.isEmpty()) tmpDirs.add(new File(System.getProperty("java.io.tmpdir")));
		SAMFileReader sfr=null;
		try
			{

			fastqCollection=SortingCollection.newInstance(
					MappedFastq.class,
					new MappedFastqCodec(),
					new MappedFastqComparator(),
					maxRecordsInRAM,
					tmpDirs
					);
			fastqCollection.setDestructiveIteration(true);
			boolean found_single=false;
			boolean found_paired=false;
			long non_primary_alignmaned_flag=0L;
			
			if(opt.getOptInd()==args.length)
				{
				info("Reading from stdin");
				sfr=new SAMFileReader(System.in);
				}
			else if(opt.getOptInd()+1==args.length)
				{
				String filename=args[opt.getOptInd()];
				sfr=new SAMFileReader(new File(filename));
				}
			else
				{
				error("Illegal parameters: check the number of arguments and the interlaced option.");
				return -1;
				}
			sfr.setValidationStringency(ValidationStringency.LENIENT);
			SAMRecordIterator iter=sfr.iterator();
			long nReads=0;
			while(iter.hasNext())
				{
				SAMRecord rec=iter.next();
				

				
				if(rec.getNotPrimaryAlignmentFlag())
					{
					if(non_primary_alignmaned_flag%maxRecordsInRAM==0)
						{
						warning("SKIPPING NON-PRIMARY "+(non_primary_alignmaned_flag+1)+" ALIGNMENTS");
						
						}
					non_primary_alignmaned_flag++;
					continue;
					}
				if(++nReads%maxRecordsInRAM==0)
					{
					info("Read "+nReads+" SAM records");
					}
				MappedFastq m=new MappedFastq();
				m.name=rec.getReadName();
				if(m.name==null)m.name="";
				m.seq=rec.getReadString();
				if(m.seq.equals(SAMRecord.NULL_SEQUENCE_STRING)) m.seq="";
				m.qual=rec.getBaseQualityString();
				if(m.qual.equals(SAMRecord.NULL_QUALS_STRING)) m.qual="";
				
				if(rec.getReadPairedFlag())
					{
					found_paired=true;
					if(found_single )
						{
						sfr.close();
						throw new PicardException("input is a mix of paired/singled reads");
						}
					m.side=(byte)(rec.getSecondOfPairFlag()?2:1);
					if(rec.getReadUnmappedFlag())
						{
						if(rec.getMateUnmappedFlag())
							{
							m.tid=Integer.MAX_VALUE-1;
							m.pos=0;
							}
						else
							{
							m.tid=rec.getMateReferenceIndex();
							m.pos=rec.getMateAlignmentStart();
							}
						}
					else
						{
						if(rec.getMateUnmappedFlag())
							{
							m.tid=rec.getReferenceIndex();
							m.pos=rec.getAlignmentStart();
							}
						else
							{
							m.tid=Math.min(rec.getReferenceIndex(),rec.getMateReferenceIndex());
							m.pos=Math.min(rec.getAlignmentStart(),rec.getMateAlignmentStart());
							}
						}
					
					}
				else
					{
					found_single=true;
					if(found_paired )
						{
						sfr.close();
						throw new PicardException("input is a mix of paired/singled reads");
						}
					m.side=(byte)0;
					if(rec.getReadUnmappedFlag())
						{
						m.tid=Integer.MAX_VALUE-1;
						m.pos=0;
						}
					else
						{
						m.tid=rec.getReferenceIndex();
						m.pos=rec.getAlignmentStart();
						}
					}
				fastqCollection.add(m);
				}
			iter.close();
			CloserUtil.close(iter);
			CloserUtil.close(sfr);
			
			fastqCollection.doneAdding();
			info("Done reading.");
			
			if(found_paired) 
				{
				FastqWriter fqw1=null;
				FastqWriter fqw2=null;
				if(forwardFile!=null)
					{
					info("Writing to "+forwardFile);
					fqw1=new BasicFastqWriter(forwardFile);
					}
				else
					{
					info("Writing to stdout");
					fqw1=new BasicFastqWriter(new PrintStream(System.out));
					}
				if(reverseFile!=null)
					{
					info("Writing to "+reverseFile);
					fqw2=new BasicFastqWriter(reverseFile);
					}
				else
					{
					info("Writing to interlaced stdout");
					fqw2=fqw1;
					}
				List<MappedFastq> row=new ArrayList<MappedFastq>();
				CloseableIterator<MappedFastq> r=fastqCollection.iterator();
				long nWrite=0;
				for(;;)
					{
					MappedFastq curr=null;
					if(r.hasNext()) curr=r.next();
					if(curr==null || (!row.isEmpty() && !row.get(0).name.equals(curr.name)))
						{
						if(!row.isEmpty())
							{
							if(row.size()>2)
								{
								warning("WTF :"+row);
								}
							boolean found_F=false;
							boolean found_R=false;
							for(MappedFastq m:row)
								{
								switch((int)m.side)
									{
									case 1:
										if(found_F) throw new PicardException("two forward reads found for "+row.get(0).name);
										found_F=true;
										echo(fqw1,m);
										break;
									case 2:
										if(found_R) throw new PicardException("two reverse reads found for "+row.get(0).name);
										found_R=true;
										echo(fqw2,m);
										break;
									default: throw new IllegalStateException("uh???");
									}
								
								}
							if(!found_F)
								{
								if(repair_missing_read)
									{
									warning("forward not found for "+row.get(0));
									MappedFastq pad=new MappedFastq();
									pad.tid=0;
									pad.pos=0;
									pad.side=(byte)1;
									pad.name=row.get(0).name;
									pad.seq="";
									pad.qual="";
									echo(fqw1,pad);
									}
								else
									{
									throw new PicardException("forward not found for "+row);
									}
								}
							if(!found_R)
								{
								if(repair_missing_read)
									{
									warning("reverse not found for "+row.get(0));
									MappedFastq pad=new MappedFastq();
									pad.tid=0;
									pad.pos=0;
									pad.side=(byte)2;
									pad.name=row.get(0).name;
									pad.seq="";
									pad.qual="";
									echo(fqw2,pad);
									}
								else
									{
									throw new PicardException("reverse not found for "+row);
									}
								}
							}
						if(curr==null) break;
						row.clear();
						}
					if(nWrite++%maxRecordsInRAM==0)
						{
						info("Wrote "+nWrite+" records");
						}
					row.add(curr);
					}
				r.close();
				fqw1.close();
				fqw2.close();
				}
			else if(found_single) 
				{
				FastqWriter fqw1=null;
				if(forwardFile!=null)
					{
					info("Writing to "+forwardFile);
					fqw1=new BasicFastqWriter(forwardFile);
					}
				else
					{
					info("Writing to stdout");
					fqw1=new BasicFastqWriter(new PrintStream(System.out));
					}
			
				CloseableIterator<MappedFastq> r=fastqCollection.iterator();
				while(r.hasNext())
					{
					echo(fqw1,r.next());
					}
				r.close();
				fqw1.close();
				}
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			if(fastqCollection!=null) fastqCollection.cleanup();
			}
		}
	public static void main(String[] args) {
		new BamToFastq().instanceMainWithExit(args);
		}
	}
