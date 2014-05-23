package com.github.lindenb.jvarkit.tools.fastq;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import com.github.lindenb.jvarkit.util.picard.PicardException;
import htsjdk.samtools.fastq.BasicFastqWriter;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.picard.SamFileReaderFactory;
import com.github.lindenb.jvarkit.util.picard.SortingCollectionFactory;

public class BamToFastq
	extends AbstractCommandLineProgram
	{
	
	@Override
	public String getProgramDescription()
		{
		return "Same as picard/SamToFastq but allow missing reads + shuffle reads using hash(name) so you can use them with bwa. Previous version was an Implementation of https://twitter.com/DNAntonie/status/402909852277932032";
		}
	
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/BamToFastq";
		}
	
	private static class MappedFastq
		{
		byte side;
		int hash;
		String name;
		String seq;
		String qual;
		@Override
		public String toString() {
			return "("+name+" "+qual+")";
			}
		}
	
	private static class MappedFastqComparator
		implements Comparator<MappedFastq>
		{
		@Override
		public int compare(MappedFastq o1, MappedFastq o2)
			{
			int i= o1.hash - o2.hash;
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
			dos.writeInt(o.hash);
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
				m.hash=dis.readInt();
			} catch (IOException e) {
				return null;
				}
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
		out.println(" -t (dir) "+getMessageBundle("illegal.number.of.arguments")+" . Optional.");
		out.println(" -F (fastq) Save fastq_R1 to file (default: stdout) . Optional.");
		out.println(" -R (fastq) Save fastq_R2 to file (default: interlaced with forward) . Optional.");
		out.println(" -r  repair: insert missing read");
		out.println(" -N (int) "+getMessageBundle("max.records.in.ram")+". Optional.");
		super.printOptions(out);
		}

	
	@Override
	public int doWork(String[] args)
		{
		boolean repair_missing_read=false;
		SortingCollectionFactory<MappedFastq> sortingFactory=new SortingCollectionFactory<MappedFastq>();
		File forwardFile=null;
		File reverseFile=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		
		sortingFactory.setComponentType(MappedFastq.class);
		sortingFactory.setCodec(new MappedFastqCodec());
		sortingFactory.setComparator(new MappedFastqComparator());
		
		while((c=opt.getopt(args,super.getGetOptDefault()+ "F:R:N:r"))!=-1)
			{
			switch(c)
				{
				case 'F': forwardFile=new File(opt.getOptArg());break;
				case 'R': reverseFile=new File(opt.getOptArg());break;
				case 't': addTmpDirectory(new File(opt.getOptArg()));break;
				case 'N': sortingFactory.setMaxRecordsInRAM(Math.max(Integer.parseInt(opt.getOptArg()),100));break;
				case 'r': repair_missing_read=true;break;
				case ':': System.err.println("Missing argument for option -"+opt.getOptOpt());return -1;
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
		SamReader sfr=null;
		SortingCollection<MappedFastq> fastqCollection=null;
		try
			{
			sortingFactory.setTmpDirs(this.getTmpDirectories());
			fastqCollection=sortingFactory.make();
			fastqCollection.setDestructiveIteration(true);
			boolean found_single=false;
			boolean found_paired=false;
			long non_primary_alignmaned_flag=0L;
			
			if(opt.getOptInd()==args.length)
				{
				sfr=SamFileReaderFactory.mewInstance().openStdin();
				}
			else if(opt.getOptInd()+1==args.length)
				{
				String filename=args[opt.getOptInd()];
				sfr=SamFileReaderFactory.mewInstance().open(filename);
				}
			else
				{
				error(getMessageBundle("illegal.number.of.arguments"));
				return -1;
				}
			SAMRecordIterator iter=sfr.iterator();
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(sfr.getFileHeader().getSequenceDictionary());
			while(iter.hasNext())
				{
				SAMRecord rec=iter.next();
				progress.watch(rec);

				
				if(rec.isSecondaryOrSupplementary())
					{
					if(non_primary_alignmaned_flag==0)
						{
						warning("SKIPPING NON-PRIMARY "+(non_primary_alignmaned_flag+1)+" ALIGNMENTS");
						}
					non_primary_alignmaned_flag++;
					continue;
					}
				
				MappedFastq m=new MappedFastq();
				m.name=rec.getReadName();
				if(m.name==null)m.name="";
				m.hash=m.name.hashCode();
				m.seq=rec.getReadString();
				
				if(m.seq.equals(SAMRecord.NULL_SEQUENCE_STRING)) m.seq="";
				m.qual=rec.getBaseQualityString();
				if(m.qual.equals(SAMRecord.NULL_QUALS_STRING)) m.qual="";
				if(!rec.getReadUnmappedFlag() && rec.getReadNegativeStrandFlag())
					{
					m.seq=AcidNucleics.reverseComplement(m.seq);
					m.qual=new StringBuilder(m.qual).reverse().toString();
					}
				if(m.seq.length()!=m.qual.length())
					{
					error("length(seq)!=length(qual) in "+m.name);
					continue;
					}
				if(m.seq.isEmpty() && m.qual.isEmpty())
					{
					m.seq="N";
					m.qual="#";
					}
				
				
				if(rec.getReadPairedFlag())
					{
					found_paired=true;
					if(found_single )
						{
						sfr.close();
						throw new PicardException("input is a mix of paired/singled reads");
						}
					m.side=(byte)(rec.getSecondOfPairFlag()?2:1);
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
					}
				fastqCollection.add(m);
				}
			iter.close();
			CloserUtil.close(iter);
			CloserUtil.close(sfr);
			progress.finish();
			
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
									pad.side=(byte)1;
									pad.name=row.get(0).name;
									pad.seq="N";
									pad.qual="#";
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
									pad.side=(byte)2;
									pad.name=row.get(0).name;
									pad.seq="N";
									pad.qual="#";
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
