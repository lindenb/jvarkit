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


*/
package com.github.lindenb.jvarkit.tools.fastq;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;

import htsjdk.samtools.fastq.BasicFastqWriter;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;

import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

public class BamToFastq
	extends AbstractBamToFastq
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(BamToFastq.class);

	
	private static class MappedFastq
		{
		byte side=0;
		//int hash;
		String name;
		String seq;
		String qual;
		@Override
		public String toString() {
			return "("+name+":"+(int)side+" "+seq+" "+qual+")";
			}
		}
	
	private static class MappedFastqComparator
		implements Comparator<MappedFastq>
		{
		@Override
		public int compare(final MappedFastq o1, final  MappedFastq o2)
			{
			//int i= o1.hash - o2.hash;
			//if(i!=0) return i;
			int i= o1.name.compareTo(o2.name);
			if(i!=0) return i;
			i= (int)o1.side-(int)o2.side;
			if(i==0) 
				{
				System.err.println(o1+" "+o2);
				}
			return i;
			}
		
		}
	
	private static class MappedFastqCodec extends AbstractDataCodec<MappedFastq>
		{
		@Override
		public void encode(final  DataOutputStream dos, final MappedFastq o)
				throws IOException
			{
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
				m.side=dis.readByte();
			} catch (IOException e) {
				return null;
				}
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
		out.println(" -t (dir) "+getMessageBundle("add.tmp.dir")+" . Optional.");
		out.println(" -F (fastq) Save fastq_R1 to file (default: stdout) . Optional.");
		out.println(" -R (fastq) Save fastq_R2 to file (default: interlaced with forward) . Optional.");
		out.println(" -r  repair: insert missing read");
		out.println(" -N (int) "+getMessageBundle("max.records.in.ram")+". Optional.");
		super.printOptions(out);
		}

	@Override
	protected Collection<Throwable> call(String inputName) throws Exception {
		SamReader sfr=null;
		SortingCollection<MappedFastq> fastqCollection=null;
		try
			{
			boolean found_single=false;
			boolean found_paired=false;
			long non_primary_alignmaned_flag=0L;
			
			sfr = openSamReader(inputName);
			
			fastqCollection = SortingCollection.newInstance(
					MappedFastq.class,
					new MappedFastqCodec(),
					new MappedFastqComparator(),
					super.maxRecordsInRam,
					getTmpDirectories()
					);
			fastqCollection.setDestructiveIteration(true);

			SAMRecordIterator iter=sfr.iterator();
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(sfr.getFileHeader().getSequenceDictionary());
			while(iter.hasNext())
				{
				final SAMRecord rec=progress.watch(iter.next());
				

				
				if(rec.isSecondaryOrSupplementary())
					{
					if(non_primary_alignmaned_flag==0)
						{
						LOG.warn("SKIPPING NON-PRIMARY "+(non_primary_alignmaned_flag+1)+" ALIGNMENTS");
						}
					non_primary_alignmaned_flag++;
					continue;
					}
				
				MappedFastq m=new MappedFastq();
				m.name=rec.getReadName();
				if(m.name==null)m.name="";
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
					LOG.error("length(seq)!=length(qual) in "+m.name);
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
						throw new RuntimeException("input is a mix of paired/singled reads");
						}
					m.side=(byte)(rec.getSecondOfPairFlag()?2:1);
					}
				else
					{
					found_single=true;
					if(found_paired )
						{
						sfr.close();
						throw new RuntimeException("input is a mix of paired/singled reads");
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
			LOG.info("Done reading.");
			
			if(found_paired) 
				{
				FastqWriter fqw1=null;
				FastqWriter fqw2=null;
				if(forwardFile!=null)
					{
					LOG.info("Writing to "+forwardFile);
					fqw1=new BasicFastqWriter(forwardFile);
					}
				else
					{
					LOG.info("Writing to stdout");
					fqw1=new BasicFastqWriter(new PrintStream(System.out));
					}
				if(reverseFile!=null)
					{
					LOG.info("Writing to "+reverseFile);
					fqw2=new BasicFastqWriter(reverseFile);
					}
				else
					{
					LOG.info("Writing to interlaced stdout");
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
								LOG.warn("WTF :"+row);
								}
							boolean found_F=false;
							boolean found_R=false;
							for(MappedFastq m:row)
								{
								switch((int)m.side)
									{
									case 1:
										if(found_F) throw new RuntimeException("two forward reads found for "+row.get(0).name);
										found_F=true;
										echo(fqw1,m);
										break;
									case 2:
										if(found_R) throw new RuntimeException("two reverse reads found for "+row.get(0).name);
										found_R=true;
										echo(fqw2,m);
										break;
									default: throw new IllegalStateException("uh???");
									}
								
								}
							if(!found_F)
								{
								if(super.repair_missing_read)
									{
									LOG.warn("forward not found for "+row.get(0));
									MappedFastq pad=new MappedFastq();
									pad.side=(byte)1;
									pad.name=row.get(0).name;
									pad.seq="N";
									pad.qual="#";
									echo(fqw1,pad);
									}
								else
									{
									throw new RuntimeException("forward not found for "+row);
									}
								}
							if(!found_R)
								{
								if(repair_missing_read)
									{
									LOG.warn("reverse not found for "+row.get(0));
									MappedFastq pad=new MappedFastq();
									pad.side=(byte)2;
									pad.name=row.get(0).name;
									pad.seq="N";
									pad.qual="#";
									echo(fqw2,pad);
									}
								else
									{
									throw new RuntimeException("reverse not found for "+row);
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
					LOG.info("Writing to "+forwardFile);
					fqw1=new BasicFastqWriter(forwardFile);
					}
				else
					{
					LOG.info("Writing to stdout");
					fqw1=new BasicFastqWriter(new PrintStream(stdout()));
					}
			
				final CloseableIterator<MappedFastq> r=fastqCollection.iterator();
				while(r.hasNext())
					{
					echo(fqw1,r.next());
					}
				r.close();
				fqw1.close();
				}
			return RETURN_OK;
			}
		catch(final Exception err)
			{
			return wrapException(err);
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
