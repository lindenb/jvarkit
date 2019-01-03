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
package com.github.lindenb.jvarkit.tools.fastq;

import java.io.DataInputStream;

import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
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

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;


/**

BEGIN_DOC


# Deprecated: 

use picard please


# Warnings

Previous version was an Implementation of https://twitter.com/DNAntonie/status/402909852277932032

	illumina  read is filtered is always "n"
	illumina control number is always 0
	Illumina index sequence is lost.

## Example
piping bwa mem

```
$ bwa mem -M  human_g1k_v37.fasta  Sample1_L001_R1_001.fastq.gz Sample2_S5_L001_R2_001.fastq.gz |\
  java -jar dist/bam2fastq.jar  -F tmpR1.fastq.gz -R tmpR2.fastq.gz

```

before:

```
$ ls -lah Sample1_L001_R1_001.fastq.gz Sample2_S5_L001_R2_001.fastq.gz
-rw-r--r-- 1 lindenb lindenb 181M Jun 14 15:20 Sample1_L001_R1_001.fastq.gz
-rw-r--r-- 1 lindenb lindenb 190M Jun 14 15:20 Sample1_L001_R2_001.fastq.gz
```

after (these are Haloplex Data, with a lot of duplicates )

```
$ ls -lah tmpR1.fastq.gz  tmpR2.fastq.gz
-rw-rw-r-- 1 lindenb lindenb  96M Nov 20 17:10 tmpR1.fastq.gz
-rw-rw-r-- 1 lindenb lindenb 106M Nov 20 17:10 tmpR2.fastq.gz
```

using BZ2:

```
$  ls -lah *.bz2
-rw-rw-r-- 1 lindenb lindenb 77M Nov 20 17:55 tmpR1.fastq.bz2
-rw-rw-r-- 1 lindenb lindenb 87M Nov 20 17:55 tmpR2.fastq.bz2
```

check the number of reads

```
$ gunzip -c Sample1_L001_R1_001.fastq.gz | wc -l
5824676
$ gunzip -c tmpR1.fastq.gz | wc -l
5824676
```

verify one read

```
$ gunzip -c Sample1_L001_R1_001.fastq.gz | cat -n | head -n 4
     1	@M00491:25:000000000-A46H3:1:1101:11697:2045 1:N:0:5
     2	AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACATTGGCAAATAGCATGCCGAGGTACGCTTAAAAAAAAAACGACGCGAGGCAGGGGGGGAGGAAGCAGGGGAGCAACAGGGGGAAGGGAAGGGAAGAGAAGAAGAACGAACGAAAG
     3	+
     4	AAAAAAAA1AC1FFGCGA0AFFBGAGHHFF2GBGHH0B2DBCF101111D211B////A11///B/1DE1E/>>E//?///</<><C////<?9-9-99A-;/---;---;-9--9=---------9:AF---9//:/9/:9---9-:-9-


$ gunzip -c tmpR1.fastq.gz | cat -n | grep  -A 3 -w "@M00491:25:000000000-A46H3:1:1101:11697:2045"
5771577	@M00491:25:000000000-A46H3:1:1101:11697:2045 1:N:0:1
5771578	AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACATTGGCAAATAGCATGCCGAGGTACGCTTAAAAAAAAAACGACGCGAGGCAGGGGGGGAGGAAGCAGGGGAGCAACAGGGGGAAGGGAAGGGAAGAGAAGAAGAACGAACGAAAG
5771579	+
5771580	AAAAAAAA1AC1FFGCGA0AFFBGAGHHFF2GBGHH0B2DBCF101111D211B////A11///B/1DE1E/>>E//?///</<><C////<?9-9-99A-;/---;---;-9--9=---------9:AF---9//:/9/:9---9-:-9-
```

## Example 2 from BAM

```
$ java -jar dist/bam2fastq.jar \
    -F tmpR1.fastq.gz -R tmpR2.fastq.gz file.bam
(...)
-rw-r--r-- 1 lindenb lindenb 565M Nov 18 10:44 Sample_S1_L001_R1_001.fastq.gz
-rw-r--r-- 1 lindenb lindenb 649M Nov 18 10:45 Sample_S1_L001_R2_001.fastq.gz
-rw-rw-r-- 1 lindenb lindenb 470M Nov 20 16:17 tmpR1.fastq.gz.fastq.gz
-rw-rw-r-- 1 lindenb lindenb 554M Nov 20 16:17 tmpR2.fastq.gz.fastq.gz
```

## Cited In:

  * "Plastomes of nine hornbeams and phylogenetic implications", Ying Li & al;  Ecology and Evolution, 2018; DOI: 10.1002/ece3.4414; https://onlinelibrary.wiley.com/doi/pdf/10.1002/ece3.4414 

END_DOC
*/
@Program(name="bam2fastq",
	description="Same as picard/SamToFastq but allow missing reads + shuffle reads using hash(name) so you can use them with bwa. ",
	deprecatedMsg="use picard",
	keywords={"fastq"}
	)
public class BamToFastq
	extends Launcher
	{

	private static final Logger LOG = Logger.build(BamToFastq.class).make();




	@Parameter(names={"-F","--forward"},description="Save fastq_R1 to file (default: stdout)")
	private File forwardFile = null;

	@Parameter(names={"-R","--reverse"},description="Save fastq_R2 to file (default: interlaced with forward)")
	private File reverseFile = null;

	@Parameter(names={"-r","--repair"},description="repair: insert missing read")
	private boolean repair_missing_read = false;
	
	@Parameter(names={"-T","--tmpDir"},description="tmp directory")
	private File tmpDir = IOUtils.getDefaultTmpDir();

	@Parameter(names={"-maxRecordsInRam","--maxRecordsInRam"},description="Max records in RAM")
	private int maxRecordsInRam =50000;

	
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
		public MappedFastq decode(final DataInputStream dis) throws IOException
			{
			final MappedFastq m=new MappedFastq();
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
	public int doWork(List<String> args) {

		SamReader sfr=null;
		SortingCollection<MappedFastq> fastqCollection=null;
		try
			{
			boolean found_single=false;
			boolean found_paired=false;
			long non_primary_alignmaned_flag=0L;
			
			sfr = super.openSamReader(oneFileOrNull(args));
			
			fastqCollection = SortingCollection.newInstance(
					MappedFastq.class,
					new MappedFastqCodec(),
					new MappedFastqComparator(),
					this.maxRecordsInRam,
					this.tmpDir.toPath()
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
					fqw1=new BasicFastqWriter(new PrintStream(stdout()));
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
								if(this.repair_missing_read)
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
			LOG.error(err);
			return -1;
			}
		finally
			{
			if(fastqCollection!=null) fastqCollection.cleanup();
			}
		}
	public static void main(final String[] args) {
		new BamToFastq().instanceMainWithExit(args);
		}
	}
