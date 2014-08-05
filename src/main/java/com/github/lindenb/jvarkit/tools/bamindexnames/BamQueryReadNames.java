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

package com.github.lindenb.jvarkit.tools.bamindexnames;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.tribble.readers.LineIterator;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.util.LinkedHashSet;
import java.util.Set;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;

public class BamQueryReadNames extends AbstractBamIndexReadNames
	{
	private RandomAccessFile raf;
	private NameIndexDef indexDef;

	private BamQueryReadNames()
		{
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return DEFAULT_WIKI_PREFIX+"BamQueryReadNames";
		}
	
	@Override
	public String getProgramDescription() {
		return "Query a Bam file indexed with BamIndexReadNames";
		}
	
	
	private NameAndPos getNameAndPosAt(long index)
		throws IOException
		{
		long fileoffset=FILE_PREFIX_SIZE//header
					+ index*(this.indexDef.sizeOfNameAndPos());
		
		ByteBuffer byteBuff= ByteBuffer.allocate( indexDef.sizeOfNameAndPos());
		this.raf.seek(fileoffset);
		this.raf.readFully(byteBuff.array());
		
		NameAndPos nap=new NameAndPos();
		
		StringBuilder b=new StringBuilder(indexDef.maxNameLengt);
		for(int i=0;i< indexDef.maxNameLengt && byteBuff.get(i)!=0;++i)
			{
			b.append((char)byteBuff.get(i));
			}
		nap.name=b.toString();
		byteBuff.position( indexDef.maxNameLengt);
		nap.tid=byteBuff.getInt();
		nap.pos=byteBuff.getInt();
		return nap;
		
		}
	
	private long lower_bound(long first, long last, final String readName) throws IOException
		{
			long len = last-first;
			while (len > 0)
			{
				long half = len /2;
				long middle = first + half;
				String middle_s = getNameAndPosAt(middle).name;
				
				if (middle_s.compareTo(readName)<0)
					{
					first = middle;
					++first;
					len = len - half - 1;
					}
				else
					{
					len = half;
					}
				

				}
			return first;
		}

	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -N (file) save unmatched names here. (Optional)");
		out.println(" -s user list of read names is sorted. (Optional)");
		out.println(" -b write binary bam (Optional)");
		super.printOptions(out);
		}

	
	@Override
	public int doWork(String[] args)
		{
		boolean binary_bam=false;
		boolean query_reads_is_sorted=false;
		PrintWriter notFoundStream=new PrintWriter(new NullOuputStream());
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"N:sb"))!=-1)
			{
			switch(c)
				{
				case 's': query_reads_is_sorted=true;break;
				case 'b': binary_bam=true;break;
				case 'N':
					{
					try
						{
						notFoundStream=new PrintWriter(new File(opt.getOptArg()));
						}
					catch(IOException err)
						{
						error(err);
						return -1;
						}
					break;
					}
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
		SAMFileWriter bamw=null;
		try
			{
			if(!(opt.getOptInd()+2==args.length || opt.getOptInd()+1==args.length))
				{
				error(getMessageBundle("illegal.number.of.arguments"));
				return -1;
				}
			File bamFile=new File(args[opt.getOptInd()+0]);
			sfr=SamReaderFactory.makeDefault().
					validationStringency(ValidationStringency.SILENT).
					open(bamFile);
			File nameIdxFile=new File(bamFile.getParentFile(), bamFile.getName()+NAME_IDX_EXTENSION);
			this.indexDef=new NameIndexDef();
			this.raf=new RandomAccessFile(nameIdxFile, "r");
			indexDef.countReads=raf.readLong();
			indexDef.maxNameLengt=raf.readInt();
			
			
			LineIterator r=null;
			if(opt.getOptInd()+2==args.length)
				{
				r=IOUtils.openURIForLineIterator(args[opt.getOptInd()+1]);
				}
			else
				{
				r=IOUtils.openStdinForLineIterator();
				}
			SAMFileHeader header=sfr.getFileHeader().clone();
			SAMProgramRecord spr=header.createProgramRecord();
			spr.setProgramName(getProgramName());
			spr.setCommandLine(getProgramCommandLine());
			spr.setProgramVersion(getVersion());
			
			SAMFileWriterFactory sfw=new SAMFileWriterFactory();
			sfw.setCreateIndex(false);
			header.setSortOrder(SortOrder.unsorted);
			if(binary_bam)
				{
				bamw=sfw.makeBAMWriter(header, false, System.out);
				}
			else
				{
				bamw=sfw.makeSAMWriter(header, false, System.out);
				}
			
			
			long iter_start = 0L;
			
			while(r.hasNext())
				{
				String line=r.next();
				String searchRead=null;
				int side=-1;
				if(line.isEmpty() || line.startsWith("#")) continue;
				
				/* forward or reverse is specified ? */
				if(line.endsWith("/1"))
					{
					side=1;
					searchRead=line.substring(0, line.length()-2);
					}
				else if(line.endsWith("/2"))
					{
					side=2;
					searchRead=line.substring(0, line.length()-2);
					}
				else
					{
					side=-1;
					searchRead=line;
					}	
				long index=lower_bound(
						iter_start,
						this.indexDef.countReads,
						searchRead
						);
				if(index>=this.indexDef.countReads)
					{
					notFoundStream.println(line);
					continue;
					}
				if(query_reads_is_sorted)
					{
					iter_start=index;
					}
				
				Set<SAMRecord> found=new LinkedHashSet<SAMRecord>();
				while(index <this.indexDef.countReads )
					{
					NameAndPos nap=getNameAndPosAt(index);
					if(nap.name.compareTo(searchRead)<0)
						{
						++index;
						continue;
						}
					else if(nap.name.compareTo(searchRead)>0)
						{
						break;
						}
					SAMRecordIterator iter;
					if(nap.tid<0)
						{
						iter=sfr.queryUnmapped();
						}
					else
						{
						iter=sfr.query(
							header.getSequence(nap.tid).getSequenceName(),
							nap.pos,
							0,
							true
							);
						}
					while(iter.hasNext())
						{
						SAMRecord rec=iter.next();
						if(nap.tid>=0)
							{
							if(nap.tid!=rec.getReferenceIndex())throw new IllegalStateException();

							if(rec.getAlignmentStart()< nap.pos)
								{
								continue;
								}
							if(rec.getAlignmentStart()> nap.pos)
								{
								break;
								}
							}
						if(rec.getReadName().equals(searchRead))
							{
							if(side==1 && !(rec.getReadPairedFlag() && rec.getFirstOfPairFlag()))
								{
								continue;
								}
							else if(side==2 && !(rec.getReadPairedFlag() && rec.getSecondOfPairFlag()))
								{
								continue;
								}
							found.add(rec);
							}
						
						}
					iter.close();
					
					++index;
					}
				if(found.isEmpty())
					{
					notFoundStream.println(line);
					}
				else
					{
					for(SAMRecord rec:found)
						{
						bamw.addAlignment(rec);
						}
					}
				
				}
			CloserUtil.close(r);
			
			notFoundStream.flush();
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(notFoundStream);
			CloserUtil.close(raf);
			CloserUtil.close(sfr);
			CloserUtil.close(bamw);
			}
		}
	public static void main(String[] args)
		{
		new BamQueryReadNames().instanceMain(args);
		}
}
