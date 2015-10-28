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
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloserUtil;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.util.Collection;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.util.command.Command;

public class BamQueryReadNames extends AbstractBamQueryReadNames
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(BamQueryReadNames.class);

    public BamQueryReadNames()
    	{
    	}
  
    @Override
	public Command createCommand()
		{
		return new MyCommand();
		}
	
	static private class MyCommand extends AbstractBamQueryReadNames.AbstractBamQueryReadNamesCommand
		{

	
	private RandomAccessFile raf;
	private NameIndexDef indexDef;

	
	
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
		protected Collection<Throwable> call(String inputName) throws Exception {
			throw new IllegalStateException("should not be called");
		}
	
	@Override
			public Collection<Throwable> call() throws Exception {
				{
					
					
				
				PrintWriter notFoundStream;
				
				if(getNotFoundStreamFile()==null)
					{
					notFoundStream =new PrintWriter(new NullOuputStream());
					}
				else
					{
					notFoundStream = new PrintWriter(getNotFoundStreamFile());
					}
				
				
				
				final List<String> args = getInputFiles();
				SamReader sfr=null;
				SAMFileWriter bamw=null;
				try
					{
					if(!(2==args.size() ||1==args.size()))
						{
						return wrapException(getMessageBundle("illegal.number.of.arguments"));
						}
					final File bamFile=new File(args.get(0));
					sfr= openSamReader(args.get(0));
					File nameIdxFile=new File(bamFile.getParentFile(), bamFile.getName()+NAME_IDX_EXTENSION);
					this.indexDef=new NameIndexDef();
					this.raf=new RandomAccessFile(nameIdxFile, "r");
					indexDef.countReads=raf.readLong();
					indexDef.maxNameLengt=raf.readInt();
					
					
					BufferedReader r=null;
					if(2==args.size())
						{
						r=IOUtils.openURIForBufferedReading(args.get(1));
						}
					else
						{
						r=IOUtils.openStreamForBufferedReader(stdin());
						}
					SAMFileHeader header=sfr.getFileHeader().clone();
					header.setSortOrder(SortOrder.unsorted);
					
					bamw= openSAMFileWriter(header, true);
					
					
					
					long iter_start = 0L;
					String line;
					while((line=r.readLine())!=null)
						{
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
					LOG.info("done");
					return RETURN_OK;
					}
				catch(Exception err)
					{
					return wrapException(err);
					}
				finally
					{
					CloserUtil.close(notFoundStream);
					CloserUtil.close(raf);
					CloserUtil.close(sfr);
					CloserUtil.close(bamw);
					}
				}
			}
		}
	
	public static void main(String[] args)
		{
		new BamQueryReadNames().instanceMain(args);
		}
}
