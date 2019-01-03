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


History:
* 2014 creation

*/

package com.github.lindenb.jvarkit.tools.bamindexnames;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
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
import java.util.List;
import java.util.Set;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

/**
 BEGIN_DOC
 
 
## Example



```bash
$cat read.names
HWI-1KL149:18:C0RNBACXX:1:1101:12800:7069
HWI-1KL149:18:C0RNBACXX:1:1101:15533:71685
HWI-1KL149:18:C0RNBACXX:1:1103:6001:91243
HWI-1KL149:18:C0RNBACXX:1:1107:2088:3461
HWI-1KL149:18:C0RNBACXX:1:1108:2098:26795
HWI-1KL149:18:C0RNBACXX:1:1110:10318:73043
HWI-1KL149:18:C0RNBACXX:1:1112:18688:6422
HWI-1KL149:18:C0RNBACXX:1:1116:8824:38450/1
HWI-1KL149:18:C0RNBACXX:1:1202:13982:33444/2
ZZZZ:X



$ java -jar dist/bamqueryreadnames.jar -b -s -N list.notfound input.bam read.names | samtools view |head

HWI-1KL149:18:C0RNBACXX:1:1101:12800:7069	77	*	0	0	*	*	0	0	NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN	##################################################	RG:Z:p1294	AS:i:0	XS:i:0
HWI-1KL149:18:C0RNBACXX:1:1101:12800:7069	141	*	0	0	*	*	0	0	NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN	##################################################	RG:Z:p1294	AS:i:0	XS:i:0
HWI-1KL149:18:C0RNBACXX:1:1101:15533:71685	99	X	238583	60	100M	=	23858972	274	(...)
(...)

$ cat list.notfound
ZZZZ:X
```


 
 END_DOC
 */
@Program(description="Query a Bam file indexed with BamIndexReadNames")
public class BamQueryReadNames extends BaseBamIndexReadNames
	{
	private static final Logger LOG=Logger.build(BamQueryReadNames.class).make();
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile=null;

	@Parameter(names={"-s"},description="user list of read names is sorted")
	private boolean query_reads_is_sorted=false;
	
	@Parameter(names={"-N"},description=" save unmatched names here")
	private File notFoundFile=null;

	@ParametersDelegate
	private WritingBamArgs writingBamArgs=new WritingBamArgs();
	
	private RandomAccessFile raf;
	private NameIndexDef indexDef;

	private BamQueryReadNames()
		{
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
	public int doWork(final List<String> args) {
		PrintWriter notFoundStream=new PrintWriter(new NullOuputStream());
		SamReader sfr=null;
		SAMFileWriter bamw=null;
		try
			{
			if(!(2==args.size() ||1==args.size()))
				{
				LOG.error("illegal.number.of.arguments");
				return -1;
				}
			
			if(this.notFoundFile==null)
				{
				notFoundStream.close();
				notFoundStream=openFileOrStdoutAsPrintWriter(notFoundFile);
				}
			
			File bamFile=new File(args.get(0));
			sfr=SamReaderFactory.makeDefault().
					validationStringency(ValidationStringency.SILENT).
					open(bamFile);
			File nameIdxFile=new File(bamFile.getParentFile(), bamFile.getName()+NAME_IDX_EXTENSION);
			this.indexDef=new NameIndexDef();
			this.raf=new RandomAccessFile(nameIdxFile, "r");
			indexDef.countReads=raf.readLong();
			indexDef.maxNameLengt=raf.readInt();
			
			
			LineIterator r=null;
			if(args.size()==2)
				{
				r=IOUtils.openURIForLineIterator(args.get(1));
				}
			else
				{
				r=IOUtils.openStdinForLineIterator();
				}
			SAMFileHeader header=sfr.getFileHeader().clone();
			
			bamw=writingBamArgs.openSAMFileWriter(this.outputFile, header, true);
			
			
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
			notFoundStream.close();notFoundStream=null;
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
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
		new BamQueryReadNames().instanceMainWithExit(args);
		}
}
