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
package com.github.lindenb.jvarkit.tools.mem;


import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;

/**
BEGIN_DOC
package/bwa/bwa-0.7.4/bwa mem \
 -M  ~/tmp/DATASANGER/hg18/chr1.fa  \
 ~/tmp/DATASANGER/4368_1_1.fastq.gz  ~/tmp/DATASANGER/4368_1_2.fastq.gz 2> /dev/null | 
 java -jar dist/bwamemdigest.jar -B <(curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gap.txt.gz" | gunzip -c | cut -f2,3,4 ) -x 500  | tee /dev/tty | gzip --best > /tmp/jeter.mem.bed.gz


END_DOC
 */
@Program(name="bwamemdigest",description="")
public class BWAMemDigest extends Launcher
	{
    private static final Logger LOG = Logger.build(BWAMemDigest.class).make();


    @Parameter(names={"-B"}, description="BED of Regions to ignore.")
    public File IGNORE_BED = null;
 
    @Parameter(names={"-x"}, description=" Extend BED Regions to ignore by 'x' bases. ")
    public int IGNORE_EXTEND = 0;
    
    
    private abstract class AbstractMemOuput
    	implements Closeable
    	{
    	
    	}
    
    
    private class DefaultMemOuput extends AbstractMemOuput
		{
    	PrintWriter out=new PrintWriter(stdout());
    	
    	public void insertion(SAMRecord record,long readNum,float countS,float countM)
	    	{
	    	stdout().println(
					bed(record)+
					"\t"+readNum+
					"\tINSERT\t"+
					record.getCigarString()+"\t"+(countS/countM)
					);
	    	}
    	
    	public void xp(final SAMRecord record,long readNum,final SAMRecord xp)
    		{
    		stdout().println(
					bed(record)+
					"\t"+readNum+
					"\tXP"+
				    "\t"+record.getCigarString()+
				    "\t"+xp.getReferenceName()+
				    "\t"+xp.getAlignmentStart()+
				    "\t"+(xp.getReadNegativeStrandFlag()?"-":"+")+
				    "\t"+xp.getCigarString()
					);
    		}
    	public void orphan(SAMRecord record,long readNum)
	    	{
    		stdout().println(
					bed(record)+
					"\t"+readNum+
				    "\tORPHAN"
					);
	    	}
    	public void deletion(SAMRecord record,long readNum)
	    	{
    		stdout().println(
					bed(record)+
					"\t"+readNum+
				    "\tDEL"+
				    "\t"+mate(record)
					);
	    	}
    	public void transloc(SAMRecord record,long readNum)
	    	{
    		stdout().println(
					bed(record)+
					"\t"+readNum+
				    "\tTRANSLOC"+
				    "\t"+mate(record)
					);
	    	}
    	@Override
    	public void close() throws IOException {
    		out.flush();
    		out.close();
    		}
    	}
    
	public BWAMemDigest()
		{
		
		}

	
	private static String bed(SAMRecord record)
		{
		return record.getReferenceName()+"\t"+
			    (record.getAlignmentStart()-1)+"\t"+
			    record.getAlignmentEnd()+"\t"+
			    (record.getReadNegativeStrandFlag()?"-":"+")
			    ;
		}
	
	private static String mate(SAMRecord record)
		{
		return record.getMateReferenceName()+"\t"+
			    (record.getMateAlignmentStart()-1)+"\t"+
			    (record.getMateNegativeStrandFlag()?"-":"+")
			    ;
		}
	
	@Override
	public int doWork(List<String> args) {
		IntervalTreeMap<Boolean> ignore=null;
		DefaultMemOuput output=new DefaultMemOuput();
		
		final float limitcigar=0.15f;
		SamReader r=null;
		try {
			r = super.openSamReader(oneFileOrNull(args));
			
		
		if(IGNORE_BED!=null)
			{
			LOG.info("open "+IGNORE_BED);
			ignore=new IntervalTreeMap<>();
			BufferedReader in=IOUtils.openFileForBufferedReading(IGNORE_BED);
			String line;
			while((line=in.readLine())!=null)
				{
				if(line.isEmpty() || line.startsWith("#")) continue;
				String tokens[]=line.split("[\t]");
				if(tokens.length<3) continue;
				if(ignore.put(new Interval(tokens[0],
						Math.max(1,Integer.parseInt(tokens[1])-(1+IGNORE_EXTEND)),
						Integer.parseInt(tokens[2])+IGNORE_EXTEND),
						Boolean.TRUE
						))
					{
					LOG.warn("BED:ignoring "+line);
					}
				}	
			in.close();
			}
		
		SAMRecordIterator iter=r.iterator();	
		long readNum=0L;
		while(iter.hasNext())
			{
			SAMRecord record=iter.next();
			++readNum;
			if(!record.getReadPairedFlag()) continue;
			if(record.getProperPairFlag()) continue;
			if(record.getReadFailsVendorQualityCheckFlag()) continue;
			if(record.getDuplicateReadFlag()) continue;
			if(record.getReadUnmappedFlag()) continue;
			if(ignore!=null && ignore.containsOverlapping(new Interval(
					record.getReferenceName(),
					record.getAlignmentStart(),
					record.getAlignmentEnd())
					))
				{
				LOG.info("ignore "+record);
				continue;
				}
			
			float countM=0f;
			float countS=0f;
			for(CigarElement c:record.getCigar().getCigarElements())
				{
				switch(c.getOperator())
					{
					case M: countM+=c.getLength(); break;
					case S: countS+=c.getLength(); break;
					default: break;
					}
				}
			
			if(countM>10 && ((countS/countM)>(1f-limitcigar) && (countS/countM)< (1f+limitcigar) ))
				{
				output.insertion(record, readNum, countS, countM);
				}
			
			for(final SAMRecord xp: SAMUtils.getOtherCanonicalAlignments(record))
				{
				if(ignore!=null &&
						ignore.containsOverlapping(new Interval(
						xp.getReferenceName(),
						xp.getAlignmentStart(),
						xp.getAlignmentStart()
						)))
					{
					LOG.info("ignore "+record);
					continue;
					}
				
				
				output.xp(record, readNum, xp);
				}
			
			
			if(record.getMateUnmappedFlag())
				{
				output.orphan(record, readNum);
				continue;
				}
			
			
			if(ignore!=null && ignore.containsOverlapping(new Interval(
					record.getMateReferenceName(),
					record.getMateAlignmentStart(),
					record.getMateAlignmentStart())
					))
				{
				LOG.info("ignore "+record);
				continue;
				}
			if(record.getReferenceIndex()==record.getMateReferenceIndex())
				{
				output.deletion(record, readNum);
				}
			else
				{
				output.transloc(record, readNum);
				}
			
			
			}
		}
		
		catch (Exception e)
			{
			e.printStackTrace();
			LOG.error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(r);
			CloserUtil.close(output);
			}
		return 0;
		}
	
	public static void main(String[] args) {
		new BWAMemDigest().instanceMainWithExit(args);
	}
}
