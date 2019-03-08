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
package com.github.lindenb.jvarkit.tools.sam4weblogo;

import java.io.File;
import java.io.PrintWriter;
import java.util.List;
import java.util.function.Function;
import java.util.function.Predicate;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.bio.IntervalParser;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.SamRecordJEXLFilter;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.StringUtil;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.fastq.FastqConstants;
import htsjdk.samtools.filter.SamRecordFilter;
/**
BEGIN_DOC

## Motivation

"Sequence logo ( http://weblogo.berkeley.edu/logo.cgi ) for different alleles or generated from SAM/BAM" http://www.biostars.org/p/73021

![ScreenShot](https://raw.github.com/lindenb/jvarkit/master/doc/sam2weblogo.png)


## Example

```bash
$ java -jar dist/sam4weblogo.jar -r seq1:80-110  sorted.bam  2> /dev/null | head -n 50
>B7_593:4:106:316:452/1
TGTTG--------------------------
>B7_593:4:106:316:452a/1
TGTTG--------------------------
>B7_593:4:106:316:452b/1
TGTTG--------------------------
>B7_589:8:113:968:19/2
TGGGG--------------------------
>B7_589:8:113:968:19a/2
TGGGG--------------------------
>B7_589:8:113:968:19b/2
TGGGG--------------------------
>EAS54_65:3:321:311:983/1
TGTGGG-------------------------
>EAS54_65:3:321:311:983a/1
TGTGGG-------------------------
>EAS54_65:3:321:311:983b/1
TGTGGG-------------------------
>B7_591:6:155:12:674/2
TGTGGGGG-----------------------
>B7_591:6:155:12:674a/2
TGTGGGGG-----------------------
>B7_591:6:155:12:674b/2
TGTGGGGG-----------------------
>EAS219_FC30151:7:51:1429:1043/2
TGTGGGGGGCGCCG-----------------
>EAS219_FC30151:7:51:1429:1043a/2
TGTGGGGGGCGCCG-----------------
>EAS219_FC30151:7:51:1429:1043b/2
TGTGGGGGGCGCCG-----------------
>B7_591:5:42:540:501/1
TGTGGGGGCCGCAGTG---------------
>EAS192_3:5:223:142:410/1
TGGGGGGGGCGCAGT----------------
>B7_591:5:42:540:501a/1
TGTGGGGGCCGCAGTG---------------
>EAS192_3:5:223:142:410a/1
TGGGGGGGGCGCAGT----------------
>B7_591:5:42:540:501b/1
TGTGGGGGCCGCAGTG---------------
>EAS192_3:5:223:142:410b/1
TGGGGGGGGCGCAGT----------------
```

### fastq-like output

```
$ java -jar dist/sam4weblogo.jar -r 'RF01:100-130' src/test/resources/S1.bam -q -c

@RF01_44_622_1:0:0_1:0:0_3a/1
TATTCTTCCAATAG-----------------
+
22222222222222                 
@RF01_44_499_0:0:0_3:0:0_7b/2
TATTCTTCCAATAG-----------------
+
22222222222222                 
@RF01_67_565_0:0:0_2:0:0_67/2
TATTCTTCCAATAGTGAATTAGAGAATAGAT
+
2222222222222222222222222222222
@RF01_94_620_1:0:0_2:0:0_15/2
TATTCTTCCAATAGTGAATTAGAGAATAGAT
+
2222222222222222222222222222222
@RF01_102_665_1:0:0_1:0:0_71/1
--TTCTTCCAATAGTGAATTAGAGAATAGAT
+
  22222222222222222222222222222
@RF01_110_504_2:0:0_1:0:0_5d/2
----------ATAGTGAATTAGATAATAGAT
+
          222222222222222222222
@RF01_121_598_1:0:0_3:0:0_6e/2
---------------------GAGAATAGAT
+
                     2222222222
```


## See also

* https://www.biostars.org/p/103052/
* http://www.sciencedirect.com/science/article/pii/S1874778715300210

END_DOC

 */
@Program(name="sam4weblogo",
	description="Sequence logo for different alleles or generated from SAM/BAM ",
	biostars= {73021,368200},
	keywords={"sam","bam","visualization","logo"}
	)
public class SAM4WebLogo extends Launcher
	{
	private static final Logger LOG = Logger.build(SAM4WebLogo.class).make();
	
	@Parameter(names={"-c","--clipped","--clip"},description="Use Clipped Bases")
	private boolean useClip = false;
	
	@Parameter(names={"-r","--region","--interval"},description="Region to observe: chrom:start-end",required=true)
	private String regionStr = null;

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@Parameter(names={"-readFilter","--readFilter"},description="[20171201](moved to jexl)"+SamRecordJEXLFilter.FILTER_DESCRIPTION)
	private SamRecordFilter SamRecordFilter = SamRecordJEXLFilter.buildAcceptAll();
	
	@Parameter(names={"-q","--fastq"},description="[20180812]print fastq-like format. Was : https://github.com/lindenb/jvarkit/issues/109")
	private boolean print_fastq = false;
	@Parameter(names={"-fqu","--fqu"},description="[20180813] fastq unknown quality character")
	private String fastq_quality_unknown_str = "!";
	@Parameter(names={"-fqp","--fqp"},description="[20180813] fastq padding quality character")
	private String fastq_quality_padding_str = "-";

	
	private final Function<SAMRecord,Integer> readStart = rec -> 
		 useClip ? rec.getUnclippedStart() : rec.getAlignmentStart() ;
		
	private final Function<SAMRecord,Integer> readEnd = rec -> 
			 useClip ? rec.getUnclippedEnd() : rec.getAlignmentEnd() ;


	private void printRead(
			final PrintWriter out,
			final SAMRecord rec,
			final Interval interval
			)
		{
		
        final Cigar cigar=rec.getCigar();
        final String recQualString = rec.getBaseQualityString();
        final Function<Integer,Character> read2qual;
        if(recQualString==null || SAMRecord.NULL_QUALS_STRING.equals(recQualString)) {
        	read2qual = IDX -> '~';
        	}
        else
        	{
        	read2qual = IDX -> {
	        	if(IDX<0 || IDX>=recQualString.length()) return '~';
	        	return recQualString.charAt(IDX);
	        	};
        	}
        final byte rec_bases[] = rec.getReadBases();
        final Function<Integer,Character> read2base;
        if(SAMRecord.NULL_SEQUENCE.equals(rec_bases)) {
        	read2base = IDX-> '?'; 
        	}
        else
        	{
        	read2base = 
		        IDX -> {
			        	if(IDX<0 || IDX>=rec_bases.length) return '?';
			        	return (char)rec_bases[IDX];
			        	};
        	}    
        	
    	final Predicate<Integer> inInterval = IDX -> IDX>=interval.getStart() && IDX<=interval.getEnd();
        final StringBuilder seq  = new StringBuilder(interval.length());
        final StringBuilder qual = new StringBuilder(interval.length());
         
         int refPos = Math.min(
        		 interval.getStart(),
        		 rec.getUnclippedStart()
        		 );
         
    	 while(refPos < rec.getUnclippedStart())
    	        	{
    		 		if(inInterval.test(refPos))
    		 			{
    		 			seq.append('-');
    		 			qual.append(this.fastq_quality_padding_str);
    		 			}
    	        	++refPos;
    	        	}
        
        int readPos=0;
        for(int i=0;i< cigar.numCigarElements();++i)
        	{
        	final CigarElement ce=cigar.getCigarElement(i);
        	final CigarOperator op=ce.getOperator();
        	switch(op)
        		{
        		case P:break;
        		case I:
        			{
    				readPos+=ce.getLength();
        			break;
        			}
        		case D: case N:
        			{
        			for(int j=0;j< ce.getLength() && refPos<= interval.getEnd();++j)
        				{
        				if(inInterval.test(refPos))
        					{
        					seq.append('-');
        					qual.append(this.fastq_quality_padding_str);
        					}
        				refPos++;
        				}
        			break;
        			}
        		case H:
        			{
        			for(int j=0;j< ce.getLength() && refPos<= interval.getEnd();++j)
        				{
        				if(inInterval.test(refPos) )
        					{
        					seq.append(this.useClip?'n':'-');
        					qual.append(this.useClip?this.fastq_quality_unknown_str:this.fastq_quality_padding_str);
        					}
        				refPos++;
        				}
        			break;
        			}
        		case S:
        			{
        			for(int j=0;j< ce.getLength() && refPos<= interval.getEnd();++j)
        				{
        				if(inInterval.test(refPos) )
        					{
        					if(this.useClip)
        						{
        						seq.append(Character.toLowerCase(read2base.apply(readPos)));
        						qual.append(read2qual.apply(readPos));
        						}
        					else
        						{
        						seq.append('-');
        						qual.append(this.fastq_quality_padding_str);
        						}
        					}
        				readPos++;
        				refPos++;
        				}
        			break;
        			}
        		case M:case X: case EQ:
        			{
        			for(int j=0;j< ce.getLength() && refPos<= interval.getEnd();++j)
        				{
        				if(inInterval.test(refPos))
        					{
        					seq.append(read2base.apply(readPos));
        					qual.append(read2qual.apply(readPos));
        					}
        				readPos++;
        				refPos++;
        				}
        			break;
        			}
        		default:throw new IllegalStateException("Not handled. op:"+op);
        		}	
        	}
       
        while(refPos<= interval.getEnd())
        	{
        	seq.append('-');
        	qual.append(this.fastq_quality_padding_str);
        	refPos++;
        	}
    	out.print(this.print_fastq?FastqConstants.SEQUENCE_HEADER:">");
    	out.print(rec.getReadName());
    	if(rec.getReadPairedFlag())
        	{
        	if(rec.getFirstOfPairFlag()) out.print(FastqConstants.FIRST_OF_PAIR);
        	if(rec.getSecondOfPairFlag()) out.print(FastqConstants.SECOND_OF_PAIR);
        	}
    	out.println();
    	out.println(seq);
	    if(this.print_fastq)
	    	{
	    	out.println(FastqConstants.QUALITY_HEADER);
	    	out.println(qual);
	    	}
		}
	
	@Override
	public int doWork(final List<String> args) {
		
		if(this.fastq_quality_padding_str.length()!=1) {
			LOG.error("Bad fastq padding character (length!=1)");
			return -1;
		}
		
		if(this.fastq_quality_unknown_str.length()!=1) {
			LOG.error("Bad fastq unknown character (length!=1)");
			return -1;
		}
		
    	if(StringUtil.isBlank(this.regionStr))
			{
			LOG.error("Undefined interval.");
			return -1;
			}
    	final IntervalParser intervalParser =new IntervalParser().setFixContigName(true);
		PrintWriter out=null;
		SamReader samReader=null;
		SAMRecordIterator iter=null;
		try {
			out = super.openFileOrStdoutAsPrintWriter(outputFile);
			for(final String inputName: IOUtils.unrollFiles(args)) {
				samReader = openSamReader(inputName);
				
				final Interval interval = intervalParser.
						setDictionary(samReader.getFileHeader().getSequenceDictionary()).
						parse(this.regionStr);
				if(interval==null)
					{
					LOG.error("Bad interval "+this.regionStr);
					return -1;
					}
				
				iter=samReader.query(
						interval.getContig(),
						interval.getStart(),
						interval.getEnd(),
						false);
				while(iter.hasNext())
	                {
	                final SAMRecord rec=iter.next();
	                if(rec.getReadUnmappedFlag()) continue;
	                if(this.SamRecordFilter.filterOut(rec)) continue;
	                final Cigar cigar =rec.getCigar();
	                if(cigar==null || cigar.isEmpty()) continue;
	               
	                if(!rec.getReferenceName().equals(interval.getContig())) continue;
	               
                    if(this.readEnd.apply(rec) < interval.getStart() ) continue;
                    if(this.readStart.apply(rec) > interval.getEnd() ) continue;
	                printRead(out,rec,interval);
	                }
				iter.close();
				iter=null;
				samReader.close();
				samReader=null;
				}
			out.flush();
	        out.close();out=null;
	        return 0;
			} 
		catch (final Exception err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(samReader);
			CloserUtil.close(out);
			}
		}

public static void main(final String[] args)
	{
	new SAM4WebLogo().instanceMainWithExit(args);
	}
}
