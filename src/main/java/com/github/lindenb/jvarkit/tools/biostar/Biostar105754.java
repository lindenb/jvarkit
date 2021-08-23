/*
The MIT License (MIT)
Copyright (c) 2021 Pierre Lindenbaum

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

package com.github.lindenb.jvarkit.tools.biostar;

import htsjdk.samtools.util.CloserUtil;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.List;

import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

/**
BEGIN_DOC

## Example

```
$  echo -e "1\t1000\t20000\n3\t100\t200\nUn\t10\t11"  |\
  java -jar dist/biostar105754.jar -B path/to//All_hg19_RS_noprefix.b


#no data found for	Un	10	11
1	1000	1001	0.0	1	1000	20000
3	100	101	0.0	3	100	200
```

END_DOC
 */
@Program(name="biostar105754",
	description="bigwig : peak distance from specific genomic region",
	biostars=105754,
	keywords={"wig","bigwig"},
	modificationDate="20210818",
	creationDate="20140708"
	)
public class Biostar105754 extends Launcher {

	private static final Logger LOG = Logger.build(Biostar105754.class).make();
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	
	@Parameter(names={"-B","--bigwig"},description="Big Wig file",required=true)
	private String bigWigFile = null;

	
	private static int distance(int start1,int end1, int start2,int end2)
		{
		if(end1< start2)
			{
			return start2-end1;
			}
		else if(end2< start1)
			{
			return start1-end2;
			}
		int d=Math.abs(start1-start2);
		d=Math.min(d,Math.abs(start1-end2));
		d=Math.min(d,Math.abs(end1-start2));
		d=Math.min(d,Math.abs(end1-end2));
		return d;
		}
	
	
		private org.broad.igv.bbfile.BBFileReader bbFileReader=null;
		private final long EXTEND_SHIFT=1000000;//
		private final long MAX_CHROM_END=Integer.MAX_VALUE-EXTEND_SHIFT;
	
		
		private void run(final BedLineReader r,final PrintWriter out)
			throws IOException
			{
			while(r.hasNext() && !out.checkError())
				{
				final BedLine bedLine = r.next();
				if(bedLine==null) {
					continue;
					}
				final String chrom = bedLine.getContig() ;
				final int chromStart0 = bedLine.getStart()-1;
				final int chromEnd0 = bedLine.getEnd();
				if(chrom.isEmpty() || chromStart0<0L || chromEnd0<chromStart0)
					{
					LOG.warn("Bad BED line: "+bedLine.join());
					continue;
					}
				
				//extends bed area until something was found
				int chromStart=chromStart0;
				int chromEnd=chromEnd0;
				
				for(;;)
					{
					BigWigIterator iter=this.bbFileReader.getBigWigIterator(
							chrom,
							chromStart,
							chrom,
							chromEnd,
							false);
					if(iter!=null)
						{
						WigItem best=null;
						while(iter.hasNext())
							{
							WigItem wigItem=iter.next();
							if(best==null || 
									distance(chromStart,chromEnd,best.getStartBase(),best.getEndBase()) >
									distance(chromStart,chromEnd,wigItem.getStartBase(),wigItem.getEndBase())
									)
								{
								best=wigItem;
								}
							}
						if(best!=null)
							{
							out.print(best.getChromosome());
							out.print("\t");
							out.print(best.getStartBase());
							out.print("\t");
							out.print(best.getEndBase());
							out.print("\t");
							out.print(best.getWigValue());
							out.print("\t");
							out.print(bedLine.join("\t"));
							out.println();
							break;
							}
						}
					//extend bed area
					long start2=chromStart-EXTEND_SHIFT;
					long end2=chromEnd+EXTEND_SHIFT;
					if(start2<0) start2=0;
					if(end2>MAX_CHROM_END) end2=MAX_CHROM_END;
					//too wide, break loop
					if(start2==0 && end2==MAX_CHROM_END)
						{
						LOG.warn("no data found for\t"+ bedLine);
						break;
						}
					chromStart=(int)start2;
					chromEnd=(int)end2;
					}
				}
			}
		
		@Override
	public int doWork(final List<String> args) {
			if(this.bigWigFile==null)
				{
				LOG.error("Big wig file undefined option");
				return -1;
				}
			
			try
				{
				LOG.info("Opening "+this.bigWigFile);
				this.bbFileReader=new BBFileReader(this.bigWigFile);
				if(!this.bbFileReader.isBigWigFile())
					{
					LOG.error("File "+this.bigWigFile+" is not a bigwig file");
					return -1;
					}
				try(PrintWriter out = super.openPathOrStdoutAsPrintWriter(outputFile)) {
					if(args.isEmpty())
						{
						try(BedLineReader r= new BedLineReader(IOUtils.openStdinForBufferedReader(),"stdin")) {
							run(r,out);
							}
						}
					else
						{
						for(final String filename : args)
							{
							try(BedLineReader r= new BedLineReader(IOUtils.openURIForBufferedReading(filename),filename)) {
								run(r,out);
								}
							}
						}
					out.flush();
					}
				return 0;
				}
			catch(final Exception err)
				{
				LOG.error(err);
				return -1;
				}
			finally
				{
				CloserUtil.close(bbFileReader);
				bbFileReader=null;
				}
			}
			
	public static void main(String[] args) {
		new Biostar105754().instanceMainWithExit(args);
		}
	}
