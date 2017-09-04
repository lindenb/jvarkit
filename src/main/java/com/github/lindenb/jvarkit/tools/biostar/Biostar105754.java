package com.github.lindenb.jvarkit.tools.biostar;

import htsjdk.samtools.util.CloserUtil;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.List;
import java.util.regex.Pattern;

import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.semontology.Term;

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
	terms=Term.ID_0000015
	)
public class Biostar105754 extends Launcher
	{

	private static final Logger LOG = Logger.build(Biostar105754.class).make();
	
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	
	
	@Parameter(names={"-B","--bigwig"},description="Big Wig file")
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
	
	
			private PrintWriter out=null;
			private org.broad.igv.bbfile.BBFileReader bbFileReader=null;
			private final long EXTEND_SHIFT=1000000;//
			private final long MAX_CHROM_END=Integer.MAX_VALUE-EXTEND_SHIFT;
	
		
		private void run(BufferedReader r)
			throws IOException
			{
			String line;
			final Pattern tab=Pattern.compile("[\t]");
			while((line=r.readLine())!=null && !this.out.checkError())
				{
				if(line.startsWith("#"))
					{
					continue;
					}
				String tokens[]=tab.split(line,4);
				if(tokens.length<3)
					{
					System.err.println("Bad BED line: "+line);
					continue;
					}
				String chrom=tokens[0];
				int chromStart0=Integer.parseInt(tokens[1]);
				int chromEnd0=Integer.parseInt(tokens[2]);
				if(chrom.isEmpty() || chromStart0<0L || chromEnd0<chromStart0)
					{
					System.err.println("Bad BED line: "+line);
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
							this.out.print(best.getChromosome());
							this.out.print("\t");
							this.out.print(best.getStartBase());
							this.out.print("\t");
							this.out.print(best.getEndBase());
							this.out.print("\t");
							this.out.print(best.getWigValue());
							this.out.print("\t");
							this.out.print(line);
							this.out.println();
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
						LOG.warn("no data found for\t"+line);
						break;
						}
					chromStart=(int)start2;
					chromEnd=(int)end2;
					}
				}
			}
		
		@Override
		public int doWork(List<String> args) {
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
					this.out = super.openFileOrStdoutAsPrintWriter(outputFile);
					if(args.isEmpty())
						{
						BufferedReader r= new BufferedReader(new InputStreamReader(stdin()));
						run(r);
						CloserUtil.close(r);
						}
					else
						{
						for(final String filename : args)
							{
							LOG.info("Reading BED from "+filename);
							final BufferedReader r= IOUtils.openURIForBufferedReading(filename);
							run(r);
							CloserUtil.close(r);
							}
						}
					this.out.flush();
					this.out.close();
					return RETURN_OK;
					}
				catch(final Exception err)
					{
					LOG.error(err);
					return -1;
					}
				finally
					{
					CloserUtil.close(bbFileReader);
					CloserUtil.close(this.out);
					bbFileReader=null;
					this.out=null;
					}
				}
			
	public static void main(String[] args) {
		new Biostar105754().instanceMainWithExit(args);
		}
	}
