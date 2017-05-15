package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.PrintStream;
import java.util.List;

import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.FastqReader;
import com.github.lindenb.jvarkit.util.picard.FourLinesFastqReader;
/**
BEGIN_DOC


## Example

```bash
$ java -jar dist/fastq2fasta.jar  -N 60 -b file.fastq.gz

>HWI-1KL149:61:D2C11TCXX:2:1213:4591:29626
GAGTTGCTTTGTTTGAATATAGGTTGACTATACGAAGTGTGCGAGGACCTGCACCACGCA
GTAGGCCAAGATCAACTGAAACAGTGCTATCTGCACGACAA
>HWI-1KL149:61:D2C11TCXX:2:1213:4525:29650
CCTAGTAGTTCGTGGCCCCGGGCCCCTACTTAAACTCCTAGAACCACTCCTAGAAAGGGG
TGTTGCAGTTCGGCTCAGTCCCCGTGGTCGACTACTGTTTC
>HWI-1KL149:61:D2C11TCXX:2:1213:4569:29706
GCGCAGAGTTGTTTTAGCTATGCTGTGTTTGCATGGTTAGGTGGTGTACCTAGTGGTTTT
CTGAGACTTCTCTGAGGTTCTTGAGTAGATTAATACATCCC
>HWI-1KL149:61:D2C11TCXX:2:1213:4594:29713

```


END_DOC

 */
@Program(name="fastq2fasta",description="fastq -> fasta",
	deprecatedMsg="use awk, samtools...",
	keywords={"fastq","fasta"})
public class FastqToFasta
	extends Launcher
	{
	private static final Logger LOG = Logger.build(FastqToFasta.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@Parameter(names="-N",description="fasta line length")
	private int fastaLineLen=50;
	@Parameter(names="-b",description="trim fasta header after space")
	private boolean trim_after_space=false;
	
	private FastqToFasta()
		{
		}
	
	
	
	private void run(FastqReader r,PrintStream out)
		{
		int wsp=0;
		long nRec=0L;
		r.setValidationStringency(ValidationStringency.LENIENT);
		while(r.hasNext())
			{
			if(++nRec%1E6==0)
				{
				LOG.info("N-Reads:"+nRec);
				}
			FastqRecord fastq=r.next();
			out.print(">");
			if(!trim_after_space || (wsp=fastq.getReadHeader().indexOf(' '))==-1)
				{
				out.println(fastq.getReadHeader());
				}
			else
				{
				out.println(fastq.getReadHeader().substring(0, wsp));
				}
			
			int readLen=fastq.getReadString().length();
			int i=0;
			while(i< readLen)
				{
				int end=Math.min(i+fastaLineLen,readLen);
				out.println(fastq.getReadString().substring(i, end));
				i=end;
				}
			
			if(out.checkError()) break;
			}
		out.flush();
		LOG.info("Done. N-Reads:"+nRec);
		}
	
	@Override
	public int doWork(final List<String> args) {
		PrintStream out=null;
		try
			{
			out= super.openFileOrStdoutAsPrintStream(outputFile);
			
			if(args.isEmpty())
				{
				LOG.info("Reading from stdin");
				FastqReader fqR=new FourLinesFastqReader(stdin());
				run(fqR,out);
				fqR.close();
				}
			else for(String fname:args)
				{
				File f=new File(fname);
				LOG.info("Reading from "+f);
				FastqReader fqR=new FourLinesFastqReader(f);
				run(fqR,out);
				fqR.close();
				}
			out.flush();
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(out);
			}
		return 0;
		}
	
	public static void main(String[] args) {
		new FastqToFasta().instanceMainWithExit(args);

	}

}
