/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.PrintStream;
import java.nio.file.Path;
import java.util.List;

import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.ValidationStringency;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.log.Logger;
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
	private static final Logger LOG = Logger.of(FastqToFasta.class);

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
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
			if(!trim_after_space || (wsp=fastq.getReadName().indexOf(' '))==-1)
				{
				out.println(fastq.getReadName());
				}
			else
				{
				out.println(fastq.getReadName().substring(0, wsp));
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
		}
	
	@Override
	public int doWork(final List<String> args) {
		try
			{
			try(PrintStream out= super.openPathOrStdoutAsPrintStream(outputFile)) {
				
				if(args.isEmpty())
					{
					try(FastqReader fqR=new FourLinesFastqReader(stdin())){
						run(fqR,out);
						}
					}
				else for(String fname:args)
					{
					final File f=new File(fname);
					try(FastqReader fqR=new FourLinesFastqReader(f)) {
						run(fqR,out);
						}
					}
				out.flush();
				}
			}
		catch(Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		return 0;
		}
	
	public static void main(String[] args) {
		new FastqToFasta().instanceMainWithExit(args);

	}

}
