/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum PhD.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.PrintStream;
import java.util.List;

import htsjdk.samtools.fastq.FastqConstants;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.FastqReader;
import com.github.lindenb.jvarkit.util.picard.FourLinesFastqReader;
/*

BEGIN_DOC
 
## Example


without revcomp
```bash
$ gunzip -c mate.fastq.gz | head -n 8

@M00785:3:000000000-A60G6:1:1101:15339:1356 1:N:0:1
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCTTTTTTCTTTCTTTTCTTTTTTTTTTTTTTCTTTTTTTTTTTTTTCTTTTTTTTTTTTCTTTTTTTTTTTTTTTTTTTCTTTTTTTTCTTTCTTTCTTTTTTTTTTCTCTCTCTTTTTTTTTCTTTCTTTTTTTCCTTTTTCCTTTTTCTTTCTTCTTTTTTCTTCTTTTTTTTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTATTTT
+
111>>>11>000//A//A>>>/>/>><//<<<//<-<@--/00000-/00:0::000000:9---------//9///--9------//9/99-9----9//////-----9-----9--///9/99-;/99//////99/;/-99--/////////9;9---9//9///;;//-9////////////9///://99//9///-//9/999//99-/9/;:/9--;---9-9---@--------;---//9/
@M00785:3:000000000-A60G6:1:1101:15206:1568 1:N:0:1
GCTGGTGTTCCTCAGCCACGGGGGTAGGGAACAGGCGTTACCACTTACATTCCCAGGACACCATGGCCTGTCTCTTATACACATCTAGATGTGTATAAGAGACAGGGTAACTACAATGGACCCCTTGCAGCCTGGAAGGGCCAGCAGTTCACTTTTCCAAGAGCAGCCGTGCATTCTGCACCTGAGTGTTGGCCTCTCCTGGCCATAGAGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATC
+
11>11BFAA@@DFB1BB1FECE?CE//AC0CFHFGFCGGGGF1BFHGHFHHHFHFHGEFBGFHFBCFGHGHFGEFGBGFGFHFHHHFE2F2F1FFHFFFHGEGB0AEHHEHHHFFFGBFFGHGGGHHFFHFFG0/0CHGGGGGHGHHHGHGHHHHHHHHHHGHGHGHHGGCGHHHHHHHHHHHHHHGHHHAEHGGGGGGGGFGGGGFGFFGGGGBBEGEGGGFFF?FEFFFFFFFFFFFFFFFF?@FFFFB
```
with revcomp
```bash
 gunzip -c mate.fastq.gz | head -n 8 | java -jar dist/fastqrevcomp.jar

@M00785:3:000000000-A60G6:1:1101:15339:1356 1:N:0:1
AAAATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAAAAAAAAGAAGAAAAAAGAAGAAAGAAAAAGGAAAAAGGAAAAAAAGAAAGAAAAAAAAAGAGAGAGAAAAAAAAAAGAAAGAAAGAAAAAAAAGAAAAAAAAAAAAAAAAAAAGAAAAAAAAAAAAGAAAAAAAAAAAAAAGAAAAAAAAAAAAAAGAAAAGAAAGAAAAAAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
+
/9//---;--------@---9-9---;--9/:;/9/-99//999/9//-///9//99//:///9////////////9-//;;///9//9---9;9/////////--99-/;/99//////99/;-99/9///--9-----9-----//////9----9-99/9//------9--///9//---------9:000000::0:00/-00000/--@<-<//<<<//<>>/>/>>>A//A//000>11>>>111
@M00785:3:000000000-A60G6:1:1101:15206:1568 1:N:0:1
GATACATCGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCTATGGCCAGGAGAGGCCAACACTCAGGTGCAGAATGCACGGCTGCTCTTGGAAAAGTGAACTGCTGGCCCTTCCAGGCTGCAAGGGGTCCATTGTAGTTACCCTGTCTCTTATACACATCTAGATGTGTATAAGAGACAGGCCATGGTGTCCTGGGAATGTAAGTGGTAACGCCTGTTCCCTACCCCCGTGGCTGAGGAACACCAGC
+
BFFFF@?FFFFFFFFFFFFFFFFEF?FFFGGGEGEBBGGGGFFGFGGGGFGGGGGGGGHEAHHHGHHHHHHHHHHHHHHGCGGHHGHGHGHHHHHHHHHHGHGHHHGHGGGGGHC0/0GFFHFFHHGGGHGFFBGFFFHHHEHHEA0BGEGHFFFHFF1F2F2EFHHHFHFGFGBGFEGFHGHGFCBFHFGBFEGHFHFHHHFHGHFB1FGGGGCFGFHFC0CA//EC?ECEF1BB1BFD@@AAFB11>11
```

END_DOC
 */
@Program(name="fastqrevcomp",
	description="produces a reverse-complement fastq (for mate pair alignment see http://seqanswers.com/forums/showthread.php?t=5085 )",
	keywords={"fastq"}
	)
public class FastqRevComp extends Launcher
	{
	private static final Logger LOG = Logger.build(FastqRevComp.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	
	@Parameter(names="-1",description=" interleaced input : only reverse complement R1")
	private boolean only_R1=false;
	@Parameter(names="-2",description=" interleaced input : only reverse complement R2")

	private boolean only_R2=false;

	private FastqRevComp()
		{
		}
	

	
	
	
	private void run(FastqReader r,PrintStream out)
		{
		String s;
		long nRec=0L;
		r.setValidationStringency(ValidationStringency.LENIENT);
		while(r.hasNext())
			{
			if(++nRec%1E6==0)
				{
				LOG.info("N-Reads:"+nRec);
				}
			FastqRecord fastq=r.next();
			
			
			out.print(FastqConstants.SEQUENCE_HEADER);
			out.println(fastq.getReadName());
			s=fastq.getReadString();

			if((this.only_R2 && nRec%2==1) || (this.only_R1 && nRec%2==0) ) //interleaced
				{
				out.print(s);
				}
			else
				{
				for(int i=s.length()-1;i>=0;i--)
					{
					out.print(AcidNucleics.complement(s.charAt(i)));
					}
				}
			out.println();
			
			out.print(FastqConstants.QUALITY_HEADER);
			s=fastq.getBaseQualityHeader();
			if(s!=null) out.print(s);
			out.println();
			
			s=fastq.getBaseQualityString();
			if((this.only_R2 && nRec%2==1) || (this.only_R1 && nRec%2==0) ) //interleaced
				{
				out.print(s);
				}
			else
				{
				for(int i=s.length()-1;i>=0;i--)
					{
					out.print(s.charAt(i));
					}
				}
			out.println();
			if(out.checkError()) break;
			}
		out.flush();
		LOG.info("Done. N-Reads:"+nRec);
		}
	
	
	@Override
	public int doWork(List<String> args) {
		
		if(only_R1 && only_R2)
			{
			LOG.error("Both options -1 && -2 used.");
			return -1;
			}
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
			else for(final String fn : args)
				{
				File f=new File(fn);
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
		new FastqRevComp().instanceMainWithExit(args);

	}

}
