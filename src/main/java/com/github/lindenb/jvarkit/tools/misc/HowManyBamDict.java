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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

/**
BEGIN_DOC

## Example

```
$ find PROJ/ -name "accepted_hits.bam" |\
 java -jar dist/howmanybamdict.jar

DICT	208e057b92b5c45262c59b40209b5fde	22	2725537669	chr1=195471971;chr10=130694993;chr11=122082543;chr12=120129022;chr13=120421639;chr14=124902244;chr15=104043685;chr16=98207768;chr17=94987271;chr18=90702639;chr19=61431566;chr2=182113224;chr3=160039680;chr4=156508116;chr5=151834684;chr6=149736546;chr7=145441459;chr8=129401213;chr9=124595110;chrM=16299;chrX=171031299;chrY=91744698	PROJ/S1/accepted_hits.bams.bam
BAM	PROJ/S1/accepted_hits.bam	208e057b92b5c45262c59b40209b5fde
DICT	86b39aef44f740da797b89baf3f505d8	66	2730871774	chr1=195471971;chr10=130694993;chr11=122082543;chr12=120129022;chr13=120421639;chr14=124902244;chr15=104043685;chr16=98207768;chr17=94987271;chr18=90702639;chr19=61431566;chr1_GL456210_random=169725;chr1_GL456211_random=241735;chr1_GL456212_random=153618;chr1_GL456213_random=39340;chr1_GL456221_random=206961;chr2=182113224;chr3=160039680;chr4=156508116;chr4_GL456216_random=66673;chr4_GL456350_random=227966;chr4_JH584292_random=14945;chr4_JH584293_random=207968;chr4_JH584294_random=191905;chr4_JH584295_random=1976;chr5=151834684;chr5_GL456354_random=195993;chr5_JH584296_random=199368;chr5_JH584297_random=205776;chr5_JH584298_random=184189;chr5_JH584299_random=953012;chr6=149736546;chr7=145441459;chr7_GL456219_random=175968;chr8=129401213;chr9=124595110;chrM=16299;chrUn_GL456239=40056;chrUn_GL456359=22974;chrUn_GL456360=31704;chrUn_GL456366=47073;chrUn_GL456367=42057;chrUn_GL456368=20208;chrUn_GL456370=26764;chrUn_GL456372=28664;chrUn_GL456378=31602;chrUn_GL456379=72385;chrUn_GL456381=25871;chrUn_GL456382=23158;chrUn_GL456383=38659;chrUn_GL456385=35240;chrUn_GL456387=24685;chrUn_GL456389=28772;chrUn_GL456390=24668;chrUn_GL456392=23629;chrUn_GL456393=55711;chrUn_GL456394=24323;chrUn_GL456396=21240;chrUn_JH584304=114452;chrX=171031299;chrX_GL456233_random=336933;chrY=91744698;chrY_JH584300_random=182347;chrY_JH584301_random=259875;chrY_JH584302_random=155838;chrY_JH584303_random=158099	PROJ/S2/accepted_hits.bam	ROJ/S2/accepted_hits.bam
BAM	PROJ/S2/accepted_hits.bam	86b39aef44f740da797b89baf3f505d8
BAM	PROJ/S3/accepted_hits.bam	86b39aef44f740da797b89baf3f505d8
```

END_DOC
 */


@Program(name="howmanybamdict",description="finds if there's are some differences in the sequence dictionaries.",
		keywords={"sam","bam","dict"}
		)
public class HowManyBamDict extends Launcher {
	private static final Logger LOG = Logger.build(HowManyBamDict.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	public HowManyBamDict()
		{
		 
		}
	
	
	private class Dict
		{
		SAMSequenceDictionary ssd;
		String hash;
		File representative;
		Dict(SAMSequenceDictionary ssd,File representative)
			{
			this.ssd=ssd;
			this.representative=representative;
	    	this.hash =ssd.md5();
			}
		
		@Override
		public boolean equals(Object obj)
			{
			if(this==obj) return true;
			if(obj==null) return false;
			return this.ssd.equals(Dict.class.cast(obj).ssd);
			}
		
		@Override
		public int hashCode() {
			return ssd.hashCode();
			}
		void print()
			{
			System.out.print("DICT");
			System.out.print("\t");
			System.out.print(this.hash);
			System.out.print("\t");
			System.out.print(ssd.size());
			System.out.print("\t");
			System.out.print(ssd.getReferenceLength());
			System.out.print("\t");
			boolean first=true;
			for(SAMSequenceRecord ssr:ssd.getSequences())
				{
				if(!first) System.out.print(";");
				first=false;
				System.out.print(ssr.getSequenceName());
				System.out.print('=');
				System.out.print(ssr.getSequenceLength());
				}
			System.out.print("\t");
			System.out.print(this.representative);
			System.out.println();
			}
		}

	private Dict empty=null;
	private Set<Dict>  allditcs=new LinkedHashSet<Dict>();
	
	
 	
 	private void handle(PrintWriter out,File f) throws IOException
 		{
 		SamReader sfr=null;
 		try {
 			LOG.info(f.getPath());
			sfr= super.openSamReader(f.getPath());
			SAMFileHeader header=sfr.getFileHeader();
			if(header==null || header.getSequenceDictionary()==null)
				{
				if(this.empty==null)
					{
					this.empty=new Dict(new SAMSequenceDictionary(),f);
					allditcs.add(this.empty);
					this.empty.print();
					}
				out.print("BAM\t");
				out.print(f.getPath());
				out.print("\t");
				out.print(this.empty.hash);
				out.println();
				}
			else
				{
				Dict d=new Dict(header.getSequenceDictionary(), f);
				if(this.allditcs.add(d))
					{
					d.print();
					}
				out.print("BAM\t");
				out.print(f.getPath());
				out.print("\t");
				out.print(d.hash);
				out.println();
				}
 			} 
 		catch (Exception e)
			{
			LOG.error(e.getMessage(),e);
			throw new IOException(e);
			}
 		finally
 			{
 			CloserUtil.close(sfr);
 			}
 		}
	
 	@Override
 	public int doWork(List<String> args) {
		PrintWriter out=null;
		try
			{
			out = super.openFileOrStdoutAsPrintWriter(this.outputFile);
			if(args.isEmpty())
				{
				LOG.info("Reading from stdin");
				String line;
				
					BufferedReader in=new BufferedReader(new InputStreamReader(System.in));
					while((line=in.readLine())!=null)
						{
						if(line.isEmpty() || line.endsWith(File.separator) || line.startsWith("#")) continue;
						handle(out,new File(line));
						}
					in.close();
					
				}
			else
				{
				for(String filename:args)
					{
					handle(out,new File(filename));
					}
				}
			out.flush();
			return RETURN_OK;
			}
		catch(IOException err)
			{
			LOG.error(err);
			return -1;
			}
		finally {
			CloserUtil.close(out);
			}
		}

 	public static void main(String[] args)
		{
		new HowManyBamDict().instanceMainWithExit(args);
		}

}
