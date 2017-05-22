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


*/
package com.github.lindenb.jvarkit.tools.structvar;

import java.io.File;

import java.util.List;

import htsjdk.samtools.fastq.BasicFastqWriter;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.util.CloserUtil;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

/**

BEGIN_DOC






### Example




```

$ curl -L -s "https://raw.githubusercontent.com/samtools/samtools/develop/test/dat/mpileup.1.sam" | java -jar dist/samextractclip.jar 2> /dev/null 
@ERR013140.3521432/1:0:17:1:99
AGAGGTCCCCAACTTCTTTGCA
+
@AEDGBHIIIIIFJGIKHGHIJ
@ERR156632.12704932/2:0:17:1:163
TGGAGAAGGGGACAAGAGGTCCCCAACTTCTTTGCA
+
BFAFGFEIGFEFHHEIDKJGHHHJIIE=@KKGGKJG
@ERR156632.9601178/1:0:17:1:99
CTATGACAGGGAGGTCATGTGCAGGCTGGAGAAGGGGACAAGAGGTCCCCAACTTCTTTGCA
+
DEEEIIHHKIJILKHLHIKEKHHMKLKKJGKKKKLKLFIHEKIKL=KLJLKIILHKMH9LJJ
@ERR162872.21706338/1:0:17:1:99
CTTCTTTGCA
+
BHBFH<EIFG
@ERR243039.1049231/2:0:17:1:163
TGCAGGCTGGAGAAGGGGACAAGAGGTCCCCAACTTCTTTGCA
+
AEEFIFHIJDGIGIJHHIAGGGLGJIEJHJHHFIJGJJDFJIG
@ERR013140.20277196/2:1:17:97:163
AAACTGTCCAGCGCATACCCGCATCCCCCCAAAAGAAGCCACCGCCCCAACACACACCCCCCACCCGCATAACC
+
00($,+3(*+..,%%+6%*#%2,/001)%%$2%%/$.%$00(,%+,1'*.%7(%&$&#'$$$#%#%#($+%+)"
@ERR013140.19887184/1:0:17:99:113
GTGTGTGTCGGGGGTGTCTGGGGTCTCACCCACGACCAAC
+
%$($&$*+#%%#1'$$%2-'0&3$/$/$-73/69:7=1%2
@ERR013140.4461044/1:0:17:114:113
ACTCCCTGGGCCTGGCA
+
/=1:/=44-348<0(91
@ERR013140.3521432/2:0:17:226:147
CACCCCTAGAAGTGACGGC
+
71%??A9A792/7-2%(&:
@ERR013140.11659627/1:0:17:645:83
TTAGCAACAAAAAGGAC
+
%5?-$)89<=;9>(.14
@ERR013140.7259970/1:0:17:660:83
ACGCCTGGTACA
+
40&/81&8:/<<
@ERR013140.29762488/1:0:17:716:83
GGACTCA
+
)/4/142
@ERR013140.11567710/1:0:17:984:83
TGCTTGA
+
/36>+5/
```





### See also


 *  https://www.biostars.org/p/125874




### Cited In


 *  Perlman syndrome nuclease DIS3L2 controls cytoplasmic non-coding RNAs and provides surveillance pathway for maturing snRNAs : http://nar.oxfordjournals.org/content/44/21/10437.full







END_DOC
*/


@Program(name="samextractclip",
	description="Extract Clipped Sequences from a SAM. Ouput is a FASTQ")
public class SamExtractClip extends Launcher
	{

	private static final Logger LOG = Logger.build(SamExtractClip.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@Parameter(names={"-m","--minsize"},description="Min size of clipped read")
	private int min_clip_length = 5 ;

	@Parameter(names={"-c","--clipped"},description="Print the original Read where the clipped regions have been removed")
	private boolean print_clipped_read = false;

	@Parameter(names={"-p","--original"},description="Print Original whole Read that contained a clipped region.")
	private boolean print_original_read = false;

	
	@Override
	public int doWork(List<String> args) {
		SamReader r=null;
		BasicFastqWriter out=null;
		try
				{
				if(this.outputFile!=null)
					{
					LOG.info("writing to "+this.outputFile);
					out=new BasicFastqWriter(this.outputFile);
					}
				else
					{
					LOG.info("writing to stdout");
					out=new BasicFastqWriter(stdout());
					}
				if(args.isEmpty())
					{
					LOG.info("Reading from stdin");
					r= createSamReaderFactory().open(SamInputResource.of(stdin()));
					run(r,out);
					r.close();
					}
				else 
					{
					for(final String filename:args)
						{
						LOG.info("Reading from "+filename);
						r= createSamReaderFactory().open(SamInputResource.of(filename));
						run(r,out);
						r.close();r=null;
						}
					}
				out.flush();
				return RETURN_OK;
				}
			catch(final Exception err)
				{
				LOG.error(err);
				return -1;
				}
			finally
				{
				CloserUtil.close(r);
				CloserUtil.close(out);
				}
		}
		
		private void run(final SamReader r,final FastqWriter out)
			{
			int startend[]=new int[2];
			final SAMFileHeader header=r.getFileHeader();
			//w=swf.make(header, System.out);
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header);
			SAMRecordIterator it= r.iterator();
			while(it.hasNext())
				{
				final SAMRecord rec=progress.watch(it.next());
				
				if(rec.getReadUnmappedFlag()) continue;
				if(rec.isSecondaryOrSupplementary()) continue;
				if(rec.getDuplicateReadFlag()) continue;
				
				final Cigar cigar=rec.getCigar();
				if(cigar==null || cigar.isEmpty()) continue;
				
				String suffix="";
				if(rec.getReadPairedFlag())
					{
					suffix=(rec.getFirstOfPairFlag()?"/1":"/2");
					}
				
			
				
				
				startend[0]=0;
				startend[1]=rec.getReadLength();
				boolean found=false;
				for(int side=0;side<2;++side)
					{
					final CigarElement ce=cigar.getCigarElement(side==0?0:cigar.numCigarElements()-1);
					if(!ce.getOperator().equals(CigarOperator.S)) continue;
					if(ce.getLength() < min_clip_length) continue;
					
					found=true;
					final String clippedSeq;
					final String clippedQual;
					if(side==0)
						{
						startend[0]=ce.getLength();
						clippedSeq= rec.getReadString().substring(0, startend[0]);
						clippedQual=rec.getBaseQualityString().substring(0, startend[0]);
						}
					else
						{
						startend[1]=rec.getReadLength()-ce.getLength();
						clippedSeq= rec.getReadString().substring(startend[1]);
						clippedQual=rec.getBaseQualityString().substring(startend[1]);
						}
					
					out.write(new FastqRecord(
							rec.getReadName()+suffix+";"+side+";"+rec.getReferenceName()+";"+rec.getAlignmentStart()+";"+rec.getFlags()+";"+rec.getCigarString()+";"+(side==0?"5'":"3'"),
							clippedSeq,
							"",
							clippedQual
							));
					}
				if(!found) continue;
				
				String bases=rec.getReadString();
				String qual=rec.getBaseQualityString();
				if( rec.getReadNegativeStrandFlag())
					{
					bases=AcidNucleics.reverseComplement(bases);
					qual=new StringBuilder(qual).reverse().toString();
					}
				
				if(print_original_read)
					{
					out.write(new FastqRecord(
							rec.getReadName()+suffix,
							bases,
							"",
							qual
							));
					}
				
				if(print_clipped_read)
					{
					out.write(new FastqRecord(
							rec.getReadName()+suffix+":clipped",
							bases.substring(startend[0], startend[1]),
							"",
							qual.substring(startend[0], startend[1])
							));
					}
				}
			
			it.close();
			progress.finish();
			}

		
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new SamExtractClip().instanceMainWithExit(args);
		}
	}
