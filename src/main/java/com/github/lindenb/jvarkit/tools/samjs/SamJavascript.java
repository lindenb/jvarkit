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
package com.github.lindenb.jvarkit.tools.samjs;



import java.io.File;

/**
 * Author: Pierre Lindenbaum PhD. @yokofakun
 * Motivation http://www.biostars.org/p/66319/ 
 */

import java.util.List;

import javax.script.Bindings;
import javax.script.CompiledScript;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloserUtil;

/**
BEGIN_DOC

## Motivation

Filters a BAM using javascript( java nashorn engine).

For each read the script injects in the context the following values:


* **'record'** a SamRecord  [https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMRecord.html](https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMRecord.html)
* **'header'** a SAMFileHeader  [https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMFileHeader.html](https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMFileHeader.html)


the script should return a boolean : true accept the read, false: discard the read.

## Example

### Example 1


get a SAM where the  read OR the mate is unmapped

```bash
java -jar dist/samjs.jar  \
	-e "record.readUnmappedFlag || record.mateUnmappedFlag;" \
	ex1.bam

@HD	VN:1.4	SO:unsorted
@SQ	SN:seq1	LN:1575
@SQ	SN:seq2	LN:1584
B7_591:4:96:693:509	73	seq1	1	99	36M	*	0	0	CACTAGTGGCTCATTGTAAATGTGTGGTTTAACTCG	<<<<<<<<<<<<<<<;<<<<<<<<<5<<<<<;:<;7	H0:i:1	H1:i:0	MF:i:18	NM:i:0	UQ:i:0	Aq:i:73
EAS54_65:7:152:368:113	73	seq1	3	99	35M	*	0	0	CTAGTGGCTCATTGTAAATGTGTGGTTTAACTCGT	<<<<<<<<<<0<<<<655<<7<<<:9<<3/:<6):H0:i:1	H1:i:0	MF:i:18	NM:i:0	UQ:i:0	Aq:i:66
EAS51_64:8:5:734:57	137	seq1	5	99	35M	*	0	0	AGTGGCTCATTGTAAATGTGTGGTTTAACTCGTCC	<<<<<<<<<<<7;71<<;<;;<7;<<3;);3*8/5H0:i:1	H1:i:0	MF:i:18	NM:i:0	UQ:i:0	Aq:i:66
B7_591:1:289:587:906	137	seq1	6	63	36M	*	0	0	GTGGCTCATTGTAATTTTTTGTTTTAACTCTTCTCT	(-&----,----)-)-),'--)---',+-,),''*,	H0:i:0	H1:i:0	MF:i:130	NM:i:5	UQ:i:38	Aq:i:63
EAS56_59:8:38:671:758	137	seq1	9	99	35M	*	0	0	GCTCATTGTAAATGTGTGGTTTAACTCGTCCATGG	<<<<<<<<<<<<<<<;<;7<<<<<<<<7<<;:<5%H0:i:1	H1:i:0	MF:i:18	NM:i:0	UQ:i:0	Aq:i:72
EAS56_61:6:18:467:281	73	seq1	13	99	35M	*	0	0	ATTGTAAATGTGTGGTTTAACTCGTCCCTGGCCCA	<<<<<<<<;<<<8<<<<<;8:;6/686&;(16666H0:i:0	H1:i:1	MF:i:18	NM:i:1	UQ:i:5	Aq:i:39
EAS114_28:5:296:340:699	137	seq1	13	99	36M	*	0	0	ATTGTAAATGTGTGGTTTAACTCGTCCATGGCCCAG	<<<<<;<<<;<;<<<<<<<<<<<8<8<3<8;<;<0;	H0:i:1	H1:i:0	MF:i:18	NM:i:0	UQ:i:0	Aq:i:73
B7_597:6:194:894:408	73	seq1	15	99	35M	*	0	0	TGTAAATGTGTGGTTTAACTCGTCCATTGCCCAGC	<<<<<<<<<7<<;<<<<;<<<7;;<<<*,;;572<H0:i:0	H1:i:1	MF:i:18	NM:i:1	UQ:i:9	Aq:i:43
EAS188_4:8:12:628:973	89	seq1	18	75	35M	*	0	0	AAATGTGTGGTTTAACTCGTCCATGGCCCAGCATT	==;=:;:;;:====;=;===:=======;==;===H0:i:1	H1:i:0	MF:i:64	NM:i:0	UQ:i:0	Aq:i:0
(...)
```

### Example 2

remove reads with indels:

```
java -jar dist/samjs.jar -e 'function accept(r) { if(r.getReadUnmappedFlag()) return false; var cigar=r.getCigar();if(cigar==null) return false; for(var i=0;i< cigar.numCigarElements();++i) {if(cigar.getCigarElement(i).getOperator().isIndelOrSkippedRegion()) return false; } return true;} accept(record);' input.bam
```


END_DOC
*/
@Program(name="samjs",
	description="Filters a BAM using a javascript expression ( java nashorn engine  ).",
	keywords={"sam","bam","nashorn","javascript","filter"},
	biostars={75168,81750,75354,77802,103052,106900,150530,253774,256615},
	references="\"bioalcidae, samjs and vcffilterjs: object-oriented formatters and filters for bioinformatics files\" . Bioinformatics, 2017. Pierre Lindenbaum & Richard Redon  [https://doi.org/10.1093/bioinformatics/btx734](https://doi.org/10.1093/bioinformatics/btx734)."
	)
public class SamJavascript
	extends Launcher
	{
	private static final Logger LOG = Logger.build(SamJavascript.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;


	@Parameter(names={"-X","--fail"},description="Save dicarded reads in that file")
	private File failingReadsFile = null;

	@Parameter(names={"-N","--limit"},description="limit to 'N' records (for debugging).")
	private long LIMIT = -1L ;

	@ParametersDelegate
	private WritingBamArgs writingBamArgs = new WritingBamArgs();

	@Parameter(names={"-e","--expression"},description="javascript expression")
	private String jsExpression=null;
	@Parameter(names={"-f","--file"},description="javascript file")
	private File jsFile =null;
	private CompiledScript  script=null;
	private SAMFileWriter failingReadsWriter=null;

	public SamJavascript()
		{
		}

	/* open failing bam if it was not already open */
	private void openFailing(final SAMFileHeader h)
		{
		if(this.failingReadsFile==null) return;
		if(this.failingReadsWriter==null)
			{
			LOG.info("Writing failings to "+ this.failingReadsFile);
			final SAMFileHeader h2= h.clone();
			this.failingReadsWriter=this.writingBamArgs.openSAMFileWriter(failingReadsFile, h2,true);
			}
		}

	private void failing(final SAMRecord rec,final SAMFileHeader h)
		{
		openFailing(h);
		if(failingReadsWriter!=null) failingReadsWriter.addAlignment(rec);
		}

	@Override
	public int doWork(final List<String> args) {
		SAMRecordIterator iter=null;
		SamReader samFileReader=null;
		SAMFileWriter sw=null;
		try
			{
			this.script  = super.compileJavascript(this.jsExpression,this.jsFile);
			samFileReader= openSamReader(oneFileOrNull(args));
			final SAMFileHeader header=samFileReader.getFileHeader();
			sw = writingBamArgs.openSAMFileWriter(outputFile,header, true);
			long count=0L;
	        	final Bindings bindings = this.script.getEngine().createBindings();
		        bindings.put("header", samFileReader.getFileHeader());
		        final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header).logger(LOG);
		        iter = samFileReader.iterator();
			while(iter.hasNext())
				{
				final SAMRecord record=iter.next();
				progress.watch(record);
				bindings.put("record", record);
				if(super.evalJavaScriptBoolean(this.script, bindings))
					{
					++count;
					sw.addAlignment(record);
					if(this.LIMIT>0L && count>=this.LIMIT) break;
					}
				else
					{
					failing(record,header);
					}
				}
			sw.close();
			/* create empty if never called */
			openFailing(header);
			return RETURN_OK;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(samFileReader);
			CloserUtil.close(sw);
			CloserUtil.close(failingReadsWriter);
			}
		}

	public static void main(String[] args) throws Exception
		{
		new SamJavascript().instanceMainWithExit(args);
		}
	}
