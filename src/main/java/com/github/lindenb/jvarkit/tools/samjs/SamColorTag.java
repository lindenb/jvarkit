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



import java.awt.Color;
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
import com.github.lindenb.jvarkit.util.swing.ColorUtils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloserUtil;

/**
BEGIN_DOC

## How it works

Fill the  UCSC 'YC' color tag in a BAM using a javascript expression using
the nashorn javascript engine. The YC is used by IGV or the UCSC to colorize the reads

The script engine inject 'header' ( a htsjdk.samtools.SAMFileHeader ) and 'record' (a htsjdk.samtools.SAMRecord).
The script should return 

* a java.awt.color
* a String for a named color ('blue', 'red'...)
* a hexa color #FFFFF
* a rgb color 'rgb(100,200,100)'
* null or empty string (no YC tag, the tag is cleared)

## Example

The script:

```
function rnd(n) {
    return  Math.floor(Math.random() *n) ;
	}

function randomColor()
	{
	switch(rnd(10))
		{
		case 1: return "beige";
		case 2: return "#FF00AA";
		default: return "rgb("+rnd(255)+","+rnd(255)+","+rnd(255)+")";
		}
	return true;
	}

randomColor();
```

usage:

```
$ java -jar dist/samcolortag.jar -f script.js -o in.bam  out.bam

[lindenb@kaamelot-master01 jvarkit-git]$ samtools view out.bam | head
rotavirus_1_317_5:0:0_7:0:0_2de	99	rotavirus	1	60	70M	=	248	317	GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAATATGGCGTCAACTCAGCAGATGGTCAGCTCTAATATT	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	YC:Z:35,107,5	MD:Z:33G4A3T14A7T4	RG:Z:S1	NM:i:5	AS:i:45	XS:i:0
rotavirus_1_535_4:0:0_4:0:0_1a6	163	rotavirus	1	60	70M	=	466	535	GGCTTTTACTGCTTTTCAGTGGTTGCTTCTCAAGATGGAGTGTACTCATCAGATGGTAAGCTCTATTATT	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	YC:Z:66,172,98	MD:Z:8A18G13C6G21	RG:Z:S1	NM:i:4	AS:i:50	XS:i:0
rotavirus_1_543_5:0:0_11:0:0_390	163	rotavirus	1	60	70M	=	487	530	GGCTTTTAATGCTTTTCATTTGATGCTGCTCAAGATGGAGTCTACACAGCAGATGGTCAGCTCTATTATT	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	YC:Z:245,245,220	MD:Z:18G1G1T22T11A12	RG:Z:S1	NM:i:5	AS:i:45	XS:i:0
rotavirus_1_578_3:0:0_7:0:0_7c	99	rotavirus	1	60	70M	=	509	578	GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAAGATGGAGTCTCCTGAGCAGCTGGTAAGCTCTATTATT	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	YC:Z:67,135,235	MD:Z:43A2C5A17	RG:Z:S1	NM:i:3	AS:i:55	XS:i:0
rotavirus_1_497_4:0:0_5:0:0_2f6	163	rotavirus	1	60	70M	=	432	497	GGCATTTAATGCTTAACAGTGGTTGCTGCTCAAGATGGAGTCTACTCAGCAGATGGTAAGCTCTCTTATT	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	YC:Z:82,150,221	MD:Z:3T10T0T48A5	RG:Z:S1	NM:i:4	AS:i:51	XS:i:0
```

Another example, using the read mapping quality

```
var c=record.getMappingQuality();
"rgb("+c+","+c+","+c+")";
```


## See also:
  * com.github.lindenb.jvarkit.util.swing.ColorUtils
  * See http://software.broadinstitute.org/software/igv/book/export/html/6
  * http://genome.ucsc.edu/goldenPath/help/hgBamTrackHelp.html

## Screenshot

<img src="https://pbs.twimg.com/media/C_eefreXUAAO-rX.jpg"/>

END_DOC
*/
@Program(
	name="samcolortag",
	description="Add the UCSC 'YC' color tag in a BAM. See http://software.broadinstitute.org/software/igv/book/export/html/6 and http://genome.ucsc.edu/goldenPath/help/hgBamTrackHelp.html",
	keywords={"sam","bam","metadata","javascript","igv","visualization"}
	)
public class SamColorTag
	extends Launcher
	{
	private static final Logger LOG = Logger.build(SamColorTag.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@ParametersDelegate
	private WritingBamArgs writingBamArgs = new WritingBamArgs();
	
	@Parameter(names={"-e","--expression"},description="javascript expression")
	private String jsExpression=null;
	@Parameter(names={"-f","--file"},description="javascript file")
	private File jsFile =null;
	@Parameter(names={"-E","--ignoreErrors"},description="Ignore javascript/color errors")
	private boolean ignoreErrors=false;
	
	
	
	@Override
	public int doWork(final List<String> args) {
		SAMRecordIterator iter=null;
		SamReader samFileReader=null;
		SAMFileWriter sw=null;
		final ColorUtils colorUtils=new ColorUtils();
		try
			{
			final CompiledScript  script  = super.compileJavascript(this.jsExpression,this.jsFile);
			samFileReader= openSamReader(oneFileOrNull(args));
			final SAMFileHeader srcheader=samFileReader.getFileHeader();
			final SAMFileHeader header = srcheader.clone();
			header.addComment(ColorUtils.YC_TAG+" attribute added with "+getProgramName()+" "+getProgramCommandLine());
			sw = this.writingBamArgs.openSAMFileWriter(outputFile,header, true);
			
	        final Bindings bindings = script.getEngine().createBindings();
	        bindings.put("header", samFileReader.getFileHeader());
	        SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header);
	        iter = samFileReader.iterator();
			while(iter.hasNext())
				{
				final SAMRecord record= progress.watch(iter.next());
				bindings.put("record", record);
				
				
				final Color color;
				
				Object result;
				
				try 
					{
					result = script.eval(bindings);
					}
				catch(final Exception err)
					{
					if(!ignoreErrors) {
						LOG.error(err);
						return -1;
						}
					result=null;
					}
				
				if(result == null)
					{
					color=null;
					}
				else if(result instanceof Color)
					{
					color=Color.class.cast(result);
					}
				else if(result instanceof String)
					{
					final String s=(String)result;
					Color c2=null;
					try
						{
						if(s.trim().isEmpty()) {
							c2=null;
							}
						else
							{
							c2 = colorUtils.parse(s);
							}
						}
					catch(final Exception err)
						{
						if(!ignoreErrors) {
							LOG.error(err);
							return -1;
						}
						c2 = null;
						}
					color = c2;
					}
				else
					{
					if(!ignoreErrors) {
						LOG.error("Cannot cast to color a "+result.getClass());
						return -1;
						}
					color =null;
					}
				if(color!=null)
					{
					record.setAttribute(
						ColorUtils.YC_TAG,
						ColorUtils.colorToSamAttribute(color)
						);
					}
				else
					{
					//clear attribute
					record.setAttribute(ColorUtils.YC_TAG,null);
					}
				sw.addAlignment(record);
				}
			sw.close();
			sw=null;
			
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
			}
		}
	
		
	public static void main(String[] args) throws Exception
		{
		new SamColorTag().instanceMainWithExit(args);
		}
	
	

	
	}
