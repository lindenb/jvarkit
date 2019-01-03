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
package com.github.lindenb.jvarkit.tools.samgrep;


import java.io.BufferedReader;
import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.github.lindenb.jvarkit.io.IOUtils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.StringUtil;

/**

BEGIN_DOC


### Examples


#### Example 1


```

java -jar  dist/samgrep.jar -R r001  -- samtools-0.1.18/examples/toy.sam 

@HD     VN:1.4  SO:unsorted
@SQ     SN:ref  LN:45
@SQ     SN:ref2 LN:40
@PG     ID:0    PN:com.github.lindenb.jvarkit.tools.samgrep.SamGrep     VN:dac03b80e9fd88a15648b22550e57d10c9bed725     CL:-R r001 samtools-0.1.18/examples/toy.sam
r001    163     ref     7       30      8M4I4M1D3M      =       37      39      TTAGATAAAGAGGATACTG     *       XX:B:S,12561,2,20,112
r001    83      ref     37      30      9M      =       7       -39     CAGCGCCAT       *

```


#### Example 4


```

java -jar  dist/samgrep.jar -R r001 -- -n 1 samtools-0.1.18/examples/toy.sam 

@HD     VN:1.4  SO:unsorted
@SQ     SN:ref  LN:45
@SQ     SN:ref2 LN:40
@PG     ID:0    PN:com.github.lindenb.jvarkit.tools.samgrep.SamGrep     VN:dac03b80e9fd88a15648b22550e57d10c9bed725     CL:-R r001 -n 1 samtools-0.1.18/examples/toy.sam
r001    163     ref     7       30      8M4I4M1D3M      =       37      39      TTAGATAAAGAGGATACTG     *       XX:B:S,12561,2,20,112

```






END_DOC
*/


import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;


/**

BEGIN_DOC

### Examples

#### Example 1

```

java -jar  dist/samgrep.jar -R r001  -- samtools-0.1.18/examples/toy.sam 

@HD	VN:1.4	SO:unsorted
@SQ	SN:ref	LN:45
@SQ	SN:ref2	LN:40
@PG	ID:0	PN:com.github.lindenb.jvarkit.tools.samgrep.SamGrep	VN:dac03b80e9fd88a15648b22550e57d10c9bed725	CL:-R r001 samtools-0.1.18/examples/toy.sam
r001	163	ref	7	30	8M4I4M1D3M	=	37	39	TTAGATAAAGAGGATACTG	*	XX:B:S,12561,2,20,112
r001	83	ref	37	30	9M	=	7	-39	CAGCGCCAT	*

```





#### Example 4



```

java -jar  dist/samgrep.jar -R r001 -- -n 1 samtools-0.1.18/examples/toy.sam 

@HD	VN:1.4	SO:unsorted
@SQ	SN:ref	LN:45
@SQ	SN:ref2	LN:40
@PG	ID:0	PN:com.github.lindenb.jvarkit.tools.samgrep.SamGrep	VN:dac03b80e9fd88a15648b22550e57d10c9bed725	CL:-R r001 -n 1 samtools-0.1.18/examples/toy.sam
r001	163	ref	7	30	8M4I4M1D3M	=	37	39	TTAGATAAAGAGGATACTG	*	XX:B:S,12561,2,20,112

```




END_DOC
*/


@Program(name="samgrep",description="grep read-names in a bam file",
		keywords={"sam","bam"})
public class SamGrep extends Launcher
	{
	private static final Logger LOG = Logger.build(SamGrep.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	
	@Parameter(names={"-R","--readname"},description="add the read name")
	private Set<String> nameStrings = new HashSet<>();
	
	
	@Parameter(names={"-f","--readfile"},description="file containing a list of read names")
	private File namefile = null;
	
	@Parameter(names={"-x","--tee"},description="if output fileame specified, continue to output original input to stdout.")
	private boolean divertToStdout = false;
	
	@Parameter(names={"-n","--stopafter"},description="when found, remove the read from the list of names when found more that 'n' time (increase speed)")
	private int n_before_remove = -1 ;
	
	@Parameter(names={"-V","--invert"},description="invert")
	private boolean inverse = false;
	
	@ParametersDelegate
	private WritingBamArgs writingBamArgs=new WritingBamArgs();
	

	private final Map<String,Integer> readNames=new HashMap<String,Integer>(); 
   
    
    @Override
    public int doWork(final List<String> args) {
    	
    	readNames.clear();
    	
    	if(namefile!=null) {
	    	BufferedReader in=null;
			try
				{
				in=IOUtils.openFileForBufferedReading(this.namefile);
				in.lines().
						filter(L->!StringUtil.isBlank(L)).
						forEach(L->readNames.put(L.trim(),0));
				}
			catch(Exception err)
				{
				LOG.error(err);
				return -1;
				}
			finally
				{
				CloserUtil.close(in);
				}
	    	}
    	for(final String line: this.nameStrings) {
    		readNames.put(line,0);
    		}
    	if(readNames.isEmpty())
			{
			LOG.warn("no read found.");
			}
    	
    	
		SAMFileWriter sfw=null;
		SAMFileWriter samStdout=null;
		SamReader sfr=null;
		try {
			sfr = super.openSamReader(oneFileOrNull(args));
			final SAMFileHeader header=sfr.getFileHeader().clone();
			final ProgressFactory.Watcher<SAMRecord> progress= ProgressFactory.
					newInstance().
					logger(LOG).
					dictionary(header).
					build();

			
			
			if(this.outputFile==null)
				{
				sfw = this.writingBamArgs.openSAMFileWriter(null, header, true);
				}
			else
				{
				samStdout= this.writingBamArgs.openSAMFileWriter(null, header, true);
				
				sfw= this.writingBamArgs.openSAMFileWriter(outputFile, header, true);
				}
		
			SAMRecordIterator iter=sfr.iterator();
			while(iter.hasNext())
				{
				boolean keep=false;
				final SAMRecord rec=progress.apply(iter.next());
				if(samStdout!=null) samStdout.addAlignment(rec);
				Integer count = readNames.get(rec.getReadName());
				if(count!=null)
					{
					keep=true;
					}
				if(this.inverse) keep=!keep;
				if(keep)
					{
					sfw.addAlignment(rec);
					}
				
				if(n_before_remove!=-1 && !inverse && keep)
					{
					count++;
					if(count>=n_before_remove)
						{
						readNames.remove(rec.getReadName());
						if(samStdout==null && readNames.isEmpty()) break;
						}
					else
						{
						readNames.put(rec.getReadName(),count);
						}
					}
				}
			progress.close();
			return RETURN_OK;
		} catch (final Exception err) {
			LOG.error(err);
			return -1;
		} finally {
			CloserUtil.close(samStdout);
			CloserUtil.close(sfw);
			CloserUtil.close(sfr);
			}
		}
    
    public static void main(final String[] argv)
		{
	    new SamGrep().instanceMainWithExit(argv);
		}	

	}
