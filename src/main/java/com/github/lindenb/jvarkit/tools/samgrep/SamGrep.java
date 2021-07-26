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
package com.github.lindenb.jvarkit.tools.samgrep;


import java.io.BufferedReader;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

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
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.OnePassBamLauncher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.StringUtil;


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
		keywords={"sam","bam"},
		modificationDate="20210726",
		creationDate="20130506"
		)
public class SamGrep extends OnePassBamLauncher {
	private static final Logger LOG = Logger.build(SamGrep.class).make();
	
	@Parameter(names={"-R","--readname"},description="add the read name")
	private Set<String> nameStrings = new HashSet<>();
	
	
	@Parameter(names={"-f","--readfile"},description="file containing a list of read names")
	private Path namefile = null;
	
	@Parameter(names={"-n","--stopafter"},description="when found, remove the read from the list of names when found more that 'n' time (increase speed)")
	private int n_before_remove = -1 ;
	
	@Parameter(names={"-V","--invert"},description="invert")
	private boolean inverse = false;
		

	private final Map<String,Integer> readNames=new HashMap<String,Integer>(); 
   
    
	@Override
	protected Logger getLogger() 	{
		return LOG;
		}
	
	@Override
	protected int beforeSam()
		{
    	if(namefile!=null) {
			try(BufferedReader in=IOUtils.openPathForBufferedReading(this.namefile)) {
				in.lines().
						filter(L->!StringUtil.isBlank(L)).
						forEach(L->readNames.put(L.trim(),0));
				}
			catch(Throwable err)
				{
				LOG.error(err);
				return -1;
				}
	    	}
    	for(final String line: this.nameStrings) {
    		readNames.put(line,0);
    		}
    	if(readNames.isEmpty())
			{
			LOG.warn("no read found.");
			}
    	return super.beforeSam();
    	}
    
	
	@Override
	protected void scanIterator(
			SAMFileHeader headerIn,
			CloseableIterator<SAMRecord> iter,
			SAMFileWriter sfw)
			{		
			while(iter.hasNext())
				{
				boolean keep=false;
				final SAMRecord rec= iter.next();
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
						}
					else
						{
						readNames.put(rec.getReadName(),count);
						}
					}
				}
			}
    
    public static void main(final String[] argv)
		{
	    new SamGrep().instanceMainWithExit(argv);
		}	

	}
