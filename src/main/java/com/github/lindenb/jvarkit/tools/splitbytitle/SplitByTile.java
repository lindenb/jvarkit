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
package com.github.lindenb.jvarkit.tools.splitbytitle;


/**
BEGIN_DOC


END_DOC
 */

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;

@Program(name="splitbytile",description="Split Bam By tile",keywords={"sam","bam"})
public class SplitByTile  extends Launcher
	{
	private static final String TILEWORD="__TILE__";
    private static final Logger LOG = Logger.build(SplitByTile.class).make();
    @Parameter(names={"-o","--output"},description="Output file. Must contain "+TILEWORD,required=true)
    private String OUTPUT = null;
    @ParametersDelegate
    private WritingBamArgs writingBamArgs =new WritingBamArgs();
    
    
    public static void main(final String[] argv)
		{
	    new SplitByTile().instanceMainWithExit(argv);
		}	
    
    @Override
    public int doWork(final List<String> args) {
    	if(OUTPUT==null || !OUTPUT.contains(TILEWORD) || !OUTPUT.endsWith(".bam"))
    		{
    		LOG.error("Bad OUPUT name "+OUTPUT+". must contain "+TILEWORD+" and ends with .bam");
    		return -1;
    		}
    	
        SamReader samReader = null;
        final  Map<Integer, SAMFileWriter> tile2writer=new HashMap<Integer, SAMFileWriter>();
        final  Pattern colon=Pattern.compile("[\\:]");
        SAMRecordIterator iter=null;
     
        
        try
	        {
	        samReader= super.openSamReader(oneFileOrNull(args));
	        final SAMFileHeader header=samReader.getFileHeader();
			iter=samReader.iterator();
			while(iter.hasNext())
				{
				final SAMRecord rec=iter.next();
				final String tokens[]=colon.split(rec.getReadName(),6);
				if(tokens.length<5)
					{
		    		LOG.error("Cannot get the 6th field in "+rec.getReadName());
		    		samReader.close();samReader=null;
		    		return -1;
					}
				int tile=-1;
				try {
					tile=Integer.parseInt(tokens[4]);
					}
				catch (Exception e)
					{
					tile=-1;
					}
				if(tile<0)
					{
					LOG.error("Bad tile in read: "+rec.getReadName());
					samReader.close();samReader=null;
					return -1;
					}
				
				SAMFileWriter sfw=tile2writer.get(tile);
				if(sfw==null)
					{
					final File outFile=new File(OUTPUT.replaceAll(TILEWORD, String.valueOf(tile)));
					LOG.info("create output for "+outFile );
					if(outFile.getParentFile()!=null)
						{
						outFile.getParentFile().mkdirs();
						}
					sfw= this.writingBamArgs.openSAMFileWriter(outFile,header,true);
					tile2writer.put(tile, sfw);
					}
				sfw.addAlignment(rec);
				}
			
			for(final SAMFileWriter sfw:tile2writer.values())
				{
				sfw.close();
				}
			} 
    	catch (final Exception e) {
    		LOG.error(e);
    		return -1;
			}
        finally
	    	{
	    	if(iter!=null) iter.close();
	    	CloserUtil.close(samReader);
	    	}
    	return 0;
    	}
    
	}
