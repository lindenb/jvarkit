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
package com.github.lindenb.jvarkit.tools.splitbytitle;


/**
BEGIN_DOC


END_DOC
 */

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Map;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.jcommander.MultiBamLauncher;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;

@Program(name="splitbytile",description="Split Bam By tile",
		keywords={"sam","bam"},
		modificationDate="20210304",
		creationDate="20130406"
		)
public class SplitByTile  extends MultiBamLauncher {
	private static final String TILEWORD="__TILE__";
    private static final Logger LOG = Logger.build(SplitByTile.class).make();
    @Parameter(names={"-o","--output"},description="Output file. Must contain "+TILEWORD,required=true)
    private Path OUTPUT = null;
    @ParametersDelegate
    private WritingBamArgs writingBamArgs =new WritingBamArgs();
    
    
    @Override
    protected Logger getLogger() {
    	return LOG;
    	}
    @Override
    protected int beforeSam() {
    	if(OUTPUT==null || !OUTPUT.toString().contains(TILEWORD) ||
    			!(OUTPUT.toString().endsWith(FileExtensions.BAM) || OUTPUT.toString().endsWith(FileExtensions.CRAM)))
			{
			LOG.error("Bad OUPUT name "+OUTPUT+". must contain "+TILEWORD+" and ends with .bam|.cram");
			return -1;
			}    
    	return super.beforeSam();
    	}
    
    
    @Override
    protected int processInput(SAMFileHeader header, CloseableIterator<SAMRecord> iter) {
        final  Map<Integer, SAMFileWriter> tile2writer=new HashMap<Integer, SAMFileWriter>();
        final CharSplitter  colon= CharSplitter.COLON;
     
        try
	        {
			while(iter.hasNext())
				{
				final SAMRecord rec=iter.next();
				final String tokens[]=colon.split(rec.getReadName(),6);
				if(tokens.length<5)
					{
		    		LOG.error("Cannot get the 6th field in "+rec.getReadName());
		    		return -1;
					}
				int tile=-1;
				try {
					tile=Integer.parseInt(tokens[4]);
					}
				catch(Throwable e)
					{
					tile=-1;
					}
				if(tile<0)
					{
					LOG.error("Bad tile in read: "+rec.getReadName());
					return -1;
					}
				
				SAMFileWriter sfw=tile2writer.get(tile);
				if(sfw==null)
					{
					final Path outFile= Paths.get(OUTPUT.toString().replaceAll(TILEWORD, String.valueOf(tile)));
					LOG.info("create output for "+outFile );
					if(outFile.getParent()!=null && !Files.exists(outFile.getParent()))
						{
						Files.createDirectories(outFile.getParent());
						}
					sfw= this.writingBamArgs.
							setReferencePath(super.faidxPath).
							openSamWriter(outFile,header,true);
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
	    	}
    	return 0;
    	}
    
    public static void main(final String[] argv)
		{
	    new SplitByTile().instanceMainWithExit(argv);
		}

	}
