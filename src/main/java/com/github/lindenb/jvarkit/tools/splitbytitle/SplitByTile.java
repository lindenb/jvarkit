package com.github.lindenb.jvarkit.tools.splitbytitle;

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
    public int doWork(List<String> args) {
    	if(OUTPUT==null || !OUTPUT.contains(TILEWORD) || !OUTPUT.endsWith(".bam"))
    		{
    		LOG.error("Bad OUPUT name "+OUTPUT+". must contain "+TILEWORD+" and ends with .bam");
    		return -1;
    		}
    	
        SamReader samReader = null;
        Map<Integer, SAMFileWriter> tile2writer=new HashMap<Integer, SAMFileWriter>();
        Pattern colon=Pattern.compile("[\\:]");
        SAMRecordIterator iter=null;
     
        
        try
	        {
	        samReader= super.openSamReader(oneFileOrNull(args));
	        final SAMFileHeader header=samReader.getFileHeader();
			iter=samReader.iterator();
			while(iter.hasNext())
				{
				SAMRecord rec=iter.next();
				String tokens[]=colon.split(rec.getReadName(),6);
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
			
			for(SAMFileWriter sfw:tile2writer.values())
				{
				sfw.close();
				}
			} 
    	catch (Exception e) {
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
