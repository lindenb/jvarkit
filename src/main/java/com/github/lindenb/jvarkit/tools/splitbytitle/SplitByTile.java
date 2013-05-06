package com.github.lindenb.jvarkit.tools.splitbytitle;

import java.io.File;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Pattern;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class SplitByTile  extends CommandLineProgram
	{
	private static final String TILEWORD="__TILE__";
	private static final Log log = Log.getInstance(SplitByTile.class);
    @Usage
    public String USAGE = getStandardUsagePreamble() + "Split Bam By tile";
    
    
    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="A BAM file to process.")
    public File INPUT=null;
    @Option(shortName= "O", doc="The output file. must ends with .bam and must contains "+TILEWORD,optional=false)
    public String OUTPUT=null;
    
    public static void main(final String[] argv)
		{
	    new SplitByTile().instanceMainWithExit(argv);
		}	

    
    @Override
    protected int doWork() {
    	IoUtil.assertFileIsReadable(INPUT);
    	if(OUTPUT==null || !OUTPUT.contains(TILEWORD) || !OUTPUT.endsWith(".bam"))
    		{
    		log.error("Bad OUPUT name "+OUTPUT+". must contain "+TILEWORD+" and ends with .bam");
    		return -1;
    		}
    	
        SAMFileReader samReader = null;
        Map<Integer, SAMFileWriter> tile2writer=new HashMap<Integer, SAMFileWriter>();
        Pattern colon=Pattern.compile("[\\:]");
        SAMRecordIterator iter=null;
        SAMFileWriterFactory sfrFactory=new SAMFileWriterFactory();
        sfrFactory.setCreateIndex(super.CREATE_INDEX);
        sfrFactory.setCreateMd5File(super.CREATE_MD5_FILE);
        
        try
	        {
	        samReader=new SAMFileReader(INPUT);
	        final SAMFileHeader header=samReader.getFileHeader();
	        samReader.setValidationStringency(super.VALIDATION_STRINGENCY);
			iter=samReader.iterator();
			while(iter.hasNext())
				{
				SAMRecord rec=iter.next();
				String tokens[]=colon.split(rec.getReadName(),6);
				if(tokens.length<5)
					{
		    		log.error("Cannot get the 6th field in "+rec.getReadName());
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
					log.error("Bad tile in read: "+rec.getReadName());
					return -1;
					}
				
				SAMFileWriter sfw=tile2writer.get(tile);
				if(sfw==null)
					{
					
					File outFile=new File(OUTPUT.replaceAll(TILEWORD, String.valueOf(tile)));
					log.info("create output for "+outFile );
					if(outFile.getParentFile()!=null)
						{
						outFile.getParentFile().mkdirs();
						}
					sfw=sfrFactory.makeBAMWriter(header, header.getSortOrder()==SortOrder.coordinate,outFile);
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
    		log.error(e);
    		return -1;
			}
        finally
	    	{
	    	if(iter!=null) iter.close();
	    	if(samReader!=null) samReader.close();
	    	}
    	return 0;
    	}
    
	}
