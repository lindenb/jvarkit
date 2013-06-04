package com.github.lindenb.jvarkit.tools.samgrep;


import java.io.BufferedReader;
import java.io.File;
import java.util.HashSet;
import java.util.Set;


import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMFileHeader.SortOrder;

public class SamGrep extends CommandLineProgram
	{
	private static final Log log = Log.getInstance(SamGrep.class);
	@Usage
	public String USAGE=getStandardUsagePreamble()+". grep read-names in a bam file";
    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="A BAM file to process.")
    public File INPUT=null;
    @Option(shortName= StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The output file. default: stdout",optional=true)
    public File OUTPUT=null;
    @Option(shortName="N", doc="file containing a list of read names",optional=false)
    public File NAMES=null;
    @Option(shortName="V", doc="inverse",optional=false)
    public boolean inverse=false;
    @Option(shortName="ONE", doc="remove read from the list.",optional=false)
    public boolean removeRead=false;
    
    
	@Override
	protected int doWork()
		{
    	IoUtil.assertFileIsReadable(INPUT);
		Set<String> readNames=new HashSet<String>(); 
    	if(OUTPUT!=null && !OUTPUT.getName().endsWith(".bam"))
    		{
    		log.error("Bad OUPUT name "+OUTPUT+".");
    		return -1;
    		}
    	
    	
        SAMFileReader samReader = null;
        SAMRecordIterator iter=null;
        SAMFileWriterFactory sfwFactory=new SAMFileWriterFactory();
        if(OUTPUT!=null)
        	{
        	sfwFactory.setCreateIndex(super.CREATE_INDEX);
        	sfwFactory.setCreateMd5File(super.CREATE_MD5_FILE);
        	}
        SAMFileWriter sfw=null;
        try
	        {
	    	BufferedReader in=IoUtil.openFileForBufferedReading(NAMES);
	    	String line;
	    	while((line=in.readLine())!=null)
	    		{
	    		if(line.isEmpty()) continue;
	    		readNames.add(line);
	    		}
	    	in.close();
	    	
	    	if(readNames.isEmpty())
	    		{
	    		
	    		}

	        
	        samReader=new SAMFileReader(INPUT);
	        final SAMFileHeader header=samReader.getFileHeader();
	        samReader.setValidationStringency(super.VALIDATION_STRINGENCY);
			iter=samReader.iterator();
			if(OUTPUT==null)
				{
				sfw=sfwFactory.makeSAMWriter(header, header.getSortOrder()==SortOrder.coordinate, System.out);
				}
			else
				{
				sfw=sfwFactory.makeBAMWriter(header, header.getSortOrder()==SortOrder.coordinate, OUTPUT);
				}
			while(iter.hasNext())
				{
				SAMRecord rec=iter.next();
				boolean keep=false;
				if(readNames.contains(rec.getReadName()))
					{
					keep=true;
					
					}
				if(inverse) keep=!keep;
				if(keep)
					{
					sfw.addAlignment(rec);
					}
				
				if(removeRead && !inverse && keep)
					{
					readNames.remove(rec.getReadName());
					if(readNames.isEmpty()) break;
					}
				}
			sfw.close();	
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
    public static void main(final String[] argv)
		{
	    new SamGrep().instanceMainWithExit(argv);
		}	

	}
