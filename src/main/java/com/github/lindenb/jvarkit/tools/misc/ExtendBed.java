package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.PrintWriter;
import java.util.regex.Pattern;

import org.broad.tribble.readers.AsciiLineReader;
import org.broad.tribble.readers.LineReader;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryFactory;

import net.sf.picard.PicardException;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.util.Log;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

/**
 * Extends a BED by 'X' bases.
 * @author lindenb
 *
 */
public class ExtendBed extends AbstractCommandLineProgram
	{
	private static Log LOG=Log.getInstance(ExtendBed.class);

    @Usage(programVersion="1.0")
    public String USAGE = getStandardUsagePreamble() + " extends a BED file by 'X' bases. ";

	
    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME,doc="BED Input URI/file. default: stdin",optional=true)
    public String IN=null;

    @Option(shortName=StandardOptionDefinitions.REFERENCE_SHORT_NAME,doc="Reference",optional=false)
    public File REF=null;
	private SAMSequenceDictionary samSequenceDictionary;
	
    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME,doc="output name (default: stdout)",optional=true)
    public File OUT=null;
    
    @Option(shortName="X",doc="extend by 'X' bases.",optional=true)
    public int EXTEND=0;

    private void message(String s)
    	{
    	switch(super.VALIDATION_STRINGENCY)
    		{
    		case SILENT:return;
    		case LENIENT: LOG.warn(s);break;
    		default: LOG.error(s); throw new PicardException(s);
    		}
    	}
    
	@Override
	protected int doWork()
		{
		if(EXTEND<0)
			{
			LOG.error("negative value for X");
			return -1;
			}
		LineReader r=null;
		PrintWriter pw=new PrintWriter(System.out);
		try {
			this.samSequenceDictionary=new SAMSequenceDictionaryFactory().load(REF);
			if(IN==null)
				{
				LOG.info("reading from stdin");
				r=new AsciiLineReader(System.in);
				}
			else
				{
				r=new AsciiLineReader(IOUtils.openURIForReading(IN));
				}
			if(OUT!=null)
				{
				pw=new PrintWriter(OUT);
				}
			Pattern tab=Pattern.compile("[\t]");
			String line;
			while((line=r.readLine())!=null)
				{
				String tokens[]=tab.split(line);
				if(tokens.length<3)
					{
					message("Not enough cols in "+line);
					continue;
					}
				SAMSequenceRecord rec;
				if((rec=this.samSequenceDictionary.getSequence(tokens[0]))==null)
					{
					message("Chromosome "+tokens[0]+" not in REF dictionary. ignoring");
					continue;
					}
				int start=Integer.parseInt(tokens[1]);
				if(start<0)
					{
					message("start<0 in "+line );
					start=0;
					}
				int end=Integer.parseInt(tokens[2]);
				if(start>end)
					{
					message("ignoring : end<start in "+line );
					continue;
					}
				start=Math.max(0, start-EXTEND);
				end=Math.min(end+EXTEND,rec.getSequenceLength());
				for(int i=0;i< tokens.length;++i)
					{
					if(i>0) pw.print('\t');
					switch(i)
							{
							case 1: pw.print(start); break;
							case 2: pw.print(end); break;
							default: pw.print(tokens[i]); break;
							}
					}
				pw.println();
				}	
			} 
		catch (Exception e)
			{
			LOG.error(e);
			super.testRemoteGit();
			return -1;
			}
		finally
			{
			if(r!=null) r.close();
			pw.flush();
			pw.close();
			}
		return 0;
	}
	public static void main(String[] args) {
		new ExtendBed().instanceMainWithExit(args);
	}
	}
