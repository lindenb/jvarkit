package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Reader;
import java.util.ArrayList;
import java.util.List;

import org.broadinstitute.variant.vcf.VCFCodec;

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.util.Log;

import com.github.lindenb.jvarkit.util.picard.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.IOUtils;

public class NoEmptyVCF extends AbstractCommandLineProgram
	{
	private static Log LOG=Log.getInstance(NoEmptyVCF.class);

    @Usage(programVersion="1.0")
    public String USAGE = getStandardUsagePreamble() + "  If VCF is empty or doesn't exists, create a dummy one. ";

    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME,doc="VCF input (or stdin)",optional=true)
    public File IN=null;
    
	
    @Option(shortName="S",doc="add this sample",minElements=0)
    public List<String> SAMPLES=new ArrayList<String>();
    
    private void write(File f) throws IOException
    	{
    	LOG.info("fixing "+f);
		PrintWriter w=new PrintWriter(IOUtils.openFileForBufferedWriting(IN));
		write(w);
		w.close();
    	}
    
    private void write(PrintWriter out) throws IOException
    	{
		out.println(VCFCodec.VCF4_MAGIC_HEADER);
		out.println("##source=noEmptyVCF");
		out.print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
		if(!this.SAMPLES.isEmpty())
			{
			System.out.print("\tFORMAT");
			for(String S:this.SAMPLES)
				{
				System.out.print("\t");
				System.out.print(S);
				}
			}
		out.println();
		out.flush();
    	}
    
	@Override
	protected int doWork()
		{
		try
			{
			if(IN!=null)
				{
				if(IN.exists())
					{
					LOG.info("opening "+IN+" ... ");
					Reader reader=IOUtils.openFileForBufferedReading(IN);
					int c=reader.read();
					reader.close();
					reader=null;
					if(c!=-1 && c!='#')
						{
						LOG.error("File "+IN+" doesn't start with # but with ascii("+c+")");
						return -1;
						}
					if(c==-1)
						{
						write(this.IN);
						}
					}
				else
					{
					LOG.error("File "+IN+" doesn't exists");
					write(this.IN);
					}
				}
			else
				{
				LOG.info("reading from stdin ... ");
				int c=System.in.read();
				if(c!=-1 && c!='#')
					{
					LOG.warn("VCF starts with # but with ascii("+c+")");
					IOUtils.copyTo(System.in, System.out);
					return 0;
					}
				if(c==-1)
					{
					LOG.warn("writing empty VCF");
					write(new PrintWriter(System.out));
					}
				else
					{
					IOUtils.copyTo(System.in, System.out);
					}
				}
			
		
			}
		catch(IOException err)
			{
			LOG.error(err,"Boum");
			return -1;
			}
		return 0;
		}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new NoEmptyVCF().instanceMainWithExit(args);

	}

}
