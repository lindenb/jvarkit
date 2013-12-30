package com.github.lindenb.jvarkit.tools.vcfconcat;

import java.io.InputStream;

import net.sf.samtools.util.CloserUtil;

import org.broad.tribble.readers.LineIterator;
import org.broad.tribble.readers.LineIteratorImpl;
import org.broad.tribble.readers.LineReader;
import org.broad.tribble.readers.LineReaderUtil;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFCodec;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.cli.GetOpt;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

public class VcfConcat extends AbstractCommandLineProgram
	{
	private VcfConcat()
		{
		}
	

	@Override
	protected String getOnlineDocUrl()
		{
		return "https://github.com/lindenb/jvarkit/wiki/VcfConcat";
		}
	
	@Override
	public String getProgramDescription()
		{
		return "Concatenate VCF files.";
		}
	
	private int compareHeader(VCFHeader h1,VCFHeader h2)
		{
		if(!h1.getSampleNameToOffset().equals(h2.getSampleNameToOffset()))
			{
			error("Not the same sample list to concatenate:"+
					h1.getSampleNameToOffset()+" vs "+
					h2.getSampleNameToOffset()
					);
			return -1;
			}
		
		return 0;
		}
	
	private int fromFiles(VariantContextWriter w,int optind,String args[])
		{
		VCFHeader header0=null;
		for(int i=optind;i< args.length;++i)
			{
			
			InputStream in=null;
			try
				{
				VCFCodec codec=new VCFCodec();
				info("Opening "+args[i]);
				in=IOUtils.openURIForReading(args[i]);
				LineReader lr=LineReaderUtil.fromBufferedStream(in);
				LineIterator li=new LineIteratorImpl(lr);
				final VCFHeader header1=(VCFHeader)codec.readActualHeader(li);
				if(header0==null)
					{
					header0=header1;
					w.writeHeader(header0);
					}
				else if(compareHeader(header0,header1)!=0)
					{
					return -1;
					}
				
				while(li.hasNext())
					{
					w.add(codec.decode(li.next()));
					}
				info("Closing "+args[i]);
				}
			catch(Exception err)
				{
				error(err);
				return -1;
				}
			finally
				{
				CloserUtil.close(in);
				}
			}
		return 0;
		}
	private int fromStdin(VariantContextWriter w)
		{
		VCFCodec codec=new VCFCodec();
		LineReader lr=LineReaderUtil.fromBufferedStream(System.in);
		LineIterator li=new LineIteratorImpl(lr);
		final VCFHeader header0=(VCFHeader)codec.readActualHeader(li);
		header0.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		header0.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
		w.writeHeader(header0);
		
		while(li.hasNext())
			{
			String line=li.peek();
			if(line.startsWith("#"))
				{
				info("Found a new header");
				codec=new VCFCodec();
				VCFHeader header1=(VCFHeader)codec.readActualHeader(li);
				if(compareHeader(header0,header1)!=0)
					{
					return -1;
					}
				}
			else
				{
				w.add(codec.decode(li.next()));
				}
			}
		return 0;
		}
	
	@Override
	public int doWork(String[] args)
		{
		GetOpt getopt=new GetOpt();
		int c;
		while((c=getopt.getopt(args, super.getGetOptDefault()+""))!=-1)
			{
			switch(c)
				{
				default: switch(super.handleOtherOptions(c, getopt))
					{
					case EXIT_FAILURE: return -1;
					case EXIT_SUCCESS: return 0;
					default:break;
					}
				}
			}
		VariantContextWriter w=null;
		try
			{
			w=VCFUtils.createVariantContextWriter(null);
			if(getopt.getOptInd()==args.length)
				{
				info("reading from stdin");
				return fromStdin(w);
				}
			else
				{
				return fromFiles(w,getopt.getOptInd(),args);
				}
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(w);
			}
		}

	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new VcfConcat().instanceMainWithExit(args);
		}

	}
