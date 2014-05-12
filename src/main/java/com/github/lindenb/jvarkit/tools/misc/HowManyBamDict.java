package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.math.BigInteger;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.LinkedHashSet;
import java.util.Set;
import java.util.logging.Level;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMFileReader.ValidationStringency;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.cli.GetOpt;

public class HowManyBamDict extends AbstractCommandLineProgram {
	private MessageDigest md5;
	private HowManyBamDict()
		{
		 
		}
	
	
	private class Dict
		{
		SAMSequenceDictionary ssd;
		String hash;
		File representative;
		Dict(SAMSequenceDictionary ssd,File representative)
			{
			this.ssd=ssd;
			this.representative=representative;
		    
	    	
	    	 md5.reset();
	    	 md5.update(String.valueOf(ssd.size()).getBytes());
	    	 for(SAMSequenceRecord ssr:ssd.getSequences())
	    	 	{
	    		 md5.update(ssr.getSequenceName().getBytes());
	    		 md5.update(String.valueOf(ssr.getSequenceLength()).getBytes());
	    	 	}
	         
	         hash = new BigInteger(1, md5.digest()).toString(16);
	         if (hash.length() != 32) {
	             final String zeros = "00000000000000000000000000000000";
	             hash= zeros.substring(0, 32 - hash.length()) + hash;
	         }
	       

			
			}
		
		@Override
		public boolean equals(Object obj)
			{
			if(this==obj) return true;
			if(obj==null) return false;
			return this.ssd.equals(Dict.class.cast(obj).ssd);
			}
		
		@Override
		public int hashCode() {
			return ssd.hashCode();
			}
		void print()
			{
			System.out.print("DICT");
			System.out.print("\t");
			System.out.print(this.hash);
			System.out.print("\t");
			System.out.print(ssd.size());
			System.out.print("\t");
			System.out.print(ssd.getReferenceLength());
			System.out.print("\t");
			boolean first=true;
			for(SAMSequenceRecord ssr:ssd.getSequences())
				{
				if(!first) System.out.print(";");
				first=false;
				System.out.print(ssr.getSequenceName());
				System.out.print('=');
				System.out.print(ssr.getSequenceLength());
				}
			System.out.print("\t");
			System.out.print(this.representative);
			System.out.println();
			}
		}

	private Dict empty=null;
	private Set<Dict>  allditcs=new LinkedHashSet<Dict>();
	
	@Override
	public String getProgramDescription() {
		return "finds if there's are some differences in the sequence dictionaries.";
		}
	
 	@Override
	public void printOptions(PrintStream out)
 		{
		out.println(" -h get help (this screen)");
		out.println(" -v print version and exit.");
		out.println(" -L (level) log level. One of "+Level.class.getName()+" currently:"+getLogger().getLevel());
		}
 	
 	private void handle(File f) throws IOException
 		{
 		SAMFileReader sfr=null;
 		try {
 			info(f);
			sfr=new SAMFileReader(f);
			sfr.setValidationStringency(ValidationStringency.SILENT);
			SAMFileHeader header=sfr.getFileHeader();
			if(header==null || header.getSequenceDictionary()==null)
				{
				if(this.empty==null)
					{
					this.empty=new Dict(new SAMSequenceDictionary(),f);
					allditcs.add(this.empty);
					this.empty.print();
					}
				System.out.print("BAM\t");
				System.out.print(f.getPath());
				System.out.print("\t");
				System.out.print(this.empty.hash);
				System.out.println();
				}
			else
				{
				Dict d=new Dict(header.getSequenceDictionary(), f);
				if(this.allditcs.add(d))
					{
					d.print();
					}
				System.out.print("BAM\t");
				System.out.print(f.getPath());
				System.out.print("\t");
				System.out.print(d.hash);
				System.out.println();
				}
 			} 
 		catch (Exception e)
			{
			error(e, e.getMessage());
			throw new IOException(e);
			}
 		finally
 			{
 			CloserUtil.close(sfr);
 			}
 		}
	
	@Override
	public int doWork(String[] args)
		{
		GetOpt getopt=new GetOpt();
		int c;
		while((c=getopt.getopt(args, "hvL:"))!=-1)
			{
			switch(c)
				{
				case 'h': printUsage();return 0;
				case 'v': System.out.println(getVersion());return 0;
				case 'L': getLogger().setLevel(Level.parse(getopt.getOptArg()));break;
				case ':': System.err.println("Missing argument for option -"+getopt.getOptOpt());return -1;
				default: System.err.println("Unknown option -"+getopt.getOptOpt());return -1;
				}
			}
		
		 try {
   		  this.md5 = MessageDigest.getInstance("MD5");
         } catch (NoSuchAlgorithmException e) {
            error(e,"MD5 algorithm not found");
            return -1;
         }
		
		try
			{
			if(getopt.getOptInd()==args.length)
				{
				info("Reading from stdin");
				String line;
				
					BufferedReader in=new BufferedReader(new InputStreamReader(System.in));
					while((line=in.readLine())!=null)
						{
						if(line.isEmpty() || line.endsWith(File.separator) || line.startsWith("#")) continue;
						handle(new File(line));
						}
					in.close();
					
				}
			else
				{
				for(int i=getopt.getOptInd();i< args.length;++i)
					{
					handle(new File(args[i]));
					}
				}
			}
		catch(IOException err)
			{
			error(err);
			return -1;
			}
		
		return 0;
		}

	
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new HowManyBamDict().instanceMainWithExit(args);
		}

}
