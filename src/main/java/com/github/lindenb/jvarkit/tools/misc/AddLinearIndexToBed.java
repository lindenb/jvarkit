package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import java.util.regex.Pattern;

import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.CloserUtil;

import org.broad.tribble.readers.LineIterator;
import org.broad.tribble.readers.LineIteratorImpl;
import org.broad.tribble.readers.LineReaderUtil;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryFactory;

public class AddLinearIndexToBed extends AbstractCommandLineProgram {
	private SAMSequenceDictionary dictionary=null;
	private Set<String> unmappedChromosomes=new HashSet<String>();
	private long tid2offset[]=null;
	
	private AddLinearIndexToBed()
		{
		
		}
	

	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/AddLinearIndexToBed";
		}
	
	@Override
	public String getProgramDescription() {
		return "Use a Sequence dictionary to create a linear index for a BED file. Can be used as a X-Axis for a chart.";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -R (file) Fasta reference indexed with samtools, or sequence dictionary. Required");
		super.printOptions(out);
		}
	
	@SuppressWarnings("resource")
	protected int doWork(InputStream in,PrintStream out)
			throws IOException
		{
		Pattern tab=Pattern.compile("[\t]");
		LineIterator lr=new LineIteratorImpl(LineReaderUtil.fromBufferedStream(in));
		while(lr.hasNext())
			{	
			String line=lr.next();
			if(line.isEmpty() || line.startsWith("#")) continue;
			String tokens[]=tab.split(line,3);
			if(tokens.length<2)
				{
				warning("Bad chrom/pos line:"+line);
				continue;
				}
			SAMSequenceRecord ssr=this.dictionary.getSequence(tokens[0]);
			if(ssr==null)
				{
				if(!unmappedChromosomes.contains(tokens[0]))
					{
					warning("chromosome is undefined in dictionary: "+tokens[0]);
					unmappedChromosomes.add(tokens[0]);
					}
				continue;
				}
			int pos0=Integer.parseInt(tokens[1]);
			if(pos0<0 || pos0>=ssr.getSequenceLength())
				{
				warning("position is out of range for : "+line+" length("+tokens[0]+")="+ssr.getSequenceLength());
				}
			out.print(this.tid2offset[ssr.getSequenceIndex()]+pos0);
			out.print('\t');
			out.print(line);
			out.println();
			if(out.checkError()) break;
			}
		
		return 0;
		}

	
	@Override
	public int doWork(String[] args)
		{
		File refFile=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"R:"))!=-1)
			{
			switch(c)
				{
				case 'R':
						{
						refFile=new File(opt.getOptArg());
						break;
						}
				default:
					{
					switch(handleOtherOptions(c, opt, null))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		if(refFile==null)
			{
			error("Reference file undefined");
			return -1;
			}
		try
			{
			this.dictionary=new SAMSequenceDictionaryFactory().load(refFile);
			this.tid2offset=new long[this.dictionary.size()];
			Arrays.fill(this.tid2offset, 0L);
			for(int i=1;i< this.dictionary.size();++i )
				{
				this.tid2offset[i] = this.tid2offset[i-1]+
							this.dictionary.getSequence(i-1).getSequenceLength();
				}
			
			
			if(opt.getOptInd()==args.length)
				{
				info("reading stdin");
				doWork(System.in, System.out);
				}
			else
				{
				for(int i=opt.getOptInd();i< args.length;++i)
				
					{
					info("opening "+args[i]);
					InputStream in=IOUtils.openURIForReading(args[i]);
					doWork(in, System.out);
					CloserUtil.close(in);
					}
				}
			if(!unmappedChromosomes.isEmpty())
				{
				warning("Unmapped chromosomes:"+unmappedChromosomes);
				}
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		}
	

	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new AddLinearIndexToBed().instanceMainWithExit(args);
		}
	}
