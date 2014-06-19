package com.github.lindenb.jvarkit.tools.jaspar;

import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.LineReader;
import htsjdk.tribble.readers.LineReaderUtil;



import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.SubSequence;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.bio.RevCompCharSequence;

public class GenomicJaspar extends AbstractCommandLineProgram
	{
	private List<Matrix> jasparDb=new ArrayList<Matrix>();
	private double fraction_of_max=0.95;
	private GenomicJaspar()
		{
		}

	
	@Override
	public String getProgramDescription()
		{
		return "Find jaspar patterns in FASTA sequences. Reports a BED file.";
		}
	
	@Override
	public String getProgramName()
		{
		return "GenomicJaspar";
		}
	@Override
	protected String getOnlineDocUrl()
		{
		return "https://github.com/lindenb/jvarkit/wiki/GenomicJaspar";
		}
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -J (uri) jaspar PFM uri. required. example: http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/nonredundant/pfm_vertebrates.txt ");
		out.println(" -f (0<ratio<1) fraction of best score. default:"+this.fraction_of_max);
		super.printOptions(out);
		}
	
	private void digest(
			String seqName,
			int position0,
			final StringBuilder sequence
			)
		{
		
		for(Matrix matrix:this.jasparDb)
			{
			if(matrix.length()>sequence.length()) continue;
			
			CharSequence forward=new SubSequence(sequence,0,matrix.length());
			CharSequence revcomp=new RevCompCharSequence(forward);
			

			
			for(int strand=0;strand<2;++strand)
				{
				double score= matrix.score(strand==0?forward:revcomp);
				if(score<=0) continue;
				
				if(score>= matrix.max()*this.fraction_of_max)
					{
					System.out.print(seqName);
					System.out.print('\t');
					System.out.print(position0);
					System.out.print('\t');
					System.out.print(position0+matrix.length());
					System.out.print('\t');
					System.out.print(matrix.getName());
					System.out.print('\t');
					System.out.print((int)(1000.0*(score/matrix.max())));
					System.out.print('\t');
					System.out.print(strand==1?'-':'+');
					System.out.print('\t');
					System.out.print(matrix.length());
					System.out.print('\t');
					System.out.print(matrix.getArchetype());
					System.out.print('\t');
					System.out.print(strand==0?forward:revcomp);
					System.out.println();
					
					break;
					}
				}
			}
		}
	
	private void run(Reader in) throws IOException
		{
		int longest=0;
		for(Matrix m:this.jasparDb)
			{
			longest=Math.max(m.length(), longest);
			}
		info("longest:"+longest);
		String seqName="";
		int position0=0;
		StringBuilder sequences=new StringBuilder(longest);
		for(;;)
			{
			int c=in.read();
			if(c==-1 || c=='>')
				{
				while(sequences.length()!=0)
					{
					digest(seqName,position0,sequences);
					++position0;
					sequences.deleteCharAt(0);
					}
				if(c==-1) break;
				StringBuilder b=new StringBuilder();
				while((c=in.read())!=-1 && c!='\n')
					{
					b.append((char)c);
					}
				seqName=b.toString();
				position0=0;
				}
			else if(!Character.isWhitespace(c))
				{
				sequences.append((char)Character.toUpperCase(c));
				if(sequences.length()==longest)
					{
					digest(seqName,position0,sequences);
					if(System.out.checkError())  return ;
					++position0;
					sequences.deleteCharAt(0);
					if(position0%1000000==0)
						{
						info(seqName+" "+position0);
						}
					}
				}
			}
		}

	
	@Override
	public int doWork(String[] args)
		{
		String jasparUri=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "J:f:"))!=-1)
			{
			switch(c)
				{
				case 'J': jasparUri=opt.getOptArg(); break;
				case 'f': fraction_of_max=Double.parseDouble(opt.getOptArg()); break;
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
		if(jasparUri==null)
			{
			error("Undefined jaspar-uri");
			return -1;
			}
		
		
		
		try
			{
			info("Reading "+jasparUri);
			LineReader lr=LineReaderUtil.fromBufferedStream(IOUtils.openURIForReading(jasparUri));
			LineIterator liter=new LineIteratorImpl(lr);
			Iterator<Matrix> miter=Matrix.iterator(liter);
			while(miter.hasNext())
				{
				Matrix matrix = miter.next();
				this.jasparDb.add(matrix.convertToPWM());
				}
			lr.close();
			info("JASPAR size: "+this.jasparDb.size());
			
			
			if(opt.getOptInd()==args.length)
				{
				info("Reading from stdin");
				run(new InputStreamReader(System.in));
				}
			else
				{
				for(int i=opt.getOptInd(); i<args.length;++i)
					{
					info("Opening "+args[i]);
					Reader in=IOUtils.openURIForBufferedReading(args[i]);
					run(in);
					in.close();
					}
					
				}
			return 0;
			}
		catch(Throwable err)
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
		new GenomicJaspar().instanceMainWithExit(args);
		}

	}
