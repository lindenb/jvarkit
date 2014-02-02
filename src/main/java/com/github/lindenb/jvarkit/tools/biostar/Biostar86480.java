package com.github.lindenb.jvarkit.tools.biostar;

import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.Reader;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;



import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.Rebase;

public class Biostar86480 extends AbstractCommandLineProgram
	{
	private Rebase rebase=Rebase.createDefaultRebase();
	
	private Biostar86480()
		{
		}

	
	@Override
	public String getProgramDescription()
		{
		return "Genomic restriction finder. See also http://www.biostars.org/p/86480/";
		}
	
	@Override
	public String getProgramName()
		{
		return "Biostar86480";
		}
	@Override
	protected String getOnlineDocUrl()
		{
		return "https://github.com/lindenb/jvarkit/wiki/Biostar86480";
		}
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -E (name) restrict to that enzyme. Can be called multiple times. Optional.");
		super.printOptions(out);
		}
	
	private void digest(
			String seqName,
			int position0,
			final List<Character> sequence
			)
		{
		for(Rebase.Enzyme enzyme:this.rebase)
			{
			if(enzyme.size()>sequence.size()) continue;
			for(int strand=0;strand<2;++strand)
				{
				int x=0;
				for(x=0;x< enzyme.size();++x)
					{
					char c=(strand==0?
							enzyme.at(x):
							AcidNucleics.complement(enzyme.at((enzyme.size()-1)-x))
							);
					if(!Rebase.compatible(sequence.get(x),c)) break;
					}
				if(x==enzyme.size())
					{
					System.out.print(seqName);
					System.out.print('\t');
					System.out.print(position0);
					System.out.print('\t');
					System.out.print(position0+enzyme.size());
					System.out.print('\t');
					for(int y=0;y< enzyme.size();++y)
						{
						System.out.print(sequence.get(y));
						}
					System.out.print('\t');
					System.out.print(1000);
					System.out.print('\t');
					System.out.print(strand==1?'-':'+');
					System.out.print('\t');
					System.out.print(enzyme.getName());
					System.out.print('\t');
					System.out.print(enzyme.getDecl());
					System.out.println();
					break;
					}
				if(enzyme.isPalindromic()) break;
				}
			}
		}
	
	private void run(Reader in) throws IOException
		{
		int longest=0;
		for(Rebase.Enzyme E:this.rebase)
			{
			longest=Math.max(E.size(), longest);
			}
		String seqName="";
		int position0=0;
		ArrayList<Character> sequences=new ArrayList<Character>(longest);
		for(;;)
			{
			int c=in.read();
			if(c==-1 || c=='>')
				{
				while(!sequences.isEmpty())
					{
					digest(seqName,position0,sequences);
					++position0;
					sequences.remove(0);
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
				sequences.add((char)Character.toUpperCase(c));
				if(sequences.size()==longest)
					{
					digest(seqName,position0,sequences);
					++position0;
					sequences.remove(0);
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
		Set<String> onlyEnz=new HashSet<String>();
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "E:"))!=-1)
			{
			switch(c)
				{
				case 'E': onlyEnz.add(opt.getOptArg()); break;
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
		if(!onlyEnz.isEmpty())
			{
			Rebase rebase2=new Rebase();
			for(String e:onlyEnz)
				{
				Rebase.Enzyme enz=this.rebase.getEnzymeByName(e);
				if(enz==null)
					{
					System.err.println("Cannot find enzyme "+enz +" in RE list.");
					System.err.println("Current list is:");
					for(Rebase.Enzyme E: this.rebase)
						{
						System.err.println("\t"+E);
						}
					return -1;
					}
				rebase2.getEnzymes().add(enz);
				}
			this.rebase=rebase2;
			}
		try
			{
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
		new Biostar86480().instanceMainWithExit(args);
		}

	}
