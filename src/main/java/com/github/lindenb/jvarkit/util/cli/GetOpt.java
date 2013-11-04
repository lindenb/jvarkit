package com.github.lindenb.jvarkit.util.cli;


/**
 * Simple and Stupid implementation of C:getopt
 *
 */
public class GetOpt
	{	

	private int optind=0;
	private int idxoptind=1;
	private char optopt='\0';
	private String optarg=null;
	public static final int EOF=-1;
	public static final int ERROR='?';


	public GetOpt()
		{
		
		}
	

	
	public String getOptArg()
		{
		return this.optarg;
		}
	
	public int getOptInd()
		{
		return this.optind;
		}
	
	public int getOptOpt()
		{
		return this.optopt;
		}
	
	
	public int getopt(String args[],String pattern)
		{
		this.optopt='\0';
		this.optarg=null;
		if(this.optind>=args.length) return EOF;
		
		
		/*"--" no more arguments */
		if(args[this.optind].equals("--"))
			{
			++optind;
			return EOF;
			}
		/* long arguments, not handled */
		else if(args[this.optind].startsWith("--"))
			{
			return ERROR;
			}
		else if(args[this.optind].startsWith("-") &&
				args[this.optind].length()>1)
			{
			this.optopt=args[this.optind].charAt(idxoptind);
			if(optopt==':') return ERROR;
			int argpos=pattern.indexOf(this.optopt);
			if(argpos==-1) return -1;
			boolean hasArg=(argpos+1< pattern.length() && pattern.charAt(argpos+1)==':');
			if(argpos==-1)
				{
				return ERROR;
				}
			if(!hasArg)//no argument
				{
				idxoptind++;
				if(idxoptind>=args[optind].length())
					{
					this.idxoptind=1;
					this.optind++;
					}
				return this.optopt;
				}
			// argument is next token
			else if(
					idxoptind+1==args[optind].length() &&
					optind+1 < args.length
					)
				{
				this.optind++;					
				this.optarg=args[this.optind];
				this.idxoptind=1;
				this.optind++;
				return this.optopt;
				}
			// argument the current token
			else if(
					idxoptind+1<args[optind].length()
					)
				{
				this.optarg=args[this.optind].substring(idxoptind+1);
				this.idxoptind=1;
				this.optind++;
				return this.optopt;
				}
			else
				{
				return ERROR;
				}
				
			}
		else
			{
			return EOF;
			}
		}
	}
