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
	private String longopt=null;
	public static final int EOF=-1;
	public static final int ERROR='?';
	public static final int MISSING_ARG=':';
	public static final int LONG_OPT=1;


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
	
	public String increaseOptind(String args[])
		{
		String s= args[optind++];
		this.idxoptind=1;
		this.optarg=null;
		return s;
		}
	
	public void setOptInd(int optind)
		{
		if(optind<=this.optind)
			{
			throw new IllegalStateException("Cannot set optind "+optind+" <="+this.optind);
			}
		this.optind = optind;
		this.idxoptind=1;
		this.optarg=null;
		}
	
	public char getOptOpt()
		{
		return this.optopt;
		}
	
	public String getLongOpt()
		{
		if(this.longopt==null)
			{
			throw new IllegalStateException("currently no a long opt.");
			}
		return this.longopt;
		}
	
	public int getopt(String args[],String pattern)
		{
		this.optopt='\0';
		this.optarg=null;
		this.longopt=null;
		if(this.optind>=args.length) return EOF;
		
		
		/*"--" no more arguments */
		if(args[this.optind].equals("--"))
			{
			++optind;
			return EOF;
			}
		/* long arguments */
		else if(args[this.optind].startsWith("--"))
			{
			this.longopt=args[this.optind].substring(2);
			this.optind++;
			this.idxoptind=1;
			this.optopt='\0';
			return LONG_OPT;
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
				return MISSING_ARG;
				}
				
			}
		else
			{
			return EOF;
			}
		}
	}
