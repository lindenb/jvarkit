package com.github.lindenb.jvarkit.util;

import java.io.PrintStream;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Simple and Stupid implementation of C:getopt
 *
 */
public class GetOpt
	{
	public class Option
		{
		private char optopt;
		private String longopt;
		private boolean takesArgument;
		private String description;
		private Option() {}
		public char getOptopt()
			{
			return optopt;
			}
		public Option setDescription(String description)
			{
			this.description=description;
			return this;
			}
		public Option setLongOption(String longopt)
			{
			while(longopt.startsWith("-")) longopt=longopt.substring(1);
			if(!longopt.isEmpty()) this.longopt=longopt;
			return this;
			}
		public String getDescription()
			{
			return description;
			}
		@Override
		public String toString()
			{
			return "("+optopt+"):"+description;
			}
		}
	private Map<Character, Option> options=new LinkedHashMap<Character, GetOpt.Option>();
	private int optind=0;
	private int idxoptind=0;
	private Option optopt=null;
	private String optarg=null;
	public static final int EOF=-1;
	public static final int ERROR='?';
	
	public GetOpt()
		{
		
		}
	
	public void usage(PrintStream out)
		{
		for(Character c:options.keySet())
			{
			Option opt=this.options.get(c);
			String s="-"+c;
			if(opt.longopt!=null)
				{
				s+="  | --"+opt.longopt;
				}
			if(opt.takesArgument)
				{
				s+=" (arg)";
				}
			s+=" : ";
			s+=opt.description;
			out.println(s);
			}
		}
	
	public String getOptArg()
		{
		return this.optarg;
		}
	
	public int getOptInd()
		{
		return this.optind;
		}
	
	public Option getOptOpt()
		{
		return this.optopt;
		}
	
	public Option createOption(char c,boolean takeArgument)
		{
		Option opt=this.options.get(c);
		if(opt==null)
			{
			throw new IllegalStateException("Option already exists "+c);
			}
		opt=new Option();
		opt.optopt=c;
		opt.takesArgument=takeArgument;
		opt.description="Option "+c;
		this.options.put(c,opt);
		return opt;
		}
	
	public int getopt(String args[])
		{
		this.optopt=null;
		this.optarg=null;
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
			for(Option opt:this.options.values())
				{
				if(opt.longopt!=null && opt.longopt.equals(args[this.optind].substring(2)))
					{
					this.optopt=opt;
					break;
					}
				}
			if(this.optopt==null) return ERROR;
			if(!optopt.takesArgument)
				{
				this.idxoptind=1;
				this.optind++;
				return this.optopt.optopt;
				}
			else if(this.optind+1 < args.length)
				{
				this.idxoptind=1;
				this.optind++;
				this.optarg=args[optind];
				this.optind++;
				return this.optopt.optopt;
				}
			else
				{
				return ERROR;
				}
			}
		else if(args[this.optind].startsWith("-") &&
				args[this.optind].length()>1)
			{
			this.optopt = this.options.get(
					args[optind].charAt(this.idxoptind)
					);
			if(this.optopt==null)
				{
				return ERROR;
				}
			if(!optopt.takesArgument)//no argument
				{
				idxoptind++;
				if(idxoptind>=args[optind].length())
					{
					this.idxoptind=1;
					this.optind++;
					}
				return this.optopt.optopt;
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
				return this.optopt.optopt;
				}
			// argument the current token
			else if(
					idxoptind+1<args[optind].length()
					)
				{
				this.optarg=args[this.optind].substring(idxoptind+1);
				this.idxoptind=1;
				this.optind++;
				return this.optopt.optopt;
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
