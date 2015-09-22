/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
package com.github.lindenb.jvarkit.util.command;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import javax.swing.Icon;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.PatternOptionBuilder;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

public abstract class AbstractCommandFactory implements CommandFactory
	{
	protected enum Status {OK,EXIT_SUCCESS,EXIT_FAILURE};
	private static final Log LOG=LogFactory.getLog(AbstractCommandFactory.class);
	
	/** error stream */
	private java.io.PrintStream _errStream = System.err;
	/** stdout stream */
	private java.io.PrintStream _outStream = System.out;	
	/** stdin stream */
	private java.io.InputStream _inStream  = System.in;	

	
	public AbstractCommandFactory()
		{
		
		}
	
	public java.io.PrintStream stderr()
		{
		return this._errStream;
		}
	
	public void stderr(final java.io.PrintStream stream)
		{
		this._errStream = stream;
		}
	
	public java.io.PrintStream stdout()
		{
		return this._outStream;
		}
	
	public void stdout(final java.io.PrintStream stream)
		{
		this._outStream = stream;
		}
		
	
	public java.io.InputStream stdin()
		{
		return this._inStream;
		}
	
	public void stdin(final java.io.InputStream stream)
		{
		this._inStream = stream;
		}	

	
	@Override
	public Icon getIcon() {
		return null;
		}
	@Override
	public String getLabel() {
		return getName();
		}
	@Override
	public String getDescription() {
		return getLabel();
		}
	
	public final List<Option> getOptionList()
		{
		List<Option> L= new ArrayList<>();
		fillOptions(L);
		return L;
		}
	
	protected void fillOptions(final List<Option> opts)
		{
		PatternOptionBuilder.parsePattern("");
		}
	
	/** creates a new Command, fill the options ready to execute */
	protected Status parseArgs(final Command cmd,final CommandLine cli)
		{
		
		return Status.OK;
		}
		
	protected Status visit(final Command cmd,final Option opt)
		{
		LOG.fatal("Unknown option "+opt.getOpt());
		return Status.EXIT_FAILURE;
		}	
		
	
	
	@Override
	public int instanceMain(final String args[])
		{
		/* create a new  empty Command  */
		final Command cmd = createCommand();
		
		try{
			/* create new options  */
			final Options options = new Options();
			/* get all options and add it to options  */
			for(Option opt:getOptionList())
				{
				options.addOption(opt);
				}
			/* create a new CLI parser */
			final CommandLineParser parser = new DefaultParser();
			CommandLine cli= parser.parse(options, args);
			/* loop over the processed options */
			for(final Option opt: cli.getOptions())
				{
				final Status status = visit(cmd,opt);
				switch(status)
					{
					case EXIT_FAILURE: return -1;
					case EXIT_SUCCESS: return 0;
					default:break;
					}
				}
			cmd.setInput(cli.getArgList());
			Collection<Throwable> errors = cmd.call(); 
			if(!(errors==null || errors.isEmpty()))
				{
				return -1;
				}

			
			return 0;
			}
		catch(Throwable err)
			{
			LOG.fatal("cannot build command", err);
			return -1;
			}
		}
	@Override
	public void instanceMainWithExit(final String args[])
		{
		int ret = instanceMain(args);
		System.exit(ret);
		}
	
	
	}
