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

import java.awt.Image;
import java.awt.image.BufferedImage;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.Collection;

import javax.swing.Icon;
import javax.swing.ImageIcon;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionGroup;
import org.apache.commons.cli.Options;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

public abstract class CommandFactory 
	{
	protected enum Status {OK,EXIT_SUCCESS,EXIT_FAILURE};
	private static final Log LOG=LogFactory.getLog(CommandFactory.class);
	
	/** error stream */
	private java.io.PrintStream _errStream = System.err;
	/** stdout stream */
	private java.io.PrintStream _outStream = System.out;	
	/** stdin stream */
	private java.io.InputStream _inStream  = System.in;	

	
	
	public CommandFactory()
		{
		
		}
	
	public String getMessageBundle(String s)
		{
		return s==null?"":s;
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

	protected static final int ICON_SIZE=64;
	private Icon _icon=null;
	protected Image getImageIcon()
		{
		BufferedImage img=new BufferedImage(ICON_SIZE, ICON_SIZE, BufferedImage.TYPE_INT_ARGB);
		
		return img;
		}

	public Icon getIcon() {
		if(_icon==null)
			{
			_icon = new ImageIcon(getImageIcon());
			}
		return _icon;
		}

	public String getName() {
		return getClass().getName();
		}
	
	public String getVersion()
		{
		return "undefined";
		}
	
	public String getLabel() {
		return getName();
		}
	
	public String getDescription() {
		return getLabel();
		}
	
	protected Option createHelpOption()
		{
		return Option.builder("h").longOpt("help").hasArg(false).desc("show help (this screen)").build();
		}
	protected Option createVersionOption()
		{
		return Option.builder("version").longOpt("version").hasArg(false).desc("show version and exit").build();
		}

	
	protected void fillOptions(final Options opts)
		{		
		final OptionGroup g = new OptionGroup();
		g.addOption(createHelpOption());
		g.addOption(createVersionOption());
		opts.addOptionGroup(g);
		}
	
	
	public abstract Command createCommand();
	
		
	protected  Status visit(final Option opt)
		{
		if(opt.getOpt().equals("version"))
			{
			stdout().println(getVersion());
			return Status.EXIT_SUCCESS;
			}
		else if(opt.getOpt().equals("h"))
			{
			usage(stdout());
			return Status.EXIT_SUCCESS;
			}
		LOG.fatal("Unknown option "+opt.getOpt());
		return Status.EXIT_FAILURE;
		}	
		
	protected void usage(final PrintStream out)
		{
		final Options opts =new Options();
		fillOptions(opts);
		final HelpFormatter fmt = new HelpFormatter();
		final PrintWriter w = new PrintWriter(out);
		fmt.setSyntaxPrefix("");
		fmt.printOptions(w, 80, opts, 4,5);
		w.flush();
		}
	
	public <X extends CommandFactory> int instanceMain(final String args[])
		{
		/* create a new  empty Command  */
		Command cmd = null;
		
		try{
			/* create new options  */
			final Options options = new Options();
			/* get all options and add it to options  */
			fillOptions(options);
			/* create a new CLI parser */
			final CommandLineParser parser = new DefaultParser();
			CommandLine cli= parser.parse(options, args);
			/* loop over the processed options */
			for(final Option opt: cli.getOptions())
				{
				final Status status = visit(opt);
				switch(status)
					{
					case EXIT_FAILURE: return -1;
					case EXIT_SUCCESS: return 0;
					default:break;
					}
				}
			cmd = createCommand();
			cmd.copyFrom(this);
			cmd.setInputFiles(cli.getArgList());
			LOG.debug("calling "+cmd);
			final Collection<Throwable> errors = cmd.call(); 
			if(!(errors==null || errors.isEmpty()))
				{
				for(Throwable err:errors)
					{
					LOG.error(err);
					}
				LOG.fatal("Command failed ");
				return -1;
				}
			return 0;
			}
		catch(org.apache.commons.cli.UnrecognizedOptionException err)
			{
			LOG.fatal("Unknown command ("+ err.getOption()+ ") , please check your input.",err);
			return -1;
			}
		catch(Throwable err)
			{
			LOG.fatal("cannot build command", err);
			return -1;
			}
		finally
			{
			if(cmd!=null) cmd.cleanup();
			}
		}

	public void instanceMainWithExit(final String args[])
		{
		int ret = instanceMain(args);
		System.exit(ret);
		}
	
	
	}
