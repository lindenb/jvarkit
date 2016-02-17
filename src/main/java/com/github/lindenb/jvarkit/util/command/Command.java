/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.util.command;

import java.awt.Image;
import java.awt.image.BufferedImage;
import java.io.InputStream;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.net.InetAddress;
import java.net.URL;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Date;
import java.util.Enumeration;
import java.util.List;
import java.util.Locale;
import java.util.ResourceBundle;
import java.util.concurrent.Callable;
import java.util.jar.Manifest;

import javax.swing.Icon;
import javax.swing.ImageIcon;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionGroup;
import org.apache.commons.cli.Options;
import org.slf4j.Logger;

import com.github.lindenb.jvarkit.util.htsjdk.HtsjdkVersion;



public abstract class Command
	implements Callable<Collection<Throwable>>
	{
	protected enum Status {OK,EXIT_SUCCESS,EXIT_FAILURE};
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(Command.class);

	protected static final Collection<Throwable> RETURN_OK = Collections.emptyList();
	/** error stream */
	private java.io.PrintStream _errStream = System.err;
	/** stdout stream */
	private java.io.PrintStream _outStream = System.out;	
	/** stdin stream */
	private java.io.InputStream _inStream  = System.in;
	/** program command line as string */
	private String programCommandLine="";
	/** git hash in the manifest */
	private String gitHash = null;
	/** compile date in the manifest */
	private String compileDate = null;
	/** apache command line */
	private CommandLine _cliCommand = null;
	private List<String> inputFiles = Collections.emptyList();
	/** messages */
	private ResourceBundle messagesBundle=null;
	protected static final int ICON_SIZE=64;
	private Icon _icon=null;
	
	public Command()
		{
		try {
			this.messagesBundle=ResourceBundle.getBundle("messages");
			} 
		catch (Exception e)
			{
			LOG.warn("Cannot get messages bundle ",e);
			}
		/** set locale http://seqanswers.com/forums/showthread.php?p=174020#post174020 */
		Locale.setDefault(Locale.US);
		System.setProperty("file.encoding", "UTF-8");
		}
	
	public Collection<Throwable> initializeKnime()	
		{
		return RETURN_OK;
		}
	
	public void disposeKnime()	
		{
		
		}
	
	public CommandLine getCommandLine()
		{
		return _cliCommand;
		}
	
	/** used by swing */
	public void setInputFiles(List<String> inputFiles)
		{
		this.inputFiles = inputFiles;
		}
	
	public List<String> getInputFiles() {
		return inputFiles;
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
	
	public String getName()
		{
		return getClass().getName();
		}
	public String getLabel()
		{
		return getName();
		}
	public String getDescription()
		{
		return getLabel();
		}
	private void loadManifest()
		{
		try
			{
			Enumeration<URL> resources = getClass().getClassLoader()
					  .getResources("META-INF/MANIFEST.MF");//not '/META-INF'
			while (resources.hasMoreElements())
				{
				URL url=resources.nextElement();
				InputStream in=url.openStream();
				if(in==null)
					{
					continue;
					}
				
				Manifest m=new Manifest(in);
				in.close();
				in=null;
				java.util.jar.Attributes attrs=m.getMainAttributes();
				if(attrs==null)
					{
					continue;
					}
				String s =attrs.getValue("Git-Hash");
				if(s!=null && !s.isEmpty() && !s.contains("$")) //ant failed
					{
					this.gitHash=s;
					}
				
				s =attrs.getValue("Compile-Date");
				if(s!=null && !s.isEmpty()) //ant failed
					{
					this.compileDate=s;
					}
				}
			}	
		catch(Exception err)
			{
			
			}
		
		}
	
	public String getCompileDate()
		{
		if(this.compileDate==null)
			{
			this.compileDate="undefined";
			loadManifest();
			}
		return compileDate;
		}

	public String getGitHash()
		{
		if(this.gitHash==null)
			{
			this.gitHash="1.0";
			loadManifest();
			}
		return this.gitHash;
		}

	public String getVersion()
		{
		return getGitHash();
		}
	
	public Logger getLog()
		{
		return LOG;
		}
	
	protected Collection<Throwable> wrapException(final String err)
		{
		return wrapException(new RuntimeException(err));
		}
	protected Collection<Throwable> wrapException(final Throwable err)
		{
		err.printStackTrace(stderr());
		return new ArrayList<>(Collections.singletonList(err));
		}
	@Deprecated
	public void info(final Object o) { getLog().info(String.valueOf(o));}
	@Deprecated
	public void info(final Object o,final  Throwable err) { getLog().info(String.valueOf(o),err);}
	@Deprecated
	public void warning(final Object o) { getLog().warn(String.valueOf(o));}
	@Deprecated
	public void warning(final Object o,final  Throwable err) { getLog().warn(String.valueOf(o),err);}
	@Deprecated
	public void error(final Object o) { getLog().error(String.valueOf(o));}
	@Deprecated
	public void error(final Object o,final  Throwable err) { getLog().error(String.valueOf(o),err);}
	
	public String getMessageBundle(String key)
		{
		if(key==null) return "(null)";
		if(this.messagesBundle==null) return key;
		try
			{
			return this.messagesBundle.getString(key);
			}
		catch(Exception err)
			{
			return key;
			}
		}
	
	
	public String getProgramCommandLine() {return this.programCommandLine==null?"":programCommandLine;}
	@Override
	public String toString() {
		return getName();
		}
	
	@Deprecated
	public void cleanup() {}
	
	
	//

	
	protected String getIconAsBase64()
		{
		return null;
		}
	
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

	
	protected Option createHelpOption()
		{
		return Option.builder("h").
				longOpt("help").
				hasArg(false).
				desc("show help (this screen)").
				build()
				;
		}
	protected Option createVersionOption()
		{
		return Option.builder("version").
				longOpt("version").
				hasArg(false).
				desc("show version and exit").
				build();
		}

	
	protected void fillOptions(final Options opts)
		{		
		final OptionGroup g = new OptionGroup();
		g.addOption(createHelpOption());
		g.addOption(createVersionOption());
		opts.addOptionGroup(g);
		}
	
	
	
		
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
		
		LOG.error("Unknown option "+opt.getOpt());
		return Status.EXIT_FAILURE;
		}	
		
	protected void printOptions(final PrintStream out)
		{
		final Options opts =new Options();
		fillOptions(opts);
		final HelpFormatter fmt = new HelpFormatter();
		final PrintWriter w = new PrintWriter(out);
		fmt.setSyntaxPrefix("");
		fmt.printOptions(w, 80, opts, 4,5);
		w.flush();
		}
	
	protected void printStandardPreamble(PrintStream out)
		{
		String prg=getName();
		out.println(prg);
		for(int i=0;i< prg.length();++i) out.print('=');
		out.println();
		out.print("\nDescription: ");
		out.println(getDescription());
		out.println();
		out.println("Author         : "+getAuthorName());
		out.println("Mail           : "+getAuthorMail());
		out.println("WWW            : "+getOnlineDocUrl());
		out.println("Compilation    : "+getCompileDate());
		out.println("Git-Hash       : "+getVersion());
		out.println("Htsjdk-version : "+HtsjdkVersion.getVersion());
		out.println("Htsjdk-home    : "+HtsjdkVersion.getHome());
		}



	
	public String getOnlineDocUrl()
		{
		return "https://github.com/lindenb/jvarkit/wiki/"+getName();
		}
	
	public String getOnlineSrcUrl()
		{
		return "https://github.com/lindenb/jvarkit";
		}
	
	public String getAuthorName()
		{
		return "Pierre Lindenbaum PhD.";
		}
	public String getAuthorMail()
		{
		return "plinden"+"baum"+
				'@'+
				"yahoo"+
				'.'+
				"fr";
		}
	
	
	

	protected void printSynopsis(PrintStream out)
		{
		out.println("\tjava -cp jar1.jar:jar2.jar:jar3.jar "+getClass().getName()+" [options] (files)");
		}
	public void usage(PrintStream out)
		{
		out.println();
		printStandardPreamble(out);
		out.println();
		out.println("Usage:");
		printSynopsis(out);
		out.println();
		out.println("Options:");
		printOptions(out);
		out.println();
		}

	protected void writeHtmlDoc(final XMLStreamWriter w)
		throws XMLStreamException
		{
		w.writeStartElement("div");
		w.writeComment("no doc available for "+getName());
		w.writeEndElement();
		}
	
	public  int instanceMainWithExceptions(final String args[],final List<Throwable> putExceptionsHere)
		{
		try{
			/* cmd line */
			final StringBuilder sb=new StringBuilder();
			for(String s:args)
				{
				if(sb.length()!=0) sb.append(" ");
				sb.append(s);
				}
			this.programCommandLine = sb.toString().
					replace('\n', ' ').
					replace('\t', ' ').
					replace('\"', ' ').
					replace('\'', ' ')
					;
			
			
			/* create new options  */
			final Options options = new Options();
			/* get all options and add it to options  */
			LOG.debug("creating options");
			fillOptions(options);
			/* create a new CLI parser */
			LOG.debug("parsing options");
			final CommandLineParser parser = new DefaultParser();
			this._cliCommand = parser.parse(options, args);
			this.setInputFiles( this._cliCommand.getArgList() );
			
			/* loop over the processed options */
			for(final Option opt: this._cliCommand.getOptions())
				{
				if(opt.hasArg() && !opt.hasOptionalArg() && opt.getValue()==null)
					{
					LOG.info("WARNING OPTION ####"+opt+" with null value ??");
					}
				final Status status = visit(opt);
				switch(status)
					{
					case EXIT_FAILURE: return -1;
					case EXIT_SUCCESS: return 0;
					default:break;
					}
				}
			
			final Date startDate=new Date();
			LOG.info("Starting JOB at "+startDate+" "+getClass().getName()+
					" version="+getVersion()+" "+
					" built="+getCompileDate());
			LOG.info("Command Line args : "+(
					this.getProgramCommandLine().isEmpty()?
						"(EMPTY/NO-ARGS)":
						this.getProgramCommandLine()
						));
			String hostname="";
			try
				{
				hostname= InetAddress.getLocalHost().getHostName();
				}
			catch(Exception err)
				{
				hostname="host";
				}
			LOG.info("Executing as " +
	                System.getProperty("user.name") + "@" + hostname +
	                " on " + System.getProperty("os.name") + " " + System.getProperty("os.version") +
	                " " + System.getProperty("os.arch") + "; " + System.getProperty("java.vm.name") +
	                " " + System.getProperty("java.runtime.version") 
	                );

			LOG.debug("initialize");
			final Collection<Throwable> initErrors = this.initializeKnime(); 
			if(!(initErrors==null || initErrors.isEmpty()))
				{
				putExceptionsHere.addAll(initErrors);
				for(Throwable err:initErrors)
					{
					LOG.error(err.getMessage());
					}
				LOG.error("Command initialization failed ");
				return -1;
				}
			LOG.debug("calling");
			final Collection<Throwable> errors = this.call(); 
			if(!(errors==null || errors.isEmpty()))
				{
				putExceptionsHere.addAll(initErrors);
				for(Throwable err:errors)
					{
					LOG.error(err.getMessage());
					}
				LOG.error("Command failed ");
				return -1;
				}
			LOG.debug("success ");
			this.disposeKnime();
			final Date endDate = new Date();
			final double elapsedMinutes = (endDate.getTime() - startDate.getTime()) / (1000d * 60d);
	        final String elapsedString  = new DecimalFormat("#,##0.00").format(elapsedMinutes);
	    	LOG.info("End JOB  [" + endDate + "] " + getName() + " done. Elapsed time: " + elapsedString + " minutes.");
			return 0;
			}
		catch(org.apache.commons.cli.UnrecognizedOptionException err)
			{
			LOG.error("Unknown command ("+ err.getOption()+ ") , please check your input.",err);
			return -1;
			}
		catch(Throwable err)
			{
			LOG.error("cannot execute command", err);
			return -1;
			}
		finally
			{
			this.disposeKnime();
			}
		}
	
	public  int instanceMain(final String args[])
		{
		return instanceMainWithExceptions(args,new ArrayList<Throwable>());
		}
		
	
	public void instanceMainWithExit(final String args[])
		{
		int ret = instanceMain(args);
		System.exit(ret);
		}
	
	
	}
