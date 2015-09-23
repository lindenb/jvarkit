package com.github.lindenb.jvarkit.util.command;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.Callable;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;



public abstract class Command
	implements Callable<Collection<Throwable>>
	{
	private static final Log LOG=LogFactory.getLog(Command.class);

	/** error stream */
	private java.io.PrintStream _errStream = System.err;
	/** stdout stream */
	private java.io.PrintStream _outStream = System.out;	
	/** stdin stream */
	private java.io.InputStream _inStream  = System.in;	
	/** input files */
	private List<String> inputFiles = Collections.emptyList();
	/** factory */
	private CommandFactory factory;
	/** log */
	private Log _log = LOG;
	
	public Command()
		{
		}
	
	public void copyFrom(final CommandFactory factory)
		{
		if(this.factory!=null) LOG.warn("copyFrom called twice ?");
		this.factory = factory;
		this.stdin(factory.stdin());
		this.stderr(factory.stderr());
		this.stdout(factory.stdout());
		}
	
	/** returns owner factory */
	public CommandFactory getFactory()
		{	
		return this.factory;
		}
	
	public void setInputFiles(final List<String> inputFiles) {
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
		
	public void setFactory(CommandFactory factory) {
		this.factory = factory;
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
		return getFactory().getName();
		}
	public String getLabel()
		{
		return getFactory().getLabel();
		}
	public String getDescription()
		{
		return getFactory().getLabel();
		}
	
	public void setLog(final Log l)
		{
		this._log = (l==null?LOG:l);
		}
	
	public Log getLog()
		{
		return _log;
		}
	protected Collection<Throwable> wrapException(final String err)
		{
		return wrapException(new RuntimeException(err));
		}
	protected Collection<Throwable> wrapException(final Throwable err)
		{
		return new ArrayList<>(Collections.singletonList(err));
		}
	@Deprecated
	public void info(final Object o) { getLog().info(o);}
	@Deprecated
	public void info(final Object o,final  Throwable err) { getLog().info(o,err);}
	@Deprecated
	public void warning(final Object o) { getLog().warn(o);}
	@Deprecated
	public void warning(final Object o,final  Throwable err) { getLog().warn(o,err);}
	@Deprecated
	public void error(final Object o) { getLog().error(o);}
	@Deprecated
	public void error(final Object o,final  Throwable err) { getLog().error(o,err);}
	public String getMessageBundle(String s) { return getFactory().getMessageBundle(s);}
	public boolean checkError() { return stdout().checkError();}
	public String getProgramCommandLine() {return "";}
	public String getVersion() { return getFactory().getVersion();}
	@Override
	public String toString() {
		return getName();
		}
	public void cleanup() {}
	}
