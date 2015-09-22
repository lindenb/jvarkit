package com.github.lindenb.jvarkit.util.command;

import java.util.Collection;
import java.util.Collections;

public abstract class AbstractCommand implements Command
	{
	private CommandFactory factory;
	/** error stream */
	protected java.io.PrintStream _errStream = System.err;
	/** stdout stream */
	protected java.io.PrintStream _outStream = System.out;	
	/** stdin stream */
	protected java.io.InputStream _inStream  = System.in;	
	
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
	public CommandFactory getFactory() {
		return this.factory;
		}
	@Override
	public String getDescription() {
		return getFactory().getDescription();
		}
	@Override
	public String getLabel() {
		return getFactory().getLabel();
		}
	
	@Override
	public String getName() {
		return getFactory().getName();
		}
	
	protected Collection<Throwable> convertMessageToErrors(final String err)
		{
		return Collections.singletonList(new Throwable(err));
		}
	
	}
