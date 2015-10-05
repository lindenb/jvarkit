/**
 * 
 */
package com.github.lindenb.jvarkit.knime;


import java.io.File;
import java.util.ArrayList;
import java.util.List;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.command.CommandFactory;

/**
 * @author lindenb
 *
 */
public abstract class AbstractKnimeApplication
		extends AbstractCommandLineProgram
		implements KnimeApplication
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(AbstractKnimeApplication.class);
	private File outputfile=null;
	protected static final String INPUT_SOURCE_STDIN="<STDIN>";
	
	protected static abstract class AbstractKnimeCommand extends Command
		{
		
		}
	public void setOutputFile(String out) {
		setOutputFile(new File(out));
		}
	@Override
	public void setOutputFile(File out) {
		this.outputfile=out;
		}
	public File getOutputFile() {
		return outputfile;
		}
	
	@Override
	public int initializeKnime() {
		return 0;
		}

	@Override
	public void disposeKnime() {
		}

	@Override
	public void checkKnimeCancelled() {
		// TODO Auto-generated method stub
		
	}
	protected boolean checkOutputError()
	{
	if(getOutputFile()!=null) return false;
	return System.out.checkError();
	}
	
	protected int mainWork(int optind,String args[])
	{
	if( initializeKnime() != 0)
		{
		error("Initialization of "+getProgramName()+" failed.");
		return -1;
		}
	try
		{
		List<String> arglist=new ArrayList<>();
		while(optind<args.length)
			{
			String arg = args[optind];
			arglist.add(arg);
			++optind;
			}
		return executeKnime(arglist);
		}
	catch(Exception err)
		{
		error(err);
		return -1;
		}
	finally
		{
		disposeKnime();
		}
	}	
}
