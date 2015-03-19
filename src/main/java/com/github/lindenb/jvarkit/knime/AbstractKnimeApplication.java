/**
 * 
 */
package com.github.lindenb.jvarkit.knime;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;

/**
 * @author lindenb
 *
 */
public abstract class AbstractKnimeApplication extends AbstractCommandLineProgram
		implements KnimeApplication
	{
	/** output File or null=stdout */
	private File outputFile=null;
	/** default name for stdin as input source */
	protected static final String INPUT_SOURCE_STDIN="<STDIN>";

	@Override
	public void setOutputFile(File out) {
		this.outputFile=out;
		}

	public void setOutputFile(String out)
		{
		setOutputFile(new File(out));
		}

	
	public File getOutputFile() {
		return outputFile;
		}
	
	
	/** 
	 * tell wether output was interrupted (pipeline)
	 * return false of getOutputFile() return not null 
	 *  or else, it returns the value of System.out.checkError();
	 */
	protected boolean checkOutputError()
		{
		if(getOutputFile()!=null) return false;
		return System.out.checkError();
		}
	
	/* (non-Javadoc)
	 * @see com.github.lindenb.jvarkit.knime.KnimeApplication#initializeKnime()
	 */
	@Override
	public int initializeKnime() {
		return 0;
	}

	/* (non-Javadoc)
	 * @see com.github.lindenb.jvarkit.knime.KnimeApplication#executeKnime(java.util.List)
	 */
	@Override
	public abstract int executeKnime(List<String> args);

	/* (non-Javadoc)
	 * @see com.github.lindenb.jvarkit.knime.KnimeApplication#disposeKnime()
	 */
	@Override
	public void disposeKnime() {

	}

	/* (non-Javadoc)
	 * @see com.github.lindenb.jvarkit.knime.KnimeApplication#checkKnimeCancelled()
	 */
	@Override
	public void checkKnimeCancelled() {
	}

	/** called by 'main(String args[])'.
	 * Calls initializeKnime,executeKnime,disposeKnime
	 * @return
	 */
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
