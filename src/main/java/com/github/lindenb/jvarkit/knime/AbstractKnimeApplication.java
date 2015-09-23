/**
 * 
 */
package com.github.lindenb.jvarkit.knime;


import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.command.CommandFactory;

/**
 * @author lindenb
 *
 */
public abstract class AbstractKnimeApplication
		extends CommandFactory
		//implements KnimeApplication
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(AbstractKnimeApplication.class);

	protected static abstract class AbstractKnimeCommand extends Command
		{
		
		}

	

}
