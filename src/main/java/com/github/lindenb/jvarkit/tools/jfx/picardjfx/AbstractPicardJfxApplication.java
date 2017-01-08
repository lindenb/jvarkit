package com.github.lindenb.jvarkit.tools.jfx.picardjfx;

import java.util.List;

import com.github.lindenb.jvarkit.tools.jfx.AbstractJfxApplication;

import javafx.stage.Stage;
import picard.cmdline.CommandLineProgram;

public abstract class AbstractPicardJfxApplication
	extends AbstractJfxApplication
	{
	
	private final Class<? extends CommandLineProgram> commandLine;
	protected AbstractPicardJfxApplication(final Class<? extends CommandLineProgram> commandLine)
		{
		this.commandLine = commandLine;
		}
	
	protected abstract List<String> buildArgs() throws JFXException;

	
	@Override
	protected Runnable createRunnable() throws JFXException
		{
		final CommandLineProgram picardApp ;
		try
			{
			picardApp =  this.commandLine.newInstance();
			}
		catch (final Exception err)
			{
			throw new RuntimeException(err);
			}
		
		final List<String> args = buildArgs();
			
		realStderr.println("ARGS="+args);
		if(args==null) return null;
		final String array[] = args.toArray(new String[args.size()]);
		return new Runnable()
				{
				@Override
				public void run()
					{
					picardApp.instanceMain(array);
					}
				};
		}

	@Override
	public void start(Stage stage) throws Exception {
		stage.setTitle(commandLine.getSimpleName());
		super.start(stage);
		}
	}
