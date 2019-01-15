/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.jfx.gatkjfx;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.github.lindenb.jvarkit.jfx.components.FileChooserPane;
import com.github.lindenb.jvarkit.tools.jfx.AbstractJfxApplication;
import javafx.application.Platform;
import javafx.fxml.FXML;
import javafx.scene.control.Alert;
import javafx.scene.control.TextArea;
import javafx.scene.control.TextField;
import javafx.scene.control.Alert.AlertType;
import javafx.scene.control.CheckBox;
import javafx.stage.Stage;

import org.apache.log4j.AppenderSkeleton;
import org.apache.log4j.Logger;
import org.apache.log4j.PatternLayout;
import org.apache.log4j.spi.LoggingEvent;
import org.apache.log4j.spi.ThrowableInformation;
import org.broadinstitute.gatk.utils.commandline.CommandLineUtils;


public abstract class AbstractGatkJfxApplication
	extends AbstractJfxApplication
	{
	private Runnable runningThread=null; 
	protected static final Logger LOG = CommandLineUtils.getStingLogger();

	@FXML
	private FileChooserPane referenceFile;
	/* intervals */
	@FXML
	private FileChooserPane intervalFile;
	@FXML
	private TextField intervalStr;
	@FXML
	private CheckBox inverseInterval;

	
	
	protected AbstractGatkJfxApplication()
		{
		}
	
	/** get option -T */
	protected abstract String getAnalysisType();
	/** get window name */
	protected String getName() {
		return getClass().getSimpleName();
		}
	protected List<String> buildArgs() throws JFXException
		{
		if(referenceFile.getSelectedFile()==null) {
			throw new JFXException("Reference file was not set");
			}
		final ArrayList<String> args = new ArrayList<>();
		args.add("-T");
		args.add(getAnalysisType());
		args.add("-R");
		args.add(referenceFile.getSelectedFile().getPath());
		
		final String intervalOpt = (this.inverseInterval.isSelected()?"--excludeIntervals":"--intervals");
		new OptionBuilder(intervalFile,intervalOpt).fill(args);
		new OptionBuilder(intervalStr,intervalOpt).split().fill(args);

		
		return args;
		}

	
	@Override
	protected Runnable createRunnable() throws JFXException
		{
		final List<String> args = buildArgs();
		realStderr.println("ARGS="+args);
		if(args==null) return null;
		this.runningThread = new GATKRunner(this,args);
		return this.runningThread;
		}

	@Override
	public void start(Stage stage) throws Exception {
		stage.setTitle(getName());
		
		
		super.start(stage);
		
		
		//final MyAppender appender = new MyAppender(super.console);
		//Logger.getRootLogger().addAppender(appender);
		}
	
	private static class MyAppender extends AppenderSkeleton
		{
		private TextArea output;
		private MyAppender(final TextArea output)
			{
			this.output=output;
			setLayout(new PatternLayout(PatternLayout.TTCC_CONVERSION_PATTERN));
			}
		@Override
		protected void append(final LoggingEvent e) {
			try {
				Platform.runLater(new Runnable()
					{
					@Override
					public void run() {
					String logString = getLayout().format(e);
					if(e.getThrowableInformation()!=null)
						{	
						ThrowableInformation ti = e.getThrowableInformation();
						if(ti.getThrowable()!=null)
							{
							logString+=""+ti.getThrowable().getMessage()+"\n";
							}
						}	
	
					if(!output.isVisible())
						{
						//System.err.println(logString);
						return;
						}
					output.appendText(logString);
					if(output.getLength()>50000)
						{
						output.deleteText(0, output.getLength()-50000);
						}
					}
				});
			} catch (Exception e2) {
				e2.printStackTrace(AbstractJfxApplication.realStderr);
				}
			}
	
		@Override
		public void close() {
			
		}
		
		@Override
		public boolean requiresLayout() {
			return true;
			}
		}
	
	private static class GATKRunner implements Runnable
		{
		private AbstractGatkJfxApplication owner;
		private final String args[];
		public GATKRunner(final AbstractGatkJfxApplication ui,final List<String> args)
			{
			this.owner=ui;
			this.args=args.toArray(new String[args.size()]);
			}
		@Override
		public void run()
			{
			LOG.info("starting "+Arrays.toString(args));
			AbstractGatkJfxApplication.realStderr.println(args);
			org.broadinstitute.gatk.engine.CommandLineGATK instance= new org.broadinstitute.gatk.engine.CommandLineGATK();
	
			try
				{
				org.broadinstitute.gatk.engine.CommandLineGATK.start(instance, this.args);
				
				if(org.broadinstitute.gatk.engine.CommandLineGATK.result == 0)
					{
					try {
						Platform.runLater(new Runnable() {
							@Override
							public void run() {
								//if(GATKRunner.this != owner.runningThread) return;
								final Alert alert = new Alert(AlertType.CONFIRMATION);
								alert.setContentText("Completed:"+Arrays.toString(GATKRunner.this.args));
								
								alert.showAndWait();
								owner.runningThread=null;
							}
						});
					} catch (Exception e) {
						LOG.warn(e);
						e.printStackTrace(AbstractGatkJfxApplication.realStderr);
						}
					}
				else
					{
					try {
						Platform.runLater(new Runnable() {
							@Override
							public void run() {
								//if(GATKRunner.this != owner.runningThread) return;
								final Alert alert = new Alert(AlertType.ERROR);
								alert.setContentText("Failure:"+Arrays.toString(GATKRunner.this.args));

								owner.runningThread=null;
							}
						});
					} catch (Exception e) {
						LOG.warn(e);
						e.printStackTrace(AbstractGatkJfxApplication.realStderr);
						}
					}
				}
			catch(final Exception err)
				{
				try {
					Platform.runLater(new Runnable() {
						@Override
						public void run() {
							//if(GATKRunner.this != owner.runningThread) return;
							final Alert alert = new Alert(AlertType.ERROR);
							alert.setContentText("Failure:"+Arrays.toString(GATKRunner.this.args));

							owner.runningThread=null;
						}
					});
				} catch (Exception e) {
					LOG.warn(e);
					e.printStackTrace(AbstractGatkJfxApplication.realStderr);
					}
				}
			}
		}
	protected static class GatkResourceOptionBuilder
		{
		private boolean name_is_required=true;
		private final TextField name;
		private final FileChooserPane fileChooser;
		private final String opt;
		public GatkResourceOptionBuilder(final TextField name,final FileChooserPane fileChooser,String opt) throws JFXException
			{
			this.name=name;
			this.fileChooser=fileChooser;
			this.opt=opt;
			if(this.name==null) throw new JFXException("name is null for option "+this.opt);
			if(this.fileChooser==null) throw new JFXException("file is null for option "+this.opt);
			}
		public void fill(final List<String> args)  throws JFXException {
			final String prefix=name.getText().trim();
			final File f=this.fileChooser.getSelectedFile();
			if(f==null && prefix.isEmpty()) return;
			if(!prefix.isEmpty() && f==null) {
				 throw new JFXException("prefix is defined but no file for option "+this.opt);
				}
			else if(prefix.isEmpty() && name_is_required && f!=null)
				{
				 throw new JFXException("prefix is not defined but file defined for option "+this.opt);
				}
			args.add(this.opt+(prefix.isEmpty()?"":":"+prefix));
			args.add(f.getPath());
			}
		}
	
	}
