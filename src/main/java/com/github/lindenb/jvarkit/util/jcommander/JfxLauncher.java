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
package com.github.lindenb.jvarkit.util.jcommander;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Supplier;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.log.Logger;

import javafx.application.Application;
import javafx.application.Platform;
import javafx.event.ActionEvent;
import javafx.scene.control.Alert;
import javafx.scene.control.Alert.AlertType;
import javafx.scene.control.ScrollPane;
import javafx.scene.control.TextArea;
import javafx.scene.layout.BorderPane;
import javafx.stage.Stage;

public abstract class JfxLauncher
	extends Application
	{
	private static final Logger LOG = Logger.build(JfxLauncher.class).make();
	private static int exitStatus = 0;

	@ParametersDelegate
	private Launcher.UsageBuider usageBuider = null;
	@Parameter(description = "Files")
	private List<String> args = new ArrayList<>();
	@Parameter(names="--testng",description = "testng",hidden=true)
	private boolean this_is_a_unit_test = false;

	/** when testing, give a chance to set the parameters */
	private Supplier<List<String>> argsargvsupplier =()-> 
		 this.getParameters().getUnnamed();
		;
	
	protected JCommander jCommander = null;
	
	protected JfxLauncher() {
		JfxLauncher.exitStatus = 0;
	}

	 public void init() throws Exception {
	 super.init();
		this.usageBuider = new Launcher.UsageBuider(getClass());
		this.jCommander = new JCommander(this);
	 }
	
	protected void doMenuAbout(final ActionEvent event) {
		final Alert alert = new Alert(AlertType.INFORMATION);
		alert.setHeaderText("About...");
		alert.setContentText("Pierre Lindenbaum PhD. Institut du Thorax. "
				+ "Nantes. France.");
		alert.showAndWait();
	}

	protected void doMenuQuit(final ActionEvent event) {
		// http://stackoverflow.com/questions/12153622
		Platform.exit();
	}

	public static void setExitStatus(int exitStatus) {
		JfxLauncher.exitStatus = exitStatus;
	}

	protected boolean isUnitText() {
		return this.this_is_a_unit_test;
	}
	
	public static int getExitStatus() {
		return exitStatus;
	}


	public void setArcArgvSupplier(final Supplier<List<String>> args) {
		this.argsargvsupplier = args;
		}
	

	protected void displayAlert(final Throwable err) {
		final Alert alert = new Alert(AlertType.ERROR);
		alert.setHeaderText("Cannot create Command.");
		alert.setContentText(String.valueOf(err.getMessage()));

		// Create expandable Exception.
		final StringWriter sw = new StringWriter();
		final PrintWriter pw = new PrintWriter(sw);
		err.printStackTrace(pw);

		final TextArea textArea = new TextArea(sw.toString());
		textArea.setEditable(false);
		textArea.setWrapText(true);

		final BorderPane pane = new BorderPane(new ScrollPane(textArea));
		alert.getDialogPane().setExpandableContent(pane);

		alert.showAndWait();
		}

	/** return 0 on successful initialization */
	protected abstract int doWork(final Stage primaryStage,final List<String> args);
	
	@Override
	public void start(final Stage primaryStage) throws Exception {
		// reset exit status
		JfxLauncher.exitStatus = 0;
		final List<String> params = this.argsargvsupplier.get();
		this.jCommander.parse(params.toArray(new String[params.size()]));
		if(this.usageBuider.shouldPrintUsage())
			{
			this.usageBuider.usage(jCommander);
			Platform.exit();
			return;
			}
		int ret  = 0;
		try
			{
			ret = doWork(primaryStage,this.args);
			}
		catch(final Throwable err)
			{
			ret= -1;
			LOG.error(err);
			}
		if(ret!=0) {
			setExitStatus(ret);
			Platform.exit();
			}
		}
	}
