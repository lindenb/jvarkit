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
package com.github.lindenb.jvarkit.tools.jfx.snpeffjfx;

import java.lang.reflect.Constructor;
import java.util.ArrayList;
import java.util.List;

import com.github.lindenb.jvarkit.jfx.components.FileChooserPane;

import javafx.fxml.FXML;
import javafx.scene.Scene;
import javafx.scene.control.CheckBox;
import javafx.scene.control.TextArea;
import javafx.scene.control.TextField;
import javafx.stage.Stage;

public class SnpSiftFilterJfx extends AbstractSnpEfffJfxApplication {
	@FXML private FileChooserPane inputvcf;
	@FXML private FileChooserPane expressionFile;
	@FXML private FileChooserPane outputvcf;

	@FXML private TextArea expression;
	@FXML private TextField addFilter;
	@FXML private TextField filterId;
	@FXML private CheckBox inverse;
	@FXML private CheckBox pass;
	@FXML private TextField rmFilter;
	@FXML private CheckBox errMissing;

	
	public SnpSiftFilterJfx()
		{
		
		}
	
	@Override
	public void start(Stage stage) throws Exception {
		final Scene scene = new Scene(fxmlLoad("SnpSiftFilterJfx.fxml"));
		stage.setScene(scene);
		super.start(stage);
		}
	
	private  List<String> buildArgs() throws JFXException {
		if(outputvcf.getSelectedFile()==null) throw new JFXException("file out missing");
		final  List<String> args = new ArrayList<>();
		args.add("filter");
		new OptionBuilder(addFilter, "--addFilter").fill(args);
		new OptionBuilder(filterId, "--filterId").fill(args);
		new OptionBuilder(inverse, "--inverse").fill(args);
		new OptionBuilder(pass, "--pass").fill(args);
		new OptionBuilder(rmFilter, "--rmFilter").fill(args);
		new OptionBuilder(errMissing, "--errMissing").fill(args);
		new OptionBuilder(inputvcf, "--file").fill(args);

		final String expressionstr = this.expression.getText().trim();
		if(expressionstr.isEmpty() && expressionFile.getSelectedFile()==null)  {
			throw new JFXException("expression is empty and no script selected");
			}
		else if(!expressionstr.isEmpty() && expressionFile.getSelectedFile()!=null)  {
			throw new JFXException("expression string and file both declared.");
			}
		else if(expressionFile.getSelectedFile()==null)
			{
			args.add(expressionstr);
			}
		else
			{
			new OptionBuilder(expressionFile, "--exprFile").fill(args);
			}	
		
		


		return args;
		}
	
	@Override
	protected Runnable createRunnable() throws JFXException {
		if(outputvcf.getSelectedFile()==null) throw new JFXException("file out missing");

		return new SnpSiftEffRunner();
		}
		
	public static void main(String[] args) {
		launch(args);
	}
	
	
	private class SnpSiftEffRunner extends AbstractSnpEfffJfxApplication.AbstractSnpEffRunner
	{
		public SnpSiftEffRunner() throws JFXException {
			super("ca.mcgill.mcb.pcingola.snpSift.SnpSift",
					buildArgs(),
					outputvcf.getSelectedFile()
					);
			}
		
		@Override
		protected boolean isExitSuccess(Object is_ok) throws JFXException {
		return is_ok==null;
		}

		
		@Override
		protected Object createInstance() throws JFXException
			{
			final Constructor<?> constructor;
			try {
				constructor = super.snpEffClass.getConstructor(super.args.getClass(),String.class);
			} catch (NoSuchMethodException err) {
				throw new JFXException(err);
			} catch (SecurityException err) {
				throw new JFXException(err);
			}
			Object instance;
			
			try
				{
				instance = constructor.newInstance((Object)super.args,String.class.cast(null));
				}
			catch(Exception err)
				{
				throw new JFXException(err);
				}
			return instance;
			}
	}


}
