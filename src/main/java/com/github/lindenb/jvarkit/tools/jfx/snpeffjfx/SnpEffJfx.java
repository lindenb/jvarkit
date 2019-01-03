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
import javafx.scene.control.TextField;
import javafx.stage.Stage;

public class SnpEffJfx extends AbstractSnpEfffJfxApplication {
	@FXML private FileChooserPane inputvcf;
	@FXML private FileChooserPane configFile;
	@FXML private FileChooserPane dataDir;
	@FXML private FileChooserPane outputvcf;
	@FXML private CheckBox chrPrefix;
	@FXML private CheckBox download;
	@FXML private FileChooserPane htmlSummary;
	@FXML private TextField genome_version;

	
	public SnpEffJfx()
		{
		
		}
	
	@Override
	public void start(Stage stage) throws Exception {
		final Scene scene = new Scene(fxmlLoad("SnpEffJfx.fxml"));
		stage.setScene(scene);
		super.start(stage);
		}
	
	private  List<String> buildArgs() throws JFXException {
		if(outputvcf.getSelectedFile()==null) throw new JFXException("file out missing");
		final  List<String> args = new ArrayList<>();
		args.add("eff");
		new OptionBuilder(configFile, "-c").fill(args);
		new OptionBuilder(dataDir, "-dataDir").fill(args);
		new OptionBuilder(chrPrefix, "-chr").fill(args);
		
		args.add(this.download.isSelected()?"-download":"-nodownload");
			
		args.add("-i");
		args.add("vcf");
		args.add("-o");
		args.add("vcf");
		if(htmlSummary.getSelectedFile()==null)
			{
			args.add("-noStats");
			}
		else
			{
			new OptionBuilder(htmlSummary, "-s").fill(args);
			}
		String genomeVersion = this.genome_version.getText().trim();
		if(genomeVersion.isEmpty())  {
			throw new JFXException("no genome version available");
		}
		args.add(genomeVersion);
		
		if(inputvcf.getSelectedFile()==null) {
			throw new JFXException("file in missing");
		}
		args.add(inputvcf.getSelectedFile().getPath());

		return args;
		}
	
	@Override
	protected Runnable createRunnable() throws JFXException {
		if(outputvcf.getSelectedFile()==null) throw new JFXException("file out missing");

		return new SnpEffRunner();
		}
		
	public static void main(String[] args) {
		launch(args);
	}

	private class SnpEffRunner extends AbstractSnpEfffJfxApplication.AbstractSnpEffRunner
		{
			public SnpEffRunner() throws JFXException {
				super("ca.mcgill.mcb.pcingola.snpEffect.commandLine.SnpEff",
						buildArgs(),
						outputvcf.getSelectedFile()
						);
				}
			
			@Override
			protected boolean isExitSuccess(Object is_ok) throws JFXException {
			return Boolean.TRUE.equals(is_ok);
			}
			
			@Override
			protected Object createInstance() throws JFXException
				{
				final Constructor<?> constructor;
				try {
					constructor = super.snpEffClass.getConstructor(super.args.getClass());
				} catch (NoSuchMethodException err) {
					throw new JFXException(err);
				} catch (SecurityException err) {
					throw new JFXException(err);
				}
				Object instance;
				
				try
					{
					instance = constructor.newInstance((Object)super.args);
					}
				catch(Exception err)
					{
					throw new JFXException(err);
					}
				return instance;
				}
		}
}
