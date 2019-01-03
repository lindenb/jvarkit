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

import javafx.scene.Scene;
import javafx.scene.control.ComboBox;
import javafx.scene.control.Spinner;
import javafx.scene.control.TextField;
import javafx.scene.control.CheckBox;
import javafx.stage.Stage;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import com.github.lindenb.jvarkit.jfx.components.FileChooserPane;
import com.github.lindenb.jvarkit.jfx.components.FilesChooserPane;

import javafx.fxml.*;


public class CombineVariantsJfx extends AbstractGatkJfxApplication
	{
	@FXML
	private FilesChooserPane inputvcfs;
	@FXML
	private FileChooserPane outputvcf;
	@FXML
	private ComboBox<String> genotypeMergeOptions;
	@FXML
	private ComboBox<String> filteredRecordsMergeType;
 	
 	@FXML
 	private CheckBox assumeIdenticalSamples;

 	@FXML
 	private CheckBox excludeNonVariants;

 	@FXML
 	private CheckBox filteredAreUncalled;

 	@FXML
 	private CheckBox mergeInfoWithMaxAC;

 	@FXML
 	private CheckBox minimalVCF;

 	@FXML
 	private CheckBox printComplexMerges;

 	@FXML
 	private CheckBox suppressCommandLineHeader;
	
 	@FXML private Spinner<Integer> minimumN;
 	@FXML private TextField setKey;

	
	public CombineVariantsJfx() {
	}
	
	@Override
	protected String getAnalysisType() {
		return "CombineVariants";
		}
	
	@Override
	public void start(Stage stage) throws Exception {
		final Scene scene = new Scene(fxmlLoad("CombineVariantsJfx.fxml"));
        stage.setScene(scene);
        super.start(stage);
    	}
	
	
	public static void main(String[] args)
		{
		launch(args);
		}
	
	@Override
	protected  List<String> buildArgs() throws JFXException {
		final List<String> args= super.buildArgs();
		
		List<File> selectFiles= this.inputvcfs.getSelectedFiles();
		List<String> priorityList=new ArrayList<>();
		int idx=0;
		for(final File f:selectFiles)
			{
			String p="v"+(idx);
			priorityList.add(p);
			args.add("--variant:"+p);
			args.add(f.getPath());
			idx++;
			}
		
		new OptionBuilder(outputvcf,"-o").fill(args);		
		new OptionBuilder(genotypeMergeOptions,"-genotypeMergeOptions").fill(args);		
		
		if("PRIORITIZE".equals(this.genotypeMergeOptions.getValue()))
			{
			args.add("-priority");
			args.add(String.join(",", priorityList));
			
			}
		new OptionBuilder(filteredRecordsMergeType,"-filteredRecordsMergeType").fill(args);	
	 	new OptionBuilder(assumeIdenticalSamples,"--assumeIdenticalSamples").fill(args);		
		new OptionBuilder(excludeNonVariants,"--excludeNonVariants").fill(args);		
		new OptionBuilder(filteredAreUncalled,"--filteredAreUncalled").fill(args);		
		new OptionBuilder(mergeInfoWithMaxAC,"--mergeInfoWithMaxAC").fill(args);		
		new OptionBuilder(minimalVCF,"--minimalVCF").fill(args);
		new OptionBuilder(printComplexMerges,"--printComplexMerges").fill(args);		
		new OptionBuilder(suppressCommandLineHeader,"--suppressCommandLineHeader").fill(args);		
		new OptionBuilder(minimumN,"--minimumN").fill(args);		
		new OptionBuilder(setKey,"--setKey").fill(args);		
		return args;
	}
	
	
	

	
	}
