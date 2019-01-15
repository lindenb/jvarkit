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
import javafx.scene.control.CheckBox;
import javafx.scene.control.TextField;
import javafx.stage.Stage;
import java.util.List;

import com.github.lindenb.jvarkit.jfx.components.FileChooserPane;

import javafx.fxml.*;


public class VariantsToTableJfx extends AbstractGatkJfxApplication
	{
	@FXML
	private FileChooserPane inputvcf;
	@FXML
	private FileChooserPane outputtable;
	
	
 	@FXML private TextField fields;
 	@FXML private TextField genotypeFields;
	
 	@FXML private CheckBox allowMissingData;

 	@FXML private TextField maxRecords;

 	@FXML private CheckBox moltenize;

 	@FXML private CheckBox showFiltered;

 	@FXML private CheckBox splitMultiAllelic;
 	
 	
	public VariantsToTableJfx() {
	}
	
	@Override
	protected String getAnalysisType() {
		return "VariantsToTable";
		}
	
	@Override
	public void start(Stage stage) throws Exception {
		final Scene scene = new Scene(fxmlLoad("VariantsToTableJfx.fxml"));
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
		new OptionBuilder(inputvcf,"--variant").fill(args);
		new OptionBuilder(outputtable,"--out").fill(args);

		new OptionBuilder(this.allowMissingData,"--allowMissingData").fill(args);		
		

			new OptionBuilder(this.fields,"--fields").split().fill(args);		

			new OptionBuilder(this.genotypeFields,"--genotypeFields").split().fill(args);		

			new OptionBuilder(this.maxRecords,"--maxRecords").itemClass(Integer.class).fill(args);		

			new OptionBuilder(this.moltenize,"--moltenize").fill(args);		

			new OptionBuilder(this.showFiltered,"--showFiltered").fill(args);		

			new OptionBuilder(this.splitMultiAllelic,"--splitMultiAllelic").fill(args);		

		return args;
	}
	
	
	

	
	}
