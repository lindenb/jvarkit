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
import javafx.scene.control.TextArea;
import javafx.stage.Stage;
import java.util.List;

import com.github.lindenb.jvarkit.jfx.components.FileChooserPane;

import javafx.fxml.*;


public class GenotypeConcordanceJfx extends AbstractGatkJfxApplication
	{
	@FXML
	private FileChooserPane eval;
	@FXML
	private FileChooserPane comp;
	@FXML
	private FileChooserPane output;
	
	
 	@FXML private CheckBox ignoreFilters;
 	@FXML private CheckBox moltenize;
 	@FXML private FileChooserPane printInterestingSites;
 	@FXML private TextArea genotypeFilterExpressionComp;
 	@FXML private TextArea genotypeFilterExpressionEval;

 	
	public GenotypeConcordanceJfx() {
	}
	
	@Override
	protected String getAnalysisType() {
		return "GenotypeConcordance";
		}
	
	@Override
	public void start(Stage stage) throws Exception {
		final Scene scene = new Scene(fxmlLoad("GenotypeConcordanceJfx.fxml"));
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
		new OptionBuilder(eval,"--eval").fill(args);
		new OptionBuilder(comp,"-comp").fill(args);
		new OptionBuilder(output,"--out").fill(args);
		
		
		

		new OptionBuilder(this.genotypeFilterExpressionComp,"--genotypeFilterExpressionComp").fill(args);		
		new OptionBuilder(this.genotypeFilterExpressionEval,"--genotypeFilterExpressionEval").fill(args);		
		new OptionBuilder(this.ignoreFilters,"--ignoreFilters").fill(args);		
		new OptionBuilder(this.moltenize,"--moltenize").fill(args);		
		new OptionBuilder(this.printInterestingSites,"--printInterestingSites").fill(args);		

		return args;
	}

	
	}
