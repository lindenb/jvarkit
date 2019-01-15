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


public class SelectHeadersJfx extends AbstractGatkJfxApplication
	{
	@FXML
	private FileChooserPane inputvcf;
	@FXML
	private FileChooserPane outputvcf;
	
	
	@FXML private TextArea exclude_header_name;
	@FXML private TextArea header_expression;
	@FXML private TextArea header_name;
	@FXML private CheckBox include_interval_names;
 	
	public SelectHeadersJfx() {
	}
	
	@Override
	protected String getAnalysisType() {
		return "SelectHeaders";
		}
	
	@Override
	public void start(Stage stage) throws Exception {
		final Scene scene = new Scene(fxmlLoad("SelectHeadersJfx.fxml"));
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
		new OptionBuilder(outputvcf,"-o").fill(args);
		
		
			new OptionBuilder(this.exclude_header_name,"--exclude_header_name").fill(args);		

			new OptionBuilder(this.header_expression,"--header_expression").fill(args);		

			new OptionBuilder(this.header_name,"--header_name").fill(args);		

			new OptionBuilder(this.include_interval_names,"--include_interval_names").fill(args);		



		return args;
	}

	
	}
