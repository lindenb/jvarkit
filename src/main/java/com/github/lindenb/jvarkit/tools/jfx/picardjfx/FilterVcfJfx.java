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
package com.github.lindenb.jvarkit.tools.jfx.picardjfx;

import javafx.scene.Scene;
import javafx.scene.control.Spinner;
import javafx.stage.Stage;
import picard.vcf.filter.FilterVcf;

import java.util.ArrayList;
import java.util.List;

import com.github.lindenb.jvarkit.jfx.components.FileChooserPane;

import javafx.fxml.*;


public class FilterVcfJfx extends AbstractPicardJfxApplication
	{
	@FXML
	private FileChooserPane inputvcf;
	@FXML
	private FileChooserPane outputvcf;
	@FXML
	private Spinner<Double> min_ab;
	@FXML
	private Spinner<Integer> min_dp;
	@FXML
	private Spinner<Integer> min_gq;
	@FXML
	private Spinner<Double> max_fs;
	@FXML
	private Spinner<Double> min_qd;
	@FXML
	private FileChooserPane javascript;
	
	public FilterVcfJfx() {
		super(FilterVcf.class);
	}
	
	@Override
	public void start(Stage stage) throws Exception {
		final Scene scene = new Scene(fxmlLoad("FilterVcfJfx.fxml"));
        stage.setScene(scene);
        super.start(stage);
    	}
	
	
	public static void main(String[] args)
		{
		launch(args);
		}
	
	@Override
	protected  List<String> buildArgs() throws JFXException {
		final List<String> args= new ArrayList<>();
		new OptionBuilder(inputvcf,"I=").fill(args);
		new OptionBuilder(outputvcf,"O=").fill(args);
		new OptionBuilder(min_ab,"MIN_AB=").fill(args);
		new OptionBuilder(min_dp,"MIN_DP=").fill(args);
		new OptionBuilder(min_gq,"MIN_GQ=").fill(args);
		new OptionBuilder(max_fs,"MAX_FS=").fill(args);
		new OptionBuilder(min_qd,"MIN_QD=").fill(args);
		new OptionBuilder(javascript,"JS=").fill(args);
		return args;
	}
	
	
	

	
	}
