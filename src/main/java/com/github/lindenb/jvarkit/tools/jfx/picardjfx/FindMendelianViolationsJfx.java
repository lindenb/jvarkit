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
import javafx.scene.control.TextArea;
import javafx.scene.control.TextField;
import javafx.stage.Stage;
import picard.vcf.MendelianViolations.FindMendelianViolations;

import java.util.ArrayList;
import java.util.List;

import com.github.lindenb.jvarkit.jfx.components.FileChooserPane;

import javafx.fxml.*;


public class FindMendelianViolationsJfx extends AbstractPicardJfxApplication
	{
	@FXML
	private FileChooserPane inputvcf;
	@FXML
	private FileChooserPane pedigree;
	@FXML
	private FileChooserPane outmetrics;
	@FXML
	private Spinner<Integer> min_gq;
	@FXML
	private Spinner<Integer> min_dp;
	@FXML
	private Spinner<Double> min_het_fraction;
	@FXML
	private FileChooserPane vcfdir;
	@FXML
	private TextField skipchroms;
	@FXML
	private TextField malechroms;
	@FXML
	private TextField femalechroms;
	@FXML
	private TextArea par_regions;
	
	
	public FindMendelianViolationsJfx() {
		super(FindMendelianViolations.class);
	}
	
	@Override
	public void start(Stage stage) throws Exception {
		final Scene scene = new Scene(fxmlLoad("FindMendelianViolationsJfx.fxml"));
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
		new OptionBuilder(pedigree,"PED=").fill(args);
		new OptionBuilder(outmetrics,"O=").fill(args);
		new OptionBuilder(min_gq,"MIN_GQ=").fill(args);
		new OptionBuilder(min_dp,"MIN_DP=").fill(args);
		new OptionBuilder(min_het_fraction,"MINHET=").fill(args);
		new OptionBuilder(vcfdir,"VCF_DIR=").fill(args);
		new OptionBuilder(skipchroms,"SKIP_CHROMS=").fill(args);
		new OptionBuilder(malechroms,"MALE_CHROMS=").fill(args);
		new OptionBuilder(femalechroms,"FEMALE_CHROMS=").fill(args);
		new OptionBuilder(par_regions,"PSEUDO_AUTOSOMAL_REGIONS=").fill(args);
		
		return args;
	}
	
	
	

	
	}
