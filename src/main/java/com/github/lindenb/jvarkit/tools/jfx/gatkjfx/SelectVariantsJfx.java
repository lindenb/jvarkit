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
import javafx.scene.control.Spinner;
import javafx.scene.control.TextArea;
import javafx.scene.control.TextField;
import javafx.stage.Stage;

import java.util.List;

import com.github.lindenb.jvarkit.jfx.components.FileChooserPane;

import javafx.fxml.*;


public class SelectVariantsJfx extends AbstractGatkJfxApplication
	{
	@FXML
	private FileChooserPane inputvcf;
	@FXML
	private FileChooserPane outputvcf;
	
	/* samples */
	@FXML
	private TextArea sampleNamesStr;
	@FXML
	private FileChooserPane sampleNamesFile;
	@FXML
	private TextArea sampleExpr;
	@FXML
	private CheckBox inverseSamples;
	
	@FXML private FileChooserPane concordance;
	@FXML private FileChooserPane discordance;
 	@FXML private CheckBox forceValidOutput;
 	
	@FXML private FileChooserPane pedigreeFile;
 	@FXML private CheckBox noPedigreeValidation;
 	@FXML private CheckBox mendelianViolation;
	@FXML private CheckBox invertMendelianViolation;
	@FXML private TextField mendelianViolationQualThreshold;

	@FXML private FileChooserPane excludeID;
	@FXML private FileChooserPane keepIDs;

	
	
 	@FXML private CheckBox excludeFiltered;
 	@FXML private CheckBox excludeNonVariants;
 	@FXML private CheckBox keepOriginalAC;
 	@FXML private CheckBox keepOriginalDP;
	@FXML private CheckBox setFilteredGtToNocall;
	@FXML private CheckBox preserveAlleles;
	
 	@FXML private TextArea selectexpressions;
 	@FXML private CheckBox invertselect;
	@FXML private TextField selectTypeToInclude;
	@FXML private TextField selectTypeToExclude;
	@FXML private TextField select_random_fraction;
	@FXML private TextField remove_fraction_genotypes;
	@FXML private Spinner<Integer> minIndelSize;
	
	@FXML private TextField minFilteredGenotypes;
	@FXML private TextField maxFilteredGenotypes;
	@FXML private TextField minFractionFilteredGenotypes;
	@FXML private TextField maxFractionFilteredGenotypes;
	@FXML private TextField maxNOCALLfraction;
	@FXML private TextField maxNOCALLnumber;
	
	public SelectVariantsJfx() {
	}
	
	@Override
	protected String getAnalysisType() {
		return "SelectVariants";
		}
	
	@Override
	public void start(Stage stage) throws Exception {
		final Scene scene = new Scene(fxmlLoad("SelectVariantsJfx.fxml"));
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
		final String prefix= (this.inverseSamples.isSelected()?"--exclude_sample":"--sample");
		
		new OptionBuilder(sampleNamesFile,prefix+"file").fill(args);
		new OptionBuilder(sampleNamesStr,prefix+"name").fill(args);
		new OptionBuilder(sampleExpr,prefix+"expressions").fill(args);
		
		new OptionBuilder(concordance,"--concordance").fill(args);
		new OptionBuilder(discordance,"--discordance").fill(args);
		new OptionBuilder(this.forceValidOutput,"--forceValidOutput").fill(args);		
		
		if(this.pedigreeFile.getSelectedFile()!=null)
			{
			new OptionBuilder(this.pedigreeFile,"--pedigree").fill(args);		
	
			args.add("--pedigreeValidationType");
			args.add(noPedigreeValidation.isSelected()?"SILENT":"STRICT");
			
			
			
			new OptionBuilder(this.mendelianViolation,"--mendelianViolation").fill(args);		
			new OptionBuilder(this.invertMendelianViolation,"--invertMendelianViolation").fill(args);		
			new OptionBuilder(this.mendelianViolationQualThreshold,"--mendelianViolationQualThreshold").itemClass(Double.class).fill(args);		
			}
		
		new OptionBuilder(this.excludeID,"--excludeID").fill(args);		
		new OptionBuilder(this.keepIDs,"--keepIDs").fill(args);		

		new OptionBuilder(this.selectexpressions,"--selectexpressions").fill(args);		
		new OptionBuilder(this.invertselect,"--invertselect").fill(args);		
		new OptionBuilder(this.selectTypeToExclude,"--selectTypeToExclude").split().fill(args);		
		new OptionBuilder(this.selectTypeToInclude,"--selectTypeToInclude").split().fill(args);		
		new OptionBuilder(this.select_random_fraction,"--select_random_fraction").itemClass(Double.class).fill(args);		
		new OptionBuilder(this.remove_fraction_genotypes,"--remove_fraction_genotypes").itemClass(Double.class).fill(args);		
		
		
		
		new OptionBuilder(this.excludeFiltered,"--excludeFiltered").fill(args);		
		new OptionBuilder(this.excludeNonVariants,"--excludeNonVariants").fill(args);		
		new OptionBuilder(this.keepOriginalAC,"--keepOriginalAC").fill(args);		
		new OptionBuilder(this.keepOriginalDP,"--keepOriginalDP").fill(args);		
		new OptionBuilder(this.setFilteredGtToNocall,"--setFilteredGtToNocall").fill(args);		
		new OptionBuilder(this.preserveAlleles,"--preserveAlleles").fill(args);		
		new OptionBuilder(this.minIndelSize,"--minIndelSize").fill(args);		
		
		
		
		
		new OptionBuilder(this.minFractionFilteredGenotypes,"--minFractionFilteredGenotypes").itemClass(Double.class).fill(args);		
		new OptionBuilder(this.minFilteredGenotypes,"--minFilteredGenotypes").itemClass(Integer.class).fill(args);		
		new OptionBuilder(this.maxFractionFilteredGenotypes,"--maxFractionFilteredGenotypes").itemClass(Double.class).fill(args);		
		new OptionBuilder(this.maxFilteredGenotypes,"--maxFilteredGenotypes").itemClass(Integer.class).fill(args);		
		
		new OptionBuilder(this.maxNOCALLfraction,"--maxNOCALLfraction").itemClass(Double.class).fill(args);		
		new OptionBuilder(this.maxNOCALLnumber,"--maxNOCALLnumber").itemClass(Integer.class).fill(args);		
		
		
		return args;
	}
	
	
	

	
	}
