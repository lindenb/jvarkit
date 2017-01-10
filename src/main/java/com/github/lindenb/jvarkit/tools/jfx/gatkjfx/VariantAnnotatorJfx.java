package com.github.lindenb.jvarkit.tools.jfx.gatkjfx;

import java.util.List;

import com.github.lindenb.jvarkit.jfx.components.FileChooserPane;

import javafx.fxml.FXML;
import javafx.scene.Scene;
import javafx.scene.control.CheckBox;
import javafx.scene.control.TextArea;
import javafx.scene.control.TextField;
import javafx.stage.Stage;

public class VariantAnnotatorJfx extends AbstractGatkJfxApplication {
	@FXML
	private FileChooserPane inputvcf;
	@FXML
	private FileChooserPane outputvcf;
	
	@FXML private FileChooserPane dbsnp;
 	@FXML private CheckBox alwaysAppendDbsnpId;
 	@FXML private CheckBox resourceAlleleConcordance;
 	
	@FXML private FileChooserPane pedigreeFile;
 	@FXML private CheckBox noPedigreeValidation;
 	@FXML private TextField mendelianViolationQualThreshold;


	
	@FXML private TextField resource_name1;
	@FXML private FileChooserPane resource_file1;
	@FXML private TextField resource_name2;
	@FXML private FileChooserPane resource_file2;
	@FXML private TextField resource_name3;
	@FXML private FileChooserPane resource_file3;
	@FXML private TextArea expression;
	
	@FXML private TextField comp_name1;
	@FXML private FileChooserPane comp_file1;
	@FXML private TextField comp_name2;
	@FXML private FileChooserPane comp_file2;
	@FXML private TextField comp_name3;
	@FXML private FileChooserPane comp_file3;
	
 	@FXML private CheckBox homopolymerRun;
 	@FXML private CheckBox gcContent;

	
	
	public VariantAnnotatorJfx() {
	}
	
	
	@Override
	public void start(Stage stage) throws Exception {
		final Scene scene = new Scene(fxmlLoad("VariantAnnotatorJfx.fxml"));
        stage.setScene(scene);
        super.start(stage);
    	}

	
	@Override
	protected String getAnalysisType() {
		return "VariantAnnotator";
	}

	public static void main(String[] args) {
		launch(args);
	}

	
	@Override
	protected  List<String> buildArgs() throws JFXException {
		final List<String> args= super.buildArgs();
		new OptionBuilder(inputvcf,"--variant").fill(args);
		new OptionBuilder(outputvcf,"-o").fill(args);
		new OptionBuilder(dbsnp,"--dbsnp").fill(args);
		new OptionBuilder(this.alwaysAppendDbsnpId,"--alwaysAppendDbsnpId").fill(args);
		new OptionBuilder(this.resourceAlleleConcordance,"--resourceAlleleConcordance").fill(args);		

		
		if(this.pedigreeFile.getSelectedFile()!=null)
			{
			new OptionBuilder(this.pedigreeFile,"--pedigree").fill(args);		
	
			args.add("--pedigreeValidationType");
			args.add(noPedigreeValidation.isSelected()?"SILENT":"STRICT");
			
			
			new OptionBuilder(this.mendelianViolationQualThreshold,"--MendelViolationGenotypeQualityThreshold").itemClass(Double.class).fill(args);		
			}
		
		
		new GatkResourceOptionBuilder(resource_name1,resource_file1,"--resource").fill(args);
		new GatkResourceOptionBuilder(resource_name2,resource_file2,"--resource").fill(args);
		new GatkResourceOptionBuilder(resource_name3,resource_file3,"--resource").fill(args);
		new OptionBuilder(this.expression,"--expression");
		
		new GatkResourceOptionBuilder(comp_name1,comp_file1,"--comp").fill(args);
		new GatkResourceOptionBuilder(comp_name2,comp_file2,"--comp").fill(args);
		new GatkResourceOptionBuilder(comp_name3,comp_file3,"--comp").fill(args);
		
		if( homopolymerRun.isSelected()) {
			args.add("-A");
			args.add("HomopolymerRun");
		}
		if( gcContent.isSelected()) {
			args.add("-A");
			args.add("GCContent");
		}
		
		return args;
		}
	
}
