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

import java.util.List;

import com.github.lindenb.jvarkit.jfx.components.FileChooserPane;
import com.github.lindenb.jvarkit.jfx.components.FilesChooserPane;

import javafx.fxml.FXML;
import javafx.scene.Scene;
import javafx.scene.control.CheckBox;
import javafx.scene.control.ComboBox;
import javafx.scene.control.TextField;
import javafx.stage.Stage;

public class DepthOfCoverageJfx extends AbstractGatkJfxApplication {
	
	
	@FXML private FilesChooserPane inputbams;
	@FXML private FileChooserPane outputdepthofcoverage;
	@FXML private CheckBox ignoreDeletionSites;
	@FXML private CheckBox includeDeletions;
	@FXML private CheckBox includeRefNSites;
	@FXML private TextField maxMappingQuality;
	@FXML private TextField minMappingQuality;
	@FXML private TextField maxBaseQuality;
	@FXML private TextField minBaseQuality;
	@FXML private TextField summaryCoverageThreshold;
	@FXML private CheckBox omitDepthOutputAtEachBase;
	@FXML private CheckBox omitIntervalStatistics;
	@FXML private CheckBox omitLocusTable;
	@FXML private CheckBox omitPerSampleStats;
	@FXML private ComboBox<?> countType;
	@FXML private ComboBox<?> partitionType;

	public DepthOfCoverageJfx() {
	}
	
	@Override
	public void start(Stage stage) throws Exception {
		final Scene scene = new Scene(fxmlLoad("DepthOfCoverageJfx.fxml"));
        stage.setScene(scene);
        super.start(stage);
    	}
	
	@Override
	protected String getAnalysisType() {
		return "DepthOfCoverage";
		}
	protected  List<String> buildArgs() throws JFXException {
		final List<String> args= super.buildArgs();

		new OptionBuilder(this.inputbams,"-I").fill(args);		
		new OptionBuilder(this.outputdepthofcoverage,"--out").fill(args);	
		
		new OptionBuilder(this.countType,"--countType").fill(args);		
		new OptionBuilder(this.ignoreDeletionSites,"--ignoreDeletionSites").fill(args);		
		new OptionBuilder(this.includeDeletions,"--includeDeletions").fill(args);		
		new OptionBuilder(this.includeRefNSites,"--includeRefNSites").fill(args);		
		new OptionBuilder(this.maxBaseQuality,"--maxBaseQuality").fill(args);		

		new OptionBuilder(this.maxMappingQuality,"--maxMappingQuality").fill(args);		

		new OptionBuilder(this.minBaseQuality,"--minBaseQuality").fill(args);		

		new OptionBuilder(this.minMappingQuality,"--minMappingQuality").fill(args);		


		new OptionBuilder(this.omitDepthOutputAtEachBase,"--omitDepthOutputAtEachBase").fill(args);		
		new OptionBuilder(this.omitIntervalStatistics,"--omitIntervalStatistics").fill(args);		
		new OptionBuilder(this.omitLocusTable,"--omitLocusTable").fill(args);		
		new OptionBuilder(this.omitPerSampleStats,"--omitPerSampleStats").fill(args);		
		new OptionBuilder(this.partitionType,"--partitionType").fill(args);		
		new OptionBuilder(this.summaryCoverageThreshold,"--summaryCoverageThreshold").split().fill(args);		
		return args;
		}
	
	public static void main(String[] args)
		{
		launch(args);
		}

}
