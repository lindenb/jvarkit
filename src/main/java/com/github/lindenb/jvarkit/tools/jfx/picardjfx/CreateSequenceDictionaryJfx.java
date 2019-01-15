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
import javafx.scene.control.CheckBox;
import javafx.scene.control.TextField;
import javafx.stage.Stage;
import picard.sam.CreateSequenceDictionary;

import java.util.ArrayList;
import java.util.List;

import com.github.lindenb.jvarkit.jfx.components.FileChooserPane;

import javafx.fxml.*;


public class CreateSequenceDictionaryJfx extends AbstractPicardJfxApplication
	{
	@FXML
	private FileChooserPane inputref;
	@FXML
	private FileChooserPane outputdict;
	
	@FXML private TextField genomeAssembly;
	@FXML private TextField uri;
	@FXML private TextField species;
	@FXML private CheckBox truncate;
	@FXML private TextField maxseq;
	
	public CreateSequenceDictionaryJfx() {
		super(CreateSequenceDictionary.class);
	}
	
	@Override
	public void start(Stage stage) throws Exception {
		final Scene scene = new Scene(fxmlLoad("CreateSequenceDictionaryJfx.fxml"));
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
		new OptionBuilder(inputref,"R=").fill(args);
		new OptionBuilder(outputdict,"O=").fill(args);
		new OptionBuilder(genomeAssembly,"GENOME_ASSEMBLY=").fill(args);
		new OptionBuilder(uri,"URI=").fill(args);
		new OptionBuilder(species,"SPECIES=").fill(args);
		new OptionBuilder(truncate,"TRUNCATE_NAMES_AT_WHITESPACE=").fill(args);
		new OptionBuilder(maxseq,"NUM_SEQUENCES=").itemClass(Integer.class).fill(args);
		return args;
	}
	
	
	

	
	}
