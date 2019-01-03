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
package com.github.lindenb.jvarkit.jfx.components;


import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import javafx.beans.property.IntegerProperty;
import javafx.beans.property.SimpleIntegerProperty;
import javafx.beans.value.ChangeListener;
import javafx.beans.value.ObservableValue;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.fxml.FXML;
import javafx.fxml.FXMLLoader;
import javafx.scene.control.Button;
import javafx.scene.control.ListView;
import javafx.scene.control.SelectionMode;

public class FilesChooserPane extends AbstractFileChooserPane {
	@FXML
	private ListView<File> files =null;
	
	@FXML
	private Button eraseButton=null;

	
    private  IntegerProperty _minCardinality = new SimpleIntegerProperty(0);
    private  IntegerProperty _maxCardinality = new SimpleIntegerProperty(-1);


	public FilesChooserPane()	{
		final FXMLLoader fxmlLoader = new FXMLLoader(getClass().getResource("FilesChooserPane.fxml"));
		fxmlLoader.setRoot(this);
		fxmlLoader.setController(this);
		try {
		    fxmlLoader.load();
		} catch (IOException exception) {
		    throw new RuntimeException(exception);
			}
		this.files.getSelectionModel().setSelectionMode(SelectionMode.MULTIPLE);
		this.eraseButton.setDisable(true);
        this.files.getSelectionModel().selectedItemProperty().addListener( new ChangeListener<File>() {
                @Override
                public void changed(ObservableValue<? extends File> observable, File oldValue, File newValue) {
                eraseButton.setDisable(files.getSelectionModel().isEmpty());
                }
        	});
        eraseButton.setOnAction(new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent event) {
				doClear(event);
			}
		});
        }
	
    public IntegerProperty minCardinalityProperty() {return this._minCardinality;}
    public void setMinCardinality(int b) {  this.minCardinalityProperty().set(b);  }
    public int getMinCardinality() { return this.minCardinalityProperty().get(); }

    public IntegerProperty maxCardinalityProperty() {return this._maxCardinality;}
    public void setMaxCardinality(int b) {  this.maxCardinalityProperty().set(b);  }
    public int getMaxCardinality() { return this.maxCardinalityProperty().get(); }
	
	
	/** get the selected files. Never null */
	public List<File> getSelectedFiles(){
		final List<File> f=this.files.getItems();
		if(f==null || f.isEmpty()) return Collections.emptyList();
		return new ArrayList<>(f);
	}
	
	@FXML
	protected void doSelectFile(ActionEvent event) {  
            if(!files.getItems().isEmpty()) {
            	fc.setInitialDirectory( files.getItems().get(files.getItems().size()-1).getParentFile()
                    );
            	}
            else
            	{
            	File prev = this.getLastSaved();
            	if(prev!=null && prev.exists() ) {
            		fc.setInitialDirectory(prev.isDirectory()?prev:prev.getParentFile());
            		}
            	}
            updateExtensionFilter();        
            List<File> fs=fc.showOpenMultipleDialog(this.getScene().getWindow());
            if(fs==null || fs.isEmpty()) return;
            for(File f:fs){
                    if(files.getItems().contains(f)) continue;
                    files.getItems().add(f);
            	}
            setLastSaved(fs.get(fs.size()-1));
			}
	
	@FXML
	protected void doClear(ActionEvent event) {  
            final List<Integer> indices=new ArrayList<>(files.getSelectionModel().getSelectedIndices());
            for(int i=indices.size()-1;i>=0;--i){
            	files.getItems().remove(indices.get(i).intValue());
            	}
            eraseButton.setDisable(files.getSelectionModel().isEmpty());
			}

	}

