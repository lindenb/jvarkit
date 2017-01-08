package com.github.lindenb.jvarkit.jfx.components;


import java.io.File;
import java.io.IOException;
import java.util.List;

import javafx.beans.value.ChangeListener;
import javafx.beans.value.ObservableValue;
import javafx.collections.ObservableList;
import javafx.fxml.FXML;
import javafx.fxml.FXMLLoader;
import javafx.scene.control.Button;
import javafx.scene.control.ListView;

public class FilesChooserPane extends AbstractFileChooserPane {
	@FXML
	final ListView<File> files =null;
	
	@FXML
	private Button eraseButton=null;


	public FilesChooserPane()	{
		final FXMLLoader fxmlLoader = new FXMLLoader(getClass().getResource("FilesChooserPane.fxml"));
		fxmlLoader.setRoot(this);
		fxmlLoader.setController(this);
		try {
		    fxmlLoader.load();
		} catch (IOException exception) {
		    throw new RuntimeException(exception);
			}
		eraseButton.setDisable(true);
        this.files.getSelectionModel().selectedItemProperty().addListener( new ChangeListener<File>() {
                @Override
                public void changed(ObservableValue<? extends File> observable, File oldValue, File newValue) {
                eraseButton.setDisable(files.getSelectionModel().isEmpty());
                }
        });
        }
	
	@FXML
	protected void doSelect() {  
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
	protected void doRemove() {  
            ObservableList<Integer> indices=files.getSelectionModel().getSelectedIndices();
            for(int i=indices.size()-1;i>=0;--i){
            files.getItems().remove(indices.get(i));
            }
            eraseButton.setDisable(files.getSelectionModel().isEmpty());
			}

		}

