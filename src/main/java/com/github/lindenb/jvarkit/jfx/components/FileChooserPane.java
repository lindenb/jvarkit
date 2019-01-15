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

import javafx.beans.property.BooleanProperty;
import javafx.beans.property.ObjectProperty;
import javafx.beans.property.SimpleBooleanProperty;
import javafx.beans.property.SimpleObjectProperty;
import javafx.beans.value.ChangeListener;
import javafx.beans.value.ObservableValue;
import javafx.fxml.FXML;
import javafx.fxml.FXMLLoader;
import javafx.scene.control.TextField;
import javafx.scene.control.Tooltip;
import javafx.stage.DirectoryChooser;
import javafx.stage.Window;

public class FileChooserPane extends AbstractFileChooserPane {
        @FXML
        private TextField textField = null;
        private final ObjectProperty<File> selectedFile = new SimpleObjectProperty<File>();        
        private BooleanProperty forReading = new SimpleBooleanProperty(true);
        private BooleanProperty required = new SimpleBooleanProperty(false);
        private BooleanProperty directory = new SimpleBooleanProperty(false);
        private final DirectoryChooser directoryChooser = new DirectoryChooser();
        private BooleanProperty remember = new SimpleBooleanProperty(false);

        
        
        public FileChooserPane() {
                final FXMLLoader fxmlLoader = new FXMLLoader(getClass().getResource("FileChooserPane.fxml"));
                fxmlLoader.setRoot(this);
                fxmlLoader.setController(this);
                try {
                    fxmlLoader.load();
                } catch (IOException exception) {
                        throw new RuntimeException(exception);
                	}
                this.required.addListener(new ChangeListener<Boolean>()
					{
					@Override
					public void changed(
							ObservableValue<? extends Boolean> observable,
							Boolean oldValue, Boolean newValue)
						{
						updateTextFieldColor();
						}
					});
                
                this.selectedFile.addListener(new ChangeListener<File>() {
                        @Override
                        public void changed(final ObservableValue<? extends File> observable, File oldValue, File newValue) {
                                if( newValue == null)
                                        {
                                        textField.setText("");
                                        textField.setTooltip(null);
                                        textField.selectPositionCaret(0);
                                        
                                        }
                                else
                                        {
                                        textField.setText(newValue.getPath());
                                        textField.setTooltip(new Tooltip(newValue.getPath()));
                                        textField.selectPositionCaret(0);
                                        }
                                updateTextFieldColor();
                                }
                        });
               this.remember.addListener(new ChangeListener<Boolean>()
					{
					@Override
					public void changed(
							ObservableValue<? extends Boolean> observable,
							Boolean oldValue, Boolean newValue)
						{
						if(!Boolean.TRUE.equals(newValue)) return;
						File f= getLastSaved();
						if(f==null) return;
						selectedFile.set(f);
						}
					});
                updateTextFieldColor();
        		}

    private void updateTextFieldColor()
    	{
    	this.textField.setStyle("-fx-background-color: "+
    			(getSelectedFile()==null && isRequired()? "red" : "white" )+
    			";"
    			);
    	}
        
	    public final File getSelectedFile() {  return selectedFileProperty().get();  }
	    public final void setSelectedFile(File f) { this.selectedFileProperty().set(f);  }
	    public ObjectProperty<File> selectedFileProperty() { return selectedFile ;  }

    
    
        
        
        public BooleanProperty openProperty() {return this.forReading;}
        public void setOpen(boolean b) {  this.openProperty().set(b);  }
        public boolean isOpen() { return this.openProperty().get(); }
 
        
        public BooleanProperty requiredProperty() { return this.required;}
		public void setRequired(boolean b) {  this.requiredProperty().set(b); }
		public boolean isRequired() { return this.requiredProperty().get(); }
	
        public BooleanProperty directoryProperty() { return this.directory;}
		public void setDirectory(boolean b) {  this.directoryProperty().set(b); }
		public boolean isDirectory() { return this.directoryProperty().get(); }
	        
        public BooleanProperty rememberProperty() { return this.remember;}
		public void setRemember(boolean b) {  this.rememberProperty().set(b); }
		public boolean isRemember() { return this.rememberProperty().get(); }
        

        @FXML
        protected void doSelectFile() {
        	
        		if(isDirectory())
        			{
        			 if (getSelectedFile() != null && getSelectedFile().isDirectory())
	                     {
	                     directoryChooser.setInitialDirectory(getSelectedFile());
	                     }
        			 else
        			 	{
        				File prev = super.getLastSaved();
 	                	if(prev!=null && prev.exists()) {
 	                		directoryChooser.setInitialDirectory(prev.isDirectory()?prev:prev.getParentFile());
 	                		
 	                		}
        			 	}
        			}
        		else
	        		{
	                if (getSelectedFile() != null)
	                        {
	                        fc.setInitialDirectory(
	                        		getSelectedFile().isDirectory()?
	                        		getSelectedFile():
	                        		getSelectedFile().getParentFile()
	                        		);
	                        }
	                else 
	                	{
	                	File prev = super.getLastSaved();
	                	if(prev!=null && prev.exists() ) {
	                		fc.setInitialDirectory(prev.isDirectory()?prev:prev.getParentFile());
	                		}
	                	}
	                updateExtensionFilter();
	        		}
                
                
                
                final Window win = this.getScene().getWindow();
                final File f ;
                if(isDirectory()) {
                	f= this.directoryChooser.showDialog(win);
                	}
                else  if(isOpen()) {
                        f = this.fc.showOpenDialog(win);
                } else
                        {
                        f = this.fc.showSaveDialog(win);
                        /* no, save does this 
                        if(f!=null && f.exists()) {
                                final Alert alert = new Alert(AlertType.WARNING);
                                alert.setTitle("Confirmation Dialog");
                                alert.setHeaderText("Delete file ?");
                                alert.setContentText("File \""+f.getName()+"\" exists. Overwrite ?");
                                final Optional<ButtonType> result = alert.showAndWait();
                                if (result.get() != ButtonType.OK){
                                        {
                                        return;
                                        }
                                        }
                                }
                        */
                        }
                if( f == null) return;
                setLastSaved(f);
                setSelectedFile(f);
                }

        @FXML
        protected void doClear() {
                setSelectedFile(null);
        }

}

