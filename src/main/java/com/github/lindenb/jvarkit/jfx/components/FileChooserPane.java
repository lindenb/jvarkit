package com.github.lindenb.jvarkit.jfx.components;


import java.io.File;
import java.io.IOException;
import java.util.Optional;

import javafx.beans.property.BooleanProperty;
import javafx.beans.property.ObjectProperty;
import javafx.beans.property.SimpleBooleanProperty;
import javafx.beans.property.SimpleObjectProperty;
import javafx.beans.value.ChangeListener;
import javafx.beans.value.ObservableValue;
import javafx.fxml.FXML;
import javafx.fxml.FXMLLoader;
import javafx.scene.control.Alert;
import javafx.scene.control.Alert.AlertType;
import javafx.scene.control.ButtonType;
import javafx.scene.control.TextField;
import javafx.scene.control.Tooltip;
import javafx.stage.Window;

public class FileChooserPane extends AbstractFileChooserPane {
        @FXML
        private TextField textField = null;
        private final ObjectProperty<File> selectedFile = new SimpleObjectProperty<File>();        
        private BooleanProperty forReading = new SimpleBooleanProperty(true);
        private BooleanProperty required = new SimpleBooleanProperty(false);

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
	
	        
        

        @FXML
        protected void doSelectFile() {
                if (getSelectedFile() != null)
                        {
                        fc.setInitialDirectory(getSelectedFile().getParentFile());
                        }
                else 	{
                	File prev = super.getLastSaved();
                	if(prev!=null && prev.exists() ) {
                		fc.setInitialDirectory(prev.isDirectory()?prev:prev.getParentFile());
                		}
                	}
              
                
                updateExtensionFilter();
                
                final Window win = this.getScene().getWindow();
                final File f ;
                if(isOpen()) {
                        f = this.fc.showOpenDialog(win);
                } else
                        {
                        f = this.fc.showSaveDialog(win);
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

