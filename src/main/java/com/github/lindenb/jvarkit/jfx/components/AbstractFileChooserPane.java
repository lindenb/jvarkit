package com.github.lindenb.jvarkit.jfx.components;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.prefs.Preferences;

import javafx.beans.property.SimpleStringProperty;
import javafx.beans.property.StringProperty;
import javafx.fxml.FXML;
import javafx.scene.layout.HBox;
import javafx.stage.FileChooser;

public abstract class AbstractFileChooserPane extends HBox 
	{
    protected final FileChooser fc = new FileChooser();
    private StringProperty saveKeyProperty = new SimpleStringProperty("");
    private StringProperty filterProperty = new SimpleStringProperty("");

	
	protected AbstractFileChooserPane() {
		}
	
    public StringProperty filterProperty() { return this.filterProperty;  }
    public void setFilter(String b) { this.filterProperty().set(b);  }
    public String getFilter() { return this.filterProperty().get(); }
    
    
    public StringProperty saveKeyProperty() { return this.saveKeyProperty;  }
    public void setSaveKey(String b) { this.saveKeyProperty().set(b);  }
    public String getSaveKey() { return this.saveKeyProperty().get(); }
   

    protected FileChooser.ExtensionFilter getExtensionFilter()
	    {
	    final String filterStr = this.getFilter();
	    if(filterStr==null || filterStr.trim().isEmpty()) return null;
	    int colon = filterStr.indexOf(":");
	    if(colon==-1) return null;
	    final String key=filterStr.substring(0, colon).trim();
	    final List<String> exts=new ArrayList<>();
	    for(final String ext:filterStr.substring(colon+1).trim().split("[,; /]+"))
	    	{
	    	if(ext.isEmpty()) continue;
	    	exts.add("*."+ext);
	    	}
	    if(exts.isEmpty()) return null;
 	    return new FileChooser.ExtensionFilter(key,exts);
		}
    
    
    protected void setLastSaved(File f) {
    	String k = getSaveKey();
		if(k==null || k.trim().isEmpty() || f==null) return;
		try {
			Preferences prefs = Preferences.userNodeForPackage(getClass());
			prefs.put(k, f.getPath());
			prefs.sync();
			}
		catch(Exception err){
			return;
			}
		
    	}
    protected File getLastSaved() {
  		String k = getSaveKey();
  		if(k==null || k.trim().isEmpty()) return null;
  		try {
			Preferences prefs = Preferences.userNodeForPackage(getClass());
			String f = prefs.get(k, null);
			if(f==null || f.isEmpty()) return null;
			return new File(f);
			}
		catch(Exception err){
			return null;
			}
    	
    	}
    
    protected void updateExtensionFilter() {
	    this.fc.getExtensionFilters().clear();
	    final FileChooser.ExtensionFilter extf = getExtensionFilter();
	    if(extf!=null)
	    	{
	    	this.fc.getExtensionFilters().add(extf);
	    	}
	    }
	}
