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
import java.util.List;
import java.util.prefs.Preferences;

import javafx.beans.property.SimpleStringProperty;
import javafx.beans.property.StringProperty;
import javafx.fxml.FXMLLoader;
import javafx.scene.layout.HBox;
import javafx.stage.FileChooser;

public abstract class AbstractFileChooserPane extends HBox 
	{
    protected final FileChooser fc = new FileChooser();
    private StringProperty saveKeyProperty = new SimpleStringProperty("");
    private StringProperty filterProperty = new SimpleStringProperty("");

	
	protected AbstractFileChooserPane() {
        final FXMLLoader fxmlLoader = new FXMLLoader(getClass().getResource("FileChooserPane.fxml"));
        fxmlLoader.setRoot(this);
        fxmlLoader.setController(this);
        try {
            fxmlLoader.load();
        } catch (IOException exception) {
                throw new RuntimeException(exception);
        	}
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
    
    
    protected void setLastSaved(final File f) {
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
  		final String k = getSaveKey();
  		if(k==null || k.trim().isEmpty()) return null;
  		try {
  			final Preferences prefs = Preferences.userNodeForPackage(getClass());
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
