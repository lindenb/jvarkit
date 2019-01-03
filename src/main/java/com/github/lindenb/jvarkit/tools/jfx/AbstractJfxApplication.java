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
package com.github.lindenb.jvarkit.tools.jfx;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

import com.github.lindenb.jvarkit.jfx.components.FileChooserPane;
import com.github.lindenb.jvarkit.jfx.components.FilesChooserPane;

import javafx.application.Application;
import javafx.application.Platform;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.fxml.FXML;
import javafx.fxml.FXMLLoader;
import javafx.geometry.Rectangle2D;
import javafx.scene.Parent;
import javafx.scene.control.Alert;
import javafx.scene.control.Alert.AlertType;
import javafx.scene.control.Button;
import javafx.scene.control.CheckBox;
import javafx.scene.control.ComboBox;
import javafx.scene.control.ScrollPane;
import javafx.scene.control.Spinner;
import javafx.scene.control.TextArea;
import javafx.scene.control.TextField;
import javafx.scene.layout.BorderPane;
import javafx.stage.Screen;
import javafx.stage.Stage;

public abstract class AbstractJfxApplication
	extends Application
	{
	protected static final PrintStream realStderr=System.err;
	protected static final PrintStream realStdout=System.out;
	@FXML
	private Button runCommandButton;
	@FXML
	private Button cancelCommandButton;
	@FXML
	protected TextArea console;
	protected PrintStream printToConsole;
	
	protected Thread commandThread=null;
	protected AbstractJfxApplication()
		{
		}
	
	
    @FXML protected void doMenuAbout(final ActionEvent event) {
    	final Alert alert = new Alert(AlertType.INFORMATION);
    	alert.setHeaderText("About...");
    	alert.setContentText("Pierre Lindenbaum PhD. Institut du Thorax. Nantes. France.");
    	
    	alert.showAndWait();
    	}
    
    @FXML protected void doMenuQuit(final ActionEvent event) {
    	// http://stackoverflow.com/questions/12153622
    	Platform.exit();
    }
	
	@Override
	public void start(Stage stage) throws Exception {
	     
	     
	    Rectangle2D primaryScreenBounds = Screen.getPrimary().getVisualBounds();

        stage.setX(50);
        stage.setY(50);
        stage.setWidth(primaryScreenBounds.getWidth()-100);
        stage.setHeight(primaryScreenBounds.getHeight()-100);
        stage.show();
        this.runCommandButton.setOnAction(new EventHandler<ActionEvent>()
			{
			@Override
			public void handle(ActionEvent event)
				{
				doCommandStart(event);
				}
			});
        this.cancelCommandButton.setOnAction(new EventHandler<ActionEvent>()
			{
			@Override
			public void handle(ActionEvent event)
				{
				doCommandEnd(event);
				}
			});
        this.cancelCommandButton.setDisable(true);
		
        this.printToConsole= new PrintStream(new Console(console),true);
		}
	
	protected abstract Runnable createRunnable() throws JFXException;
	
	protected void displayAlert(final Throwable err) {
		final Alert alert = new Alert(AlertType.ERROR);
		alert.setHeaderText("Cannot create Command.");
		alert.setContentText(String.valueOf(err.getMessage()));
		
		// Create expandable Exception.
		StringWriter sw = new StringWriter();
		PrintWriter pw = new PrintWriter(sw);
		err.printStackTrace(pw);
		
		TextArea textArea = new TextArea(sw.toString());
		textArea.setEditable(false);
		textArea.setWrapText(true);

		BorderPane pane=new BorderPane(new ScrollPane(textArea));
		alert.getDialogPane().setExpandableContent(pane);
		
		alert.showAndWait();
		
	}
	
	private void doCommandStart(final ActionEvent event) {
		doCommandEnd(event);
		final Runnable target;
		try {
			 target = createRunnable();
			if( target ==null ) return;
			
		} catch(final Throwable err)
			{
			err.printStackTrace(realStderr);
			displayAlert(err);
			return;
			}
		
		synchronized(AbstractJfxApplication.class) {
			try {
				this.commandThread=null;
				this.runCommandButton.setDisable(true);
				this.cancelCommandButton.setDisable(false);
				this.commandThread = new Thread(new RunnerDelegate(target));
				this.commandThread.start();
				
			} catch(final Throwable err)
				{
				err.printStackTrace(realStderr);
				System.setErr(realStderr);
				System.setOut(realStdout);
				}
			}
		}
	
	private void doCommandEnd(final  ActionEvent event) {
		synchronized(AbstractJfxApplication.class) {
			this.runCommandButton.setDisable(false);
			this.cancelCommandButton.setDisable(true);

			if(this.commandThread==null) return;
			try {
				this.commandThread.interrupt();
			} catch(Throwable err)
				{
				
				}
			this.commandThread=null;
			System.setErr(realStderr);
			System.setOut(realStdout);
			}
		}
	
	protected Parent fxmlLoad(final String resource) throws Exception
		{
		try
			{
			java.net.URL url= getClass().getResource(resource);
			if(url==null) throw new java.io.IOException("cannot get resource \""+resource+"\" for Class:"+this.getClass());
			final FXMLLoader loader = new FXMLLoader(url);
			loader.setController(this);
			return loader.load();
			}
		catch(final Exception err)
			{
			err.printStackTrace();
			throw err;
			}
		}
	
	private class RunnerDelegate implements Runnable
		{
		final Runnable delegate;
		RunnerDelegate(final Runnable delegate) {
			this.delegate = delegate;
			}
		@Override
		public void run()
			{
			try {
				System.setErr(AbstractJfxApplication.this.printToConsole);
				//not stdout, things like snpeff write to stdout
				this.delegate.run();
				}
			catch(final Throwable err)
				{
				err.printStackTrace(realStderr);
				}
			Platform.runLater(new Runnable()
				{
				@Override
				public void run()
					{
					doCommandEnd(null);
					}
				});
			}
		}
	
	
	
	
	protected class OptionBuilder
		{
		protected final Parent component;
		protected final String option;
		protected int _minCardinality=0;
		protected int _maxCardinality=-1;
		protected Pattern splitPattern=null;
		protected Class<?> itemClass=null;
		
		public OptionBuilder(final Parent component,final String option) throws JFXException {
			this.component=component;
			this.option=option;
			if(this.component==null) throw new JFXException("component is null in ctor ("+option+")");
			if(this.option==null) throw new JFXException("opt is null in ctor");
			}
		
		
		
		
		protected void fill(final List<String> args,String s) {
			if(this.option.startsWith("-")) {
				args.add(this.option);
				args.add(s);
				}
			else if(this.option.endsWith("="))//picard like option
				{
				args.add(this.option+s);
				}
			}
		
		protected String validateType(String s)  throws JFXException {
			if(itemClass!=null) {
				try {
					this.itemClass.getConstructor(String.class).newInstance(s);
				} catch (Exception e) {
					throw new  JFXException("Cannot cast "+s +" to "+itemClass,e);
					}
				}
			return s;
		}
		
		public OptionBuilder minCardinality(int v) { this._minCardinality=v; return this;}
		public OptionBuilder maxCardinality(int v) { this._maxCardinality=v; return this;}
		public OptionBuilder split(final Pattern pat) { this.splitPattern=pat; return this;}
		public OptionBuilder split() {return this.split(Pattern.compile("[\\s]+"));}
		public OptionBuilder itemClass(Class<?> C) { this.itemClass=C;return this;}
		
		protected void validateCardinality(List<?> list)  throws JFXException
			{
			if(list.size()<this._minCardinality) throw new JFXException("Expected at least "+this._minCardinality+" items for "+this.option);
			if(this._maxCardinality!=-1 && list.size()>this._maxCardinality) throw new JFXException("Expected less than "+this._maxCardinality+" items for "+this.option);
			}
		
		public void fill(final List<String> args) throws JFXException {
			if( component==null)  {
				throw new JFXException("component is null ("+this.option+")");
				}
			else if(component instanceof FileChooserPane){
				final FileChooserPane comp = FileChooserPane.class.cast(component);
				final File f= comp.getSelectedFile();
				if(comp.isRequired() && f==null)  throw new JFXException("missing file for "+this.option);
				if(f==null) return;
				
				fill(args,f.getPath());
				}
			else if(component instanceof FilesChooserPane){
				final FilesChooserPane comp = FilesChooserPane.class.cast(component);
				final List<File> list= comp.getSelectedFiles();
				if(list.size()<comp.getMinCardinality()) throw new JFXException("Expected at least "+comp.getMinCardinality()+" items for "+this.option);
				if(comp.getMaxCardinality()!=-1 && list.size()>comp.getMaxCardinality()) throw new JFXException("Expected less than "+comp.getMaxCardinality()+" items for "+this.option);
				for(final File f: list) {
					fill(args,f.getPath());
					}
				}
			else if(component instanceof CheckBox) {
				final CheckBox comp = CheckBox.class.cast(component);
				if(this.option.endsWith("=")) //picard
					{
					fill(args,Boolean.toString(comp.isSelected()));
					}
				else if(comp.isSelected())
					{
					args.add(this.option);
					}
				}
			else if(component instanceof Spinner) {
				Spinner<?> comp = Spinner.class.cast(component);
				Object v=comp.getValue();
				fill(args,String.valueOf(v));
				}
			else if(component instanceof ComboBox)
				{
				ComboBox<?> comp = ComboBox.class.cast(component);
				Object o=comp.getValue();
				if(o==null) return;
				String s=o.toString();
				if(s.isEmpty()) return;
				fill(args,validateType(s));
				}
			else if(component instanceof TextArea) {
				TextArea comp = TextArea.class.cast(component);
				final List<String> list=new ArrayList<>();
				for(final String s: comp.getText().split("[\n]"))
					{
					if(s.trim().isEmpty() || s.startsWith("#")) continue;
					list.add(s);
					}
				validateCardinality(list);
				for(final String s: list) {
					fill(args,validateType(s));
					}
				}
			else if(component instanceof TextField) {
				final TextField comp = TextField.class.cast(component);
				if(this.splitPattern==null)
					{
					String s= comp.getText().trim();
					if(s.isEmpty()) return;
					fill(args,validateType(s));
					}
				else
					{
					List<String> list=new ArrayList<>();
					for(final String s: this.splitPattern.split(comp.getText()))
						{
						if(s.trim().isEmpty()) continue;
						list.add(s);
						}
					validateCardinality(list);
					for(final String s: list) {
						fill(args,validateType(s));
						}
					}
				}
			else
				{
				throw new JFXException("undefined Class of component ("+this.option+") "+this.component.getClass());
				}
			}
		}
	
	
    private class Console extends OutputStream {
	
	    private final TextArea output;
	
	    public Console(TextArea ta) {
	        this.output = ta;
	    }
	    @Override
	    public void write(byte[] b, int off, int len) throws IOException
	    	{
	    	final String str=new String(b, off, len);
	    	if(str.isEmpty()) return;
	    	Platform.runLater(new Runnable()
					{
					@Override
					public void run()
						{
						output.appendText(str);
						if(output.getLength()>50000)
							{
							output.deleteText(0, output.getLength()-1000);
							}
						}
					});
	    	}
	    @Override
	    public void write(byte[] b) throws IOException
	        {
	        super.write(b,0,b.length);
	        }
	    @Override
	    public void write( int c) throws IOException {
	        if(c==-1) return;
	        byte a[]=new byte[]{(byte)c};
	        a[0]=(byte)c;
	        write(a);
	        }
		}
    
    public static class JFXException extends Exception
    	{
		private static final long serialVersionUID = 1L;
		public JFXException() { super();}
    	public JFXException(final String msg) { super(msg);}
    	public JFXException(final Exception err) { super(err);}
    	public JFXException(final String msg,final Exception err) { super(msg,err);}
    	}	
	}
