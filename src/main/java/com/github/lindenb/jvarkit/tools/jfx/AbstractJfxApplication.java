package com.github.lindenb.jvarkit.tools.jfx;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.List;

import com.github.lindenb.jvarkit.jfx.components.FileChooserPane;

import javafx.application.Application;
import javafx.application.Platform;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.fxml.FXML;
import javafx.geometry.Rectangle2D;
import javafx.scene.Parent;
import javafx.scene.control.Button;
import javafx.scene.control.Spinner;
import javafx.scene.control.TextArea;
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
	private TextArea console;

	protected Thread commandThread=null;
	protected AbstractJfxApplication()
		{
		}
	
	@Override
	public void start(Stage stage) throws Exception {
		 final PrintStream pr=new PrintStream(new Console(console),true);
	     System.setErr(pr);
	     System.setOut(pr);
	     
	     
	    Rectangle2D primaryScreenBounds = Screen.getPrimary().getVisualBounds();

        stage.setX(50);
        stage.setY(50);
        stage.setWidth(primaryScreenBounds.getWidth()-100);
        stage.setHeight(primaryScreenBounds.getHeight()-100);
        stage.show();
        runCommandButton.setOnAction(new EventHandler<ActionEvent>()
			{
			@Override
			public void handle(ActionEvent event)
				{
				doCommandStart(event);
				}
			});
        cancelCommandButton.setOnAction(new EventHandler<ActionEvent>()
			{
			@Override
			public void handle(ActionEvent event)
				{
				doCommandEnd(event);
				}
			});
        cancelCommandButton.setDisable(true);
		}
	
	protected abstract Runnable createRunnable() throws JFXException;
	
	private void doCommandStart(ActionEvent event) {
		doCommandEnd(event);
		synchronized(AbstractJfxApplication.class) {
			try {
				commandThread=null;
				final Runnable target = createRunnable();
				if( target ==null ) return;
				this.runCommandButton.setDisable(true);
				this.cancelCommandButton.setDisable(false);
				this.commandThread = new Thread(new RunnerDelegate(target));
				this.commandThread.start();
				
			} catch(Throwable err)
				{
				err.printStackTrace(realStderr);
				System.setErr(realStderr);
				System.setOut(realStdout);
				}
			}
		}
	
	private void doCommandEnd(ActionEvent event) {
		synchronized(AbstractJfxApplication.class) {
			runCommandButton.setDisable(false);
			cancelCommandButton.setDisable(true);

			if(this.commandThread==null) return;
			try {
				this.commandThread.interrupt();
			} catch(Throwable err)
				{
				
				}
			commandThread=null;
			System.setErr(realStderr);
			System.setOut(realStdout);
			}
		}
	
	
	private class RunnerDelegate implements Runnable
		{
		final Runnable delegate;
		RunnerDelegate(Runnable delegate) {
			this.delegate = delegate;
			}
		@Override
		public void run()
			{
			try {
				this.delegate.run();
				}
			catch(Throwable err)
				{
				err.printStackTrace(realStderr);
				}
			doCommandEnd(null);
			}
		}
	
	protected class OptionBuilder
		{
		final Parent component;
		final String option;
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
		public void fill(final List<String> args) throws JFXException {
			if( component==null)  {
				throw new JFXException("component is null ("+this.option+")");
				}
			else if(component instanceof FileChooserPane){
				FileChooserPane comp = FileChooserPane.class.cast(component);
				File f= comp.getSelectedFile();
				if(comp.isRequired() && f==null)  throw new JFXException("missing file for "+this.option);
				if(f==null) return;
				fill(args,f.getPath());
				}
			else if(component instanceof Spinner) {
				Spinner<?> comp = Spinner.class.cast(component);
				Object v=comp.getValue();
				fill(args,String.valueOf(v));
				}
			else
				{
				throw new JFXException("undefined component ("+this.option+")");
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
