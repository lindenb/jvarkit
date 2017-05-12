package com.github.lindenb.jvarkit.tools.jfx;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;


import com.github.lindenb.jvarkit.tools.bioalcidae.BioAlcidae;
import com.github.lindenb.jvarkit.tools.vcffilterjs.VCFFilterJS;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import javafx.application.Application;
import javafx.application.Platform;
import javafx.beans.property.SimpleStringProperty;
import javafx.geometry.Insets;
import javafx.scene.Scene;
import javafx.scene.control.Alert;
import javafx.scene.control.Button;
import javafx.scene.control.Menu;
import javafx.scene.control.MenuBar;
import javafx.scene.control.MenuItem;
import javafx.scene.control.Tab;
import javafx.scene.control.TabPane;
import javafx.scene.control.TableColumn;
import javafx.scene.control.TableView;
import javafx.scene.control.TextArea;
import javafx.scene.control.Alert.AlertType;
import javafx.scene.control.TabPane.TabClosingPolicy;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.FlowPane;
import javafx.scene.layout.StackPane;
import javafx.scene.layout.VBox;
import javafx.stage.FileChooser;
import javafx.stage.Stage;

public class JvarkitCentral extends Application {
	private static final Logger LOG = Logger.build(JvarkitCentral.class).make();

	@Override
	public void start(Stage primaryStage) throws Exception {
		primaryStage.setTitle("Jvarkit-Central");
		StackPane root = new StackPane();
		root.setPadding(new Insets(2));
		final TableView<Class<?>> tableView =new TableView<>();
		
        final TableColumn<Class<?>, String> nameCol = new TableColumn<>("Name");
        nameCol.setCellValueFactory(CB-> {
        	final String value;
        	Class<?> clazz = CB.getValue();
        	final Program program=clazz.getAnnotation(Program.class);
        	if(!Launcher.class.isAssignableFrom(clazz) || program==null)
        		{
        		value=null;
        		}
        	else
        		{
        		value=program.name();
        		}
	        return new SimpleStringProperty(value);
	        });
        final TableColumn<Class<?>, String> descCol = new TableColumn<>("Description");
        descCol.setCellValueFactory(CB-> {
        	final String value;
        	Class<?> clazz = CB.getValue();
        	final Program program=clazz.getAnnotation(Program.class);
        	if(!Launcher.class.isAssignableFrom(clazz) || program==null)
        		{
        		value=null;
        		}
        	else
        		{
        		value=program.description();
        		}
	        return new SimpleStringProperty(value);
	        });
        
        
        tableView.getColumns().addAll(nameCol,descCol);
		final BorderPane borderPane1 = new BorderPane(tableView);
		borderPane1.setPadding(new Insets(10));
		final Button but=new Button("New Instance...");
		but.setOnAction(AE->{
			final Class<?> clazz=tableView.getSelectionModel().getSelectedItem();
			if(clazz==null) return;
			final Program program=clazz.getAnnotation(Program.class);
			if(!Launcher.class.isAssignableFrom(clazz) || program==null) return;
			createNewInstanceOf(clazz);
		});
		FlowPane bottom=new FlowPane(but);
		borderPane1.setBottom(bottom);
		
		tableView.getItems().addAll(String.class,Integer.class,
				VCFFilterJS.class,
				BioAlcidae.class
				);
		
		
		root.getChildren().add(borderPane1);
		
        //root.getChildren().add(btn);
        primaryStage.setScene(new Scene(root, 300, 250));
        primaryStage.show();
		}
	
	private void createNewInstanceOf(final Class<?> clazz)
		{
		final LauncherStage stage = new LauncherStage(clazz);
		
		stage.show();
		}
	
	private static class LaunchAndRun extends Thread
		{
		private LauncherStage owner=null;
		private  Launcher instance=null;
		private  String args[]=null;
		private int returnStatus=0;
		@Override
		public void run() {
			
			final PrintStream redirectStream=new PrintStream(new ConsolePrintStream()
				{
				@Override
				public void writeToDevice(final String textToPrint) {
						Platform.runLater(()->
							{
							if(owner.runner != LaunchAndRun.this) return;
							int L = owner.outputArea.getText().length();
							while(L>10000) {
								owner.outputArea.deleteText(0,10000);
								L = owner.outputArea.getText().length();
								}
							owner.outputArea.appendText(textToPrint);
							}
						);
					}	
				});
			
			final InputStream no_stdin= new ByteArrayInputStream(new byte[0]);

			
			try
				{
				
				instance.stdin(no_stdin);
				instance.stdout(redirectStream);
				instance.stdout().flush();
				instance.stderr(redirectStream);
				int ret= instance.instanceMain(args);
				this.returnStatus=ret;
				}
			catch(Exception err)
				{
				LOG.error(err);
				Platform.runLater(()->{
					if(this.owner.runner != LaunchAndRun.this) return;
					this.owner.runner=null;
					Alert alert=new Alert(AlertType.ERROR,"An error occured "+err.getMessage());
				 	alert.showAndWait();
				 	});
				returnStatus=-1;
				return;
				}
		
		 Platform.runLater(()->{
			if(this.owner.runner != LaunchAndRun.this) return;
			this.owner.runner=null;
			Alert alert=new Alert(returnStatus==0?AlertType.INFORMATION:AlertType.ERROR,
					"Program exited with status "+returnStatus
					);
			alert.showAndWait();
			});
		}
		}
	private static  class LauncherStage extends Stage
		{
		private File lastDir=null;
		final TabPane tabPane ;
		private final TextArea outputArea;
		private final TextArea textArea;
		private final Class<?> clazz;
		private volatile LaunchAndRun runner=null;
		
		LauncherStage(Class<?> clazz) {
			this.clazz = clazz;
			this.setTitle(clazz.getSimpleName());
			VBox root = new VBox();
			root.setPadding(new Insets(2));

			final Menu menuFile=new Menu("File");	
			
			MenuItem menuItem=new MenuItem("Close");
			menuItem.setOnAction(AE->{doMenuStop();LauncherStage.this.close();});
			menuFile.getItems().add(menuItem);
			menuItem=new MenuItem("Insert Path to Open");
			menuItem.setOnAction(AE->{doMenuInsertPath(true);});
			menuFile.getItems().add(menuItem);
			menuItem=new MenuItem("Insert Path to Save");
			menuItem.setOnAction(AE->{doMenuInsertPath(false);});
			menuFile.getItems().add(menuItem);

			
			final MenuBar menuBar=new MenuBar(menuFile);
			root.getChildren().add(menuBar);
			
			this.tabPane = new TabPane();
			this.tabPane.setTabClosingPolicy(TabClosingPolicy.UNAVAILABLE);
			this.tabPane.setPadding(new Insets(10));
			
			this.textArea=new TextArea("");
			this.textArea.setPrefColumnCount(80);
			this.textArea.setPrefRowCount(25);
			BorderPane borderPane1 = new BorderPane(this.textArea);
			
			final Button runButton=new Button("Run");
			runButton.setOnAction(AE->doMenuStart());
			FlowPane bottom=new FlowPane(runButton);
			borderPane1.setBottom(bottom);
			
			
			final Tab writeCmdTab=new Tab("Command",borderPane1);
			tabPane.getTabs().add(writeCmdTab);
			

			this.outputArea=new TextArea("");
			this.outputArea.setPrefColumnCount(80);
			this.outputArea.setPrefRowCount(25);
			this.outputArea.setEditable(false);
			borderPane1 = new BorderPane(this.outputArea);
			
			final Button stopBut=new Button("Stop");
			stopBut.setOnAction(AE->doMenuStop());
			bottom=new FlowPane(stopBut);
			borderPane1.setBottom(bottom);
			final Tab outputTab=new Tab("Output",borderPane1);			
			this.tabPane.getTabs().add(outputTab);
			
			
			root.getChildren().add(tabPane);
			this.setScene(new Scene(root, 300, 250));
			this.setOnCloseRequest(AE->doMenuStop());
			}
		
		protected void doMenuInsertPath(boolean openDialog) {
			final FileChooser fc=new FileChooser();
			fc.setInitialDirectory(this.lastDir);
			File f=null;
			if(openDialog)
				{
				f=fc.showOpenDialog(LauncherStage.this);
				}
			else
				{
				f=fc.showSaveDialog(LauncherStage.this);
				}
			if(f==null) return;
			this.lastDir=f.getParentFile();
			this.textArea.insertText(this.textArea.getCaretPosition(),f.getPath());
			}

		
		private synchronized void doMenuStop()
			{
			if(this.runner!=null)
				{
				try {this.runner.interrupt(); } catch(Exception err) {err.printStackTrace();}
				this.runner=null;
				}
			this.tabPane.getSelectionModel().select(0);
			}
		private synchronized void doMenuStart()
			{
			if(this.runner!=null)
				{	
				final Alert alert=new Alert(AlertType.WARNING,  "Program is already running....");
				alert.showAndWait();
				return;
				}
			this.tabPane.getSelectionModel().select(1);
			
			
			final List<String> args;
			try {
				args = this.getArguments();
				}
			catch(final IllegalArgumentException err)
				{
				LOG.error(err);
				final Alert alert=new Alert(AlertType.ERROR, 
						String.valueOf(err.getMessage())
						);
				alert.showAndWait();
				return;
				}
			LaunchAndRun run =new LaunchAndRun();
			//final Launcher instance;
			try {
				run.instance = (Launcher)this.clazz.newInstance();
				}
			catch(final Throwable err)
				{
				LOG.error(err);
				final Alert alert=new Alert(AlertType.ERROR, 
						"Error creating a new instance of "+clazz.getName()
						);
				alert.showAndWait();
				return;
				}
			run.args = args.toArray(new String[args.size()]);
			this.runner=run;
			this.runner.owner = this;
			this.runner.start();
			}
		
		protected List<String> getArguments() {
			final List<String> args = new ArrayList<>();
			final String s=this.textArea.getText();
			int i=0;
			while(i<s.length())
				{
				if(Character.isWhitespace(s.charAt(i))) {++i;continue;}
				if(s.charAt(i)=='\"' || s.charAt(i)=='\'')
					{
					char quote=s.charAt(i);
					i++;
					final StringBuilder sb=new StringBuilder();
					while(i< s.length())
						{
						char c = s.charAt(i);
						++i;
						if(c==quote) break;
						if(c=='\\')
							{
							if(i+1>=s.length())
								{	
								throw new IllegalArgumentException("Unclosed string after "+sb.toString());
								}
							c= s.charAt(i);
							switch(c)
								{
								case '\'': sb.append('\'');break;
								case '\"': sb.append('\"');break;
								case 'n': sb.append('\n');break;
								case 't': sb.append('\t');break;
								case '\\': sb.append('\\');break;
								default: throw new IllegalArgumentException("Unknown escape sequence after: "+sb.toString());
								}
							}
						else
							{
							sb.append(c);
							}
						}
					args.add(sb.toString());
					}
				else
					{
					final StringBuilder sb=new StringBuilder();
					while(i< s.length() && !Character.isWhitespace(s.charAt(i)))
						{
						sb.append(s.charAt(i));
						i++;
						}
					args.add(sb.toString());
					}
				}
			LOG.info(this.clazz.getName()+" args are "+args);
			return args;
			}
		
		}
	
	private static class ConsolePrintStream extends OutputStream
		{
		private StringBuilder buffer= new StringBuilder();
		@Override
		public void write(int b) throws IOException {
			if(b==-1) flush();
			buffer.append((char)b);
			if(b=='\n' || buffer.length()>1000) flush();
			}
		@Override
		public void flush() throws IOException {
			if(buffer.length()==0) return;
			writeToDevice(buffer.toString());
			buffer.setLength(0);
			}
		public void writeToDevice(String s)
			{
			
			}
		}

	
    public static void main(String[] args) {
        launch(args);
    }
}
