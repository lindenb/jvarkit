package com.github.lindenb.jvarkit.tools.vcfviewgui;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Vector;
import java.util.function.Function;
import java.util.logging.Logger;

import javax.script.Bindings;
import javax.script.CompiledScript;
import javax.script.ScriptException;
import javax.xml.namespace.QName;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;

import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.ChartFactory;
import com.github.lindenb.jvarkit.util.igv.IgvSocket;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import javafx.application.Platform;
import javafx.beans.property.ReadOnlyObjectWrapper;
import javafx.beans.value.ObservableValue;
import javafx.collections.FXCollections;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.scene.Scene;
import javafx.scene.chart.Chart;
import javafx.scene.control.Label;
import javafx.scene.control.Menu;
import javafx.scene.control.MenuBar;
import javafx.scene.control.MenuItem;
import javafx.scene.control.Spinner;
import javafx.scene.control.Tab;
import javafx.scene.control.TabPane;
import javafx.scene.control.TableColumn;
import javafx.scene.control.TableView;
import javafx.scene.control.TextArea;
import javafx.scene.control.TextField;
import javafx.scene.control.Alert.AlertType;
import javafx.scene.control.TableColumn.CellDataFeatures;
import javafx.scene.layout.VBox;
import javafx.scene.paint.Color;
import javafx.stage.FileChooser;
import javafx.stage.Stage;
import javafx.stage.WindowEvent;
import javafx.util.Callback;

/**
 * Abstract Stage for Bam and VCF file
 * @author lindenb
 *
 */
public abstract class NgsStage extends Stage {
    protected static final Logger LOG= Logger.getLogger("NgsStage");
	protected static final String JAVASCRIPT_TAB_KEY="JS";
    /** owner Application */
    protected JfxNgs owner;
    /** src file */
	protected final Object urlOrFile;
	/** javascript filtering */
	protected final TextArea javascriptArea=new TextArea();
	/** message stuff */
	protected final Label messageLabel=new Label();
	/** message stuff */
	protected final TextField gotoField=new TextField();
	/** menuBar */
	protected final MenuBar menuBar=new MenuBar();
	/** File Menu */
	protected final Menu fileMenu=new Menu("File");
    /** limit number of items */
	protected final Spinner<Integer> maxItemsLimitSpinner=
			new Spinner<>(0, 10000, 1);

	
	protected abstract class AbstractQualityStage<T>
		extends Stage
		{
		protected abstract class ScanThread 
			extends Thread
			{
			protected long nItems=0L;
			protected final File source;
			protected final CompiledScript compiledScript;
			protected volatile boolean kill_flag=false;
    		protected long lastRefresh =System.currentTimeMillis();
    		protected final List<ChartFactory<T>> factories = new Vector<>();
    		ScanThread(final File source,CompiledScript compiledScript)
				{
				this.source=source;
				this.compiledScript = compiledScript;
				}
            /** called by javascript filters */
            protected boolean accept(final Bindings bindings)
    			{
            	if(compiledScript==null) return true;
    			final Object result;
    			try  {
    				result = this.compiledScript.eval(bindings);
    			} catch(final ScriptException err)
    				{
    				return false;
    				}
    			
    			if(result==null) return false;;
    			if(result instanceof Boolean)
    				{
    				if(Boolean.FALSE.equals(result)) return false;
    				}
    			else if(result instanceof Number)
    				{
    				if(((Number)result).intValue()!=1) return false;
    				}
    			else
    				{
    				return false;
    				}
    			return true;
    			}
            protected void atEnd() {
            	if(kill_flag) {
					LOG.warning("Thread was killed");
					}
				else
    				{
    				repaint();
    				Platform.runLater(new Runnable() {
        				 @Override
        				public void run() {
        					 AbstractQualityStage.this.countItemsLabel.setText(
        						"Done... Number of items: "+nItems+ (kill_flag?" [KILLED]":""));
        				 }
        			 	});
    				}
            	}
            protected void repaint()
				{
				final List<Chart> L=new ArrayList<>(this.factories.size());
				for(final ChartFactory<T> chartter:this.factories) L.add(chartter.build());
				 Platform.runLater(new Runnable() {
					 @Override
					public void run() {
						 for(int i=0;i< L.size();++i)
						 	{
							AbstractQualityStage.this.tabs.get(i).setContent(L.get(i));
						 	}
						AbstractQualityStage.this.countItemsLabel.setText("Running... Number of items: "+nItems);
					 	}
				 	});
				}
            protected void update()
				{
				long now  =System.currentTimeMillis(); 
				if( kill_flag || now - lastRefresh < 5*1000) return ;//5 seconds;
				lastRefresh = now;
				repaint();
				}
			}
		
		protected final ScanThread thread;
		protected final List<Tab> tabs=new Vector<>();
		protected final Label countItemsLabel=new Label();

		protected AbstractQualityStage(File file,final CompiledScript compiledScript)
			{
			setTitle(file.getPath());
			this.thread = createThread(file,compiledScript);
	    	this.addEventHandler(
	    			WindowEvent.WINDOW_CLOSE_REQUEST ,new EventHandler<WindowEvent>() {
	                    @Override
	                    public void handle(WindowEvent event) {
	                    	kill();
	                        }
	                    });

        	this.addEventHandler(
        			WindowEvent.WINDOW_SHOWN ,new EventHandler<WindowEvent>() {
                        @Override
                        public void handle(WindowEvent event) {
                        	thread.start();
                            }
                        });
        	final TabPane tabPane=new TabPane();
        	for(int i=0;i< this.thread.factories.size();++i)
        		{
        		final Tab tab=new Tab(
        				this.thread.factories.get(i).getName(),
        				this.thread.factories.get(i).build()
        				);
        		tab.setClosable(false);
        		this.tabs.add(tab);
        		tabPane.getTabs().add(tab);
        		}
        	final VBox box1=new VBox(tabPane,this.countItemsLabel);
        	final Scene scene=new Scene(box1,1000,800);
        	setScene(scene);
			}
		
		protected abstract ScanThread createThread(File file,final CompiledScript compiledScript);
		@Override
    	protected void finalize() throws Throwable {
    		kill();
	    	super.finalize();
	    	}
    	void kill()
    		{
    		thread.kill_flag=true;
    		}
    	}
	
    public NgsStage(final JfxNgs owner,final Object urlOrFile) throws IOException
    	{
    	this.owner= owner;
    	this.urlOrFile= urlOrFile;
    	this.setTitle(this.urlOrFile.toString());
    	this.maxItemsLimitSpinner.setEditable(true);
    	
    	
    	if(this.owner.javascriptEngine==null) {
    		this.javascriptArea.setEditable(false);
    		this.javascriptArea.setPromptText("Javascript engine is not available");
    	} else
    		{
    		this.javascriptArea.setPromptText("Use this area to create a javascript-bases filter to ignore some items");
    		}
    	
    	
    	this.addEventHandler(
    			WindowEvent.WINDOW_SHOWING ,new EventHandler<WindowEvent>() {
                    @Override
                    public void handle(WindowEvent event) {
                        owner.registerStage(NgsStage.this);
                        NgsStage.this.reloadData();
                        }
                    });
    	this.addEventHandler(
    			WindowEvent.WINDOW_CLOSE_REQUEST ,new EventHandler<WindowEvent>() {
                    @Override
                    public void handle(WindowEvent event) {
                    	NgsStage.this.closeNgsResource();
                    	owner.unregisterStage(NgsStage.this);
                        }
                    });
       
        this.gotoField.setPrefColumnCount(15);
        this.gotoField.setEditable(true);
        
        this.fileMenu.getItems().add(JfxNgs.createMenuItem("Open",new Runnable() {
			@Override
			public void run() {
				owner.openNgsFiles(NgsStage.this);
			}
		}));
        this.fileMenu.getItems().add(JfxNgs.createMenuItem("Save Filtered Data as... ",new Runnable() {
			@Override
			public void run() {
				doMenuSaveAs();
			}
		}));
    	this.menuBar.getMenus().add(this.fileMenu);
        
    	final Menu statsMenu=new Menu("Stats");
    	this.menuBar.getMenus().add(statsMenu);

        statsMenu.getItems().add(JfxNgs.createMenuItem("Show for current Reads",new Runnable() {
			@Override
			public void run() {
				doMenuShowLocalStats();
			}
		}));
        statsMenu.getItems().add(JfxNgs.createMenuItem("Show for whole file",new Runnable() {
			@Override
			public void run() {
				doMenuShowWholeStats();
			}
		}));
    	}
    
    @Override
    protected void finalize() throws Throwable {
    	closeNgsResource();
    	super.finalize();
    	}
    
    /** path to snippets */
    protected String getSnippetResourcePath()
    	{
    	return null;
    	}
    
    
    protected Menu createSnippetMenu() {
    	final Menu menu=new Menu("Snippets");
    	final String rsrc = getSnippetResourcePath();
    	if(rsrc!=null) {
    		InputStream in=null;
    		XMLEventReader r=null;
    		try
    			{
    			in = getClass().getResourceAsStream(rsrc);
    			if(in!=null)
    				{
    				r=XMLInputFactory.newFactory().createXMLEventReader(in);
    				final QName labelAtt=new QName("label");
    				while(r.hasNext())
    					{
    					XMLEvent evt=r.nextEvent();
    					if(!evt.isStartElement() ) continue;
    					StartElement start=evt.asStartElement();
    					if(!start.getName().getLocalPart().equals("code")) continue;
    					Attribute attLabel=start.getAttributeByName(labelAtt);
    					if(attLabel!=null && r.hasNext() && r.peek().isCharacters())
    						{
        					final MenuItem item=new MenuItem(attLabel.getValue());
    						final String code= r.nextEvent().asCharacters().getData();

    						item.setOnAction(new EventHandler<ActionEvent>() {
								@Override
								public void handle(ActionEvent event) {
									NgsStage.this.javascriptArea.setText(code);
								}
							});
    						menu.getItems().add(item);
    						}
    					}
    				}
    			}
    		catch(Exception err)
    			{
    			LOG.warning(err.getMessage());
    			}
    		finally
    			{
    			CloserUtil.close(r);
    			CloserUtil.close(in);
    			}
    		}
    	else
    		{
    		LOG.warning("No snippets defined for "+getClass());
    		}
    	return menu;
    	}
    
    /** send those command to IGV */
    protected void openInIgv(final List<String> commands)
    	{
    	if(commands==null || commands.isEmpty()) return;
    	@SuppressWarnings("resource")
		final IgvSocket socket=new IgvSocket();
    	final Runnable r=socket.buildRunnable(commands);
    	new Thread(r).start();
    	}

    /** create a MenuItem saving a table */
    protected <T> MenuItem menuForSavingTable(final String tableName,TableView<T> table)
    	{
    	final MenuItem item=new MenuItem("Save "+tableName+" as ...");
    	item.setOnAction(new EventHandler<ActionEvent>()
			{
			@Override
			public void handle(ActionEvent event)
				{
				FileChooser fc= owner.newFileChooser();
				fc.setSelectedExtensionFilter(new FileChooser.ExtensionFilter("TSV","tsv"));
				File fout= owner.updateLastDir(fc.showSaveDialog(NgsStage.this));
				if(fout==null) return ;
				PrintWriter out=null;
				try
					{
					out = new PrintWriter(fout);
					for(int x=0;x< table.getColumns().size();++x)
						{
						out.print(x==0?"#":"\t");
						out.print(table.getColumns().get(x).getText());
						}
					out.println();
					for(int y=0;y< table.getItems().size();++y)
						{
						final T row=table.getItems().get(y);
						for(int x=0;x< table.getColumns().size();++x)
							{
							if(x>0) out.print("\t");
							out.print(table.getColumns().get(x).getCellObservableValue(row).getValue());
							}
						out.println();
						}
					out.flush();
					}
				catch(Exception err)
					{
					err.printStackTrace();
					JfxNgs.showExceptionDialog(NgsStage.this, err);
					}
				finally
					{
					CloserUtil.close(out);
					}
				
				}
			});
    	return item;
    	}
    
    /** send a goto command to IGV */
    protected void openInIgv(final Locatable feature)
    	{
    	if(feature==null) {
    		updateStatusBar(AlertType.WARNING,"No Feature was selected");
    		return;
    	}
    	openInIgv(
    			Collections.singletonList("goto "+feature.getContig()+":"+feature.getStart()+"-"+feature.getEnd())
    			);
    	}
    
    abstract void openInIgv();
    /** close the NGS resource , even if the window was not opened */
    abstract void closeNgsResource();
    /** reload all data */
    abstract void reloadData();
    
    /** show stats */
    protected abstract void doMenuShowWholeStats();
    /** show stats for whole file */
    protected abstract void doMenuShowLocalStats();
    
    protected void updateStatusBar(final AlertType type,final Object o)
    	{
    	final Color textColor;
    	switch(type)
    		{
    		case CONFIRMATION: textColor=Color.BLUE; break;
    		case ERROR: textColor=Color.RED;break;
    		case INFORMATION: textColor=Color.GREEN;break;
    		case NONE: textColor=Color.BLACK; break;
    		case WARNING: textColor=Color.ORANGE;break;
    		default: textColor=Color.BLACK; break;
    		}
    	this.messageLabel.setTextFill(textColor);
    	if(o==null) {
    		this.messageLabel.setText("");
    		}
    	else if(o instanceof Throwable )
    		{
    		this.messageLabel.setText(String.valueOf(Throwable.class.cast(o).getMessage()));
    		}
    	else
    		{
    		this.messageLabel.setText(String.valueOf(o));
    		}
    	}
    
    protected <T,R> TableColumn<T,R> makeColumn(final String tag,final Function<T,R> supplier)
    	{
        final TableColumn<T,R>  col = new TableColumn<>(tag);
	        col.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<T,R>, ObservableValue<R>>() {				
				@Override
				public ObservableValue<R> call(CellDataFeatures<T, R> param) {
					return new ReadOnlyObjectWrapper<R>(supplier.apply(param.getValue()));
					}
				});
	        return col;
    	}
    /** called by javascript filters */
    protected boolean accept(final CompiledScript script,final Bindings bindings)
		{
		final Object result;
		try  {
			result = script.eval(bindings);
		} catch(ScriptException err)
		{
			LOG.severe(err.getMessage());
			err.printStackTrace();
			updateStatusBar(AlertType.WARNING,err);
			return false;
		}
		
		if(result==null) return false;;
		if(result instanceof Boolean)
			{
			if(Boolean.FALSE.equals(result)) return false;
			}
		else if(result instanceof Number)
			{
			if(((Number)result).intValue()!=1) return false;
			}
		else
			{
			updateStatusBar(AlertType.WARNING,"Script returned something that is not a boolean or a number:"+result.getClass());
			return false;
			}
		return true;
		}
    
    /** build a table view for a Dictionary */
    protected Tab buildDictTab(final SAMSequenceDictionary dict)
        {
    	
		/* build INFO Table */
		final TableView<SAMSequenceRecord> table=new TableView<>(
				dict==null?
				FXCollections.observableArrayList():
				FXCollections.observableArrayList(dict.getSequences())
				);
		
        table.getColumns().add(makeColumn("Name",SSR->SSR.getSequenceName()));
        table.getColumns().add(makeColumn("Length",SSR->SSR.getSequenceLength()));
        table.getColumns().add(makeColumn("Assembly",SSR->SSR.getAssembly()));
        final Tab tab=new Tab("Dict", table);
        tab.setClosable(false);
        
        return tab;
        }
    /** save filtered Data As ... */
    protected abstract void doMenuSaveAs();
    
    protected abstract SAMSequenceDictionary getSAMSequenceDictionary();
    
    /** called by main stage: set location box content and call reloadData */
    protected void moveTo(final String s)
    	{
    	this.gotoField.setText(s);
    	this.reloadData();
    	}
    
    protected Interval parseInterval(final String location)
    	{
    	final SAMSequenceDictionary dict=getSAMSequenceDictionary();
		final String contig;
		int colon =location.indexOf(":");
		if(colon==-1)
			{
			contig=location;
			}
		else
			{
			contig=location.substring(0,colon);
			}
		
		SAMSequenceRecord ssr= dict.getSequence(contig);
		if(ssr==null && !contig.startsWith("chr"))
			{
			ssr= dict.getSequence("chr"+contig);
			}
		if(ssr==null && contig.startsWith("chr"))
			{
			ssr= dict.getSequence( contig.substring(3));
			}
		if(ssr==null)
			{
			updateStatusBar(AlertType.WARNING, "Cannot find contig in dictionary: "+location);
			return null;
			}
		
		if(colon!=-1)
			{
			int hyphen=location.indexOf('-');
			Integer start=null,end=null;
			if(hyphen==-1)
				{
				final String startStr=location.substring(colon+1).trim();
				try {
					start= Math.max(0, Integer.parseInt(startStr));
					end=ssr.getSequenceLength();
					}
				catch(final NumberFormatException err ) {
					start=null;
					LOG.warning(startStr);
					updateStatusBar(AlertType.WARNING, "Bad Start in : "+location);
					return null;
					}
				}
			else
				{
				try {
					start= Math.max(0, new Integer(location.substring(colon+1,hyphen).trim()));
					end= Math.min(
							Integer.parseInt(location.substring(hyphen+1).trim()),
							ssr.getSequenceLength()
							);
					}
				catch(NumberFormatException err )
					{
					start=null;end=null;
					LOG.warning(location);
					updateStatusBar(AlertType.WARNING, "Bad Start/End in : "+location);
					return null;
					}
				}
			if(start!=null && end!=null && start.compareTo(end)<=0)
				{
				return new Interval(ssr.getSequenceName(), start, end);
				}
			else
				{
				return null;
				}
			}
		else
			{
			return new Interval(ssr.getSequenceName(), 0, ssr.getSequenceLength());
			}
    	}
    
}
