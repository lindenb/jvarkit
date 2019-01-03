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
package com.github.lindenb.jvarkit.tools.vcfviewgui;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FilterOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.function.Supplier;
import java.util.stream.Collectors;

import javax.script.CompiledScript;
import javax.script.ScriptException;
import javax.script.SimpleBindings;
import javax.xml.namespace.QName;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;

import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.ChartFactory;
import com.github.lindenb.jvarkit.util.hershey.JfxHershey;
import com.github.lindenb.jvarkit.util.igv.IgvSocket;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import javafx.application.Platform;
import javafx.beans.property.ReadOnlyObjectWrapper;
import javafx.beans.value.ObservableValue;
import javafx.collections.FXCollections;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.geometry.Insets;
import javafx.geometry.Orientation;
import javafx.scene.Scene;
import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.chart.Chart;
import javafx.scene.control.Label;
import javafx.scene.control.Menu;
import javafx.scene.control.MenuBar;
import javafx.scene.control.MenuItem;
import javafx.scene.control.ScrollPane;
import javafx.scene.control.Separator;
import javafx.scene.control.SeparatorMenuItem;
import javafx.scene.control.Spinner;
import javafx.scene.control.Tab;
import javafx.scene.control.TableCell;
import javafx.scene.control.TableColumn;
import javafx.scene.control.TableView;
import javafx.scene.control.TextArea;
import javafx.scene.control.TextField;
import javafx.scene.control.Tooltip;
import javafx.scene.input.Clipboard;
import javafx.scene.input.ClipboardContent;
import javafx.scene.control.Alert.AlertType;
import javafx.scene.control.Alert;
import javafx.scene.control.Button;
import javafx.scene.control.ButtonType;
import javafx.scene.control.Hyperlink;
import javafx.scene.control.TableColumn.CellDataFeatures;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.FlowPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Pane;
import javafx.scene.layout.Priority;
import javafx.scene.layout.VBox;
import javafx.scene.paint.Color;
import javafx.scene.paint.CycleMethod;
import javafx.scene.paint.LinearGradient;
import javafx.scene.paint.Paint;
import javafx.scene.paint.Stop;
import javafx.scene.text.Font;
import javafx.scene.text.Text;
import javafx.scene.text.TextFlow;
import javafx.stage.FileChooser;
import javafx.stage.Stage;
import javafx.stage.WindowEvent;
import javafx.stage.FileChooser.ExtensionFilter;
import javafx.stage.Modality;
import javafx.util.Callback;

/**
 * Abstract Stage for Bam and VCF file
 * @author lindenb
 *
 */
public abstract class NgsStage<HEADERTYPE,ITEMTYPE extends Locatable> extends Stage {
    protected static final Logger LOG= Logger.build(NgsStage.class).make();
	protected static final String JAVASCRIPT_TAB_KEY="JS";
	protected static final String TOOL_CONTEXT_KEY="tools";
	protected static final String PEDIGREE_CONTEXT_KEY="pedigree";
	protected static final String HEADER_CONTEXT_KEY="header";
	protected static final String OUT_CONTEXT_KEY="out";
	protected static final String ITER_CONTEXT_KEY="iter";
	private static final int REFRESH_SECOND=Integer.parseInt(System.getProperty("jfxngs.refresh.seconds","5"));
    static final List<ExtensionFilter> JS_EXTENSION_FILTERS=Arrays.asList(
    		new ExtensionFilter("Javascript Files", "*.js","*.javascript")
			);

	/** owner Application */
    protected final JfxNgs owner;
    /** src file */
	private final NgsFile<HEADERTYPE,ITEMTYPE> ngsFile;
	/** javascript filtering */
	protected final TextArea javascriptArea=new TextArea();
	/** draw the karyotype */
	protected final SeqDictionaryCanvas seqDictionaryCanvas;
	/** message stuff */
	protected final Label messageLabel=new Label();
	/** message stuff */
	protected final TextField gotoField=new TextField();
	/** menuBar */
	protected final MenuBar menuBar=new MenuBar();
	/** File Menu */
	protected final Menu fileMenu=new Menu("File");
	/** Stats Menu */
	protected final Menu statsMenu=new Menu("Stats");
	/** Bioalcidae Menu */
	protected final Menu bioalcidaeMenu=new Menu("Bioalcidae");
	/** accelerator to get the starting position of the whole genome */
	private final Map<String,Long> chrom2start;


    /** limit number of items */
	protected final Spinner<Integer> maxItemsLimitSpinner=
			new Spinner<>(0, 10000, 1);
	
	/** limited writer for Bioalcidae Security */
	private static class LimitSecurityStream extends FilterOutputStream
		{
		private final long limit;
		private long written=0L;
		private boolean error_raised=false;
		LimitSecurityStream( final OutputStream delegate,long limit){
			super(delegate);
			this.limit=limit;
			}
		@Override
		public void write(byte[] cbuf, int off, int len) throws IOException {
			if(error_raised) return;
			if(this.limit!=-1L && this.written+len>=this.limit) {
				error_raised=true;
				final String msg="For security Reason I cannot write more than "+this.limit+" characters.";
				LOG.severe(msg);
				throw new IOException(msg);
				}
			super.out.write(cbuf, off, len);
			this.written+=len;
			}
		
		@Override
		public void write(int b) throws IOException {
			byte array[]=new byte[]{(byte)b};
			this.write(array,0,array.length);
			}
		
		
		}
	
	/** resizable canvas http://stackoverflow.com/questions/31761361/ */
	 protected static abstract class ResizableCanvas extends Pane {

	        private final Canvas canvas;

	        ResizableCanvas(final double width,final  double height) {
	            this.canvas = new Canvas(width, height);
	            this.getChildren().add(this.canvas);
	            this.canvas.heightProperty().addListener(CL->{repaintCanvas();});
	            this.canvas.widthProperty().addListener(CL->{repaintCanvas();});
	        }

	        public Canvas getCanvas() {
	            return canvas;
	        }
	        
	        public abstract void repaintCanvas();
	  
	        @Override
	        protected void layoutChildren() {
	            final double x = snappedLeftInset();
	            final double y = snappedTopInset();
	            final double w = snapSize(getWidth()) - x - snappedRightInset();
	            final double h = snapSize(getHeight()) - y - snappedBottomInset();
	            this.canvas.setLayoutX(x);
	            this.canvas.setLayoutY(y);
	            this.canvas.setWidth(w);
	            this.canvas.setHeight(h);
	        }
	    }

	
	
	/** Simple Wrapper for Iterator, logging things every xxx seconds */
	protected  class LogCloseableIterator
		extends AbstractIterator<ITEMTYPE>
		implements CloseableIterator<ITEMTYPE>
		{
		volatile int __stop_flag=0; 
		private boolean _closed=false;
		private final  CloseableIterator<ITEMTYPE> delegate;
		private long clock = System.currentTimeMillis();
		private long count=0L;
		private final Label messageLabel;
		LogCloseableIterator(final CloseableIterator<ITEMTYPE> iter, final Label messageLabel) {
			this.delegate = iter;
			this.messageLabel=messageLabel;
			if(iter==null) throw new NullPointerException("delegate is null");
			}
		
		@Override
		protected ITEMTYPE advance() {
			if(this._closed) return null;
			if(__stop_flag!=0) {
				this.close();
				return null;
				}
			if( this.delegate.hasNext()) {
				final ITEMTYPE t= this.delegate.next();
				++count;
				long now = System.currentTimeMillis();
				if(now - this.clock  > 1*1000) {
					final String msg;
					if(t!=null && t.getContig()!=null)
						{
						msg = "N="+count+" last:"+t.getContig()+":"+t.getStart()+"-"+t.getEnd();
						}
					else
						{
						msg = "N="+count;
						}
					if(this.messageLabel!=null)
						{
						Platform.runLater(()->{messageLabel.setText(msg);});
						}
					else
						{
						LOG.info(msg);
						}
					this.clock = now;
					}
				return t;
				} else {
					close();
					return null;
				}
			}
		
		@Override
		public void close() {
			CloserUtil.close(this.delegate);
			this._closed=true;
			}
		}

	/** simple pair chromosome/pos, used when clicking on ideogram */
	protected static class ContigPos
	implements Comparable<ContigPos>, Locatable
		{
		final String contig;
		final int position;
		ContigPos(final String contig,final int position) {
			this.contig=contig;
			this.position=position;
		}
		
		@Override
		public String getContig() {
			return this.contig;
			}
		@Override
		public int getStart() {
			return this.position;
			}
		@Override
		public int getEnd() {
			return this.position;
			}
		@Override
		public int compareTo(final ContigPos o) {
			int i=contig.compareTo(o.contig);
			if(i!=0) return i;
			return position-o.position;
			}
		@Override
		public String toString() {
			return this.contig+":"+this.position;
			}
		}
	
	/* pane used to draw the karyotype of the reference */
	protected class SeqDictionaryCanvas
		extends ResizableCanvas
		{
		private final SAMSequenceDictionary dict;
		private final double refLength;
		private ContigPos itemsStart=null;
		private ContigPos itemsEnd=null;
		private Interval selectInterval=null;
	    private final Tooltip mousePositionToolTip = new Tooltip("dddd");
    	private java.text.DecimalFormat niceIntFormat = new java.text.DecimalFormat("###,###");

		SeqDictionaryCanvas(final SAMSequenceDictionary dict) {
			super(200,100);
			this.setMinHeight(20.0);
			this.setMinWidth(500.0);
			this.dict = dict;
			this.refLength = dict.getReferenceLength();
			this.setPrefWidth(Double.MAX_VALUE);
			//this.setPrefWidth(25.0);
           
            this.getCanvas().setOnMouseClicked(AE->{
            	if(AE.getClickCount()<2 || this.getCanvas().getWidth()==0) return;
            	final ContigPos cp = pixel2base(AE.getX());
            	if(cp==null) return;
            	NgsStage.this.gotoField.setText(cp.contig+":"+cp.position);
                NgsStage.this.reloadData();		
            	});
            Tooltip.install(this.getCanvas(),this.mousePositionToolTip);
            this.getCanvas().setOnMouseMoved(AE->{
           		final ContigPos cp = pixel2base(AE.getX());
           		this.mousePositionToolTip.setText(cp==null?"":cp.contig+":"+this.niceIntFormat.format(cp.position));});
			}
		
		private ContigPos pixel2base(final double pixx)
			{
			final long x= (long)((pixx/this.getCanvas().getWidth())*this.refLength);
			return NgsStage.this.convertGenomicIndexToContigPos(x);
			}
		
		public void setItemsInterval(final ContigPos start,final ContigPos end) {
			this.itemsStart = start;
			this.itemsEnd = end;
			repaintCanvas();
			}
		public void setSelectInterval(final Interval selectInterval) {
			this.selectInterval = selectInterval;
			repaintCanvas();
			}
		
		private double base2pixel(final String contig,final int v) {
			final long n= NgsStage.this.convertContigPosToGenomicIndex(contig,v);
			if(n==-1L) return -9999.99;
			return ((double)(n)/this.refLength)*this.getCanvas().getWidth();
			}
		
		@Override
		public void repaintCanvas() {
            final double boundrec=5.0;
			final double width = getCanvas().getWidth();
			final double height = getCanvas().getHeight();
			if(width<=1.0 || height<=1.0) return;
            final GraphicsContext gc = this.getCanvas().getGraphicsContext2D();
            gc.setGlobalAlpha(1.0);
            gc.clearRect(0, 0, width, height);
            final Paint p1= new LinearGradient(0, 0, 0, 1.0, true, CycleMethod.NO_CYCLE,
            		new Stop(0.0, Color.DARKGRAY),
            		new Stop(0.5, Color.WHITE),
            		new Stop(1.0,Color.DARKGRAY)
            		);
            final Paint p2= new LinearGradient(0, 0, 0, 1.0, true, CycleMethod.NO_CYCLE,
            		new Stop(0.0, Color.DARKSLATEBLUE),
            		new Stop(0.5, Color.WHITE),
            		new Stop(1.0,Color.DARKSLATEBLUE)
            		);
            final JfxHershey hershey = new JfxHershey();
            /* paint each chrom */
            for(int i=0;i< this.dict.size();++i)
            	{
            	final SAMSequenceRecord ssr = this.dict.getSequence(i);
            	final double x0 = this.base2pixel(ssr.getSequenceName(),1);
            	final double x1 = this.base2pixel(ssr.getSequenceName(),ssr.getSequenceLength());
            	double labelh = height-5;
            	double labelw=Math.min(ssr.getSequenceName().length()*labelh,(x1-x0));
            	
            	gc.setLineWidth(0.5);
            	gc.setFill(i%2==0?p1:p2);
            	gc.fillRoundRect(x0, 1, (x1-x0), height-2, boundrec, boundrec);
            	gc.setStroke(Color.BLACK);
            	gc.strokeRoundRect(x0, 1, (x1-x0), height-2, boundrec, boundrec);
            	if(labelh>3 && labelw>3) {
	            	gc.setStroke(Color.BLACK);
	            	gc.setLineWidth(1.0);
	            	hershey.paint(gc, ssr.getSequenceName(), x0+(x1-x0)/2.0-labelw/2.0,2.0, labelw, labelh);
	            	}
            	}
            /** draw current region of items */
            if(this.itemsStart!=null && this.itemsEnd!=null) {
            	final double x0 = base2pixel(this.itemsStart.contig,this.itemsStart.position);
            	final double x1 =  base2pixel(this.itemsEnd.contig,this.itemsEnd.position);
            	gc.setGlobalAlpha(0.8);
            	gc.setFill(Color.ORANGE);
            	gc.fillRoundRect(x0, 1, (x1-x0)+1.0, height-2, 1, 1);
            	}
            /* draw current selected item */
            if(selectInterval!=null) {
            	final double x0 = base2pixel(selectInterval.getContig(),selectInterval.getStart());
            	final double x1 =  base2pixel(selectInterval.getContig(),selectInterval.getEnd());
            	gc.setGlobalAlpha(0.8);
            	gc.setFill(Color.RED);
            	gc.fillRoundRect(x0, 1, (x1-x0)+1.0, height-2, 1, 1);
            	}
			}

		}
	
	/** generate Bioalcidae-like window */
	 protected  abstract class AbstractAwkLike {
		 protected final TextArea scriptArea= new TextArea();
		 protected final Label progessLabel= new Label();
		 private volatile Runner curentrunner=null;
		 
		 private class Runner extends Thread
		 	{
			CompiledScript compiledScript = null;
			NgsFile<HEADERTYPE,ITEMTYPE> copyNgsFile=null;
			FileOutputStream fw=null;
			File saveAsFile=null;
			ByteArrayOutputStream saveToStringWriter=null;
			PrintStream out = null;
			SimpleBindings bindings=null;
			LogCloseableIterator iter=null;
			
			@Override
			public void run() {
				try 
					{
					this.bindings.put(HEADER_CONTEXT_KEY, this.copyNgsFile.getHeader());
					this.bindings.put(OUT_CONTEXT_KEY, this.out);
					this.bindings.put(ITER_CONTEXT_KEY, this.iter);
					
					this.compiledScript.eval(this.bindings);
					this.iter.close();
					this.out.flush();
					
					
					Platform.runLater(()->{
						AbstractAwkLike.this.progessLabel.setText("Done.");
						});
					
					if(this.iter.__stop_flag!=0)
						{
						//nothing
						}
					else if(this.out.checkError())
						{
						Platform.runLater(()->{
							final Alert alert = new Alert(
									AlertType.ERROR,
									"I/O Error. Check Stream limits in preferences.",
									ButtonType.OK);
							alert.showAndWait();
							});
						}
					else if(this.saveToStringWriter!=null) {
						final String output = new String(this.saveToStringWriter.toByteArray());
						Platform.runLater(()->{
							final Stage showResultStage=new Stage();
							showResultStage.setTitle("BioAlcidae");
							showResultStage.setScene(new Scene(new ScrollPane(new TextArea(output))));
							showResultStage.show();
							});
						} 
					else
						{
						Platform.runLater(()->{
							final Alert alert = new Alert(AlertType.CONFIRMATION, "Done", ButtonType.OK);
							alert.showAndWait();
							});
						}
					this.out.close();//we close after to avoid check error is user closed the stream
					this.out=null;
					this.iter=null;
					}
				catch(final Exception err)
					{
					Platform.runLater(()->{
						AbstractAwkLike.this.progessLabel.setText("Error :"+err.getMessage());
						});
					err.printStackTrace();
					}
				finally
					{
					AbstractAwkLike.this.curentrunner=null;
					dispose();
					}
				}
			
			void dispose()
				{
				CloserUtil.close(this.iter);
				CloserUtil.close(this.copyNgsFile);
				CloserUtil.close(this.fw);				
				}
			
		 	}
		 
		 AbstractAwkLike()
		 	{
			this.scriptArea.setPromptText("insert your javascript expression");
			this.scriptArea.setText("/** count even start */\nvar count=0;\nwhile(iter.hasNext()) {\n"
					+ " var item= iter.next();\n"
					+ " if(item.getStart()%2==0) count++;\n"
					+ "}\n"
					+ "out.println(count);\n"
					);
			this.scriptArea.setMaxWidth(Double.MAX_VALUE);
			this.scriptArea.setMaxHeight(Double.MAX_VALUE);
			this.scriptArea.setFont(Font.font("Courier",14));
			this.scriptArea.setPrefRowCount(100);
		 	}
		 
		 protected SimpleBindings completeBindings(final SimpleBindings sb,final HEADERTYPE h) {
			 return sb;
		 }
		 
		 protected TextFlow getHelpString() {
			 
			 return new TextFlow(new Text(
			     "This is a graphical interface to an instance of  "),createHyperLink("BioAlcidae","https://github.com/lindenb/jvarkit/wiki/BioAlcidae"),new Text(".\n"+
					 "The filter (including javascript) of the original window are **not** used.\n"+
					 "Use a your own risk. If your produce an infinite loop, the script might run forever and potentialy might feel your hard-drive.\n"+
					 "The script injects:\n* out : the output stream, a '"),javadocFor(PrintStream.class),new Text("'\n")
					 );
		 }
		
		private synchronized void stopRunner()
			{
			final Runner runner=this.curentrunner;
			if(runner!=null) {
				if(runner.iter!=null) runner.iter.__stop_flag=1;
				try { Thread.sleep(2000);} catch (InterruptedException e) {}
				try {runner.dispose();} catch (Exception e) {}
				try {runner.interrupt();} catch (Exception e) {}
				this.curentrunner=null;
				this.progessLabel.setText("Runner stopped.");
				}	
			}
		
		private void execute(final Stage ownerWindow,boolean saveToFile) 
			{
			AbstractAwkLike.this.progessLabel.setText("");
			if(this.curentrunner!=null)
				{
				JfxNgs.showExceptionDialog(ownerWindow, "Process is already running.");
				return;
				}
			
			if(!NgsStage.this.owner.javascriptCompiler.isPresent())
				{
				JfxNgs.showExceptionDialog(ownerWindow, "javascript is not supported");
				return;
				}
			final String expr=this.scriptArea.getText();
			if(expr.trim().isEmpty())
				{
				JfxNgs.showExceptionDialog(ownerWindow, "Empty expression");
				return;
				}
			
			final Runner newrunner=new Runner();

			try {
				newrunner.compiledScript = NgsStage.this.owner.javascriptCompiler.get().compile(expr);
				} 
			catch(final Throwable err)
				{
				JfxNgs.showExceptionDialog(ownerWindow,err);
				return;
				}
			if(saveToFile)
				{
				final FileChooser fc=NgsStage.this.owner.newFileChooser();
				fc.setInitialFileName("output.txt");
				fc.setSelectedExtensionFilter(new ExtensionFilter("text file", ".txt",".csv"));
				newrunner.saveAsFile =fc.showSaveDialog(ownerWindow);
				if(newrunner.saveAsFile==null) return;
				}
			else
				{
				newrunner.saveToStringWriter=new ByteArrayOutputStream();
				}
			
			try {
				final LimitSecurityStream limitSecurityStream;
				if( saveToFile)
					{
					long limit_bytes;
					try
						{
						limit_bytes=Long.parseLong(
								owner.preferences.get(
										owner.pref_bioalcidae_max_stream.key,
										owner.pref_bioalcidae_max_stream.defaultValue
										)
								);
						}
					catch(final Exception err)
						{
						err.printStackTrace();
						limit_bytes=Long.parseLong(owner.pref_bioalcidae_max_stream.defaultValue);
						}
					newrunner.fw = new FileOutputStream(newrunner.saveAsFile);
					limitSecurityStream = new LimitSecurityStream(newrunner.fw,limit_bytes);
					}
				else
					{
					long limit_bytes;
					try
						{
						limit_bytes=Long.parseLong(
								owner.preferences.get(
										owner.pref_bioalcidae_max_string.key,
										owner.pref_bioalcidae_max_string.defaultValue
										)
								);
						}
					catch(final Exception err)
						{
						err.printStackTrace();
						limit_bytes=Long.parseLong(owner.pref_bioalcidae_max_string.defaultValue);
						}
					newrunner.saveToStringWriter = new ByteArrayOutputStream();
					limitSecurityStream= new LimitSecurityStream(newrunner.saveToStringWriter,limit_bytes);
					}
				newrunner.out = new PrintStream(limitSecurityStream);
				newrunner.copyNgsFile = getNgsFile().reOpen();
				newrunner.iter = new LogCloseableIterator(newrunner.copyNgsFile.iterator(),this.progessLabel);
				newrunner.bindings= completeBindings(new  SimpleBindings(),newrunner.copyNgsFile.getHeader());

				this.curentrunner=newrunner;
				this.curentrunner.start();
				}
			catch(final Throwable err) {
				err.printStackTrace();
				newrunner.dispose();
				JfxNgs.showExceptionDialog(ownerWindow, err);
				}
			finally {
				
				}
			}
		 
		public void show() {
			final Stage dialog = new Stage();
			dialog.setResizable(true);
			dialog.setTitle("Bioalcidae for "+getNgsFile().getSource());
			
			final VBox contentPane=new VBox(5);
			VBox.setVgrow(contentPane, Priority.ALWAYS);
			
			final Menu fileMenu = new Menu("File");
			
			MenuItem menu=new MenuItem("Save script as...");
			menu.setOnAction(AE->{});
			fileMenu.getItems().add(menu);
			menu=new MenuItem("Load Script ...");
			menu.setOnAction(AE->{actionLoadScript(dialog, this.scriptArea);});
			fileMenu.getItems().add(menu);
			menu=new MenuItem("Validate Script ...");
			menu.setOnAction(AE->{actionValidateScript(dialog, this.scriptArea);});
			fileMenu.getItems().add(menu);

			
			menu=new MenuItem("Run To File");
			menu.setOnAction(AE->{execute(dialog,true);});
			fileMenu.getItems().add(menu);
			menu=new MenuItem("Run To Dialog");
			menu.setOnAction(AE->{execute(dialog,false);});
			fileMenu.getItems().add(menu);

			fileMenu.getItems().add(new SeparatorMenuItem());
			menu=new MenuItem("Close");
			menu.setOnAction(AE->{dialog.hide();});
			fileMenu.getItems().add(menu);
			
			/** snippets */
			
			final Menu snippetmenu=new Menu("Snippets");

			loadSnippets().stream().filter(C->C.isBioalcidaeScope()).forEach(
    				C->{
    					final MenuItem item=new MenuItem(C.label+
    							(C.function?"[Function]":"")
    							);
						item.setOnAction(C.handler(AbstractAwkLike.this.scriptArea) );
						snippetmenu.getItems().add(item);	
    				} );
		    		
			final MenuBar menuBar=new MenuBar(fileMenu,snippetmenu);
			contentPane.getChildren().add(menuBar);

			
			final FlowPane top=new FlowPane(5,5);
			Button button=new Button("Save Script as...");
	    	button.setOnAction(AE->{actionSaveScript(dialog, this.scriptArea);});
	    	top.getChildren().add(button);
	    	
	    	button=new Button("Open Script as ...");
	    	button.setOnAction(AE->{actionLoadScript(dialog, this.scriptArea);});
	    	top.getChildren().add(button);
	    	
	    	button=new Button("Validate");
	    	button.setTooltip(new Tooltip("Validate javascript syntax"));
	    	button.setOnAction(AE->{actionValidateScript(dialog, this.scriptArea);});
	    	top.getChildren().add(button);

	    	contentPane.getChildren().add(top);
	    	
			contentPane.getChildren().add(new Separator(Orientation.HORIZONTAL));
			final BorderPane scriptPane=new BorderPane(this.scriptArea);
			scriptPane.setPadding(new Insets(5));
			scriptPane.setMaxHeight(Double.MAX_VALUE);
			scriptPane.setMaxWidth(Double.MAX_VALUE);
			contentPane.getChildren().add(scriptPane);
			
			
			contentPane.getChildren().add(new Separator(Orientation.HORIZONTAL));			
			final FlowPane helpLabel=new FlowPane(getHelpString());
			contentPane.getChildren().add(helpLabel);
			
			final HBox buttonPane=new HBox(5);
			buttonPane.setPadding(new Insets(5));
			
			button= new Button("Run To File....");
			button.setTooltip(new Tooltip("Run the process and save the result in a new file"));
			button.setFont(Font.font(button.getFont().getFamily(),24));
			button.setTextFill(Color.GREEN);
			button.setOnAction(AE->{ execute(dialog,true);});
			buttonPane.getChildren().add(button);
			
			button= new Button("Run To Dialog....");
			button.setTextFill(Color.GREEN);
			button.setTooltip(new Tooltip("Run the process and display the result in a new window"));
			button.setFont(Font.font(button.getFont().getFamily(),24));
			button.setOnAction(AE->{ execute(dialog,false);});
			buttonPane.getChildren().add(button);
			
			
			button= new Button("Stop");
			button.setTextFill(Color.RED);
			button.setTooltip(new Tooltip("Stop current running process"));
			button.setFont(Font.font(button.getFont().getFamily(),24));
			button.setOnAction(AE->{stopRunner();});
			buttonPane.getChildren().add(button);
			
			contentPane.getChildren().add(buttonPane);
			
			contentPane.getChildren().add(this.progessLabel);
			
			dialog.setScene(new Scene(contentPane));
			
			dialog.initOwner(NgsStage.this);
			dialog.initModality(Modality.APPLICATION_MODAL); 
			dialog.showAndWait();
			}
	 	}
	 
	/** javascript filter used to filter a NgsFile*/
    protected static abstract class JavascriptFilter<HEADER,DATATYPE> {
    	/** the compiled script */
		protected final Optional<CompiledScript> compiledScript;
		/** bindings that will be injected in the script */
		protected final SimpleBindings bindings=new SimpleBindings();
		Optional<Throwable> encounteredException=Optional.empty();
		protected JavascriptFilter(
				final HEADER header,
				final Optional<CompiledScript> compiledScript)
			{
			this.bindings.put("header", header);
			this.compiledScript=compiledScript;
			}
		public abstract DATATYPE eval(DATATYPE v);
		
		 /** called by javascript filters */
	    protected boolean accept()
			{
	    	if(!this.compiledScript.isPresent()) return true;
			final Object result;
			try  {
				result = this.compiledScript.get().eval(this.bindings);
			} catch(final ScriptException err)
			{
				if(!this.encounteredException.isPresent())
					{
					LOG.severe(err.getMessage());
					err.printStackTrace();
					this.encounteredException = Optional.of(err);
					}
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
				if(!this.encounteredException.isPresent())
					{
					final String err="Script returned something that is not a boolean or a number:"+result.getClass();
					LOG.warning(err);
					this.encounteredException = Optional.of(new ScriptException(err));
					}
				return false;
				}
			return true;
			}
		}

	
	protected abstract class AbstractQualityStage
		extends Stage
		{
		protected abstract class ScanThread 
			extends Thread
			{
			/** number of items scanned so far */
			protected long nItems=0L;
			/** file source */
			protected final NgsFile<HEADERTYPE, ITEMTYPE> ngsReader;
			/** script for filtering */
			protected final Optional<CompiledScript> compiledScript;
			/** should we stop the scanning */
			protected volatile boolean kill_flag=false;
			/** time of last refresh */
    		protected long lastRefresh =System.currentTimeMillis();
    		/** active chart factory */
    		protected final ChartFactory<HEADERTYPE,ITEMTYPE> factory;
    		/**  other filters */
    		protected final Predicate<ITEMTYPE> otherFilters;
    		/** last seen exception */
    		protected Optional<Throwable> encounteredException=Optional.empty();

    		ScanThread(
    				final ChartFactory<HEADERTYPE,ITEMTYPE> factory,
    				final NgsFile<HEADERTYPE, ITEMTYPE> ngsReader,
    				final Optional<CompiledScript> compiledScript,
    				final Predicate<ITEMTYPE> otherFilters)
				{
				this.ngsReader=ngsReader;
				this.compiledScript = compiledScript;
				this.factory=factory;
				this.otherFilters= otherFilters;
				}
            
    		
            /** called at end */
            protected void atEnd() {
            	CloserUtil.close(this.ngsReader);
            	if(kill_flag) {
					LOG.warning("Thread was killed");
					}
				else
    				{
					/** repaint for the last time */
    				repaint();
    				Platform.runLater(new Runnable() {
        				 @Override
        				public void run() {
        					 AbstractQualityStage.this.countItemsLabel.setText(
        						"Done... Number of items: "+nItems+ (kill_flag?" [KILLED]":""));
        					 if(encounteredException.isPresent())
        					 	{
        						JfxNgs.showExceptionDialog(NgsStage.this, encounteredException.get());
        					 	}
        				 }
        			 	});
    				}
            	}
            protected void repaint()
				{
				final Chart chart= this.factory.build();
				 Platform.runLater(new Runnable() {
					 @Override
					public void run() {
						AbstractQualityStage.this.contentPane.setCenter(chart);
						AbstractQualityStage.this.countItemsLabel.setText("Running... Number of items: "+nItems);
					 	}
				 	});
				}
            protected void update()
				{
				long now  =System.currentTimeMillis(); 
				if( kill_flag || now - lastRefresh < REFRESH_SECOND*1000) return ;//5 seconds;
				lastRefresh = now;
				repaint();
				}
            protected void onError(final Throwable err)
            	{
            	CloserUtil.close(this.ngsReader);
            	LOG.severe(err.getMessage());
				Platform.runLater(new Runnable() {
    				 @Override
    				public void run() {
    					 AbstractQualityStage.this.countItemsLabel.setText(
    						"ERROR "+err.getMessage());
    				 	}
    			 	});
            	}
            
			}
		
		protected final ScanThread thread;
		protected final BorderPane contentPane=new BorderPane();
		protected final Label countItemsLabel=new Label();

		protected AbstractQualityStage(
				final ChartFactory<HEADERTYPE,ITEMTYPE> factory,
				final NgsFile<HEADERTYPE,ITEMTYPE> ngsReader,
				final Optional<CompiledScript> compiledScript,
				final Predicate<ITEMTYPE> otherFilters
				)
			{
			this.setTitle(factory.getName()+" : "+ngsReader.getSource());
			this.thread = createThread(factory,ngsReader,compiledScript,otherFilters);
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
        	this.contentPane.setCenter(factory.build());
        	this.contentPane.setBottom(this.countItemsLabel);
        	final Scene scene=new Scene(this.contentPane,1000,500);
        	this.setScene(scene);
			}
		/** create the thread that will scan the file in the background */
		protected abstract ScanThread createThread(
				final ChartFactory<HEADERTYPE,ITEMTYPE> factory,
				final NgsFile<HEADERTYPE, ITEMTYPE> ngsReader,
				final Optional<CompiledScript> compiledScript,
				final Predicate<ITEMTYPE> otherFilters
				);
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
	
	/** Constructor */
    protected NgsStage(
    		final JfxNgs owner,
    		final NgsFile<HEADERTYPE,ITEMTYPE> ngsFile
    		) throws IOException
    	{
    	this.owner= owner;
    	this.ngsFile= ngsFile;
    	this.setTitle(this.ngsFile.getSource());
    	this.maxItemsLimitSpinner.setEditable(true);
    	this.maxItemsLimitSpinner.setTooltip(new Tooltip(
    			"When Manually editing, press <RETURN> to commit the new value." +
    			"The whole file is NOT loaded, only a subset of data will be read."
    			));
    	if(ngsFile.getSequenceDictionary()==null)
    		{
    		throw new IOException("There is no associated dictionary for "+ngsFile);
    		}
    	this.chrom2start=new HashMap<>(ngsFile.getSequenceDictionary().size());
    	
    	// for faster result in canvas things, get the genomic index of each chromosome
    	{
    	long genomeLength=0L;
    	for(final SAMSequenceRecord ssr: ngsFile.getSequenceDictionary().getSequences())
	     	{
	     	this.chrom2start.put(ssr.getSequenceName(),genomeLength);
	     	genomeLength+=ssr.getSequenceLength();
	     	}
    	}
    	
    	this.seqDictionaryCanvas = new  SeqDictionaryCanvas(ngsFile.getSequenceDictionary());
    	
    	if(!this.owner.javascriptCompiler.isPresent()) {
    		this.javascriptArea.setEditable(false);
    		this.javascriptArea.setPromptText("Javascript engine is not available.");
    	} else
    		{
    		this.javascriptArea.setPromptText("Use this area to create a javascript-bases filter to ignore some items");
    		}
    	
    	
    	this.addEventHandler(
    			WindowEvent.WINDOW_SHOWING ,WE-> {
                        owner.registerStage(NgsStage.this);
                        NgsStage.this.reloadData();   
                    });
    	this.addEventHandler(
    			WindowEvent.WINDOW_CLOSE_REQUEST ,WE->{
                    	NgsStage.this.ngsFile.close();
                    	owner.unregisterStage(NgsStage.this);
                    });
       
        this.gotoField.setPrefColumnCount(25);
        this.gotoField.setEditable(true);
        
        
        this.fileMenu.getItems().addAll(owner.createCommonMenuItems(this));
        this.fileMenu.getItems().add(JfxNgs.createMenuItem("Save filtered data as...", ()->doMenuSaveAs()));

    	
        this.menuBar.getMenus().add(this.fileMenu);
    	this.menuBar.getMenus().add(this.statsMenu);
    	this.menuBar.getMenus().add(this.createJavascriptSnippetMenu());
    	
    	if(this.owner.javascriptCompiler.isPresent()) {
	    	this.menuBar.getMenus().add(this.bioalcidaeMenu);
	    	
	    	final MenuItem runBioalcidae = new MenuItem("Open Dialog...");
	    	runBioalcidae.setOnAction(AE->{invokeBioalcidae();});
	    	this.bioalcidaeMenu.getItems().add(runBioalcidae);
	    	this.bioalcidaeMenu.getItems().add(new SeparatorMenuItem());
	    	}
    	}
    
    @Override
    protected void finalize() throws Throwable {
    	this.ngsFile.close();
    	super.finalize();
    	}
    
    /** path to java script snippets */
    protected abstract String getSnippetResourcePath();
    
    /** generate a set of common button for handling javascript */
    protected List<Button> makeJavascriptButtons()
    	{
    	List<Button> buttons=new ArrayList<>();
    	if(!owner.javascriptCompiler.isPresent()) return buttons;
    	Button button= setTooltip( new Button("Save as..."),
    			"Save the Script below in a text file"
    			);
    	button.setOnAction(AE->{actionSaveScript(NgsStage.this, this.javascriptArea);});
    	buttons.add(button);
    	
    	button= setTooltip(new Button("Open..."),
    			"Open a javascript file  that will be used to filter the data");
    	button.setOnAction(AE->{actionLoadScript(NgsStage.this, this.javascriptArea);});
    	buttons.add(button);
    	
    	button=setTooltip(
    			new Button("Validate"),
    			"Validate the Javascript syntax of the script below"
    			);
    	button.setOnAction(AE->{actionValidateScript(NgsStage.this, this.javascriptArea);});
    	buttons.add(button);
   	
    	
    	return buttons;
    	}
    
    
    /** A simple snippet bean */
    private static class SnippetCode
    	{
    	String label="";
    	String code="";
    	boolean function=false;
    	String scope="";
    	
    	public boolean isFilterScope() { return !isBioalcidaeScope();}
    	public boolean isBioalcidaeScope() { return scope!=null && scope.equalsIgnoreCase("bioalcidae");}
    	public EventHandler<ActionEvent> handler(final TextArea textArea) {
    		return new EventHandler<ActionEvent>() {
				@Override
				public void handle(final ActionEvent event) {
					if(!function)
						{
						textArea.setText(code);
						}
					else
						{
						final int caret = textArea.getCaretPosition();
						textArea.insertText(caret,code);
						}	
				textArea.requestFocus();
				}
			};
    	}
    	}
    
    protected List<SnippetCode> loadSnippets()
    	{
    	final String rsrc = getSnippetResourcePath();
    	if(rsrc==null || rsrc.isEmpty()) return Collections.emptyList();
    	final List<SnippetCode> snippets = new ArrayList<>();
		InputStream in=null;
		XMLEventReader r=null;
		try
			{
			in = getClass().getResourceAsStream(rsrc);
			if(in!=null)
				{
				final XMLInputFactory xif=XMLInputFactory.newFactory();
				xif.setProperty(XMLInputFactory.IS_COALESCING,Boolean.TRUE);
				r=xif.createXMLEventReader(in);
				final QName isFunctionAtt=new QName("is-function");
				final QName scopeAtt=new QName("scope");
				final QName labelAtt=new QName("label");
				final QName nameAtt=new QName("name");
				while(r.hasNext())
					{
					final XMLEvent evt=r.nextEvent();
					if(!evt.isStartElement() ) continue;
					final StartElement start=evt.asStartElement();
					if(!start.getName().getLocalPart().equals("code")) continue;
					final Attribute isFunction =start.getAttributeByName(isFunctionAtt);
					final Attribute scope =start.getAttributeByName(scopeAtt);

					Attribute attLabel=start.getAttributeByName(labelAtt);
					if(attLabel==null) attLabel=start.getAttributeByName(nameAtt);
					if(attLabel!=null && r.hasNext() && r.peek().isCharacters())
						{
						final SnippetCode snippet=new SnippetCode();
						snippet.label = attLabel.getValue();
						snippet.code = r.nextEvent().asCharacters().getData();
						snippet.function = isFunction!=null && isFunction.getValue().equals("true");
						snippet.scope = (scope==null ?"":scope.getValue());
						snippets.add(snippet);
						}
					}
				}
			else
				{
				LOG.warning("Cannot read snippets "+rsrc);
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
		return snippets;
    	}
    
    /** generate the javascript Menu, containing the snippets.
     * the snippet are stored as a xml file in the jar file
     * @return
     */
    protected Menu createJavascriptSnippetMenu() {
    	final Menu menu=new Menu("Snippets");
    		loadSnippets().stream().filter(C->C.isFilterScope()).forEach(
    				C->{
    					final MenuItem item=new MenuItem(C.label+
    							(C.function?"[Function]":"")
    							);
						item.setOnAction(new EventHandler<ActionEvent>() {
							@Override
							public void handle(ActionEvent event) {
								if(!C.function)
									{
									NgsStage.this.javascriptArea.setText(C.code);
									}
								else
									{
									final int caret = NgsStage.this.javascriptArea.getCaretPosition();
									NgsStage.this.javascriptArea.insertText(caret, C.code);
									}
								
							}
						});
					menu.getItems().add(item);	
    				}
    				
    				);
    	
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
    	final MenuItem item=new MenuItem("Save table '"+tableName+"' as ...");
    	item.setOnAction(new EventHandler<ActionEvent>()
			{
			@Override
			public void handle(final ActionEvent event)
				{
				final FileChooser fc= owner.newFileChooser();
				fc.setSelectedExtensionFilter(new FileChooser.ExtensionFilter("TSV","tsv"));
				final File fout= owner.updateLastDir(fc.showSaveDialog(NgsStage.this));
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
    	updateStatusBar(AlertType.NONE,"");
    	if(feature==null) {
    		updateStatusBar(AlertType.WARNING,"No Feature was selected");
    		return;
    	}
    	openInIgv(
    			Collections.singletonList("goto "+feature.getContig()+":"+feature.getStart()+"-"+feature.getEnd())
    			);
    	}
    
    abstract void openInIgv();
    /** reload all data */
    abstract void reloadData();
    
    /** show stats */
    protected abstract void doMenuShowWholeStats(final ChartFactory<HEADERTYPE,ITEMTYPE> factory);
    /** show stats for whole file */
    protected final <T> void doMenuShowLocalStats(final ChartFactory<HEADERTYPE,T> factory,final Supplier<List<T>> data)
    	{
    	LOG.info("creating chart "+factory.getName());
    	factory.setHeader(getNgsFile().getHeader());
    	factory.setPedigree(getPedigree());
    	final List<T> L=data.get();
    	LOG.info("creating n items "+L.size());

    	for(final T o:L) factory.visit(o);
		final Chart chart=factory.build();
    	final BorderPane contentPane=new BorderPane(chart);
    	contentPane.setPadding(new Insets(10));
    	final Stage dialog = new Stage();
    	dialog.initOwner(this);
    	dialog.setTitle(factory.getName());
    	contentPane.setTop(new Label("Data for "+this.ngsFile.getSource()));
       	dialog.setScene(new Scene(contentPane));
    	contentPane.setBottom(new Label("Number of items: "+L.size()));
    	LOG.info("Showing chart");
    	dialog.show();
    	}
    
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
   
    
    /** build a table view for a Dictionary */
    protected Tab buildDictTab(final SAMSequenceDictionary dict)
        {
    	
		/* build dict Table */
		final TableView<SAMSequenceRecord> table=new TableView<>(
				dict==null?
				FXCollections.observableArrayList():
				FXCollections.observableArrayList(dict.getSequences())
				);
		
        table.getColumns().add(makeColumn("Name",SSR->SSR.getSequenceName()));
        table.getColumns().add(formatIntegerColumn(makeColumn("Length",SSR->SSR.getSequenceLength())));
        
        final Set<String> all_attributes=new HashSet<>();
        
       for(final SAMSequenceRecord ssr:dict.getSequences())
       		{
        	all_attributes.addAll(ssr.getAttributes().stream().
		        map(A->A.getKey()).
		        filter(S->!(S.equals("Name") || S.equals("Length"))).
		        collect(Collectors.toSet())
		        );
       		}
        for(final String key:all_attributes)
	        {
            if(dict.getSequences().stream().filter(S->S.getSpecies()!=null && !S.getSpecies().trim().isEmpty()).findAny().isPresent()) {
            	table.getColumns().add(makeColumn(key,SSR->SSR.getAttribute(key)));
            	}
	        }
        final Tab tab=new Tab("Dict", table);
        tab.setClosable(false);
        table.setPlaceholder(new Label("No Dictionary."));
        return tab;
        }
    /** save filtered Data As ... */
    protected abstract void doMenuSaveAs();
    
    
    /** called by main stage: set location box content and call reloadData */
    protected void moveTo(final String s)
    	{
    	this.gotoField.setText(s);
    	this.reloadData();
    	}
    
    private int parsePosition(String num) throws NumberFormatException
    	{
    	if( num==null || num.trim().isEmpty()) throw new NumberFormatException("bad number \""+num+"\"");
    	num = num.replace(",", "").trim();
    	if(num.toLowerCase().endsWith("bp"))
    		{
    		return Integer.parseInt(num.substring(0,num.length()-2))*1;
    		}
    	else if(num.toLowerCase().endsWith("kb"))
			{
			return Integer.parseInt(num.substring(0,num.length()-2))*1000;
			}
    	else if(num.toLowerCase().endsWith("mb"))
			{
			return Integer.parseInt(num.substring(0,num.length()-2))*1000000;
			}
    	else
    		{
    		return Integer.parseInt(num);
    		}
    	}
    
    protected Interval parseInterval(final String location)
    	{
    	updateStatusBar(AlertType.NONE,"");
    	final SAMSequenceDictionary dict=this.ngsFile.getSequenceDictionary();
		final String contig;
		int colon =location.indexOf(":");
		if(colon==-1)
			{
			contig=location.trim();
			}
		else
			{
			contig=location.substring(0,colon).trim();
			}
		
		SAMSequenceRecord ssr= dict.getSequence(contig);
		if(ssr==null)
			{
			ssr= dict.getSequence(JfxNgs.ContigToUCSC.apply(contig));
			}
		if(ssr==null)
			{
			ssr= dict.getSequence(JfxNgs.ContigToEnseml.apply(contig));
			}
		if(ssr==null)
			{
			updateStatusBar(AlertType.WARNING, "Cannot find contig in dictionary: "+location);
			return null;
			}
		
		if(colon!=-1)
			{
			final int hyphen = location.indexOf('-');
			final int plus = location.indexOf('+');
			
			if(hyphen!=-1 && plus!=-1) {
				updateStatusBar(AlertType.WARNING, "both '+' and '-' in "+location);
				return null;
				}
			
			Integer start=null,end=null;
			
			if(hyphen==-1 && plus==-1)
				{
				try {
					start= Math.max(0, parsePosition(location.substring(colon+1)));
					end=ssr.getSequenceLength();
					}
				catch(final NumberFormatException err ) {
					start=null;
					LOG.warning(location);
					updateStatusBar(AlertType.WARNING, "Bad Start in : "+location);
					return null;
					}
				}
			else
				{
				int delimidx=(hyphen==-1?plus:hyphen);
				try {
					start= Math.max(0, parsePosition(location.substring(colon+1,delimidx)));
					int num2= Math.min(
							parsePosition(location.substring(delimidx+1)),
							ssr.getSequenceLength()
							);
					if(plus!=-1) {
						int center = start;
						start=Math.max(center-num2,0);
						end = Math.min(center+num2,ssr.getSequenceLength());
						}
					else
						{
						end=num2;
						}
					}
				catch(final NumberFormatException err )
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
    
    protected <T extends  javafx.scene.control.Control> T setTooltip(final T control, final String text)
    	{
    	if(text!=null && !text.trim().isEmpty())
    		{
    		control.setTooltip(new Tooltip(text));
    		}
    	return control;
    	}
    /** create a button opening the current item in Broad/IGV */
    protected Button createIgvButton() {
        final Button igvButton =setTooltip(new Button("IGV"),
        		"Open the current selected item in Broad/Integrative Genome Viewer"
        		);
        igvButton.setOnAction(AE->openInIgv());
		return igvButton;
        }
    
    static <T>  TableColumn<T, Integer> formatIntegerColumn(final TableColumn<T, Integer> tc) {
    	tc.setCellFactory(tv ->new TableCell<T,Integer>() {
    	private java.text.DecimalFormat niceIntFormat = new java.text.DecimalFormat("###,###");
	    @Override
	    protected void updateItem(final Integer pos, boolean empty) {
	        super.updateItem(pos, empty);
	        setText(null);
	        if(pos==null)
	        	{
	        	setText(null);
	            setGraphic(null);
	        	}
	        else
	        	{
	        	setText(niceIntFormat.format(pos));
	        	 setGraphic(null);
	        	}
	        	}
	    	});
    	return tc;
    	}
    
    /** returns  the underlying instance of NgsFile (BamFile, VcfFile...) */
    protected NgsFile<HEADERTYPE, ITEMTYPE> getNgsFile() {
		return ngsFile;
	}
    
    /** returns the pedigree associated to this file */
    public PedFile getPedigree() {
    	return PedFile.getEmptyInstance();
    }
    
    /** returns the currently selected item , used for contextual menus */
    protected abstract Optional<ITEMTYPE> getCurrentSelectedItem(); 
    
    /** build context menu for current selected locatable */
    protected List<MenuItem> buildItemsForContextMenu() {
    	final List<MenuItem> L=new ArrayList<>();
    	MenuItem menuItem;
    	
    	menuItem=new MenuItem("Copy Interval in Clipboard");
		menuItem.setOnAction(AE->{
			Optional<ITEMTYPE> sel=getCurrentSelectedItem();
			if(!sel.isPresent() || sel.get().getContig()==null) return;
			final Clipboard clipboard = Clipboard.getSystemClipboard();
		    final ClipboardContent content = new ClipboardContent();
		    content.putString(sel.get().getContig()+":"+sel.get().getStart()+"-"+sel.get().getEnd());
		    clipboard.setContent(content);
		});

		L.add(menuItem);
    	
    	for(final String build:new String[]{"hg38","hg19","hg18"})
    		{
    		menuItem=new MenuItem("Open in UCSC "+build+" ... ");
    		menuItem.setOnAction(AE->{
    			final Optional<ITEMTYPE> sel=getCurrentSelectedItem();
    			if(!sel.isPresent() || sel.get().getContig()==null) return;
    			NgsStage.this.owner.getHostServices().showDocument(
    				"http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position="+
    						JfxNgs.ContigToUCSC.apply(sel.get().getContig())+"%3A"+sel.get().getStart()+"-"+sel.get().getEnd()
    					);
    		});
    		L.add(menuItem);
    		}
    	for(final String build:new String[]{"grch37.Homo_sapiens","www.Homo_sapiens","www.Rattus_norvegicus"})
			{
    		int dot=build.indexOf('.');
    		final String host=build.substring(0,dot);
    		final String org=build.substring(dot+1);
    		menuItem=new MenuItem("Open in Ensembl "+org+(host.equals("www")?"":"["+host+"]")+" ... ");
    		menuItem.setOnAction(AE->{
			
    		final Optional<ITEMTYPE> sel=getCurrentSelectedItem();
			if(!sel.isPresent() || sel.get().getContig()==null) return;
			
			NgsStage.this.owner.getHostServices().showDocument(
				"http://"+host +".ensembl.org/"+ org+"/Location/View?r="+
						JfxNgs.ContigToEnseml.apply(sel.get().getContig())+"%3A"+sel.get().getStart()+"-"+sel.get().getEnd()
					);
				});
			L.add(menuItem);
			}
    	
	    	for(final String database : new String[]{"Exac","gnomAD"})
	    	{
	    	menuItem=new MenuItem("Open Region in "+database+" ... ");
			menuItem.setOnAction(AE->{
				final Optional<ITEMTYPE> sel=getCurrentSelectedItem();
				if(!sel.isPresent() || sel.get().getContig()==null) return;
				NgsStage.this.owner.getHostServices().showDocument(
					"http://"+database.toLowerCase()+".broadinstitute.org/region/"+
							JfxNgs.ContigToEnseml.apply(sel.get().getContig())+"-"+sel.get().getStart()+"-"+sel.get().getEnd()
						);
			});
			L.add(menuItem);
	    	}
    	
    	return L;
    	}
    /** convert a contig pos to genomic index, return -1 on error (contig not found ) */
    protected long convertContigPosToGenomicIndex(final ContigPos cp) {
    	if(cp==null || cp.getContig()==null) return -1L;
    	return convertContigPosToGenomicIndex(cp.getContig(),cp.position);
    	}
    
    /** convert a contig pos to genomic index, return -1 on error (contig not found ) */
    protected long convertContigPosToGenomicIndex(final String contig,int pos) {
    	if(contig==null) return -1L;
    	final Long n= this.chrom2start.get(contig);
    	if(n==null) {
    		System.err.println("Cannot find "+contig);
    		return -1L;
    	}
    	final SAMSequenceRecord ssr=  getNgsFile().getSequenceDictionary().getSequence(contig);
    	if(ssr==null){
    		System.err.println("Cannot find "+contig);
    		return -1L;
    	}
    	return n.longValue() + Math.min( Math.max(1L,pos),ssr.getSequenceLength());
    	}
    
    /** convert a contig pos to genomic index, return null on error (contig not found ) */
    protected ContigPos convertGenomicIndexToContigPos(long gi) {
    	if(gi<0L) gi=0L;
    	long n=0L;
    	for(final SAMSequenceRecord ssr: getNgsFile().getSequenceDictionary().getSequences())
	    	{
	    	if(n<=gi && gi<=n+ssr.getSequenceLength())
    				{
    				return new ContigPos(ssr.getSequenceName(),(int)(gi-n));
    				}
	    	n+=ssr.getSequenceLength();
	    	}
    	return null;
    	}
    
    /** called by JfxNgs add an event handler to show a specific location when the window opens
     * @param s the location
     * @return 'this'
     */
    NgsStage<HEADERTYPE, ITEMTYPE> setLocationOnOpen(final String s)
    	{
    	if(!(s==null || s.trim().isEmpty()))
    		{
    		this.addEventHandler(
        			WindowEvent.WINDOW_SHOWING ,WE-> {
        					NgsStage.this.gotoField.setText(s);
                            NgsStage.this.reloadData();   
                        });
    		}
    	return this;
    	}
    
    /** invoke bioalcidae */
    abstract void invokeBioalcidae();
    
    
    private void actionLoadScript(final Stage dialog,final TextArea scriptArea) {
		final FileChooser fc= NgsStage.this.owner.newFileChooser();
		fc.getExtensionFilters().addAll(JS_EXTENSION_FILTERS);
		final File js=NgsStage.this.owner.updateLastDir(fc.showOpenDialog(dialog));
		if(js==null) return;
		InputStream in=null;
		try 
			{
			in= new FileInputStream(js);
			scriptArea.setText(IOUtil.readFully(in));
			in.close();
			}
		catch(final Exception err)
			{
			JfxNgs.showExceptionDialog(dialog, err);
			}
		finally
			{
			CloserUtil.close(in);
			}
    	}
    private void actionSaveScript(final Stage dialog,final TextArea scriptArea) {
		final FileChooser fc= this.owner.newFileChooser();
		fc.getExtensionFilters().addAll(JS_EXTENSION_FILTERS);
		fc.setInitialFileName("script.js");
		final File js= this.owner.updateLastDir(fc.showSaveDialog(dialog));
		if(js==null) return;
		PrintWriter pw=null;
		try 
			{
			pw=new PrintWriter(js);
			pw.write(scriptArea.getText());
			pw.flush();
			pw.close();
			}
		catch(final Exception err)
			{
			JfxNgs.showExceptionDialog(dialog, err);
			}
		finally
			{
			CloserUtil.close(pw);
			}
    	}
    
    private void actionValidateScript(final Stage dialog,final TextArea scriptArea) {
    	if(!NgsStage.this.owner.javascriptCompiler.isPresent()) return;
    	try {
			NgsStage.this.owner.javascriptCompiler.get().compile(scriptArea.getText());
			final Alert alert=new Alert(AlertType.CONFIRMATION);
			alert.setAlertType(AlertType.CONFIRMATION);
			alert.setTitle("OK");
			alert.setContentText("OK. At least, the script is compilable.");
			alert.showAndWait();
			}
		catch(final Exception err)
			{
			JfxNgs.showExceptionDialog(dialog, err);
			}
			
		}
    protected Hyperlink createHyperLink(final String label,final String  url)
    	{
    	final Hyperlink a = new Hyperlink(label==null || label.isEmpty()?url:label);
    	a.setTooltip(new Tooltip(url));
    	a.setOnAction(event -> {
    		NgsStage.this.owner.getHostServices().showDocument(url);
    	});
    	return a;
    	}
    
    /** create a hyperlink to javadoc, used in 'help' sections of javascript */
    protected Hyperlink javadocFor(final Class<?> clazz) {
		String name=clazz.getName();

		int dollar=name.indexOf("$");
		if(dollar!=-1) name=name.substring(0,dollar);
		String url;
		
		if(name.contains("htsjdk"))
			{	
			url="https://samtools.github.io/htsjdk/javadoc/htsjdk/"+name.replace(".", "/")+".html";
			}
		else if(name.startsWith("java"))
			{	
			url="https://docs.oracle.com/javase/8/docs/api/"+name.replace(".", "/")+".html";
			}
		else
			{
			url="https://github.com/lindenb/jvarkit/blob/master/src/main/java/"+name.replace(".", "/")+".java";
			}
    	return createHyperLink(clazz.getName(),url);
    	}
    
	}
