/*
The MIT License (MIT)

Copyright (c) 2017 Pierre Lindenbaum

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

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;
import java.util.function.Function;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.prefs.BackingStoreException;
import java.util.prefs.Preferences;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import javax.script.Bindings;
import javax.script.Compilable;
import javax.script.CompiledScript;
import javax.script.ScriptEngineManager;
import javax.script.ScriptException;
import javax.script.SimpleBindings;

import com.github.lindenb.jvarkit.util.Hershey;
import com.github.lindenb.jvarkit.util.igv.IgvSocket;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecord.SAMTagAndValue;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import javafx.application.Application;
import javafx.application.Platform;
import javafx.beans.property.ReadOnlyObjectWrapper;
import javafx.beans.value.ChangeListener;
import javafx.beans.value.ObservableValue;
import javafx.collections.FXCollections;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.geometry.Insets;
import javafx.geometry.Orientation;
import javafx.geometry.Rectangle2D;
import javafx.scene.Scene;
import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.control.Alert;
import javafx.scene.control.Alert.AlertType;
import javafx.scene.control.Button;
import javafx.scene.control.ButtonType;
import javafx.scene.control.CheckBox;
import javafx.scene.control.CheckMenuItem;
import javafx.scene.control.Label;
import javafx.scene.control.Menu;
import javafx.scene.control.MenuBar;
import javafx.scene.control.MenuItem;
import javafx.scene.control.ScrollBar;
import javafx.scene.control.ScrollPane;
import javafx.scene.control.Separator;
import javafx.scene.control.Spinner;
import javafx.scene.control.Tab;
import javafx.scene.control.TabPane;
import javafx.scene.control.TableColumn;
import javafx.scene.control.TableColumn.CellDataFeatures;
import javafx.scene.control.TableView;
import javafx.scene.control.TextArea;
import javafx.scene.control.TextField;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.FlowPane;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.Priority;
import javafx.scene.layout.VBox;
import javafx.scene.paint.Color;
import javafx.stage.FileChooser;
import javafx.stage.FileChooser.ExtensionFilter;
import javafx.stage.Screen;
import javafx.stage.Stage;
import javafx.stage.Window;
import javafx.stage.WindowEvent;
import javafx.util.Callback;
/**
 * GUI viewer for NGS files (BAM, VCF... )
 * @author lindenb
 *
 */
public class JfxNgs extends Application {
    private static final Logger LOG= Logger.getLogger("JfxNgs");
    private final Preferences preferences ;
    private final Compilable javascriptEngine;
    private static final int DEFAULT_BAM_RECORDS_COUNT=Integer.parseInt(System.getProperty("jxf.ngs.default.sam", "10000"));
    private static final int DEFAULT_VCF_RECORDS_COUNT=Integer.parseInt(System.getProperty("jxf.ngs.default.vcf", "1000"));
    private final List<StageContent> all_opened_stages=new ArrayList<>();
    private final Random randomMoveWindow=new Random();
    
    /** abstract base class for NGS window */
    private abstract class StageContent
        extends Stage
        {
    	/** javascript filtering */
    	protected final TextArea javascriptArea=new TextArea();
    	/** message stuff */
    	protected final Label messageLabel=new Label();
    	/** message stuff */
    	protected final TextField gotoField=new TextField();

    	
        public StageContent()
        	{
        	if(JfxNgs.this.javascriptEngine==null) {
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
                            JfxNgs.this.registerStage(StageContent.this);
                            StageContent.this.reloadData();
                            }
                        });
        	this.addEventHandler(
        			WindowEvent.WINDOW_CLOSE_REQUEST ,new EventHandler<WindowEvent>() {
                        @Override
                        public void handle(WindowEvent event) {
                        	StageContent.this.closeNgsResource();
                        	JfxNgs.this.unregisterStage(StageContent.this);
                            }
                        });
           
            this.gotoField.setPrefColumnCount(15);
            }
        @Override
        protected void finalize() throws Throwable {
        	closeNgsResource();
        	super.finalize();
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
        abstract void reloadData();
        
        
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

    /** describe the state of a SamFlag */
    private static class SamFlagRow
		{
    	private final SAMRecord record;
    	private final SAMFlag flag ;
    	SamFlagRow(final SAMRecord record,final SAMFlag flag)
    		{
    		this.record=record;
    		this.flag=flag;
    		}
		}
    
    /** describe the base of a read in the cigar context */
    private static class CigarAndBase
		{
    	private final String ref;
    	private final CigarOperator op;
    	private final Integer posInRead ;
    	private final Integer posInRef ;
    	private final Integer count;
    	private final Byte base ;
    	CigarAndBase(final String ref,final CigarOperator op,final Integer posInRead,final Integer posInRef,Integer count,Byte base)
    		{
    		this.ref=ref;
    		this.op=op;
    		this.posInRead = posInRead;
    		this.posInRef = posInRef;
    		this.base = base;
    		this.count=count;
    		}
		}

    private static class ContigPos
    	implements Comparable<ContigPos>
		{
		final String contig;
		final int position;
		ContigPos(final String contig,final int position) {
			this.contig=contig;
			this.position=position;
		}
		@Override
		public int compareTo(ContigPos o) {
			int i=contig.compareTo(o.contig);
			if(i!=0) return i;
			return position-o.position;
			}
		}
    private static class Pileup extends ContigPos
    	{
    	
    	final int count[]=new int[5];
    	final StringBuilder seq=new StringBuilder();
    	final StringBuilder qual=new StringBuilder();
    	final StringBuilder operators=new StringBuilder();
    	
    	Pileup(final String contig,final int position) {
			super(contig,position);
			Arrays.fill(count,0);
			}
    	int depth() { return count[0]+count[1]+count[2]+count[3]+count[4];}
    	}
    
    /** NGS window for BAM */
    private class BamStageContent extends StageContent
        {
        private final SamReader samReader;
        private final TableView<SAMRecord> recordTable;
        private final TableView<SamFlagRow> flagsTable;
        private final TableView<SAMTagAndValue> metaDataTable;
        private final TableView<CigarAndBase> cigarTable;
        private final TableView<Pileup> pileupTable;
        private final Map<SAMFlag,CheckMenuItem> flag2filterInMenuItem=new HashMap<>();
        private final Map<SAMFlag,CheckMenuItem> flag2filterOutMenuItem=new HashMap<>();
        private final Spinner<Integer> maxReadLimitSpinner;
        private final Canvas canvas = new Canvas(900, 600);
        private final ScrollBar canvasScrollV = new ScrollBar();
        private final CheckBox canvasShowReadName = new CheckBox("Show Read Name");
        
        BamStageContent(final String url) throws IOException
        	{
        	this.setTitle(url);
            final SamReaderFactory srf= SamReaderFactory.makeDefault();
            srf.validationStringency(Level.OFF.equals(LOG.getLevel())?
            		ValidationStringency.SILENT:
            		ValidationStringency.LENIENT
            		);
            LOG.info("Opening "+url);
            
            
            this.samReader=srf.open(SamInputResource.of(url));
            if(!this.samReader.hasIndex())
            	{
            	this.samReader.close();
            	throw new IOException("Bam without index "+url);
            	}

            /** Build menu for SAM Flags */
            for(final SAMFlag flg:SAMFlag.values())
            	{
            	flag2filterInMenuItem.put(flg,new CheckMenuItem("Filter In "+flg.name()));
            	flag2filterOutMenuItem.put(flg,new CheckMenuItem("Filter Out "+flg.name()));
            	}
            final Menu fileMenu=new Menu("File");
            fileMenu.getItems().add(createMenuItem("Open",new Runnable() {
				@Override
				public void run() {
					openNgsFiles(BamStageContent.this);
				}
			}));
            fileMenu.getItems().add(createMenuItem("Close",new Runnable() {
				@Override
				public void run() {
					BamStageContent.this.hide();
				}
			}));
           
            //selectFlagMenu.getItems().addAll(flag2filterOutMenuItem.values());
            final MenuBar menuBar=new MenuBar(fileMenu);
            final VBox vbox1 = new VBox();
            vbox1.getChildren().add(menuBar);
            
            FlowPane top1= new FlowPane(5,5);
            top1.setPadding(new Insets(10, 10, 10, 10));
            vbox1.getChildren().add(top1);
            top1.getChildren().add(new Label("GoTo:"));
            top1.getChildren().add(this.gotoField);
            final Button gotoButton=new Button("Go");
            gotoButton.setOnAction(new EventHandler<ActionEvent>()
				{
				@Override
				public void handle(ActionEvent event)
					{
					reloadData();
					}
				});
            top1.getChildren().add(gotoButton);
            top1.getChildren().add(new Separator(Orientation.VERTICAL));
            top1.getChildren().add(new Label("Limit:"));
            top1.getChildren().add(this.maxReadLimitSpinner=new Spinner<Integer>(0,100000,JfxNgs.DEFAULT_BAM_RECORDS_COUNT));
            top1.getChildren().add(new Separator(Orientation.VERTICAL));
            List<CheckMenuItem> menuFlags=new ArrayList<>(flag2filterOutMenuItem.values());
            menuFlags.addAll(flag2filterInMenuItem.values());
            CheckMenuItem tmp[]=new CheckMenuItem[menuFlags.size()];
            menuFlags.toArray(tmp);
            top1.getChildren().add(new MenuBar(new Menu("Flags",null,tmp)));
            this.maxReadLimitSpinner.setEditable(true);
            this.gotoField.setOnAction(new EventHandler<ActionEvent>()
				{
				@Override
				public void handle(ActionEvent event)
					{
					reloadData();
					}
				});
            final Button igvButton =new Button("IGV");
            top1.getChildren().add(igvButton);
            igvButton.setOnAction(new EventHandler<ActionEvent>() {
				@Override
				public void handle(ActionEvent event) {
					openInIgv();
				}
			});
            
            
            TabPane tabbedPane = new TabPane();
            tabbedPane.setPadding(new Insets(10, 10, 10, 10));
            Tab tab= new Tab("Reads");
            tab.setClosable(false);
            tabbedPane.getTabs().add(tab);
            
            this.recordTable = new TableView<>();
            /** create columns */
            
            this.recordTable.getColumns().add(makeColumn("Read-Name",REC->REC.getReadName()));
            this.recordTable.getColumns().add(makeColumn("Flag",REC->REC.getFlags()));
            this.recordTable.getColumns().add(makeColumn("Ref",REC->REC.getReferenceName()));
            this.recordTable.getColumns().add(makeColumn("Read-Pos",REC->REC.getAlignmentStart()));
            this.recordTable.getColumns().add(makeColumn("MAPQ",REC->REC.getMappingQuality()));
            this.recordTable.getColumns().add(makeColumn("CIGAR",REC->REC.getCigarString()));
            this.recordTable.getColumns().add(makeColumn("LEN",REC->REC.getInferredInsertSize()));
            this.recordTable.getColumns().add(makeColumn("Mate-Ref",REC->REC.getMateReferenceName()));
            this.recordTable.getColumns().add(makeColumn("Mate-Pos",REC->REC.getMateAlignmentStart()));
            this.recordTable.getColumns().add(makeColumn("SEQ",REC->REC.getReadString()));
            this.recordTable.getColumns().add(makeColumn("QUAL",REC->REC.getBaseQualityString()));
            
            VBox borderPane = new VBox();
            borderPane.setPadding(new Insets(10, 10, 10, 10));
           
  
            
            //ScrollPane scroll = new ScrollPane(this.recordTable);
            //scroll.setFitToHeight(true);
            //scroll.setFitToWidth(true);
            
            borderPane.getChildren().add(this.recordTable);
            
            GridPane tilePane = new GridPane();
            tilePane.setPadding(new Insets(10, 10, 10, 10));
            tilePane.setVgap(4);
            tilePane.setHgap(4);
           
            
            borderPane.getChildren().add(tilePane);
            
            /* define SAM Flag table */
            this.flagsTable= createSamFlagTable();
           
            
            //scroll.setFitToHeight(true);
            GridPane.setConstraints(  this.flagsTable,1, 1); // column=1 row=1
            tilePane.getChildren().add( this.flagsTable);
            
            
            /* define Meta Data table */
            this.metaDataTable = createMetaDataTable();

            GridPane.setConstraints( this.metaDataTable,2, 1); // column=2 row=1
            tilePane.getChildren().add(this.metaDataTable);

            
            /* build the cigar table */
            this.cigarTable = createCigarTable();
            
          
            //scroll = new ScrollPane();
            //scroll.setFitToHeight(true);
            //scroll.setFitToWidth(true);
            GridPane.setConstraints( this.cigarTable,3, 1); // column=3 row=1
            tilePane.getChildren().add(this.cigarTable);

            
            this.pileupTable = createPileupTable();
            tilePane.getChildren().add(this.cigarTable);
            
            /* when a read is selected update the flagsTable and metaDataTable */
            this.recordTable.getSelectionModel().selectedItemProperty().addListener((obs, oldSelection, newSelection) -> {
                if(newSelection==null)
                	{
                	flagsTable.getItems().clear();
                	metaDataTable.getItems().clear();
                	cigarTable.getItems().clear();
                	}
                else
                	{
                	final List<SamFlagRow> L=new ArrayList<>();
                	for(final SAMFlag flag: SAMFlag.values())
                		{
                		L.add(new SamFlagRow(newSelection,flag));
                		}
                	flagsTable.getItems().setAll(L);
                	
                	/* update meta data */
                	metaDataTable.getItems().setAll(newSelection.getAttributes());
                	
                	if(!newSelection.getReadUnmappedFlag() && newSelection.getCigar()!=null)
                		{
                		final List<CigarAndBase> M = new ArrayList<>();
                		int posInRead=0;
                		int posInRef=newSelection.getUnclippedStart();
                		final byte readString[] = newSelection.getReadBases();
                		for(CigarElement ce: newSelection.getCigar())
                			{
                			final CigarOperator op= ce.getOperator();
                			
                			switch(op)
                				{
                				case H: case D: case N:
                					{
                					M.add(new CigarAndBase(newSelection.getReferenceName(), op, null, posInRef,ce.getLength(),null));
                					posInRef+=ce.getLength();
                					break;
                					}
                				case P:
                					{
                					M.add(new CigarAndBase(newSelection.getReferenceName(), op, null, null,ce.getLength(),null));
                					break;
                					}
                				case I: 
	            					{
	            					for(int i=0;i< ce.getLength();++i)
	    	                			{
	                					M.add(new CigarAndBase(null, op, posInRead, null,1,
	                							readString==null || posInRead>=readString.length?null:readString[posInRead]
	                							));
	                					posInRead++;
	    	                			}
	            					break;
	            					}
                				case M: case X: case EQ: case S: 
                					{
                					for(int i=0;i< ce.getLength();++i)
        	                			{
	                					M.add(new CigarAndBase(newSelection.getReferenceName(), op, posInRead, posInRef,1,
	                							readString==null || posInRead>=readString.length?null:readString[posInRead]
	                							));
	                					posInRead++;
	                					posInRef++;
        	                			}
                					break;
                					}
                				}
	                		}
                			
                		cigarTable.getItems().setAll(M);
                		}
                	else
                		{
                		cigarTable.getItems().clear();
                		}
                	}
            });

            
            tab.setContent(borderPane);
           
            
            tab=new Tab("Header");
            tab.setClosable(false);
            tabbedPane.getTabs().add(tab);
            final StringWriter headerTextBuffer = new StringWriter();
            new SAMTextHeaderCodec().encode(headerTextBuffer, this.samReader.getFileHeader());
            final TextArea textAreaHeader =new TextArea(headerTextBuffer.toString());
            textAreaHeader.setEditable(false);
            
            ScrollPane scroll = new ScrollPane(textAreaHeader);
            scroll.setFitToHeight(true);
            scroll.setFitToWidth(true);
            tab.setContent(scroll);
            
            tabbedPane.getTabs().add(buildDictTab( this.samReader.getFileHeader().getSequenceDictionary()));
            tabbedPane.getTabs().add(createReadGroupPane(this.samReader.getFileHeader()));
            tabbedPane.getTabs().add(createProgramRecordPane(this.samReader.getFileHeader()));
            tabbedPane.getTabs().add(buildJavascriptPane());
            
            /* CANVAS STUFFF */
            final BorderPane canvasPane = new BorderPane(this.canvas);
            this.canvasScrollV.setOrientation(Orientation.VERTICAL);
            canvasPane.setRight(this.canvasScrollV);
            final FlowPane canvasTop=new FlowPane(this.canvasShowReadName);
            canvasPane.setTop(canvasTop);
            
            this.canvasShowReadName.setOnAction(new EventHandler<ActionEvent>() {
				@Override
				public void handle(ActionEvent event) {
					repaintCanvas();
				}
			});
            this.canvasScrollV.valueProperty().addListener(new ChangeListener<Number>() {
                public void changed(ObservableValue<? extends Number> ov, Number old_val, Number new_val) {
                	repaintCanvas();
                }
            });
            
            tab=new Tab("Canvas",canvasPane);
            tab.setClosable(false);
            tabbedPane.getTabs().add(tab);
            /* END CANVAS STUFF */
            
            vbox1.getChildren().add(tabbedPane);
            
            final FlowPane bottom=new FlowPane(super.messageLabel);
            vbox1.getChildren().add(bottom);
            
            this.setScene(new Scene(vbox1,
            		1000+randomMoveWindow.nextInt(10),
            		800+randomMoveWindow.nextInt(10)
            		));
            
            
        }

        
        @Override
        protected SAMSequenceDictionary getSAMSequenceDictionary() {
        	return this.samReader.getFileHeader().getSequenceDictionary();
        	}
        
        @Override
        void closeNgsResource()
        	{
        	LOG.info("closing"+samReader.getResourceDescription());
            CloserUtil.close(samReader);
        	}
        
        /** get a stream of read we can display on canvas */
        private Stream<SAMRecord> getDisplayableSamRecordStream()
        	{
        	return this.recordTable.getItems().
        			stream().
        			filter(R->!R.getReadUnmappedFlag() && R.getCigar()!=null)
        			;
        	}
        
        private void repaintCanvas()
        	{
        	final boolean showReadName = this.canvasShowReadName.isSelected();
        	final int baseSize=15;
        	final double canvaswidth= this.canvas.getWidth();
        	final double canvasheight= this.canvas.getHeight();
        	final GraphicsContext gc=this.canvas.getGraphicsContext2D();
        	gc.setFill(Color.WHITE);
        	gc.fillRect(0, 0, canvaswidth, canvasheight);
        	double y=baseSize*2;
        	final List<SAMRecord> records=getDisplayableSamRecordStream().collect(Collectors.toList());
        	if(records.isEmpty()) return;
     
        	final int recordStart=(int)this.canvasScrollV.getValue();
        	if(recordStart>=records.size()) return;
        	int recordIndex=recordStart;
        	final int chromStart=records.get(recordStart).getUnclippedStart();
        	final int chromLen=(int)(canvaswidth/baseSize);
        	if(chromLen==0) return;
        	
        	
        	Function<Integer,Double> position2pixel=new Function<Integer, Double>() {
				@Override
				public Double apply(Integer pos) {
					return ((pos-(double)chromStart)/(double)chromLen)*canvaswidth;
				}
			};
			
        	

			
			final Hershey hershey=new Hershey();
			
			hershey.paint(gc,records.get(recordStart).getContig(),1,0,baseSize*records.get(recordStart).getContig().length(),baseSize-2);
			
			for(int x=chromStart;x<chromStart+chromLen;++x)
				{
				double px=position2pixel.apply(x);
				gc.setStroke(x%10==0?Color.BLACK:Color.GRAY);
				gc.setLineWidth(x%10==0?5:0.5);
				gc.strokeLine(px, baseSize, px, canvasheight);
				if(x%10==0)
					{
					String s=String.valueOf(x);
					gc.setLineWidth(1.0);
					hershey.paint(gc,s,px,baseSize,baseSize*s.length(),baseSize-1);
					}
				}
			gc.setLineWidth(1);
        	while(y<canvasheight && recordIndex < records.size())
        		{
        		final SAMRecord rec=records.get(recordIndex);
        		if(!rec.getReferenceName().equals(records.get(recordStart).getReferenceName())) {
        			++recordIndex;
        			continue;
        			}
        		int baseIndex=0;
        		int refIndex=rec.getUnclippedStart();
        		final byte bases[]=rec.getReadBases();
        		
        		final Function<Integer,String> getBaseAt = new Function<Integer, String>() {
					@Override
					public String apply(Integer readPos) {
						char c;
						if(showReadName)
							{
							if(rec.getReadNameLength()<=readPos) return "";
							c= rec.getReadName().charAt(readPos);
							}
						else if(bases==null || bases.length<=readPos)
							{
							return "";
							}
						else
							{
							c=(char)bases[readPos];
							}	
						c=(rec.getReadNegativeStrandFlag()?
								Character.toLowerCase(c):
								Character.toUpperCase(c)
								);
						return String.valueOf(c);
					}
				};
				
        		final Function<Integer,Color> getColorAt = new Function<Integer, Color>() {
					
					public Color apply(Integer readPos) {
						if(bases==null || bases.length<=readPos)
							{
							return Color.BLACK;
							}
						switch(Character.toUpperCase((char)bases[readPos]))
							{
							case 'A': return Color.BLUE;
							case 'T': return Color.GREEN;
							case 'C': return Color.YELLOW;
							case 'G': return Color.RED;
							default: return Color.BLACK;
							}
					}
				};

				
				
				gc.setLineWidth(1.0);
        		//arrow end
    			{
    			double endpos=position2pixel.apply(rec.getReadNegativeStrandFlag()?rec.getUnclippedStart():rec.getUnclippedEnd()+1);
    			double radius=baseSize/4.0;
    			gc.setFill(Color.BLACK);
    			gc.fillOval(
    					endpos- radius,
        				y+baseSize/2.0 - radius,
        				radius*2 ,
        				radius*2
        				);
    			}

        		
        		final Set<Integer> referenceEvents=new HashSet<>();
        		for(CigarElement ce:rec.getCigar()) {
        			switch(ce.getOperator())
        				{
        				case P: break;
        				case I: 
        					{
        					baseIndex+=ce.getLength();
        					referenceEvents.add(refIndex);
        					break;
        					}
        				case D: case N:
        					{
        					gc.setFill(Color.RED);
        					for(int x=0;x< ce.getLength();++x)
        						{
        						gc.fillRect(position2pixel.apply(refIndex),y,baseSize,baseSize-1);
        						refIndex++;
        						}
        					break;
        					}
        				case H:
        					{
    						gc.setFill(Color.YELLOW);
        					for(int x=0;x< ce.getLength();++x)
        						{
        						gc.fillRect(position2pixel.apply(refIndex),y,baseSize,baseSize-1);
        						refIndex++;
        						}
        					break;
        					}
        				case S:
	    					{
	    					for(int x=0;x< ce.getLength();++x)
	    						{
	    						gc.setFill(Color.YELLOW);
	    						gc.fillRect(position2pixel.apply(refIndex),y,baseSize,baseSize-1);
	    						gc.setStroke(getColorAt.apply(baseIndex));
	    						hershey.paint(gc,getBaseAt.apply(baseIndex), position2pixel.apply(refIndex),y,baseSize-1,baseSize-2);
	    						refIndex++;
	    						baseIndex++;
	    						}
	    					break;
	    					}
        				case EQ:case X:case M:
    						{
	    					for(int x=0;x< ce.getLength();++x)
	    						{
        						gc.setFill(ce.getOperator()==CigarOperator.X?Color.RED:Color.LIGHTGRAY);
	    						gc.fillRect(position2pixel.apply(refIndex),y,baseSize,baseSize-1);
	    						gc.setStroke(getColorAt.apply(baseIndex));
	    						hershey.paint(gc,getBaseAt.apply(baseIndex), position2pixel.apply(refIndex),y,baseSize-1,baseSize-2);
	    						refIndex++;
	    						baseIndex++;
	    						}
        					break;
        					}
        				
        				default:break;
        				}
        			if(refIndex> chromStart+chromLen) break;
        			}
        		
        		
        		gc.setStroke(Color.BLACK);
        		gc.strokeRect(
        				position2pixel.apply(rec.getUnclippedStart()),
        				y,
        				position2pixel.apply(rec.getUnclippedEnd()+1)-position2pixel.apply(rec.getUnclippedStart()) ,
        				baseSize-1
        				);
        		
        		for(Integer pos:referenceEvents)
        			{
        			double x=position2pixel.apply(pos);
        			gc.setStroke(Color.RED);
        			gc.setLineWidth(0.5);
        			gc.strokeLine(x, y, x, y+baseSize);
        			}
        		
        		recordIndex++;
        		y+=baseSize;
        		}
        	gc.setStroke(Color.BLACK);
    		gc.rect(0,0,canvaswidth-1,canvasheight-1);
        	}
        
        private TableView<SamFlagRow> createSamFlagTable()
	        {
	        final TableView<SamFlagRow> table=new TableView<>();
	    	table.getColumns().add(makeColumn("Flag", O->O.flag.name()));
	    	table.getColumns().add(makeColumn("Status",new Function<SamFlagRow,Boolean>() {
	    		@Override
	    		public Boolean apply(final SamFlagRow param) {
	    			return param.flag.isSet(param.record.getFlags());
	    			}
				} ));
	    	
	        return table;
	        }        
        
        private TableView<SAMTagAndValue> createMetaDataTable()
	        {
	        final TableView<SAMTagAndValue> table=new TableView<>();
	    	table.getColumns().add(makeColumn("Key", O->O.tag));
	    	table.getColumns().add(makeColumn("Value", O->O.value));
	        return table;
	        }
        
        
        private TableView<Pileup> createPileupTable()
        	{
        	final TableView<Pileup> table=new TableView<>();
        	table.getColumns().add(makeColumn("REF", O->O.contig));
        	table.getColumns().add(makeColumn("POS", O->O.position));
        	table.getColumns().add(makeColumn("Depth", O->O.depth()));
        	table.getColumns().add(makeColumn("A", O->O.count[0]));
        	table.getColumns().add(makeColumn("T", O->O.count[1]));
        	table.getColumns().add(makeColumn("G", O->O.count[2]));
        	table.getColumns().add(makeColumn("C", O->O.count[3]));
        	table.getColumns().add(makeColumn("N", O->O.count[4]));
        	table.getColumns().add(makeColumn("Bases", O->O.seq.toString()));
        	table.getColumns().add(makeColumn("Qual", O->O.qual.toString()));
        	table.getColumns().add(makeColumn("Operators", O->O.operators.toString()));
        	return table;
        	}
        
        private TableView<CigarAndBase> createCigarTable() 
        	{
        	final TableView<CigarAndBase> table=new TableView<>();
        	table.getColumns().add(makeColumn("REF", O->O.ref));
        	table.getColumns().add(makeColumn("Read-Pos", O->O.posInRead));
        	table.getColumns().add(makeColumn("Ref-Pos", O->O.posInRef));
        	table.getColumns().add(makeColumn("OP",new Function<CigarAndBase,String>() {
        		@Override
        		public String apply(CigarAndBase param) {
        			return param.op==null?null:param.op.name();
        		}
				} ));
        	table.getColumns().add(makeColumn("Len", O->O.count));
        	table.getColumns().add(makeColumn("Read-Bases",new Function<CigarAndBase,String>() {
        		@Override
        		public String apply(CigarAndBase param) {
        			return param.base==null?null:String.valueOf((char)param.base.intValue());
        		}
				} ));
            return table;
        	}
        
        private Tab buildJavascriptPane()
			{
			
			final ScrollPane scroll=new ScrollPane(super.javascriptArea);
			scroll.setFitToWidth(true);
			scroll.setFitToHeight(true);
			final BorderPane pane=new BorderPane(scroll);
			
			final Label helpLabel=new Label("The script injects:\n"+
    				"* header ( htsjdk.samtools.SAMFileHeader )\n"+
					"* record ( htsjdk.samtools.SAMRecord )\n"+
    				"The script should return a boolean: true (accept record) or false (reject record)"
					);
    		helpLabel.setWrapText(true);
    		pane.setBottom(helpLabel);
			
			final Tab tab=new Tab("Javascript",pane);
			tab.setClosable(false);
			return tab;
			}
        
        private Tab createReadGroupPane(final SAMFileHeader header)
        	{
        	final TableView<SAMReadGroupRecord> table=new TableView<>(header==null?
        			FXCollections.observableArrayList():
        			FXCollections.observableArrayList(header.getReadGroups())
        			);
        	table.getColumns().add(makeColumn("ID", G->G.getId()));
        	table.getColumns().add(makeColumn("Sample", G->G.getSample()));
        	table.getColumns().add(makeColumn("Center", G->G.getSequencingCenter()));
        	table.getColumns().add(makeColumn("Platform", G->G.getPlatform()));
        	table.getColumns().add(makeColumn("PlatformModel", G->G.getPlatformModel()));
        	table.getColumns().add(makeColumn("PlatformUnit", G->G.getPlatformUnit()));
        	table.getColumns().add(makeColumn("MedianInsertSize", G->G.getPredictedMedianInsertSize()));
        	table.getColumns().add(makeColumn("Desc", G->G.getDescription()));
        	table.getColumns().add(makeColumn("PU", G->G.getPlatformUnit()));
        	table.getColumns().add(makeColumn("Lib", G->G.getLibrary()));
        	table.getColumns().add(makeColumn("Run-Date", G->G.getRunDate()));
        	
        	final Tab tab=new Tab("ReadGroups", table);
        	tab.setClosable(false);
        	return tab;
        	}
        
        private Tab createProgramRecordPane(final SAMFileHeader header)
	    	{
	    	final TableView<SAMProgramRecord> table=new TableView<>(header==null?
	    			FXCollections.observableArrayList():
	    			FXCollections.observableArrayList(header.getProgramRecords())
	    			);
        	table.getColumns().add(makeColumn("ID", G->G.getId()));
        	table.getColumns().add(makeColumn("PG-ID", G->G.getProgramGroupId()));
        	table.getColumns().add(makeColumn("Prev-PG-ID", G->G.getPreviousProgramGroupId()));
        	table.getColumns().add(makeColumn("Version", G->G.getProgramVersion()));
        	table.getColumns().add(makeColumn("Command", G->G.getCommandLine()));
	        
	    	final Tab tab=new Tab("PG", table);
	    	tab.setClosable(false);
	    	return tab;
	    	}
        
        
        @Override
        void reloadData() {
        	final int max_items= this.maxReadLimitSpinner.getValue();
        	final List<SAMRecord> L= new ArrayList<SAMRecord>(max_items);
        	final String location = this.gotoField.getText().trim();
        	final SAMRecordIterator iter;
        	java.util.function.Predicate<SAMRecord> recFilter= x -> true;
        	
        	for(final SAMFlag flag: this.flag2filterInMenuItem.keySet())
        		{
        		CheckMenuItem cbox = this.flag2filterInMenuItem.get(flag);
        		if(!cbox.isSelected()) continue;
        		recFilter=recFilter.and(R-> flag.isSet(R.getFlags()));
        		}
        	for(final SAMFlag flag: this.flag2filterOutMenuItem.keySet())
        		{
        		CheckMenuItem cbox = this.flag2filterOutMenuItem.get(flag);
        		if(!cbox.isSelected()) continue;
        		recFilter=recFilter.and(R-> !flag.isSet(R.getFlags()));
        		}
        	
        	if(location.isEmpty())
        		{
        		iter = this.samReader.iterator();
        		}
        	else if(location.equalsIgnoreCase("unmapped"))
        		{
        		iter = this.samReader.queryUnmapped();
        		}
        	else
        		{
        		final Interval interval=parseInterval(location);
        		if(interval==null)
        			{
        			iter=null;
        			}
        		else
        			{
        			iter=samReader.queryOverlapping(interval.getContig(),interval.getStart(), interval.getEnd());
        			}
        		}
        	
        	CompiledScript compiledScript=null;
        	SimpleBindings binding=null;
        	if(JfxNgs.this.javascriptEngine!=null && !this.javascriptArea.getText().trim().isEmpty())
        		{
        		try
        			{
        			binding=new SimpleBindings();
        			compiledScript = JfxNgs.this.javascriptEngine.compile(this.javascriptArea.getText());
        			binding.put("header", this.samReader.getFileHeader());
        			}
        		catch(Exception err)
        			{
        			LOG.warning(err.getMessage());
        			updateStatusBar(AlertType.ERROR, err);
        			compiledScript=null;
        			}
        		}
        	final Map<ContigPos,Pileup> pos2pileup=new TreeMap<>();
        	Function<ContigPos,Pileup> getpileup=new  Function<ContigPos, Pileup>() {
				@Override
				public Pileup apply(ContigPos t) {
					Pileup p =pos2pileup.get(t);
					if(p==null) { p=new Pileup(t.contig,t.position);pos2pileup.put(t,p);}
					return p;
				}
			};
			
        	int count_items=0;
        	while(iter!=null && iter.hasNext() && count_items<max_items)
        		{
        		final SAMRecord rec = iter.next();
        		++count_items;
        		if(compiledScript!=null)
        			{
        			binding.put("record", rec);
        			if(!super.accept(compiledScript, binding)) continue;
        			}
        		
        		if(!recFilter.test(rec)) continue;
        		L.add(rec);
        		
        		/* FILL pileup */
        		if(!rec.getReadUnmappedFlag() && rec.getCigar()!=null)
        			{
        			int refpos=rec.getUnclippedStart();
        			int readpos=0;
        			final byte bases[]=rec.getReadBases();
        			final byte quals[]=rec.getOriginalBaseQualities();
        			
            		final Function<Integer,Character> getBaseAt = new Function<Integer, Character>() {
    					@Override
    					public Character apply(Integer readPos) {
    						char c;
    						if(bases==null || bases.length<=readPos)
    							{
    							return '?';
    							}
    						else
    							{
    							c=(char)bases[readPos];
    							}	
    						c=(rec.getReadNegativeStrandFlag()?
    								Character.toLowerCase(c):
    								Character.toUpperCase(c)
    								);
    						return c;
    					}
    				};
        			
            		final Function<Integer,Character> getQualAt = new Function<Integer, Character>() {
    					@Override
    					public Character apply(Integer readPos) {
    						char c;
    						if(quals==null || quals.length<=readPos)
    							{
    							return '#';
    							}
    						else
    							{
    							c=(char)quals[readPos];
    							}	
    						return c;
    					}
    				};

    				
        			for(final CigarElement ce:rec.getCigar())
        				{
        				switch(ce.getOperator())
        					{
        					case P:break;
        					case N:case D:
        						{
        						refpos+=ce.getLength();
        						break;
        						}
        					case H:
        						{
    							for(int i=0;i< ce.getLength();++i) {
        							Pileup p = getpileup.apply(new ContigPos(rec.getContig(),refpos));
        							p.seq.append('?');
        							p.qual.append('?');
        							p.operators.append((char)CigarOperator.enumToCharacter(ce.getOperator()));
        							++refpos;
        							}
        						break;
        						}
        					case I:
        						{
    							Pileup p = getpileup.apply(new ContigPos(rec.getContig(),refpos));
    							p.seq.append(getBaseAt.apply(readpos));
    							p.qual.append(getQualAt.apply(readpos));
    							p.operators.append((char)CigarOperator.enumToCharacter(ce.getOperator()));
    							readpos+=ce.getLength();
    							break;
        						}
        					case S: case EQ: case X: case M:
        						for(int i=0;i< ce.getLength();++i) {
        							Pileup p = getpileup.apply(new ContigPos(rec.getContig(),refpos));
        							p.seq.append(getBaseAt.apply(readpos));
        							p.qual.append(getQualAt.apply(readpos));
        							p.operators.append((char)CigarOperator.enumToCharacter(ce.getOperator()));
        							++readpos;
        							++refpos;
        							}
        						break;
        					}
        				}
        			}
        		
        		
        		}
        	if(iter!=null) iter.close();
        	this.recordTable.getItems().setAll(L);
        	this.pileupTable.getItems().setAll(pos2pileup.values());
        	
        	
        	final int countDisplayable = (int)getDisplayableSamRecordStream().count();
        	this.canvasScrollV.setMin(0);
        	this.canvasScrollV.setMax(countDisplayable);
        	this.canvasScrollV.setValue(0);
        	
        	repaintCanvas();
        	}
        @Override
    	void openInIgv() {
	    	final SAMRecord ctx=this.recordTable.getSelectionModel().getSelectedItem();
	    	if(ctx==null) {
	    		updateStatusBar(AlertType.WARNING,"no variant selected");
	    		return;
	    		}
	    	if(ctx.getReadUnmappedFlag()) return;
	    	openInIgv(ctx);
	    	}
        }
    
    /** display item in VCF/INFO */
    private static class InfoTableRow
		{
		final String key;
		final Integer index;
		final Object value;
		InfoTableRow(final String key,final Integer index,final Object value)
			{
			this.key = key;
			this.index = index;
			this.value=value;
			}
		}
    private class VcfStage extends StageContent
    	{
    	private final TextField javascriptFILTERfield=new TextField("");
    	private final VCFFileReader vcfFileReader;
    	private Spinner<Integer> maxReadLimitSpinner;
    	private TableView<VariantContext> variantTable;
    	private TableView<Genotype> genotypeTable;
    	private TableView<InfoTableRow> infoTableRow;
    	private TableView<String> filterTableRow;
    	
    	VcfStage(final String path) {
    		this.setTitle(path);
    		this.vcfFileReader = new VCFFileReader(new File(path),true);
    		final VCFHeader header=this.vcfFileReader.getFileHeader();
    		
    		
            final Menu fileMenu=new Menu("File");
            fileMenu.getItems().add(createMenuItem("Open",new Runnable() {
				@Override
				public void run() {
					openNgsFiles(VcfStage.this);
				}
			}));
            fileMenu.getItems().add(createMenuItem("Save Filtered VCF as... ",new Runnable() {
				@Override
				public void run() {
					doMenuSaveAs();
				}
			}));
            fileMenu.getItems().add(createMenuItem("Close",new Runnable() {
				@Override
				public void run() {
					VcfStage.this.hide();
				}
			}));
            //selectFlagMenu.getItems().addAll(flag2filterOutMenuItem.values());
            final MenuBar menuBar=new MenuBar(fileMenu);
            final VBox vbox1 = new VBox();
            vbox1.getChildren().add(menuBar);
            
            FlowPane top1= new FlowPane();
            vbox1.getChildren().add(top1);
            top1.getChildren().add(new Label("GoTo:"));
            top1.getChildren().add(this.gotoField);
            final Button gotoButton=new Button("Go");
            gotoButton.setOnAction(new EventHandler<ActionEvent>()
				{
				@Override
				public void handle(ActionEvent event)
					{
					reloadData();
					}
				});
            top1.getChildren().add(gotoButton);
            top1.getChildren().add(new Separator(Orientation.VERTICAL));
            top1.getChildren().add(new Label("Limit:"));
            top1.getChildren().add(this.maxReadLimitSpinner=new Spinner<Integer>(0,100000,JfxNgs.DEFAULT_VCF_RECORDS_COUNT));
            top1.getChildren().add(new Separator(Orientation.VERTICAL));
    		
            final Button igvButton =new Button("IGV");
            top1.getChildren().add(igvButton);
	        igvButton.setOnAction(new EventHandler<ActionEvent>() {
					@Override
					public void handle(ActionEvent event) {
						openInIgv();
					}
				});
            this.gotoField.setOnAction(new EventHandler<ActionEvent>()
				{
				@Override
				public void handle(ActionEvent event)
					{
					reloadData();
					}
				});
    		final TabPane tabPane=new TabPane();
    		tabPane.setPadding(new Insets(10, 10, 10, 10));
    		vbox1.getChildren().add(tabPane);
    		
    		
    		GridPane gridPane = new GridPane();
    		gridPane.setPadding(new Insets(10, 10, 10, 10));
    		gridPane.setVgap(4);
    		gridPane.setHgap(4);

            
    		/* build variant table */
    		this.variantTable = this.buildVariantTable();
    		GridPane.setConstraints( this.variantTable,0, 0,5,10); // column=3 row=1
    		gridPane.getChildren().add(this.variantTable);
    		
    		/* build genotype table */
    		this.genotypeTable =this.buildGenotypeTableRow(header);
    		GridPane.setConstraints( this.genotypeTable,5, 0,5,10); 
    		gridPane.getChildren().add(this.genotypeTable);
    		
    		/* filter table */
    		this.filterTableRow = this.buildFilterTable();
    		GridPane.setConstraints( this.filterTableRow,0, 10,3,2);
    		gridPane.getChildren().add(this.filterTableRow);

    		
    		/* build info Table table */
    		this.infoTableRow = this.buildInfoTableRow();
    		GridPane.setConstraints( this.infoTableRow,3, 10,8,2); // column=3 row=1
    		gridPane.getChildren().add(this.infoTableRow);
    		
          
           
            
            //vbox1.getChildren().add(gridPane);
    		
    		final Tab tab=new Tab("Variants", gridPane);
    		tab.setClosable(false);
    		tabPane.getTabs().add(tab);
    		tabPane.getTabs().add(buildInfoHeaderTab(header));
    		tabPane.getTabs().add(buildFormatHeaderTab(header));
    		tabPane.getTabs().add(buildFilterHeaderTab(header));
            tabPane.getTabs().add(buildDictTab( header.getSequenceDictionary()));
            tabPane.getTabs().add(buildJavascriptPane());
            /* register selection handler */
            /* when a read is selected update the flagsTable and metaDataTable */
            this.variantTable.getSelectionModel().selectedItemProperty().addListener((obs, oldSelection, newSelection) -> {
                fireSelectedVariantChanged(newSelection);
            	});
            
            
            final FlowPane bottom=new FlowPane(super.messageLabel);
            vbox1.getChildren().add(bottom);
            final Rectangle2D primaryScreenBounds = Screen.getPrimary().getVisualBounds();

            final Scene scene = new Scene(vbox1,
            		primaryScreenBounds.getWidth()-200+randomMoveWindow.nextInt(10),
            		primaryScreenBounds.getHeight()-200+randomMoveWindow.nextInt(10));
           
            this.setScene(scene);
            
           

    		}
    	
        @Override
        protected SAMSequenceDictionary getSAMSequenceDictionary() {
        	return this.vcfFileReader.getFileHeader().getSequenceDictionary();
        	}

    	
    	private void doMenuSaveAs()
    		{
    		final FileChooser fc=new FileChooser();
        	fc.setSelectedExtensionFilter(new ExtensionFilter("VCF Files", "*.vcf.gz"));
    		File saveAs=fc.showSaveDialog(this);
    		if(saveAs==null) return;
    		if(!saveAs.getName().endsWith(".vcf.gz"))
    			{
    			Alert alert=new Alert(AlertType.ERROR, "Output should end with .vcf.gz", ButtonType.OK);
    			alert.showAndWait();
    			return;
    			}
    		
        	CompiledScript compiledScript=null;
        	SimpleBindings binding=null;
        	if(JfxNgs.this.javascriptEngine!=null && !this.javascriptArea.getText().trim().isEmpty())
        		{
        		try
        			{
        			binding=new SimpleBindings();
        			compiledScript = JfxNgs.this.javascriptEngine.compile(this.javascriptArea.getText());
        			binding.put("header", this.vcfFileReader.getFileHeader());
        			}
        		catch(final Exception err)
        			{
        			showExceptionDialog(this, err);
        			return;
        			}
        		}

    		
    		VariantContextWriter w=null;
    		CloseableIterator<VariantContext> iter=null;
    		try
    			{
    			final VariantContextWriterBuilder vcwb=new VariantContextWriterBuilder();
    			vcwb.setOutputFile(saveAs);
    			w= vcwb.build();
    			VCFHeader header2= new VCFHeader(this.vcfFileReader.getFileHeader());
            	final String filterName = this.javascriptFILTERfield.getText().trim();
            	if(!filterName.isEmpty())
            		{
            		header2.addMetaDataLine(new VCFFilterHeaderLine(filterName, "Set by User in JfxNgs"));
            		}
    			w.writeHeader(header2);
    			iter=this.vcfFileReader.iterator();
    			while(iter.hasNext())
    				{
    				VariantContext ctx=iter.next();
    				if(compiledScript!=null)
	        			{
	        			binding.put("variant", ctx);
	        			if(!super.accept(compiledScript,binding))
	        				{
	        				if(filterName.isEmpty()) continue;
	        				ctx=new VariantContextBuilder(ctx).filter(filterName).make();
	        				}
	        			}
    				w.add(ctx);
    				}
    			w.close();
    			}
    		catch(Exception err)
    			{
    			showExceptionDialog(this, err);
    			return;
    			}
    		finally
    			{
    			CloserUtil.close(w);
    			}    		}
    	
    	private Tab buildJavascriptPane()
    		{
    		final ScrollPane scroll=new ScrollPane(super.javascriptArea);
    		scroll.setFitToWidth(true);
    		scroll.setFitToHeight(true);
    		final BorderPane pane=new BorderPane(scroll);
    		if(JfxNgs.this.javascriptEngine!=null) {
	    		final FlowPane flowPane = new FlowPane(
	    				new Label("set following FILTER on rejection:"),
	    				javascriptFILTERfield
	    				);
	    		
	    		this.javascriptFILTERfield.setPromptText("If not empty , don't discard the variant but set the FILTER");
	    		pane.setTop(flowPane);
	    		}
    		final Label helpLabel=new Label("The script injects:\n"+
    				"* header ( htsjdk.variant.vcf.VCFHeader )\n"+
					"* variant ( htsjdk.variant.variantcontext.VariantContext )\n"+
    				"The script should return a boolean: true (accept variant) or false (reject variant)"
					);
    		helpLabel.setWrapText(true);
    		pane.setBottom(helpLabel);
    		
    		final Tab tab=new Tab("Javascript",pane);
    		tab.setClosable(false);
    		return tab;
    		}
    	
    	/** build FILTER table */
    	private TableView<String> buildFilterTable()
    		{
			final TableView<String> table=new TableView<>();
			
			TableColumn<String,String>  scol = new TableColumn<>("Filter");
			scol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<String,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<String, String> param) {
					return new ReadOnlyObjectWrapper<String>(param.getValue());
					}
				});
            table.getColumns().add(scol);
			
            
			return table;    
			}
    	
    	@Override
    	void openInIgv() {
	    	final VariantContext ctx=this.variantTable.getSelectionModel().getSelectedItem();
	    	if(ctx==null) {
	    		updateStatusBar(AlertType.WARNING,"no variant selected");
	    		return;
	    		}
	    	openInIgv(ctx);
	    	}
    	
    	
    	/** build table of variants */
    	private TableView<VariantContext> buildVariantTable()
			{
			final TableView<VariantContext> table=new TableView<>();
			table.getColumns().add(makeColumn("CHROM", V->V.getContig()));
            table.getColumns().add(makeColumn("POS", V->V.getStart()));
            table.getColumns().add(makeColumn("ID", V->V.hasID()?V.getID():null));
            table.getColumns().add(makeColumn("REF", V->V.getReference().getDisplayString()));
            table.getColumns().add(makeColumn("ALT", V->V.getAlternateAlleles().stream().map(A->A.getDisplayString()).collect(Collectors.joining(","))));
            table.getColumns().add(makeColumn("FILTER", V->V.getFilters().stream().collect(Collectors.joining(","))));
            table.getColumns().add(makeColumn("QUAL", V->V.hasLog10PError()?V.getPhredScaledQual():null));
            return table;
			}

    	/** build INFO table */
    	private TableView<InfoTableRow> buildInfoTableRow()
    		{
    		final TableView<InfoTableRow> table=new TableView<JfxNgs.InfoTableRow>();
    		table.getColumns().add(makeColumn("Key", R->R.key));
    		table.getColumns().add(makeColumn("Index", R->R.index));
    		table.getColumns().add(makeColumn("Value", R->R.value));
    		return table;
    		}
          
    	/** build Genotype table */
    	private TableView<Genotype> buildGenotypeTableRow(final VCFHeader header)
			{
			final TableView<Genotype> table=new TableView<Genotype>();
			
			/* sample */
			table.getColumns().add(makeColumn("Sample", G->G.getSampleName()));
			
			for(final VCFFormatHeaderLine h:header.getFormatHeaderLines())
				{
	            final TableColumn<Genotype,String>  newcol = new TableColumn<>(h.getID());
	            newcol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<Genotype,String>, ObservableValue<String>>() {				
					@Override
					public ObservableValue<String> call(CellDataFeatures<Genotype, String> param) {
						Object o = param.getValue().getAnyAttribute(param.getTableColumn().getText());
						if(o==null)
							{
							return new ReadOnlyObjectWrapper<String>(null);
							}
						if(o instanceof List)
							{
							List<?> L=(List<?>)o;
							String delim=(param.getTableColumn().getText().equals(VCFConstants.GENOTYPE_KEY) && param.getValue().isPhased()?"|":",");
							o = L.stream().map(S -> String.valueOf(S)).collect(Collectors.joining(delim)).toString();
							}
						return new ReadOnlyObjectWrapper<String>(String.valueOf(o));
						}
					});
	            table.getColumns().add(newcol);
				}
			/* type */
			table.getColumns().add(makeColumn("Type", G->G.getType().name()));
			return table;
			}

    	
    	
        private void fireSelectedVariantChanged(final VariantContext ctx)
        	{
        	if(ctx!=null)
        		{
        		this.genotypeTable.getItems().setAll(ctx.getGenotypes());
        		
        		this.filterTableRow.getItems().setAll(
        				ctx.getFilters()
        				);
        		
        		final List<InfoTableRow> infos=new ArrayList<>();
        		final Map<String,Object> atts = ctx.getAttributes();
        		for(final String key:atts.keySet())
        			{
        			Object v= atts.get(key);
        			final List<?> L;
        			if(v instanceof List)
        				{
        				L=(List<?>)v;
        				}
        			else if(v.getClass().isArray())
        				{
        				Object a[]=(Object[])v;
        				L=Arrays.asList(a);
        				}
        			else
        				{
        				L=Collections.singletonList(v);
        				}
        			for(int x=0;x< L.size();++x)
        				{
        				
        				infos.add(new InfoTableRow(key,(L.size()==1?null:x+1),L.get(x)));
        				}
        			}
        		this.infoTableRow.getItems().setAll(infos);
        		}
        	else
        		{
            	this.genotypeTable.getItems().clear();
            	this.infoTableRow.getItems().clear();
            	this.filterTableRow.getItems().clear();
        		}
        	}
    	
    	/** build a table describing the FORMAT column */
    	private Tab buildFormatHeaderTab(final VCFHeader header)
    		{
            final TableView<VCFFormatHeaderLine> table=new TableView<>(FXCollections.observableArrayList(header.getFormatHeaderLines()));
            table.getColumns().add(makeColumn("ID", F->F.getID()));
            table.getColumns().add(makeColumn("Type", F->F.getType()==null?null:F.getType().name()));
            table.getColumns().add(makeColumn("Count", F->F.isFixedCount()?F.getCount():null));
            table.getColumns().add(makeColumn("Description", F->F.getDescription()));
            final Tab tab=new Tab("FORMAT",table);
            tab.setClosable(false);
            return tab;
    		}
  
    	/** build a table describing the INFO column */
    	private Tab buildInfoHeaderTab(final VCFHeader header)
    		{
            final TableView<VCFInfoHeaderLine> table=new TableView<>(FXCollections.observableArrayList(header.getInfoHeaderLines()));
            table.getColumns().add(makeColumn("ID", F->F.getID()));
            table.getColumns().add(makeColumn("Type", F->F.getType()==null?null:F.getType().name()));
            table.getColumns().add(makeColumn("Count", F->F.isFixedCount()?F.getCount():null));
            table.getColumns().add(makeColumn("Description", F->F.getDescription()));
            final Tab tab=new Tab("INFO",table);
            tab.setClosable(false);
            return tab;
    		}
    	
    	/** build a table describing the INFO column */
    	private Tab buildFilterHeaderTab(final VCFHeader header)
    		{
            final TableView<VCFFilterHeaderLine> table=new TableView<>(FXCollections.observableArrayList(header.getFilterLines()));
            table.getColumns().add(makeColumn("ID", F->F.getID()));
            final Tab tab=new Tab("FILTER",table);
            tab.setClosable(false);
            return tab;
    		}
   	
    	
    	@Override
        void closeNgsResource() {
        	CloserUtil.close(this.vcfFileReader);
        	}
    	
    	@Override
        void reloadData()
    		{
        	final int max_items= this.maxReadLimitSpinner.getValue();
        	final List<VariantContext> L= new ArrayList<>(max_items);
        	final String location = this.gotoField.getText().trim();
        	final CloseableIterator<VariantContext> iter;
        	
        	
        	if(location.isEmpty())
        		{
        		iter = this.vcfFileReader.iterator();
        		}
        	else
        		{
        		final Interval interval=this.parseInterval(location);
        		if(interval==null)
        			{
        			iter=null;
        			}
        		else
        			{
    				iter=this.vcfFileReader.query(interval.getContig(),interval.getStart(),interval.getEnd());
        			}
        		}
        	CompiledScript compiledScript=null;
        	SimpleBindings binding=null;
        	if(JfxNgs.this.javascriptEngine!=null && !this.javascriptArea.getText().trim().isEmpty())
        		{
        		try
        			{
        			binding=new SimpleBindings();
        			compiledScript = JfxNgs.this.javascriptEngine.compile(this.javascriptArea.getText());
        			binding.put("header", this.vcfFileReader.getFileHeader());
        			}
        		catch(Exception err)
        			{
        			LOG.warning(err.getMessage());
        			compiledScript=null;
        			}
        		}
        	final String filterName = this.javascriptFILTERfield.getText().trim();
        	int count_items=0;
        	while(iter!=null && iter.hasNext() && count_items<max_items)
        		{
        		VariantContext rec = iter.next();
        		++count_items;
        		if(compiledScript!=null)
        			{
        			binding.put("variant", rec);
        			if(!super.accept(compiledScript,binding))
        				{
        				if(filterName.isEmpty()) continue;
        				rec=new VariantContextBuilder(rec).filter(filterName).make();
        				}
        			}
        		
        		L.add(rec);
        		}
        	if(iter!=null) iter.close();
        	this.variantTable.getItems().setAll(L);
        	}

    	}
    
    
    public JfxNgs()
		{
		this.preferences = Preferences.userNodeForPackage(JfxNgs.class);
		Compilable engine;
		try {
			final ScriptEngineManager manager = new ScriptEngineManager();
			engine = (Compilable)manager.getEngineByName("js");
			}
		catch(Exception err)
			{
			engine=null;
			LOG.warning("Cannot get Compilable JS engine "+err.getMessage());
			}
		this.javascriptEngine = engine;
		}

    @Override
    public void stop() throws Exception
    	{
    	try {
    		this.preferences.sync();
    		}
    	catch(BackingStoreException err)
    		{
    		LOG.warning(err.getMessage());
    		}
    	try {
    		LOG.info("flush preferences");
    		this.preferences.flush();
    		}
    	catch(BackingStoreException err)
    		{
    		LOG.warning(err.getMessage());
    		}
    	super.stop();
    	}
    
    /*
    private void showPreferenceDialoge(Window parentStage)
	    {
	    Stage dialog = new Stage();
	
	     
	     dialog.setTitle("Preferences");
		 dialog.initOwner(parentStage);
		 dialog.initModality(Modality.APPLICATION_MODAL); 
		 dialog.showAndWait();
	    }*/
    
    @Override
    public void start(final Stage primaryStage) throws Exception {
        final Parameters params = this.getParameters();
        
        
       
       
        primaryStage.setTitle(getClass().getSimpleName());
        Menu menu=new Menu("File");
        
        MenuItem menuItem=new MenuItem("About...");
        menuItem.setOnAction(new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent event) {
				Alert alert=new Alert(AlertType.INFORMATION);
				alert.setHeaderText("JFXNGS");
				alert.setContentText("Pierre Lindenbaum PhD. 2017.\n"
						+ "@yokofakun\n"
						+ "Institut du Thorax - Nantes - France\n"
						+ "https://github.com/lindenb/jvarkit"
						);
				alert.showAndWait();
			}
		});
        menu.getItems().add(menuItem);
        
        menuItem=new MenuItem("Open...");
        menuItem.setOnAction(new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent event) {
				openNgsFiles(primaryStage);
				
			}
		});
        menu.getItems().add(menuItem);
        
        
        
        menuItem=new MenuItem("Quit...");
        menuItem.setOnAction(new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent event) {
				doMenuQuit();
			}
		});
        menu.getItems().add(menuItem);
        
        MenuBar bar=new MenuBar(menu);
        FlowPane flow=new FlowPane(5,5);
        flow.setPadding(new Insets(10,10,10,10));
        flow.getChildren().add(new Label("Set Location of all frames to:"));
        final TextField textField=new TextField();
        textField.setPrefColumnCount(15);
        textField.setPromptText("Location");
        flow.getChildren().add(textField);
        Button button=new Button("Go");
        flow.getChildren().add(button);
        
        EventHandler<ActionEvent> handler=new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent event) {
				final String loc=textField.getText().trim();
				LOG.info("moveTo all to "+loc);
				for(final StageContent sc:all_opened_stages )
					{
					LOG.info("moveTo "+sc.getTitle()+" to "+loc);
					sc.moveTo(loc);
					}
				}
		};
        button.setOnAction(handler);
        textField.setOnAction(handler);
        
        
        BorderPane pane=new BorderPane();
        
        
        pane.setBottom(new Label("Author: Pierre Lindenbaum PhD."));
        
       
        VBox vbox1= new VBox(bar,flow,pane);
        
        
        final Rectangle2D primaryScreenBounds = Screen.getPrimary().getVisualBounds();

        final Scene scene = new Scene(vbox1,
        		primaryScreenBounds.getWidth()-200+randomMoveWindow.nextInt(10),
        		primaryScreenBounds.getHeight()-200+randomMoveWindow.nextInt(10));
       
        primaryStage.setScene(scene);
        primaryStage.setOnHidden(e -> doMenuQuit());
       
        Collection<StageContent> newStages=new ArrayList<>();
        Exception lastException=null;
        for(final String arg: params.getUnnamed())
	        {
	        try 
		        {
	        	newStages.addAll(openNgsFiles(arg));
		        }
	        catch(final Exception err)
	        	{
	        	lastException=err;
	        	}
	        }
        if(lastException!=null)
        	{
        	final Exception error=lastException;
        	primaryStage.setOnShown(new EventHandler<WindowEvent>() {
				@Override
				public void handle(WindowEvent event) {
					showExceptionDialog(primaryStage,error);
				}
			});
        	}
        primaryStage.show();
        for(StageContent sc:newStages)
	    	{
	    	sc.show();
	    	}
        }

    private void unregisterStage(StageContent s) {
    	this.all_opened_stages.remove(s);
    	LOG.info("unregister nbr.win"+this.all_opened_stages.size());

    }
    private void registerStage(StageContent s) {
    	this.all_opened_stages.add(s);  
    	LOG.info("register nbr.win"+this.all_opened_stages.size());
    }
    
    
    private static void showExceptionDialog(final Window owner,Throwable error)
    	{
    	final Alert alert = new Alert(AlertType.WARNING);
		alert.setTitle("Error");
		alert.setHeaderText("Error");
		// Create expandable Exception.
		
		StringWriter sw = new StringWriter();
		PrintWriter pw = new PrintWriter(sw);
		error.printStackTrace(pw);
		String exceptionText = sw.toString();

		Label label = new Label("The exception stacktrace was:");

		TextArea textArea = new TextArea(exceptionText);
		textArea.setEditable(false);
		textArea.setWrapText(true);

		textArea.setMaxWidth(Double.MAX_VALUE);
		textArea.setMaxHeight(Double.MAX_VALUE);
		GridPane.setVgrow(textArea, Priority.ALWAYS);
		GridPane.setHgrow(textArea, Priority.ALWAYS);

		GridPane expContent = new GridPane();
		expContent.setMaxWidth(Double.MAX_VALUE);
		expContent.add(label, 0, 0);
		expContent.add(textArea, 0, 1);

		// Set expandable Exception into the dialog pane.
		alert.getDialogPane().setExpandableContent(expContent);
		alert.showAndWait();
    	}
    
    private void doMenuQuit()
    	{
    	Platform.exit();
    	}
    
    private void openNgsFiles(final Window owner)
    	{
    	final String LAST_USED_DIR="last.used.dir";
    	final FileChooser fc=new FileChooser();
    	String lastDirStr= preferences.get(LAST_USED_DIR, null);
    	if(lastDirStr!=null && !lastDirStr.isEmpty())
    		{
    		fc.setInitialDirectory(new File(lastDirStr));
    		}
    	
    	fc.setSelectedExtensionFilter(new ExtensionFilter("NGS Files", "*.bam","*.vcf","*.vcf.gz","*.list"));
    	final List<File> selFiles = fc.showOpenMultipleDialog(owner);
    	if(selFiles==null || selFiles.isEmpty()) return ;
    	List<StageContent> stages =new ArrayList<>();
    	File parentDir=null;
    	try 
    		{
	    	for(final File f:selFiles)
	    		{
	    		stages.addAll( openNgsFiles(f.getPath()) );
	    		parentDir=f.getParentFile();
	    		}
	    	for(StageContent sc:stages)
	    		{	
	    		sc.show();
	    		}
	    	if(parentDir!=null)
	    		{
	    		preferences.put(LAST_USED_DIR, parentDir.getPath());
	    		}
    		}
    	catch(final Exception err)
    		{
    		for(StageContent sc:stages) sc.closeNgsResource();
    		showExceptionDialog(owner, err);
    		}		
    	}

    private Collection<StageContent> openNgsFiles(final String path0) throws Exception
    	{
    	final List<String> pathList;
    	if(!IOUtil.isUrl(path0) && path0.endsWith(".list"))
    		{
    		pathList = Files.lines(Paths.get(path0)).filter(L-> ! (L.isEmpty() && L.startsWith("#"))).collect(Collectors.toList());
    		}
    	else
    		{
    		pathList = Collections.singletonList(path0);
    		}
		final List<StageContent> stages =new ArrayList<>();

    	for(final String uri:pathList)
	    	{
	    	//try as BAM
	    	if(uri.endsWith(".bam") )
		    	{
	    		final StageContent sc=openBam(uri);
		    	if(sc!=null) stages.add(sc);
		    	}
	    	else if(uri.endsWith(".vcf") || uri.endsWith(".bcf")  || uri.endsWith(".vcf.gz") )
		    	{
	    		final StageContent sc=openVcf(uri);
		    	if(sc!=null) stages.add(sc);
		    	}
	    	else
	    		{
	    		LOG.warning("Don't know how to open "+uri);
	    		}	
	    	}
    	LOG.info("N opened stages = "+stages.size());
    	return stages;
    	}
    
    private BamStageContent openBam(final String uri)
    	{
    	SamReader samIn=null;
    	try
			{
			final SamReaderFactory srf = SamReaderFactory.makeDefault();
			srf.validationStringency(ValidationStringency.LENIENT);
			samIn = srf.open(SamInputResource.of(uri));
			if(!samIn.hasIndex())
				{
				LOG.warning("No index for "+uri);
				return null;
				}
			if(samIn.getFileHeader()==null)
				{
				LOG.warning("cannot get SAM header for "+uri);
				return null;
				}
			if(samIn.getFileHeader().getSequenceDictionary()==null)
				{
				LOG.warning("cannot get SAM Dictionary for "+uri);
				return null;
				}	
			samIn.close();
			LOG.info("OK for BAM "+uri);
			return new BamStageContent(uri);
			}
		catch(final Exception err)
			{
			LOG.warning("not an indexed bam file : "+uri);
			return null;
			}
		finally
			{
			CloserUtil.close(samIn);
			}
    	}
    
    private VcfStage openVcf(final String uri)
    	{
    	VCFFileReader vcfIn=null;
    	try
    		{
    		vcfIn = new VCFFileReader(new File(uri), true);
    		if(vcfIn.getFileHeader()==null) {
    			LOG.info("No VCF header in "+uri);
    			return null;
    			}
    		if(vcfIn.getFileHeader().getSequenceDictionary()==null) {
    			LOG.info("No VCF idctionary in "+uri);
    			return null;
    			}
    		if(vcfIn.getFileHeader()!=null && vcfIn.getFileHeader().getSequenceDictionary()!=null)
	    		{
	    		vcfIn.close();
	    		LOG.info("OK for VCF "+uri);
				return new VcfStage(uri);
	    		}
    		return null;
    		}
    	catch(Exception err)
    		{
    		LOG.warning("not an indexed vcf file : "+uri);
    		return null;
    		}
    	finally
    		{
    		CloserUtil.close(vcfIn);
    		}
    	}
    
    private MenuItem createMenuItem(final String label,final Runnable runner)
    	{
    	final MenuItem menu=new MenuItem(label);
    	menu.setOnAction(new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent event) {
				runner.run();
			}
		});
    	return menu;
    	}

    public static void main(String[] args) {
    	
        launch(args);
    }

}

