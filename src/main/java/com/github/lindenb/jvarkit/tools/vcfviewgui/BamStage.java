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
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.TreeMap;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.function.Supplier;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import javax.script.CompiledScript;
import javax.script.ScriptException;

import com.github.lindenb.jvarkit.tools.samjs.SamJavascript;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.BasesPerPositionChartFactory;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.ChartFactory;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.CigarOpPerPositionChartFactory;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.GCPercentChartFactory;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.MapqChartFactory;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.QualityPerPositionChartFactory;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.ReadChartFactory;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.ReadLengthChartFactory;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.ReadQualityChartFactory;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.SamFlagsChartFactory;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.VariantContextChartFactory;
import com.github.lindenb.jvarkit.util.Hershey;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.SAMRecord.SAMTagAndValue;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import javafx.application.Platform;
import javafx.beans.value.ChangeListener;
import javafx.beans.value.ObservableValue;
import javafx.collections.FXCollections;
import javafx.event.ActionEvent;
import javafx.event.Event;
import javafx.event.EventHandler;
import javafx.geometry.Insets;
import javafx.geometry.Orientation;
import javafx.geometry.Rectangle2D;
import javafx.scene.Scene;
import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.chart.Chart;
import javafx.scene.control.Alert;
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
import javafx.scene.control.SeparatorMenuItem;
import javafx.scene.control.SpinnerValueFactory;
import javafx.scene.control.SplitPane;
import javafx.scene.control.Tab;
import javafx.scene.control.TabPane;
import javafx.scene.control.TableCell;
import javafx.scene.control.TableColumn;
import javafx.scene.control.TableView;
import javafx.scene.control.TextArea;
import javafx.scene.input.ScrollEvent;
import javafx.scene.control.Alert.AlertType;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.FlowPane;
import javafx.scene.layout.VBox;
import javafx.scene.paint.Color;
import javafx.scene.text.Font;
import javafx.scene.text.Text;
import javafx.scene.text.TextFlow;
import javafx.stage.FileChooser;
import javafx.stage.Screen;
import javafx.stage.Stage;
import javafx.stage.WindowEvent;
import javafx.stage.FileChooser.ExtensionFilter;

@SuppressWarnings("unused")
public class BamStage extends NgsStage {
    private static final int DEFAULT_BAM_RECORDS_COUNT=Integer.parseInt(System.getProperty("jxf.ngs.default.sam", "1000"));
    static final ExtensionFilter EXTENSION_FILTER=new ExtensionFilter("Bam Files", ".bam");
    private static final Logger LOG= Logger.getLogger("BamStage");
    private static final String SPINNER_VALUE_KEY="bam.spinner.value";
    /** shor-Read oriented chart-factories */
    private static final List<Supplier<ChartFactory<SAMRecord>>> READ_CHART_LIST=Arrays.asList(
    		()->new BasesPerPositionChartFactory(),
    		()->new QualityPerPositionChartFactory(),
    		()->new GCPercentChartFactory(),
    		()->new ReadLengthChartFactory(),
    		()->new SamFlagsChartFactory(),
    		()->new MapqChartFactory(),
    		()->new ReadQualityChartFactory(),
    		()->new CigarOpPerPositionChartFactory()
    		);
    
    
    private static class BamJavascripFilter
    	extends JavascriptFilter<SAMFileHeader,SAMRecord>
		{
		protected BamJavascripFilter(
				final SAMFileHeader header,
				final Optional<CompiledScript> compiledScript
				)
			{
			super(header,compiledScript);
			}
		@Override public SAMRecord eval(final SAMRecord rec)
			{
			if(super.compiledScript==null) return rec;
			super.bindings.put("record", rec);
			return super.accept()?rec:null;
			}
		}

    
    private class ReadQualityStage
		extends AbstractQualityStage<SAMRecord>
			{
    		private class ScanQualThread extends ScanThread
    			{
    			private final Predicate<SAMRecord> flagFilters;
        		ScanQualThread(
        				final ChartFactory<SAMRecord> factory,
        				final JfxNgs.InputSource source,
        				final Optional<CompiledScript> compiledScript,
        				final Predicate<SAMRecord> flagFilters
        				)
    				{
    				super(factory,source,compiledScript,flagFilters);
    				this.flagFilters=flagFilters;
    				}
        		
    			@Override
    			public void run() {
    				SamReader samReader =null;
    				SAMRecordIterator samIter=null;
    				try 
    					{
	    				BamJavascripFilter bamJavascripFilter=null;
	    				final SamReaderFactory srf=SamReaderFactory.makeDefault().
	    						validationStringency(ValidationStringency.SILENT);
    					samReader = this.source.isFile()?
    							srf.open(this.source.asFile()):
    							srf.open(SamInputResource.of(this.source.asUrl()))
    							;
    					samIter = samReader.iterator();
    					
    					if(this.compiledScript!=null)
    						{
    						bamJavascripFilter=new BamJavascripFilter(samReader.getFileHeader(),
    								compiledScript);
    						}	
    					
    					while(!kill_flag && samIter.hasNext())
    						{
    						final SAMRecord rec=samIter.next();
    						
    						nItems++;
    						if(!flagFilters.test(rec)) continue;
    						if(bamJavascripFilter!=null)
    							{
    							if(bamJavascripFilter.eval(rec)==null) continue;
    							}
    						this.factory.visit(rec);
    						update();
    						}
    					samIter.close();
    					samReader.close();
    					super.atEnd();	    				
	    				}
    				catch(final Throwable err)
    					{
    					super.onError(err);
    					}
    				finally
    					{
    					CloserUtil.close(samIter);
    					CloserUtil.close(samReader);
    					CloserUtil.close(samReader);
    					}
    				}
    			
    			}
    		
	    	ReadQualityStage(
	    			final ChartFactory<SAMRecord> factory,
	    			final JfxNgs.InputSource file,
	    			final Optional<CompiledScript> compiledScript,
	    			final Predicate<SAMRecord> filters)
	    		{
	    		super(factory,file,compiledScript,filters);
	    		}
	    	
			@Override
			protected AbstractQualityStage<SAMRecord>.ScanThread createThread(
					final ChartFactory<SAMRecord> factory,
					JfxNgs.InputSource file,
					final Optional<CompiledScript> compiledScript,
					final Predicate<SAMRecord> filters
					) {
				return new ScanQualThread(factory,file,compiledScript,filters);
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
    	void watch(char base,char q,CigarOperator op)
    		{
    		seq.append(base);
    		switch(Character.toUpperCase(base))
    			{
    			case 'A': count[0]++;break;
    			case 'T': count[1]++;break;
    			case 'G': count[2]++;break;
    			case 'C': count[3]++;break;
    			case '<':case '>': case '-':break;
    			default: count[4]++;break;
    			}
    		qual.append(q);
    		this.operators.append((char)CigarOperator.enumToCharacter(op));
    		}
    	int depth() { return count[0]+count[1]+count[2]+count[3];}
    	}
	
    private final SamReader samReader;
    private final TableView<SAMRecord> recordTable;
    private final TableView<SamFlagRow> flagsTable;
    private final TableView<SAMTagAndValue> metaDataTable;
    private final TableView<CigarAndBase> cigarTable;
    private final TableView<Pileup> pileupTable;
    private final Map<SAMFlag,CheckMenuItem> flag2filterInMenuItem=new HashMap<>();
    private final Map<SAMFlag,CheckMenuItem> flag2filterOutMenuItem=new HashMap<>();
    private final Canvas canvas = new Canvas(900, 600);
    private final ScrollBar canvasScrollV = new ScrollBar();
    private final CheckBox canvasShowReadName = new CheckBox("Show Read Name");
    
    BamStage(final JfxNgs owner,final JfxNgs.InputSource urlOrFile) throws IOException
    	{
    	super(owner,urlOrFile);
        final SamReaderFactory srf= SamReaderFactory.makeDefault();
        srf.validationStringency(Level.OFF.equals(LOG.getLevel())?
        		ValidationStringency.SILENT:
        		ValidationStringency.LENIENT
        		);
        LOG.info("Opening "+urlOrFile);
        
        if(urlOrFile.isFile())
        	{
        	 this.samReader=srf.open(urlOrFile.asFile());
        	}
        else
        	{
        	 this.samReader=srf.open(SamInputResource.of(urlOrFile.asUrl()));
        	}
       
        if(!this.samReader.hasIndex())
        	{
        	this.samReader.close();
        	throw new IOException("Bam without index "+urlOrFile);
        	}
        
        if(this.samReader.getFileHeader()==null)
    		{
	    	this.samReader.close();
	    	throw new IOException("Bam without header "+urlOrFile);
	    	}
        if(this.samReader.getFileHeader().getSequenceDictionary()==null)
			{
	    	this.samReader.close();
	    	throw new IOException("Bam without dictionary "+urlOrFile);
	    	}
        
        /** Build menu for SAM Flags */
        for(final SAMFlag flg:SAMFlag.values())
        	{
        	this.flag2filterInMenuItem.put(flg,new CheckMenuItem("Filter In "+flg.name()));
        	this.flag2filterOutMenuItem.put(flg,new CheckMenuItem("Filter Out "+flg.name()));
        	}
        
       
        
        
        final VBox vbox1 = new VBox();
        vbox1.getChildren().add(super.menuBar);
        
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
        
        super.maxItemsLimitSpinner.setValueFactory(
        		new SpinnerValueFactory.IntegerSpinnerValueFactory(0,100000,
        				owner.preferences.getInt(SPINNER_VALUE_KEY, DEFAULT_BAM_RECORDS_COUNT)
        				));;
        top1.getChildren().add(gotoButton);
        top1.getChildren().add(new Separator(Orientation.VERTICAL));
        top1.getChildren().add(new Label("Limit:"));
        top1.getChildren().add(super.maxItemsLimitSpinner);
        top1.getChildren().add(new Separator(Orientation.VERTICAL));
        List<CheckMenuItem> menuFlags=new ArrayList<>(flag2filterOutMenuItem.values());
        menuFlags.addAll(flag2filterInMenuItem.values());
        CheckMenuItem tmp[]=new CheckMenuItem[menuFlags.size()];
        menuFlags.toArray(tmp);
        top1.getChildren().add(new MenuBar(new Menu("Flags",null,tmp)));
        super.gotoField.setOnAction(new EventHandler<ActionEvent>()
			{
			@Override
			public void handle(ActionEvent event)
				{
				reloadData();
				}
			});
      
        top1.getChildren().add(createIgvButton());
       
        
        
        TabPane tabbedPane = new TabPane();
        tabbedPane.setPadding(new Insets(10, 10, 10, 10));
        Tab tab= new Tab("Reads");
        tab.setClosable(false);
        tabbedPane.getTabs().add(tab);
        
        this.recordTable = makeRecordTable();
        /** create columns */
        
        
        
        

        
        //ScrollPane scroll = new ScrollPane(this.recordTable);
        //scroll.setFitToHeight(true);
        //scroll.setFitToWidth(true);
        
        BorderPane borderPane1=new BorderPane(this.recordTable);
        borderPane1.setPadding(new Insets(5));
        SplitPane split1= new SplitPane();
        split1.setOrientation(Orientation.VERTICAL);
        split1.getItems().add(borderPane1);
        
       
        SplitPane split2= new SplitPane();
        split1.getItems().add(split2);
        split2.setOrientation(Orientation.HORIZONTAL);
        
        /* define SAM Flag table */
        this.flagsTable= createSamFlagTable();
       
        
        //scroll.setFitToHeight(true);
        split2.getItems().add(this.flagsTable); 
        
        /* define Meta Data table */
        this.metaDataTable = createMetaDataTable();
        split2.getItems().add(this.metaDataTable); 

        
        /* build the cigar table */
        this.cigarTable = createCigarTable();
        split2.getItems().add(this.cigarTable); 

        
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

        split2.setDividerPositions(0.1f, 0.6f, 0.9f);
        tab.setContent(split1);
       
        this.pileupTable = createPileupTable();
        tab=new Tab("Pileup",this.pileupTable);
        tab.setClosable(false);
        tabbedPane.getTabs().add(tab);

        
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
        
        this.canvasShowReadName.setOnAction(E->repaintCanvas());
        this.canvasScrollV.valueProperty().addListener(E->repaintCanvas());
        
        tab=new Tab("Canvas",canvasPane);
        tab.setClosable(false);
        tabbedPane.getTabs().add(tab);
        /* END CANVAS STUFF */
        
        vbox1.getChildren().add(tabbedPane);
        
        final FlowPane bottom=new FlowPane(super.messageLabel);
        vbox1.getChildren().add(bottom);
        
        final Rectangle2D primaryScreenBounds = Screen.getPrimary().getVisualBounds();

        this.setScene(new Scene(vbox1,
        		primaryScreenBounds.getWidth()-200,
        		primaryScreenBounds.getHeight()-200)
        		);
        		
        
        super.fileMenu.getItems().addAll(
        		menuForSavingTable("SamRecord",this.recordTable),
        		menuForSavingTable("MetaData",this.metaDataTable),
        		menuForSavingTable("Flags",this.flagsTable),
        		menuForSavingTable("Pileup",this.pileupTable),
        		menuForSavingTable("Cigar",this.cigarTable)
        		);
        
        /* fill stats menu */
        final Supplier<List<SAMRecord>> variantsProvider=()->this.recordTable.getItems();
        
        for(final Supplier<ChartFactory<SAMRecord>> supplier: READ_CHART_LIST)
	        {
        	final ChartFactory<SAMRecord> factory = supplier.get();
        	final MenuItem menuItem=new MenuItem("Local "+factory.getName());
        	menuItem.setOnAction(new EventHandler<ActionEvent>() {
				@Override
				public void handle(ActionEvent event) {
					doMenuShowLocalStats(factory, variantsProvider);
				}
			});
        	statsMenu.getItems().add(menuItem);
	        }
        super.statsMenu.getItems().add(new SeparatorMenuItem());
        for(final Supplier<ChartFactory<SAMRecord>> supplier: READ_CHART_LIST)
	        {
	    	final ChartFactory<SAMRecord> factory = supplier.get();
	    	final MenuItem menuItem=new MenuItem("Whole"+factory.getName());
	    	menuItem.setOnAction(new EventHandler<ActionEvent>() {
				@Override
				public void handle(ActionEvent event) {
					doMenuShowWholeStats(factory);
				}
			});
	    	super.statsMenu.getItems().add(menuItem);
	        }
        
        
        this.addEventHandler(
    			WindowEvent.WINDOW_CLOSE_REQUEST ,new EventHandler<WindowEvent>() {
                    @Override
                    public void handle(final WindowEvent event) {
                    owner.preferences.putInt(SPINNER_VALUE_KEY,maxItemsLimitSpinner.getValue().intValue());
                    }
                });
        
        
    	}

    private TableView<SAMRecord> makeRecordTable() {
    	final TableView<SAMRecord>  table = new TableView<SAMRecord>();
    	table.getColumns().add(makeColumn("Read-Name",REC->REC.getReadName()));
    	table.getColumns().add(makeColumn("Flag",REC->REC.getFlags()));
    	table.getColumns().add(makeColumn("Ref",REC->REC.getReferenceName()));
    	table.getColumns().add(makeColumn("Read-Pos",REC->REC.getAlignmentStart()));
    	table.getColumns().add(makeColumn("MAPQ",REC->REC.getMappingQuality()));
    	table.getColumns().add(makeColumn("CIGAR",REC->REC.getCigarString()));
    	table.getColumns().add(makeColumn("LEN",REC->REC.getInferredInsertSize()));
    	table.getColumns().add(makeColumn("Mate-Ref",REC->REC.getMateReferenceName()));
    	table.getColumns().add(makeColumn("Mate-Pos",REC->REC.getMateAlignmentStart()));
    	
    	final Font font1=new Font("Courier", 9);
    	
    	final TableColumn<SAMRecord, String> tc=makeColumn("SEQ",REC->REC.getReadString());
    	
    	// http://stackoverflow.com/questions/42187987/
    	tc.setCellFactory(tv -> new TableCell<SAMRecord, String>() { 
    		final TextFlow textFlow = new TextFlow();
    		
    	    @Override
    	    protected void updateItem(final String item, boolean empty) {
    	        super.updateItem(item, empty);
    	        setText(null);
    	        if(item==null)
    	        	{
    	        	//setText(null);
    	            setGraphic(null);
    	        	return;
    	        	}
    	        final List<Text> L=new ArrayList<>(item.length());
    	        for(int i=0;i< item.length();++i) {
    	        	final Text txt=new Text(String.valueOf(item.charAt(i)));new Text(String.valueOf(item.charAt(i)));
    	        	txt.setFont(font1);
    	        	txt.setStroke(JfxNgs.BASE2COLOR.apply(item.charAt(i)));
    	        	L.add(txt);
    	        	}
    	        this.textFlow.setLineSpacing(0.1);
    	        this.textFlow.setMaxHeight(10);
    	        this.textFlow.setPrefHeight(10);
    	        this.textFlow.getChildren().setAll(L);
    	        this.setGraphic(textFlow);
    	    }
    	});
	    table.getColumns().add(tc);
	    table.getColumns().add(makeColumn("QUAL",REC->REC.getBaseQualityString()));
	    return table;
	    }
    
    @Override
    protected SAMSequenceDictionary getSAMSequenceDictionary() {
    	return this.samReader.getFileHeader().getSequenceDictionary();
    	}
    
    @Override
    void closeNgsResource()
    	{
    	LOG.info("closing"+super.urlOrFile);
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
    
    
    /** repaint the canvas area */
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
				public Color apply(final Integer readPos) {
					if(bases==null || bases.length<=readPos)
						{
						return Color.BLACK;
						}
					return JfxNgs.BASE2COLOR.apply((char)bases[readPos]);
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
    	table.getColumns().add(makeColumn("Status",param-> param.flag.isSet(param.record.getFlags())?"\u2611":"\u2610"));
    		
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
    		public String apply(final CigarAndBase param) {
    			return param.op==null?null:param.op.name();
    		}
			} ));
    	table.getColumns().add(makeColumn("Len", O->O.count));
    	table.getColumns().add(makeColumn("Read-Bases",new Function<CigarAndBase,String>() {
    		@Override
    		public String apply(final CigarAndBase param) {
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
		pane.setPadding(new Insets(10));
		
		final FlowPane top=new FlowPane();
		top.getChildren().addAll(super.makeJavascriptButtons());
		pane.setTop(top);
		
		final Label helpLabel=new Label("The script injects:\n"+
				"* header ( htsjdk.samtools.SAMFileHeader )\n"+
				"* record ( htsjdk.samtools.SAMRecord )\n"+
				"The script should return a boolean: true (accept record) or false (reject record)"
				);
		helpLabel.setWrapText(true);
		pane.setBottom(helpLabel);
		
		final Tab tab=new Tab(JAVASCRIPT_TAB_KEY,pane);
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
    
    /** build a Predicate for filtering on SAM FLAG using the checkboxes */
    private Predicate<SAMRecord> makeFlagPredicate()
    	{
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
		return recFilter;
    	}
    
    @Override
    void reloadData() {
    	final int max_items= super.maxItemsLimitSpinner.getValue();
    	final List<SAMRecord> L= new ArrayList<SAMRecord>(max_items);
    	final String location = this.gotoField.getText().trim();
    	final SAMRecordIterator iter;
    	final java.util.function.Predicate<SAMRecord> recFilter=makeFlagPredicate();
    	
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
    	Optional<BamJavascripFilter> bamjsfilter=Optional.empty();
    	if(this.owner.javascriptCompiler.isPresent() &&
    			!this.javascriptArea.getText().trim().isEmpty())
    		{
    		try
    			{
    			bamjsfilter=Optional.of(new BamJavascripFilter(
    					this.samReader.getFileHeader(),
    					Optional.of(this.owner.javascriptCompiler.get().compile(this.javascriptArea.getText()))
    					));
    			}
    		catch(Exception err)
    			{
    			LOG.warning(err.getMessage());
    			updateStatusBar(AlertType.ERROR, err);
    			bamjsfilter=null;
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
    		if(bamjsfilter.isPresent())
    			{
    			if(bamjsfilter.get().eval(rec)==null) continue;
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
    			
    			/** function getting the ith base */
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
    			
				/** function getting the ith base quality */
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
    							final Pileup p = getpileup.apply(new ContigPos(rec.getContig(),refpos));
    							p.watch('-','#',ce.getOperator());
    							++refpos;
    							}
    						break;
    						}
    					case I:
    						{
    						for(int i=0;i< ce.getLength();++i) {
    						final Pileup p = getpileup.apply(new ContigPos(rec.getContig(),refpos));
    							p.watch('<',getQualAt.apply(readpos),ce.getOperator());
    							readpos++;
        						}
							break;
    						}
    					case S:
    						for(int i=0;i< ce.getLength();++i) {
    							final Pileup p = getpileup.apply(new ContigPos(rec.getContig(),refpos));
    							p.watch('-',getQualAt.apply(readpos),ce.getOperator());
    							++readpos;
    							++refpos;
    							}
    						break;
    					case EQ: case X: case M:
    						for(int i=0;i< ce.getLength();++i) {
    							final Pileup p = getpileup.apply(new ContigPos(rec.getContig(),refpos));
    							p.watch(getBaseAt.apply(readpos),getQualAt.apply(readpos),ce.getOperator());
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
    		updateStatusBar(AlertType.WARNING,"no Read selected");
    		return;
    		}
    	if(ctx.getReadUnmappedFlag())
    		{
    		updateStatusBar(AlertType.WARNING,"read is not mapped");
    		return;
    		}
    	openInIgv(ctx);
    	}
    
  


	@Override
	protected void doMenuShowWholeStats(final ChartFactory<?> factory) {
    	Optional<CompiledScript> compiledScript=Optional.empty();
    	if(this.owner.javascriptCompiler.isPresent() &&
    		!this.javascriptArea.getText().trim().isEmpty())
    		{
    		try
    			{
    			compiledScript =Optional.of(this.owner.javascriptCompiler.get().compile(this.javascriptArea.getText()));
    			}
    		catch(final Exception err)
    			{
    			JfxNgs.showExceptionDialog(this, err);
    			return;
    			}
    		}

    	@SuppressWarnings("unchecked")
		final ReadQualityStage qcstage=new ReadQualityStage(
    			(ChartFactory<SAMRecord>)factory,
    			super.urlOrFile,
    			compiledScript,
    			makeFlagPredicate()
    			);
		qcstage.show();
		}


	@Override
	protected void doMenuSaveAs() {
		final FileChooser fc= owner.newFileChooser();
    	fc.setSelectedExtensionFilter(EXTENSION_FILTER);
		final File saveAs= owner.updateLastDir(fc.showSaveDialog(this));
		if(saveAs==null) return;
		if(saveAs.equals(super.urlOrFile)) return;
		if(!saveAs.getName().endsWith(".bam"))
			{
			final Alert alert=new Alert(AlertType.ERROR, "Output should end with .bam", ButtonType.OK);
			alert.showAndWait();
			return;
			}
		
    	Optional<BamJavascripFilter> bamjsfilter=Optional.empty();
    	if(this.owner.javascriptCompiler.isPresent() && 
    			!this.javascriptArea.getText().trim().isEmpty())
    		{
    		try
    			{
    			bamjsfilter=Optional.of(
    					new BamJavascripFilter(this.samReader.getFileHeader(),
    							Optional.of(this.owner.javascriptCompiler.get().compile(this.javascriptArea.getText()))
    							));
    			}
    		catch(final Exception err)
    			{
    			JfxNgs.showExceptionDialog(this, err);
    			bamjsfilter=Optional.empty();
    			return;
    			}
    		}

    	Predicate<SAMRecord> filter=makeFlagPredicate();
		final SAMFileWriterFactory swf= new SAMFileWriterFactory();
		swf.setCreateIndex(true);
		SAMRecordIterator iter=null;
		SAMFileWriter w=null;
		try
			{
			w = swf.makeBAMWriter(this.samReader.getFileHeader(), true, saveAs);
			
			iter=this.samReader.iterator();
			while(iter.hasNext())
				{
				final SAMRecord rec=iter.next();
				if(!filter.test(rec)) continue;
				if(bamjsfilter.isPresent())
        			{
        			if(bamjsfilter.get().eval(rec)==null) continue;
        			}
				
				w.addAlignment(rec);
				}
			w.close();
			}
		catch(Exception err)
			{
			JfxNgs.showExceptionDialog(this, err);
			return;
			}
		finally
			{
			CloserUtil.close(w);
			}    		
		}
    @Override
    protected String getSnippetResourcePath() {
    	return "/com/github/lindenb/jvarkit/tools/vcfviewgui/bam.snippets.xml";
    	}
    }
