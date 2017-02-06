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
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.prefs.BackingStoreException;
import java.util.prefs.Preferences;
import java.util.stream.Collectors;

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
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import javafx.application.Application;
import javafx.application.Platform;
import javafx.beans.property.ReadOnlyObjectWrapper;
import javafx.beans.value.ObservableValue;
import javafx.collections.FXCollections;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.geometry.Insets;
import javafx.geometry.Orientation;
import javafx.scene.Scene;
import javafx.scene.control.Alert;
import javafx.scene.control.Alert.AlertType;
import javafx.scene.control.Button;
import javafx.scene.control.CheckMenuItem;
import javafx.scene.control.Label;
import javafx.scene.control.Menu;
import javafx.scene.control.MenuBar;
import javafx.scene.control.MenuItem;
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
import javafx.stage.FileChooser;
import javafx.stage.FileChooser.ExtensionFilter;
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
    
    
    /** abstract base class for NGS window */
    private abstract class StageContent
        extends Stage
        {
        public StageContent() {
            this.setOnCloseRequest(new EventHandler<WindowEvent>() {
                @Override
                public void handle(WindowEvent event) {
                    unregisterStage(StageContent.this);
                    }
                });
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
        	if(feature==null) return;
        	openInIgv(
        			Collections.singletonList("goto "+feature.getContig()+":"+feature.getStart()+"-"+feature.getEnd())
        			);
        	}
        
        abstract void openInIgv();
        /** close the NGS resource , even if the window was not opened */
        abstract void closeNgsResource();
        abstract void reloadData();
        
        
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


    /** NGS window for BAM */
    private class BamStageContent extends StageContent
        {
        private final SamReader samReader;
        private final TableView<SAMRecord> recordTable;
        private final TableView<SamFlagRow> flagsTable;
        private final TableView<SAMTagAndValue> metaDataTable;
        private final TableView<CigarAndBase> cigarTable;
        private final Map<SAMFlag,CheckMenuItem> flag2filterInMenuItem=new HashMap<>();
        private final Map<SAMFlag,CheckMenuItem> flag2filterOutMenuItem=new HashMap<>();
        private final TextField gotoField;
        private final Spinner<Integer> maxReadLimitSpinner;
       
        
        BamStageContent(final String url) {
        	this.setTitle(url);
            final SamReaderFactory srf= SamReaderFactory.makeDefault();
            srf.validationStringency(Level.OFF.equals(LOG.getLevel())?
            		ValidationStringency.SILENT:
            		ValidationStringency.LENIENT
            		);
            LOG.info("Opening "+url);
            this.samReader=srf.open(SamInputResource.of(url));
            if(this.samReader.hasIndex())
            	{
            	
            	}

            /** Build menu for SAM Flags */
            for(final SAMFlag flg:SAMFlag.values())
            	{
            	flag2filterInMenuItem.put(flg,new CheckMenuItem("Filter In "+flg.name()));
            	flag2filterOutMenuItem.put(flg,new CheckMenuItem("Filter Out "+flg.name()));
            	}
            final Menu fileMenu=new Menu("File");
            //selectFlagMenu.getItems().addAll(flag2filterOutMenuItem.values());
            final MenuBar menuBar=new MenuBar(fileMenu);
            final VBox vbox1 = new VBox();
            vbox1.getChildren().add(menuBar);
            
            FlowPane top1= new FlowPane(5,5);
            top1.setPadding(new Insets(10, 10, 10, 10));
            vbox1.getChildren().add(top1);
            top1.getChildren().add(new Label("GoTo:"));
            top1.getChildren().add(this.gotoField = new TextField());
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
            top1.getChildren().add(this.maxReadLimitSpinner=new Spinner<Integer>(0,100000,1000));
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
            this.flagsTable= new TableView<>();
            /* create SamFlag columns */
            final TableColumn<SamFlagRow,String>  flagNameCol = new TableColumn<>("FLAG");
            flagNameCol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<SamFlagRow,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<SamFlagRow, String> param) {
					return new ReadOnlyObjectWrapper<String>(param.getValue().flag.name());
					}
				});
            this.flagsTable.getColumns().add(flagNameCol);
            
            /* create value set/notset for columns */
            final TableColumn<SamFlagRow,Boolean>  flagStatusCol = new TableColumn<>("Status");
            flagStatusCol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<SamFlagRow,Boolean>, ObservableValue<Boolean>>() {				
				@Override
				public ObservableValue<Boolean> call(final CellDataFeatures<SamFlagRow, Boolean> param) {
					return new ReadOnlyObjectWrapper<Boolean>(param.getValue().flag.isSet(param.getValue().record.getFlags()));
					}
				});
            this.flagsTable.getColumns().add(flagStatusCol);
            
            //scroll.setFitToHeight(true);
            GridPane.setConstraints(  this.flagsTable,1, 1); // column=1 row=1
            tilePane.getChildren().add( this.flagsTable);
            
            
            /* define Meta Data table */
            this.metaDataTable = new TableView<>();
            final TableColumn<SAMTagAndValue,String>  metaDataKey = new TableColumn<>("Key");
            metaDataKey.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<SAMTagAndValue,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<SAMTagAndValue, String> param) {
					return new ReadOnlyObjectWrapper<String>(param.getValue().tag);
					}
				});
            this.metaDataTable.getColumns().add(metaDataKey);
            final TableColumn<SAMTagAndValue,String>  metaDataValue = new TableColumn<>("Value");
            metaDataValue.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<SAMTagAndValue,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<SAMTagAndValue, String> param) {
					return new ReadOnlyObjectWrapper<String>(String.valueOf(param.getValue().value));
					}
				});
            this.metaDataTable.getColumns().add(metaDataValue);

            GridPane.setConstraints( this.metaDataTable,2, 1); // column=2 row=1
            tilePane.getChildren().add(this.metaDataTable);

            
            /* build the cigar table */
            this.cigarTable = new TableView<>();
            
            
            final TableColumn<CigarAndBase,String>  cigarRefCol = new TableColumn<>("REF");
            cigarRefCol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<CigarAndBase,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<CigarAndBase, String> param) {
					return new ReadOnlyObjectWrapper<String>(param.getValue().ref);
					}
				});
            this.cigarTable.getColumns().add(cigarRefCol);
           
            final TableColumn<CigarAndBase,Integer>  cigarReadPosCol = new TableColumn<>("Read-Pos");
            cigarReadPosCol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<CigarAndBase,Integer>, ObservableValue<Integer>>() {				
				@Override
				public ObservableValue<Integer> call(CellDataFeatures<CigarAndBase, Integer> param) {
					return new ReadOnlyObjectWrapper<Integer>(param.getValue().posInRead);
					}
				});
            this.cigarTable.getColumns().add(cigarReadPosCol);
            
            final TableColumn<CigarAndBase,Integer>  cigarRefPosCol = new TableColumn<>("Ref-Pos");
            cigarRefPosCol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<CigarAndBase,Integer>, ObservableValue<Integer>>() {				
				@Override
				public ObservableValue<Integer> call(CellDataFeatures<CigarAndBase, Integer> param) {
					return new ReadOnlyObjectWrapper<Integer>(param.getValue().posInRef);
					}
				});
            this.cigarTable.getColumns().add(cigarRefPosCol);
            
            final TableColumn<CigarAndBase,String>  cigarOpCol = new TableColumn<>("OP");
            cigarOpCol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<CigarAndBase,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<CigarAndBase, String> param) {
					return new ReadOnlyObjectWrapper<String>(param.getValue().op==null?null:param.getValue().op.name());
					}
				});
            this.cigarTable.getColumns().add(cigarOpCol);
            
            
            final TableColumn<CigarAndBase,Integer>  cigarLenCol = new TableColumn<>("Len");
            cigarLenCol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<CigarAndBase,Integer>, ObservableValue<Integer>>() {				
				@Override
				public ObservableValue<Integer> call(CellDataFeatures<CigarAndBase, Integer> param) {
					return new ReadOnlyObjectWrapper<Integer>(param.getValue().count);
					}
				});
            this.cigarTable.getColumns().add(cigarLenCol);
            
            final TableColumn<CigarAndBase,String>  cigarBaseCol = new TableColumn<>("Read-Bases");
            cigarBaseCol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<CigarAndBase,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<CigarAndBase, String> param) {
					return new ReadOnlyObjectWrapper<String>(param.getValue().base==null?null:String.valueOf((char)param.getValue().base.intValue()));
					}
				});
            this.cigarTable.getColumns().add(cigarBaseCol);
            
            
            //scroll = new ScrollPane();
            //scroll.setFitToHeight(true);
            //scroll.setFitToWidth(true);
            GridPane.setConstraints( this.cigarTable,3, 1); // column=3 row=1
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
            
            vbox1.getChildren().add(tabbedPane);
            this.setScene(new Scene(vbox1,1000,800));
            
            this.setOnShowing(new EventHandler<WindowEvent>() {
				@Override
				public void handle(WindowEvent event) {
					reloadData();
				}
			});
            
            this.setOnCloseRequest(new EventHandler<WindowEvent>() {
                @Override
                public void handle(WindowEvent event) {
                	closeNgsResource();
                    unregisterStage(BamStageContent.this);
                }
            });
        }

        @Override
        void closeNgsResource()
        	{
        	LOG.info("closing"+samReader.getResourceDescription());
            CloserUtil.close(samReader);
        	}
        
        
        private Tab createReadGroupPane(final SAMFileHeader header)
        	{
        	TableView<SAMReadGroupRecord> table=new TableView<>(header==null?
        			FXCollections.observableArrayList():
        			FXCollections.observableArrayList(header.getReadGroups())
        			);
        	
            TableColumn<SAMReadGroupRecord,String>  scol = new TableColumn<>("ID");
            scol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<SAMReadGroupRecord,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<SAMReadGroupRecord, String> param) {
					return new ReadOnlyObjectWrapper<String>(param.getValue().getId());
					}
				});
            table.getColumns().add(scol);
            
            scol = new TableColumn<>("Sample");
            scol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<SAMReadGroupRecord,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<SAMReadGroupRecord, String> param) {
					return new ReadOnlyObjectWrapper<String>(param.getValue().getSample());
					}
				});
            table.getColumns().add(scol);

            scol = new TableColumn<>("SC");
            scol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<SAMReadGroupRecord,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<SAMReadGroupRecord, String> param) {
					return new ReadOnlyObjectWrapper<String>(param.getValue().getSequencingCenter());
					}
				});
            table.getColumns().add(scol);
           
            scol = new TableColumn<>("PF");
            scol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<SAMReadGroupRecord,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<SAMReadGroupRecord, String> param) {
					return new ReadOnlyObjectWrapper<String>(param.getValue().getPlatform());
					}
				});
            table.getColumns().add(scol);
            
           	
            scol = new TableColumn<>("Desc.");
            scol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<SAMReadGroupRecord,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<SAMReadGroupRecord, String> param) {
					return new ReadOnlyObjectWrapper<String>(param.getValue().getDescription());
					}
				});
            table.getColumns().add(scol);
            
            scol = new TableColumn<>("PU");
            scol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<SAMReadGroupRecord,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<SAMReadGroupRecord, String> param) {
					return new ReadOnlyObjectWrapper<String>(param.getValue().getPlatformUnit());
					}
				});
            table.getColumns().add(scol);
        	
        	Tab tab=new Tab("RG", table);
        	tab.setClosable(false);
        	return tab;
        	}
        
        private Tab createProgramRecordPane(final SAMFileHeader header)
	    	{
	    	TableView<SAMProgramRecord> table=new TableView<>(header==null?
	    			FXCollections.observableArrayList():
	    			FXCollections.observableArrayList(header.getProgramRecords())
	    			);
	    	
	        TableColumn<SAMProgramRecord,String>  scol = new TableColumn<>("ID");
	        scol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<SAMProgramRecord,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<SAMProgramRecord, String> param) {
					return new ReadOnlyObjectWrapper<String>(param.getValue().getId());
					}
				});
	        table.getColumns().add(scol);
	        
	        scol = new TableColumn<>("GroupID");
	        scol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<SAMProgramRecord,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<SAMProgramRecord, String> param) {
					return new ReadOnlyObjectWrapper<String>(param.getValue().getProgramGroupId());
					}
				});
	        table.getColumns().add(scol);
	        
	        scol = new TableColumn<>("Prev-GroupID");
	        scol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<SAMProgramRecord,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<SAMProgramRecord, String> param) {
					return new ReadOnlyObjectWrapper<String>(param.getValue().getPreviousProgramGroupId());
					}
				});
	        table.getColumns().add(scol);

	        scol = new TableColumn<>("Version");
	        scol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<SAMProgramRecord,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<SAMProgramRecord, String> param) {
					return new ReadOnlyObjectWrapper<String>(param.getValue().getProgramVersion());
					}
				});
	        table.getColumns().add(scol);
	        
	        
	        scol = new TableColumn<>("Command");
	        scol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<SAMProgramRecord,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<SAMProgramRecord, String> param) {
					return new ReadOnlyObjectWrapper<String>(param.getValue().getCommandLine());
					}
				});
	        table.getColumns().add(scol);
	        
	    	Tab tab=new Tab("PG", table);
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
        		CheckMenuItem cbox = this.flag2filterInMenuItem.get(flag);
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
        		final String contig;
        		int colon =location.indexOf(":");
        		if(colon==-1)
        			{
        			contig=location;
        			iter= this.samReader.queryAlignmentStart(contig, 1);
        			}
        		else
        			{
        			contig=location.substring(0,colon);
        			int hyphen=location.indexOf('-');
        			Integer start=null,end=null;
        			if(hyphen==-1)
        				{
        				try { start= new Integer(location.substring(colon+1).trim());}
        				catch(NumberFormatException err ) {start=null;}
        				}
        			else
        				{
        				try {
    						start= new Integer(location.substring(colon+1,hyphen).trim());
    						end= new Integer(location.substring(hyphen+1).trim());
        					}
        				catch(NumberFormatException err ) {start=null;end=null;}
        				}
        			if(start!=null && end!=null && start.compareTo(end)<=0)
        				{
        				iter=samReader.queryOverlapping(contig, start, end);
        				}
        			else
        				{
        				iter=null;
        				}
        			}
        		}
        	int count_items=0;
        	while(iter!=null && iter.hasNext() && count_items<max_items)
        		{
        		SAMRecord rec = iter.next();
        		++count_items;
        		if(!recFilter.test(rec)) continue;
        		L.add(rec);
        		}
        	if(iter!=null) iter.close();
        	this.recordTable.getItems().setAll(L);
        	}
        @Override
    	void openInIgv() {
	    	final SAMRecord ctx=this.recordTable.getSelectionModel().getSelectedItem();
	    	if(ctx==null) {
	    		LOG.info("no variant selected");
	    		}
	    	if(ctx.getReadUnmappedFlag()) return;
	    	openInIgv(ctx);
	    	}
        }
    
    private static class InfoTableRow
		{
		final String key;
		final Object value;
		InfoTableRow(final String key,final Object value)
			{
			this.key = key;
			this.value=value;
			}
		}
    private class VcfStage extends StageContent
    	{
    	
    	private final TextField gotoField;
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
            //selectFlagMenu.getItems().addAll(flag2filterOutMenuItem.values());
            final MenuBar menuBar=new MenuBar(fileMenu);
            final VBox vbox1 = new VBox();
            vbox1.getChildren().add(menuBar);
            
            FlowPane top1= new FlowPane();
            vbox1.getChildren().add(top1);
            top1.getChildren().add(new Label("GoTo:"));
            top1.getChildren().add(this.gotoField = new TextField());
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
            top1.getChildren().add(this.maxReadLimitSpinner=new Spinner<Integer>(0,100000,1000));
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

            /* register selection handler */
            /* when a read is selected update the flagsTable and metaDataTable */
            this.variantTable.getSelectionModel().selectedItemProperty().addListener((obs, oldSelection, newSelection) -> {
                fireSelectedVariantChanged(newSelection);
            	});
            
            
            
            this.setScene(new Scene(vbox1,1000,500));
            
            this.setOnShowing(new EventHandler<WindowEvent>() {
				@Override
				public void handle(WindowEvent event) {
					reloadData();
				}
			});
            
            this.setOnCloseRequest(new EventHandler<WindowEvent>() {
                @Override
                public void handle(WindowEvent event) {
                	closeNgsResource();
                    unregisterStage(VcfStage.this);
                	}
            	});

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
	    		LOG.info("no variant selected");
	    		}
	    	openInIgv(ctx);
	    	}
    	
    	
    	/** build table of variants */
    	private TableView<VariantContext> buildVariantTable()
			{
			final TableView<VariantContext> table=new TableView<>();
			
			TableColumn<VariantContext,String>  scol = new TableColumn<>("CHROM");
			scol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<VariantContext,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<VariantContext, String> param) {
					return new ReadOnlyObjectWrapper<String>(param.getValue().getContig());
					}
				});
            table.getColumns().add(scol);
			
            TableColumn<VariantContext,Integer>  icol = new TableColumn<>("POS");
			icol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<VariantContext,Integer>, ObservableValue<Integer>>() {				
				@Override
				public ObservableValue<Integer> call(CellDataFeatures<VariantContext, Integer> param) {
					return new ReadOnlyObjectWrapper<Integer>(param.getValue().getStart());
					}
				});
            table.getColumns().add(icol);
			
            scol = new TableColumn<>("ID");
			scol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<VariantContext,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<VariantContext, String> param) {
					return new ReadOnlyObjectWrapper<String>(param.getValue().hasID()?param.getValue().getID():null);
					}
				});
            table.getColumns().add(scol);
            
            scol = new TableColumn<>("REF");
			scol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<VariantContext,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<VariantContext, String> param) {
					return new ReadOnlyObjectWrapper<String>(param.getValue().getReference().getDisplayString());
					}
				});
            table.getColumns().add(scol);

            scol = new TableColumn<>("ALT");
			scol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<VariantContext,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<VariantContext, String> param) {
					return new ReadOnlyObjectWrapper<String>(
							param.getValue().getAlternateAlleles().stream().map(A->A.getDisplayString()).collect(Collectors.joining(","))
							);
					}
				});
            table.getColumns().add(scol);
            
            
            scol = new TableColumn<>("FILTER");
			scol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<VariantContext,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<VariantContext, String> param) {
					return new ReadOnlyObjectWrapper<String>(
							param.getValue().getFilters().stream().collect(Collectors.joining(","))
							);
					}
				});
            table.getColumns().add(scol);
            
            TableColumn<VariantContext,Double> fcol = new TableColumn<>("QUAL");
			fcol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<VariantContext,Double>, ObservableValue<Double>>() {				
				@Override
				public ObservableValue<Double> call(CellDataFeatures<VariantContext, Double> param) {
					return new ReadOnlyObjectWrapper<Double>(param.getValue().hasLog10PError()?
							param.getValue().getPhredScaledQual():
							null
							);
					}
				});
            table.getColumns().add(fcol);
            
            
			return table;
			}

    	/** build INFO table */
    	private TableView<InfoTableRow> buildInfoTableRow()
    		{
    		TableView<InfoTableRow> table=new TableView<JfxNgs.InfoTableRow>();
           
    		final TableColumn<InfoTableRow,String>  keycol = new TableColumn<>("Key");
    		keycol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<InfoTableRow,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<InfoTableRow, String> param) {
					return new ReadOnlyObjectWrapper<String>(param.getValue().key);
					}
				});
            table.getColumns().add(keycol);
            
            final TableColumn<InfoTableRow,String>  valuecol = new TableColumn<>("Value");
            valuecol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<InfoTableRow,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<InfoTableRow, String> param) {
					return new ReadOnlyObjectWrapper<String>(String.valueOf(param.getValue().value));
					}
				});
            table.getColumns().add(valuecol);
    		return table;
    		}
          
    	/** build Genotype table */
    	private TableView<Genotype> buildGenotypeTableRow(final VCFHeader header)
			{
			TableView<Genotype> table=new TableView<Genotype>();
			
			/* sample */
			final TableColumn<Genotype,String>  scol = new TableColumn<>("Sample");
			scol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<Genotype,String>, ObservableValue<String>>() {				
			@Override
			public ObservableValue<String> call(CellDataFeatures<Genotype, String> param) {
				return new ReadOnlyObjectWrapper<String>(param.getValue().getSampleName());
				} });
	        table.getColumns().add(scol);
			
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
			final TableColumn<Genotype,String>  newcol = new TableColumn<>("Type");
			newcol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<Genotype,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<Genotype, String> param) {
				return new ReadOnlyObjectWrapper<String>(param.getValue().getType().name());
				}});
			table.getColumns().add(newcol);
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
        			for(final Object o2:L)
        				{
        				infos.add(new InfoTableRow(key,o2));
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
           
    		/*  ID */
            final TableColumn<VCFFormatHeaderLine,String>  formatIDcol = new TableColumn<>("ID");
            formatIDcol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<VCFFormatHeaderLine,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<VCFFormatHeaderLine, String> param) {
					return new ReadOnlyObjectWrapper<String>(param.getValue().getID());
					}
				});
            table.getColumns().add(formatIDcol);

    		/*  type */
            final TableColumn<VCFFormatHeaderLine,String>  typecol = new TableColumn<>("Type");
            typecol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<VCFFormatHeaderLine,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<VCFFormatHeaderLine, String> param) {
					final VCFHeaderLineType type=param.getValue().getType();
					return new ReadOnlyObjectWrapper<String>(type==null?null:type.name());
					}
				});
            table.getColumns().add(typecol);
            
    		/*  Count */
            final TableColumn<VCFFormatHeaderLine,Integer>  countcol = new TableColumn<>("Count");
            countcol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<VCFFormatHeaderLine,Integer>, ObservableValue<Integer>>() {				
				@Override
				public ObservableValue<Integer> call(CellDataFeatures<VCFFormatHeaderLine, Integer> param) {
					if(!param.getValue().isFixedCount()) return null;
					return new ReadOnlyObjectWrapper<Integer>(param.getValue().getCount());
					}
				});
            table.getColumns().add(countcol);
          
            
            /* description */
            final TableColumn<VCFFormatHeaderLine,String>  desccol = new TableColumn<>("Description");
            desccol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<VCFFormatHeaderLine,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<VCFFormatHeaderLine, String> param) {
					return new ReadOnlyObjectWrapper<String>(param.getValue().getDescription());
					}
				});
            table.getColumns().add(desccol);
            
            Tab tab=new Tab("FORMAT",table);
            tab.setClosable(false);
            return tab;
    		}
  
    	/** build a table describing the INFO column */
    	private Tab buildInfoHeaderTab(final VCFHeader header)
    		{
            final TableView<VCFInfoHeaderLine> table=new TableView<>(FXCollections.observableArrayList(header.getInfoHeaderLines()));
           
    		/*  ID */
            final TableColumn<VCFInfoHeaderLine,String>  formatIDcol = new TableColumn<>("ID");
            formatIDcol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<VCFInfoHeaderLine,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<VCFInfoHeaderLine, String> param) {
					return new ReadOnlyObjectWrapper<String>(param.getValue().getID());
					}
				});
            table.getColumns().add(formatIDcol);

    		/*  type */
            final TableColumn<VCFInfoHeaderLine,String>  typecol = new TableColumn<>("Type");
            typecol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<VCFInfoHeaderLine,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<VCFInfoHeaderLine, String> param) {
					final VCFHeaderLineType type=param.getValue().getType();
					return new ReadOnlyObjectWrapper<String>(type==null?null:type.name());
					}
				});
            table.getColumns().add(typecol);
            
    		/*  Count */
            final TableColumn<VCFInfoHeaderLine,Integer>  countcol = new TableColumn<>("Count");
            countcol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<VCFInfoHeaderLine,Integer>, ObservableValue<Integer>>() {				
				@Override
				public ObservableValue<Integer> call(CellDataFeatures<VCFInfoHeaderLine, Integer> param) {
					if(!param.getValue().isFixedCount()) return null;
					return new ReadOnlyObjectWrapper<Integer>(param.getValue().getCount());
					}
				});
            table.getColumns().add(countcol);
          
            
            /* description */
            final TableColumn<VCFInfoHeaderLine,String>  desccol = new TableColumn<>("Description");
            desccol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<VCFInfoHeaderLine,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<VCFInfoHeaderLine, String> param) {
					return new ReadOnlyObjectWrapper<String>(param.getValue().getDescription());
					}
				});
            table.getColumns().add(desccol);
            
            Tab tab=new Tab("INFO",table);
            tab.setClosable(false);
            return tab;
    		}
    	
    	/** build a table describing the INFO column */
    	private Tab buildFilterHeaderTab(final VCFHeader header)
    		{
            final TableView<VCFFilterHeaderLine> table=new TableView<>(FXCollections.observableArrayList(header.getFilterLines()));
           
    		/*  ID */
            final TableColumn<VCFFilterHeaderLine,String>  formatIDcol = new TableColumn<>("ID");
            formatIDcol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<VCFFilterHeaderLine,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<VCFFilterHeaderLine, String> param) {
					return new ReadOnlyObjectWrapper<String>(param.getValue().getID());
					}
				});
            table.getColumns().add(formatIDcol);

            Tab tab=new Tab("FILTER",table);
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
        		final String contig;
        		int colon =location.indexOf(":");
        		if(colon==-1)
        			{
        			
        			contig=location;
        			LOG.info("VCF query chromosome "+location);
        			iter=this.vcfFileReader.query(contig, 1, Integer.MAX_VALUE);
        			}
        		else
        			{
        			contig=location.substring(0,colon);
        			int hyphen=location.indexOf('-');
        			Integer start=null,end=null;
        			if(hyphen==-1)
        				{
        				final String startStr=location.substring(colon+1).trim();
        				try { start= new Integer(startStr);end=Integer.MAX_VALUE;}
        				catch(NumberFormatException err ) {start=null;LOG.warning(startStr);}
        				}
        			else
        				{
        				try {
    						start= new Integer(location.substring(colon+1,hyphen).trim());
    						end= new Integer(location.substring(hyphen+1).trim());
        					}
        				catch(NumberFormatException err ) {start=null;end=null;LOG.warning(location);}
        				}
        			if(start!=null && end!=null && start.compareTo(end)<=0)
        				{
        				iter=this.vcfFileReader.query(contig, start, end);
        				}
        			else
        				{
        				iter=null;
        				}
        			}
        		}
        	int count_items=0;
        	while(iter!=null && iter.hasNext() && count_items<max_items)
        		{
        		VariantContext rec = iter.next();
        		++count_items;
        		L.add(rec);
        		}
        	if(iter!=null) iter.close();
        	this.variantTable.getItems().setAll(L);
        	}

    	}
    
    
    public JfxNgs()
		{
		this.preferences = Preferences.userNodeForPackage(JfxNgs.class);
		}

    @Override
    public void stop() throws Exception
    	{
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
        MenuItem menuItem=new MenuItem("Open...");
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
        
        BorderPane pane=new BorderPane();
        
        
        
        VBox vbox1= new VBox(bar,pane);
        
      
        final Scene scene = new Scene(vbox1,500,300);
       
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
    
    private Collection<StageContent> openNgsFiles(final Window owner)
    	{
    	FileChooser fc=new FileChooser();
    	fc.setSelectedExtensionFilter(new ExtensionFilter("NGS Files","*.bam","*.vcf","*.vcf.gz","list"));
    	List<File> selFiles = fc.showOpenMultipleDialog(owner);
    	if(selFiles==null || selFiles.isEmpty()) return Collections.emptyList();
    	List<StageContent> stages =new ArrayList<>();
    	try 
    		{
	    	for(final File f:selFiles)
	    		{
	    		stages.addAll( openNgsFiles(f.getPath()) );
	    		}
	    	return stages;
    		}
    	catch(final Exception err)
    		{
    		for(StageContent sc:stages) sc.closeNgsResource();
    		showExceptionDialog(owner, err);
    		return Collections.emptyList();
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

    	for(String uri:pathList)
	    	{
	    	SamReader samIn=null;
	    	//try as BAM
	    	if(IOUtil.isUrl(uri) || uri.endsWith(".bam") )
		    	{
		    	try
		    		{
		    		SamReaderFactory srf = SamReaderFactory.makeDefault();
		    		srf.validationStringency(ValidationStringency.LENIENT);
		    		samIn = srf.open(SamInputResource.of(uri));
		    		if(samIn.hasIndex())
		    			{
		    			samIn.close();
		    			stages.add(new BamStageContent(uri));
		    			continue;
		    			}
		    		}
		    	catch(Exception err)
		    		{
		    		LOG.info("not an indexed bam file : "+uri);
		    		}
		    	finally
		    		{
		    		CloserUtil.close(samIn);
		    		}
		    	}
	    	//try as VCF
	    	
	    	if(IOUtil.isUrl(uri) || uri.endsWith(".vcf") || uri.endsWith(".bcf")  || uri.endsWith(".vcf.gz") )
		    	{
		    	VCFFileReader vcfIn=null;
		    	try
		    		{
		    		vcfIn = new VCFFileReader(new File(uri), true);
		    		vcfIn.getFileHeader();
		    		vcfIn.close();
		    		
					stages.add(new VcfStage(uri));
					continue;
		    		}
		    	catch(Exception err)
		    		{
		    		LOG.info("not an indexed vcf file : "+uri);
		    		}
		    	finally
		    		{
		    		CloserUtil.close(vcfIn);
		    		}
		    	}
	    	throw new IOException("Not a NGS file : "+uri);
	    	}
    	return stages;
    	}
    

    public static void main(String[] args) {
    	args=new String[]{
    			"/commun/data/users/lindenb/src/gatk-ui/testdata/mutations.vcf",
    			"/commun/data/users/lindenb/src/gatk-ui/testdata/S1.bam",
    			"/commun/data/users/lindenb/src/gatk-ui/testdata/S2.bam",
    			};
        launch(args);
    }

}

