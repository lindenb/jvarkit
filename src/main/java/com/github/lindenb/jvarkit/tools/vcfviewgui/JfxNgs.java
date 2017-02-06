package com.github.lindenb.jvarkit.tools.vcfviewgui;

import java.io.File;
import java.io.IOException;
import java.io.StringWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.prefs.BackingStoreException;
import java.util.prefs.Preferences;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.tools.vcfstats.VcfStats;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecord.SAMTagAndValue;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SamReaderFactory.Option;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import javafx.application.Application;
import javafx.beans.property.ReadOnlyObjectWrapper;
import javafx.beans.value.ObservableValue;
import javafx.collections.FXCollections;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.geometry.Insets;
import javafx.geometry.Orientation;
import javafx.scene.Parent;
import javafx.scene.Scene;
import javafx.scene.control.Button;
import javafx.scene.control.CheckBox;
import javafx.scene.control.CheckMenuItem;
import javafx.scene.control.Label;
import javafx.scene.control.Menu;
import javafx.scene.control.MenuBar;
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
import javafx.scene.layout.FlowPane;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.VBox;
import javafx.stage.FileChooser;
import javafx.stage.FileChooser.ExtensionFilter;
import javafx.stage.Modality;
import javafx.stage.Stage;
import javafx.stage.Window;
import javafx.stage.WindowEvent;
import javafx.util.Callback;

public class JfxNgs extends Application {
    private static final Logger LOG= Logger.getLogger("JfxNgs");
    private final Preferences preferences ;

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
        abstract void closeNgsResource();
        abstract void reloadData();
        
        protected TableView buildDictTableView(final SAMSequenceDictionary dict)
	        {
			/* build INFO Table */
			final TableView<SAMSequenceRecord> table=new TableView<>(
					dict==null?
					FXCollections.observableArrayList():
					FXCollections.observableArrayList(dict.getSequences())
					);
			
			/* Name */
	        TableColumn<SAMSequenceRecord,String>  scol = new TableColumn<>("Name");
	        scol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<SAMSequenceRecord,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<SAMSequenceRecord, String> param) {
					return new ReadOnlyObjectWrapper<String>(param.getValue().getSequenceName());
					}
				});
	        table.getColumns().add(scol);
	        
			/* Length */
	        TableColumn<SAMSequenceRecord,Integer>  lcol = new TableColumn<>("Length");
	        lcol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<SAMSequenceRecord,Integer>, ObservableValue<Integer>>() {				
				@Override
				public ObservableValue<Integer> call(CellDataFeatures<SAMSequenceRecord, Integer> param) {
					return new ReadOnlyObjectWrapper<Integer>(param.getValue().getSequenceLength());
					}
				});
	        table.getColumns().add(lcol);
			/* Name */
	        scol = new TableColumn<>("Assembly");
	        scol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<SAMSequenceRecord,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<SAMSequenceRecord, String> param) {
					return new ReadOnlyObjectWrapper<String>(param.getValue().getAssembly());
					}
				});
	        table.getColumns().add(scol);
	        return table;
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
            
            
            
            TabPane tabbedPane = new TabPane();
            Tab tab= new Tab("Reads");
            tab.setClosable(false);
            tabbedPane.getTabs().add(tab);
            
            this.recordTable = new TableView<>();
            /** create columns */
            
            /* create READ NAME columns */
            final TableColumn<SAMRecord,String>  readNameCol = new TableColumn<>("Read-Name");
            readNameCol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<SAMRecord,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<SAMRecord, String> param) {
					return new ReadOnlyObjectWrapper<String>(param.getValue().getReadName());
					}
				});
            this.recordTable.getColumns().add(readNameCol);

            /* create READ Flag columns */
            final TableColumn<SAMRecord,Integer>  readFlagCol = new TableColumn<>("Flag");
            readFlagCol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<SAMRecord,Integer>, ObservableValue<Integer>>() {				
				@Override
				public ObservableValue<Integer> call(CellDataFeatures<SAMRecord, Integer> param) {
					return new ReadOnlyObjectWrapper<Integer>(param.getValue().getFlags());
					}
				});
            this.recordTable.getColumns().add(readFlagCol);

            /* create READ Reference columns */
            final TableColumn<SAMRecord,String>  readRefCol = new TableColumn<>("Read-Ref");
            readRefCol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<SAMRecord,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<SAMRecord, String> param) {
					return new ReadOnlyObjectWrapper<String>(param.getValue().getReferenceName());
					}
				});
            this.recordTable.getColumns().add(readRefCol);
            
            /* create readPos columns */
            final TableColumn<SAMRecord,Integer>  readPos = new TableColumn<>("Read-Pos");
            readPos.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<SAMRecord,Integer>, ObservableValue<Integer>>() {				
				@Override
				public ObservableValue<Integer> call(CellDataFeatures<SAMRecord, Integer> param) {
					return new ReadOnlyObjectWrapper<Integer>(param.getValue().getAlignmentStart());
					}
				});
            this.recordTable.getColumns().add(readPos);

            
            /* create Mapq columns */
            final TableColumn<SAMRecord,Integer>  readMAPQCol = new TableColumn<>("MAPQ");
            readMAPQCol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<SAMRecord,Integer>, ObservableValue<Integer>>() {				
				@Override
				public ObservableValue<Integer> call(CellDataFeatures<SAMRecord, Integer> param) {
					return new ReadOnlyObjectWrapper<Integer>(param.getValue().getMappingQuality());
					}
				});
            this.recordTable.getColumns().add(readMAPQCol);
            
            
            /* create Mate Reference columns */
            final TableColumn<SAMRecord,String>  cigarCol = new TableColumn<>("CIGAR");
            cigarCol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<SAMRecord,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<SAMRecord, String> param) {
					return new ReadOnlyObjectWrapper<String>(param.getValue().getCigarString());
					}
				});
            this.recordTable.getColumns().add(cigarCol);

            
            /* create LEN columns */
            final TableColumn<SAMRecord,Integer>  readLEN = new TableColumn<>("LEN");
            readLEN.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<SAMRecord,Integer>, ObservableValue<Integer>>() {				
				@Override
				public ObservableValue<Integer> call(CellDataFeatures<SAMRecord, Integer> param) {
					return new ReadOnlyObjectWrapper<Integer>(param.getValue().getInferredInsertSize());
					}
				});
            this.recordTable.getColumns().add(readLEN);
            
            /* create Mate Reference columns */
            final TableColumn<SAMRecord,String>  mateRefCol = new TableColumn<>("Mate-Ref");
            mateRefCol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<SAMRecord,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<SAMRecord, String> param) {
					return new ReadOnlyObjectWrapper<String>(param.getValue().getMateReferenceName());
					}
				});
            this.recordTable.getColumns().add(mateRefCol);
            
            /* create matePos columns */
            final TableColumn<SAMRecord,Integer>  matePos = new TableColumn<>("Mate-Pos");
            matePos.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<SAMRecord,Integer>, ObservableValue<Integer>>() {				
				@Override
				public ObservableValue<Integer> call(CellDataFeatures<SAMRecord, Integer> param) {
					return new ReadOnlyObjectWrapper<Integer>(param.getValue().getMateAlignmentStart());
					}
				});
            this.recordTable.getColumns().add(matePos);

            
            /* create Sequence columns */
            final TableColumn<SAMRecord,String>  readSequenceCol = new TableColumn<>("SEQ");
            readSequenceCol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<SAMRecord,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<SAMRecord, String> param) {
					return new ReadOnlyObjectWrapper<String>(param.getValue().getReadString());
					}
				});
            this.recordTable.getColumns().add(readSequenceCol);
            
            /* create QUAL columns */
            final TableColumn<SAMRecord,String>  readQualCol = new TableColumn<>("QUAL");
            readQualCol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<SAMRecord,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<SAMRecord, String> param) {
					return new ReadOnlyObjectWrapper<String>(param.getValue().getBaseQualityString());
					}
				});
            this.recordTable.getColumns().add(readQualCol);

            
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
            
            tab=new Tab("Dict",buildDictTableView( this.samReader.getFileHeader().getSequenceDictionary()));
            tab.setClosable(false);
            tabbedPane.getTabs().add(tab);
            
            
            vbox1.getChildren().add(tabbedPane);
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
        }

    private class VcfStage extends StageContent
    	{
    	final VCFFileReader vcfFileReader;
    	VcfStage(final String path) {
    		this.vcfFileReader = new VCFFileReader(new File(path),true);
    		final VCFHeader header=this.vcfFileReader.getFileHeader();
    		
    		final TabPane tabPane=new TabPane();
    		
    		
    		
    		
    		/* build FILTER Table */
    		final TableView<VCFFilterHeaderLine> filtersTable=new TableView<>(FXCollections.observableArrayList(header.getFilterLines()));
    		Tab tab=new Tab("FILTER",filtersTable);
    		tab.setClosable(false);
    		tabPane.getTabs().add(tab);
    		
    		
            final TableColumn<VCFFilterHeaderLine,String>  filterIDcol = new TableColumn<>("ID");
            filterIDcol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<VCFFilterHeaderLine,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<VCFFilterHeaderLine, String> param) {
					return new ReadOnlyObjectWrapper<String>(param.getValue().getID());
					}
				});
            filtersTable.getColumns().add(filterIDcol);
            
    		/* build INFO Table */
    		final TableView<VCFInfoHeaderLine> infosTable=new TableView<>(FXCollections.observableArrayList(header.getInfoHeaderLines()));
    		tab=new Tab("INFO",infosTable);
    		tab.setClosable(false);
    		tabPane.getTabs().add(tab);
    		
    		/* INFO ID */
            final TableColumn<VCFInfoHeaderLine,String>  infoIDcol = new TableColumn<>("ID");
            infoIDcol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<VCFInfoHeaderLine,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<VCFInfoHeaderLine, String> param) {
					return new ReadOnlyObjectWrapper<String>(param.getValue().getID());
					}
				});
            infosTable.getColumns().add(infoIDcol);
    		
    		/* INFO Descr */
            final TableColumn<VCFInfoHeaderLine,String>  infoDescrcol = new TableColumn<>("Description");
            infoDescrcol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<VCFInfoHeaderLine,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<VCFInfoHeaderLine, String> param) {
					return new ReadOnlyObjectWrapper<String>(param.getValue().getDescription());
					}
				});
            infosTable.getColumns().add(infoDescrcol);
           
            
    		/* build FORMAT Table */
    		final TableView<VCFFormatHeaderLine> formatTable=new TableView<>(FXCollections.observableArrayList(header.getFormatHeaderLines()));
    		tab=new Tab("FORMAT",formatTable);
    		tab.setClosable(false);
    		tabPane.getTabs().add(tab);
    		
    		/* INFO ID */
            final TableColumn<VCFFormatHeaderLine,String>  formatIDcol = new TableColumn<>("ID");
            formatIDcol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<VCFFormatHeaderLine,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<VCFFormatHeaderLine, String> param) {
					return new ReadOnlyObjectWrapper<String>(param.getValue().getID());
					}
				});
            formatTable.getColumns().add(formatIDcol);
    		
    		/* INFO Descr */
            final TableColumn<VCFFormatHeaderLine,String>  formatDescrcol = new TableColumn<>("Description");
            formatDescrcol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<VCFFormatHeaderLine,String>, ObservableValue<String>>() {				
				@Override
				public ObservableValue<String> call(CellDataFeatures<VCFFormatHeaderLine, String> param) {
					return new ReadOnlyObjectWrapper<String>(param.getValue().getDescription());
					}
				});
            formatTable.getColumns().add(formatDescrcol);

            
            
            tab=new Tab("Dict",buildDictTableView( header.getSequenceDictionary()));
            tab.setClosable(false);
            tabPane.getTabs().add(tab);

            
            
            
    		}
    	@Override
        void closeNgsResource() {
        	CloserUtil.close(this.vcfFileReader);
        	}
    	@Override
        void reloadData()
    		{
        
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
    
    private void showPreferenceDialoge(Window parentStage)
	    {
	    Stage dialog = new Stage();
	
	     
	     dialog.setTitle("Preferences");
		 dialog.initOwner(parentStage);
		 dialog.initModality(Modality.APPLICATION_MODAL); 
		 dialog.showAndWait();
	    }
    
    @Override
    public void start(Stage primaryStage) throws Exception {
        final Parameters params = this.getParameters();
        List<String> optargs = params.getUnnamed();
        /*
        Parent rootNode=null;
        for(int i=0;i< optargs.size();++i)
            {
            Parent node = createNodeForParam(optargs.get(i));
            if(rootNode==null) {
                rootNode = node;
                }
            else
                {
                final Stage stage=new Stage(primaryStage.getStyle());
                stage.setScene(new Scene(node));
                stage.show();
                }
            }
        if(rootNode==null) {
            rootNode = null;//TODO
            }
        primaryStage.setScene(new Scene(rootNode));
        primaryStage.show();
        */
        BamStageContent stage= new BamStageContent("/home/lindenb/src/mod_bio/examples/rf.bam");
        stage.show();
        }

    private void unregisterStage(StageContent s) {

    }

    private Collection<StageContent> openNgsFiles(final Window owner)
    	{
    	FileChooser fc=new FileChooser();
    	fc.setSelectedExtensionFilter(new ExtensionFilter("NGS Files","bam","vcf","vcf.gz","list"));
    	List<File> selFiles = fc.showOpenMultipleDialog(owner);
    	List<StageContent> stages =new ArrayList<>();
    	try 
    		{
    		
	    	for(File f:selFiles)
	    		{
	    		stages.addAll( openNgsFiles(f.getPath()) );
	    		}
	    	return stages;
    		}
    	catch(final IOException err)
    		{
    		for(StageContent sc:stages) sc.hide();//TODO dispose stuff
    		return Collections.emptyList();
    		}		
    	}

    private Collection<StageContent> openNgsFiles(String path) throws IOException
    	{
    	if(!IOUtil.isUrl(path) && path.endsWith(".list"))
    		{
    		List<StageContent> stages =new ArrayList<>();
    		for(final String s:Files.lines(Paths.get(path)).filter(L-> ! (L.isEmpty() && L.startsWith("#"))).collect(Collectors.toList()))
    			{
    			stages.addAll(openNgsFiles(s));
    			}
    		return stages;
    		}
    	SamReader samIn=null;
    	//try as BAM
    	try
    		{
    		SamReaderFactory srf = SamReaderFactory.makeDefault();
    		srf.validationStringency(ValidationStringency.LENIENT);
    		samIn = srf.open(SamInputResource.of(path));
    		if(samIn.hasIndex())
    			{
    			samIn.close();
    			return Collections.singleton(new BamStageContent(path));
    			}
    		}
    	finally
    		{
    		CloserUtil.close(samIn);
    		}
    	//try as VCF
    	
    	VCFFileReader vcfIn=null;
    	try
    		{
    		vcfIn = new VCFFileReader(new File(path), true);
    		vcfIn.getFileHeader();
    		vcfIn.close();
    		
			return Collections.singleton(new VcfStage(path));
    		}
    	finally
    		{
    		CloserUtil.close(vcfIn);
    		}
    	//TODO throw new IOException("Cannot open "+path);
    	}
    

    public static void main(String[] args) {
        launch(args);
    }

}

