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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.function.Predicate;
import java.util.function.Supplier;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;
import java.util.stream.Collectors;

import javax.script.CompiledScript;

import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.AlleleFrequencyChartFactory;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.ChartFactory;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.GenotypeTypeChartFactory;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.TiTvChartFactory;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.VariantContextChartFactory;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.VariantDepthChartFactory;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.VariantQualChartFactory;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.VariantTypeChartFactory;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParserFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParserFactory;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder.OutputType;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import javafx.beans.property.BooleanProperty;
import javafx.beans.property.ReadOnlyObjectWrapper;
import javafx.beans.property.SimpleBooleanProperty;
import javafx.beans.value.ObservableValue;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.geometry.Insets;
import javafx.geometry.Orientation;
import javafx.geometry.Rectangle2D;
import javafx.scene.Scene;
import javafx.scene.chart.PieChart;
import javafx.scene.control.Alert;
import javafx.scene.control.Button;
import javafx.scene.control.ButtonType;
import javafx.scene.control.CheckBox;
import javafx.scene.control.ContextMenu;
import javafx.scene.control.Label;
import javafx.scene.control.MenuItem;
import javafx.scene.control.ScrollPane;
import javafx.scene.control.Separator;
import javafx.scene.control.SeparatorMenuItem;
import javafx.scene.control.SpinnerValueFactory;
import javafx.scene.control.SplitPane;
import javafx.scene.control.Tab;
import javafx.scene.control.TabPane;
import javafx.scene.control.TabPane.TabClosingPolicy;
import javafx.scene.control.TableColumn;
import javafx.scene.control.TableView;
import javafx.scene.control.TextField;
import javafx.scene.control.cell.CheckBoxTableCell;
import javafx.scene.control.cell.PropertyValueFactory;
import javafx.scene.control.Alert.AlertType;
import javafx.scene.control.TableColumn.CellDataFeatures;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.FlowPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.VBox;
import javafx.stage.FileChooser;
import javafx.stage.Screen;
import javafx.stage.WindowEvent;
import javafx.stage.FileChooser.ExtensionFilter;
import javafx.util.Callback;

public class VcfStage extends NgsStage<VCFHeader,VariantContext> {
    static final String SPINNER_VALUE_KEY="vcf.spinner.value";
	static final int DEFAULT_VCF_RECORDS_COUNT= 100;
    static final ExtensionFilter EXTENSION_FILTER=new ExtensionFilter("Variant Files", ".vcf",".vcf.gz");
    
    /** variant oriented chart-factories */
    private static final List<Supplier<ChartFactory<VariantContext>>> VARIANT_CHART_FACTORIES= Arrays.asList(
		()->new VariantTypeChartFactory(),
		()->new GenotypeTypeChartFactory(),
		()->new AlleleFrequencyChartFactory(),
		()->new VariantQualChartFactory(),
		()->new TiTvChartFactory(),
		()->new VariantDepthChartFactory()
		);
 
    /** base class to represent local file or remote tabix file */
    public static interface VariantFileReader
    	{
    	public VCFHeader getFileHeader(); 
    	public CloseableIterator<VariantContext> iterator() throws IOException; 
    	public CloseableIterator<VariantContext> iterator(final String contig,int start,int end) throws IOException; 
    	public void close();
    	}
    
    /** class used to display samples */
    public static class SampleDef
    	{
    	private final String name;
    	private final String description;
        private final BooleanProperty _displayedProperty = new SimpleBooleanProperty(true);
        public BooleanProperty displayedProperty()  { return _displayedProperty; }
        public boolean isDisplayed() { return displayedProperty().get(); }
        public void setDisplayed(boolean v) { displayedProperty().set(v); }
    
    	SampleDef(final String name,final String description)
	    	{
	    	this.name=name;
	    	this.description=description;
	    	}
    
    	public String getName() {
			return name;
			}
    	public String getDescription() {
			return description;
			}
    	}
    
    private final Map<String,SampleDef> name2sampledef;
    
    /** Misc VCF tools that will be injected in the nashorn context */
    public static class VcfTools
    	{
    	private final AnnPredictionParser annParser;
    	private final VepPredictionParser vepParser;
    	VcfTools(final VCFHeader header) {
    		this.annParser  = new AnnPredictionParserFactory(header).get();
    		this.vepParser  = new VepPredictionParserFactory(header).get();
    	}
    	
    	public List<AnnPredictionParser.AnnPrediction> getAnnPredictions(final VariantContext ctx) {
    		return this.annParser.getPredictions(ctx);
    		}
    	
    	public List<VepPredictionParser.VepPrediction> getVepPredictions(final VariantContext ctx) {
    		return this.vepParser.getPredictions(ctx);
    		}

    	
    	public boolean isMendelianIncompatibility(final Genotype child,final Genotype parent)
    		{
    		if(child==null || parent==null) return false;
    		if(child.isNoCall() || parent.isNoCall()) return false;
    		if(child.getPloidy()!=2 || parent.getPloidy()!=2) return false;
    		for(final Allele childAllele:child.getAlleles())
    			{
    			if(parent.getAlleles().contains(childAllele)) return false;
    			}
    		
    		return true;
    		}
    	public boolean isMendelianIncompatibility(final Genotype child,final Genotype father,final Genotype mother)
			{
			if(child==null || child.isNoCall()) return false;
			if(father==null || father.isNoCall()) {
				return this.isMendelianIncompatibility(child,mother);
				}
			if(mother==null || mother.isNoCall()) {
				return this.isMendelianIncompatibility(child,father);
				}
			final Allele alleles[]=new Allele[2];
			for(final Allele af:father.getAlleles())
			{
				alleles[0]=af;
				for(final Allele am:mother.getAlleles())
				{
					alleles[1]=am;
					final Genotype sim = new GenotypeBuilder(child.getSampleName()).alleles(Arrays.asList(alleles)).make();
					if(sim.sameGenotype(sim, true)) return false;
				}	
			}
			
			return true;
			}
    	}
    
    /** instance of JavascriptFilter for VCFs */
    private static class VcfJavascripFilter
	extends JavascriptFilter<VCFHeader,VariantContext>
		{
    	/** set the FILTER column instead of deleting the variant */
    	private final String filter;
		VcfJavascripFilter(
				final String filter,
				final VCFHeader header,
				final Optional<CompiledScript> compiledScript)
			{
			super(header,compiledScript);
			super.bindings.put("tools", new VcfTools(header));
			this.filter=filter;
			}
		@Override public VariantContext eval(final VariantContext ctx)
			{
			if(!super.compiledScript.isPresent()) return ctx;
			super.bindings.put("variant", ctx);
			if(super.accept())
				{
				return ctx;
				}
			else
				{
				if(filter==null || filter.trim().isEmpty()) return null;
				return new VariantContextBuilder(ctx).filter(this.filter).make();
				}
			}
		}
    
    private class VcfQualityStage
		extends AbstractQualityStage
			{
    		private class ScanVcfThread extends ScanThread
    			{
    			ScanVcfThread(final ChartFactory<VariantContext> factory,
    						final VcfFile vcfinput,
    						final Optional<CompiledScript> compiledScript,
    						final Predicate<VariantContext> otherFilters)
    				{
    				super(factory,vcfinput,compiledScript,otherFilters);
    				}
        		        		
    			@Override
    			public void run() {
    				CloseableIterator<VariantContext> iter=null;
    				Optional<VcfJavascripFilter> javascriptFilter=Optional.empty();
    				try 
	    				{
    						if(this.compiledScript.isPresent() )
    							{
    							javascriptFilter=Optional.of(
    									new VcfJavascripFilter("", super.ngsReader.getHeader(), compiledScript));
    							}
    						
	    					iter = super.ngsReader.iterator();
	    					
	    					while(!kill_flag && iter.hasNext())
	    						{
	    						final VariantContext ctx=iter.next();
	    						
	    						nItems++;
	    						if(javascriptFilter.isPresent() &&
	    							javascriptFilter.get().eval(ctx)==null) continue;
	    							
	    						if(!super.otherFilters.test(ctx)) continue;
	    						super.factory.visit(ctx);
	    						update();
	    						}
	    					iter.close();
	    					super.ngsReader.close();
	    				
	    					
	    					if(javascriptFilter.isPresent() && javascriptFilter.get().encounteredException.isPresent())
	    						{
	    						this.encounteredException=javascriptFilter.get().encounteredException;
	    						}
	    					atEnd();
	    				
	    				}
    				catch(final Throwable err)
    					{
    					super.onError(err);
    					}
    				finally
    					{
    					CloserUtil.close(iter);
    					CloserUtil.close(super.ngsReader);
    					}
    				}
    			
    			}
    		
    		VcfQualityStage(final ChartFactory<VariantContext> factory,
    				final VcfFile vcfinput,
    				final Optional<CompiledScript> compiledScript,
    				final Predicate<VariantContext> otherFilters
    				)
	    		{
    			super(factory,vcfinput,compiledScript,otherFilters);
	    		}

			@Override
			protected ScanVcfThread createThread(
					final ChartFactory<VariantContext> factory,
					NgsFile<VCFHeader,VariantContext> ngsfile,
					Optional<CompiledScript> compiledScript,
					final Predicate<VariantContext> otherFilters
					) {
				return new ScanVcfThread(factory,(VcfFile)ngsfile, compiledScript,otherFilters);
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
		boolean accept(final Pattern pat)
			{
			if(pat==null) return true;
			if(pat.matcher(key).find()) return true;
			if(value!=null && pat.matcher(String.valueOf(value)).find()) return true;
			return false;
			}
		}
   
    
    
	private final TextField javascriptFILTERfield=new TextField("");
	private final TableView<VariantContext> variantTable;
	private final TableView<Genotype> genotypeTable;
	private final TableView<InfoTableRow> infoTableRow;
	private final TableView<SampleDef> sampleTable = new TableView<>();
	/** annotation SNPEFF table ROW */
	private final TableView<AnnPredictionParser.AnnPrediction> annPredictionTable;
	private final TableView<VepPredictionParser.VepPrediction> vepPredictionTable;
	private final TableView<String> filterTableRow;
	private final TableView<Allele> allelesTable;
	private final BorderPane genotypeChartPane;
	private final CheckBox cboxShowHomRef=new CheckBox("HomRef");
	private final CheckBox cboxShowNoCall=new CheckBox("NoCall");
	private final CheckBox cboxShowFiltered=new CheckBox("Filtered");
	private final TextField tfFilterInfo=new TextField();
	private final AnnPredictionParser annPredictionParser;
	private final VepPredictionParser vepPredictionParser;
	VcfStage(final JfxNgs owner,final VcfFile vcfFile) throws IOException {
		super(owner,vcfFile);
		final VCFHeader header= vcfFile.getHeader();
		
        this.annPredictionParser=new AnnPredictionParserFactory(header).get();
        this.vepPredictionParser=new VepPredictionParserFactory(header).get();
        
        //selectFlagMenu.getItems().addAll(flag2filterOutMenuItem.values());
        final VBox vbox1 = new VBox();
        vbox1.getChildren().add(super.menuBar);
        
        FlowPane top1= new FlowPane();
        vbox1.getChildren().add(top1);
        top1.getChildren().add(new Label("GoTo:"));
        top1.getChildren().add(super.gotoField);
        final Button gotoButton=new Button("Go");
        gotoButton.setOnAction(A->reloadData());
        top1.getChildren().add(gotoButton);
        top1.getChildren().add(new Separator(Orientation.VERTICAL));
        top1.getChildren().add(new Label("Limit:"));
        
        int number_of_items_in_spinner;
        try {
        	number_of_items_in_spinner= owner.preferences.getInt(SPINNER_VALUE_KEY, DEFAULT_VCF_RECORDS_COUNT);
			if(number_of_items_in_spinner<0) number_of_items_in_spinner=0;
			if(number_of_items_in_spinner>100000) number_of_items_in_spinner=100000;
        	}
        catch(Exception err)
        	{
        	number_of_items_in_spinner = DEFAULT_VCF_RECORDS_COUNT;
        	}
        
        super.maxItemsLimitSpinner.setValueFactory(
        		new SpinnerValueFactory.IntegerSpinnerValueFactory(0,100000,
        				number_of_items_in_spinner
        				));;

        top1.getChildren().add(super.maxItemsLimitSpinner);
        top1.getChildren().add(new Separator(Orientation.VERTICAL));
		
        top1.getChildren().add(createIgvButton());
        
        this.gotoField.setOnAction(A->reloadData());
		final TabPane tabPane=new TabPane();
		tabPane.setPadding(new Insets(10, 10, 10, 10));
		
		
        vbox1.getChildren().add(super.seqDictionaryCanvas);
		vbox1.getChildren().add(tabPane);
		
		
		/*GridPane gridPane = new GridPane();
		gridPane.setPadding(new Insets(10, 10, 10, 10));
		gridPane.setVgap(4);
		gridPane.setHgap(4);*/

        
		final SplitPane split1=new SplitPane();
		split1.setOrientation(Orientation.HORIZONTAL);
		
		/* build variant table */
		this.variantTable = this.buildVariantTable();
		split1.getItems().add(this.variantTable);
		
		/* build genotype table */
		
		this.genotypeTable =this.buildGenotypeTableRow(header);
		FlowPane flow3 = new FlowPane(this.cboxShowHomRef,this.cboxShowNoCall,this.cboxShowFiltered);
		this.cboxShowHomRef.setSelected(true);
		this.cboxShowNoCall.setSelected(true);
		this.cboxShowFiltered.setSelected(true);
		EventHandler<ActionEvent> repaintGTTable=new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent event) {
				refreshGenotypeTable(variantTable.getSelectionModel().getSelectedItem());
				}
			};
		this.cboxShowHomRef.setOnAction(repaintGTTable);
		this.cboxShowNoCall.setOnAction(repaintGTTable);
		this.cboxShowFiltered.setOnAction(repaintGTTable);

		
		if(header.getNGenotypeSamples()>0)
			{
			final BorderPane pane2=new BorderPane(this.genotypeTable);
			pane2.setTop(flow3);
			//GridPane.setConstraints(pane2,5, 0,5,14); 
			split1.getItems().add(pane2);
			}
		

		/* filter table */
		this.filterTableRow = this.buildFilterTable();
		//GridPane.setConstraints( this.filterTableRow,0, 10,3,1);
		//gridPane.getChildren().add(this.filterTableRow);
		//split2.getItems().add(this.filterTableRow);
		/* genotype pane */
		
		
		this.genotypeChartPane = new BorderPane(this.makeGenotypePie(null));
		if(header.getNGenotypeSamples()>0)
			{
			//GridPane.setConstraints( this.genotypeChartPane,0, 12,3,4);
	    	//gridPane.getChildren().add(this.genotypeChartPane);
			
			this.name2sampledef=new LinkedHashMap<>(header.getNGenotypeSamples());
			for(final String sampleName:header.getSampleNamesInOrder())
				{
				String desc=null;
				for(final VCFHeaderLine h:header.getOtherHeaderLines()) {
					final String key = h.getKey();
					if(!key.equals("Sample")) continue;
					final String value =h.getValue();
					if(value.contains("ID="+sampleName+",") || value.contains("ID="+sampleName+">"))
						{
						desc = value;
						break;
						}
					}
				final SampleDef sampledef = new SampleDef(sampleName, desc);
				this.name2sampledef.put(sampleName, sampledef);
				}
		
			}
		else
			{
			this.name2sampledef=Collections.emptyMap();
			}
		
		/* build info Table table */
		this.infoTableRow = this.buildInfoTableRow();
		this.annPredictionTable = this.buildAnnTableRow(this.annPredictionParser);
		this.vepPredictionTable = this.buildVepTableRow(this.vepPredictionParser);
		flow3 = new FlowPane(new Label("INFO Filter:"),tfFilterInfo);
		tfFilterInfo.setPrefColumnCount(10);
		tfFilterInfo.setOnAction(AE->reloadInfoTable(variantTable.getSelectionModel().getSelectedItem()));
		final TabPane tabPane2=new TabPane();
		tabPane2.setTabClosingPolicy(TabClosingPolicy.UNAVAILABLE);
		Tab tab=new Tab("INFO",this.infoTableRow);
		tab.setClosable(false);
		tabPane2.getTabs().add(tab);
		if(this.annPredictionParser.isValid())
			{
			tab=new Tab("ANN",this.annPredictionTable);
			tab.setClosable(false);
			tabPane2.getTabs().add(tab);
			}
		if(this.vepPredictionParser.isValid())
			{
			tab=new Tab("VEP",this.vepPredictionTable);
			tab.setClosable(false);
			tabPane2.getTabs().add(tab);
			}
		tabPane2.getTabs().add(new Tab("FILTER",this.filterTableRow));
		if(header.getNGenotypeSamples()>0)
			{
			this.genotypeChartPane.setPadding(new Insets(5));
			tabPane2.getTabs().add(new Tab("G-Chart",this.genotypeChartPane));
			}
		this.allelesTable = buildAlleleTable();
		tabPane2.getTabs().add(new Tab("ALLELES",this.allelesTable));
		
		BorderPane pane2=new BorderPane(tabPane2);
		pane2.setTop(flow3);
		
		
		
		final SplitPane split2=new SplitPane();
		split2.setOrientation(Orientation.VERTICAL);
		split2.getItems().addAll(split1,pane2);
		//GridPane.setConstraints( this.infoTableRow,3, 10,8,5); // column=3 row=1
		//gridPane.getChildren().add(this.infoTableRow);
		
		
				
       
		
        
        //vbox1.getChildren().add(gridPane);

		tab=new Tab("Variants", split2);
		tab.setClosable(false);
		tabPane.getTabs().add(tab);
		tabPane.getTabs().add(buildInfoHeaderTab(header));
		if(header.getNGenotypeSamples()>0)
			{
			tabPane.getTabs().add(buildFormatHeaderTab(header));
			}
		tabPane.getTabs().add(buildFilterHeaderTab(header));
        tabPane.getTabs().add(buildDictTab( header.getSequenceDictionary()));
        tabPane.getTabs().add(buildJavascriptPane());
        if(header.getNGenotypeSamples()>0)
	        {
	        tabPane.getTabs().add(buildSampleTab());
	        }
        /* register selection handler */
        /* when a read is selected update the flagsTable and metaDataTable */
        this.variantTable.getSelectionModel().selectedItemProperty().addListener((obs, oldSelection, newSelection) -> {
            fireSelectedVariantChanged(newSelection);
        	});
        
        
        final FlowPane bottom=new FlowPane(super.messageLabel);
        vbox1.getChildren().add(bottom);
        final Rectangle2D primaryScreenBounds = Screen.getPrimary().getVisualBounds();

        final Scene scene = new Scene(vbox1,
        		primaryScreenBounds.getWidth()-200,
        		primaryScreenBounds.getHeight()-200);
       
        this.setScene(scene);
        
        fileMenu.getItems().addAll(
        		menuForSavingTable("Variants",this.variantTable),
        		menuForSavingTable("INFO",this.infoTableRow),
        		menuForSavingTable("FILTER",this.filterTableRow),
        		menuForSavingTable("Genotypes",this.genotypeTable)
        		);
        
        /* fill stats menu */
        final Supplier<List<VariantContext>> variantsProvider=()->this.variantTable.getItems();
        
        for(final Supplier<ChartFactory<VariantContext>> supplier: VARIANT_CHART_FACTORIES)
	        {
        	final ChartFactory<VariantContext> factory = supplier.get();
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
        for(final Supplier<ChartFactory<VariantContext>> supplier: VARIANT_CHART_FACTORIES)
	        {
	    	final ChartFactory<VariantContext> factory = supplier.get();
	    	final MenuItem menuItem=new MenuItem("Whole"+factory.getName());
	    	menuItem.setOnAction(new EventHandler<ActionEvent>() {
				@Override
				public void handle(ActionEvent event) {
					doMenuShowWholeStats(factory);
				}
			});
	    	statsMenu.getItems().add(menuItem);
	        }

        
        this.addEventHandler(
    			WindowEvent.WINDOW_CLOSE_REQUEST ,new EventHandler<WindowEvent>() {
                    @Override
                    public void handle(final WindowEvent event) {
                    owner.preferences.putInt(SPINNER_VALUE_KEY,maxItemsLimitSpinner.getValue().intValue());
                    }
                });

		}
	

	@Override
	protected void doMenuSaveAs()
		{
		final FileChooser fc= owner.newFileChooser();
    	fc.setSelectedExtensionFilter(EXTENSION_FILTER);
		final File saveAs= owner.updateLastDir(fc.showSaveDialog(this));
		if(saveAs==null) return;
		if(!saveAs.getName().endsWith(".vcf.gz"))
			{
			final Alert alert=new Alert(AlertType.ERROR, "Output should end with .vcf.gz", ButtonType.OK);
			alert.showAndWait();
			return;
			}
		
    	
    	VcfJavascripFilter javascriptFilter=null;
    	if(	this.owner.javascriptCompiler.isPresent() &&
    		!this.javascriptArea.getText().trim().isEmpty())
    		{
    		try
    			{
    			javascriptFilter=new VcfJavascripFilter(
    				this.javascriptFILTERfield.getText().trim(),
    				this.getVcfFile().getHeader(),
    				Optional.of(this.owner.javascriptCompiler.get().compile(this.javascriptArea.getText()))
    				);
    			
    			}
    		catch(final Exception err)
    			{
    			JfxNgs.showExceptionDialog(this, err);
    			return;
    			}
    		}

		
		VariantContextWriter w=null;
		CloseableIterator<VariantContext> iter=null;
		try
			{
			final VariantContextWriterBuilder vcwb=new VariantContextWriterBuilder();
			vcwb.setOutputFileType(OutputType.BLOCK_COMPRESSED_VCF);
			vcwb.setOutputFile(saveAs);
			w= vcwb.build();
			final VCFHeader header2= new VCFHeader(this.getVcfFile().getHeader());
        	if(!javascriptFilter.filter.isEmpty())
        		{
        		header2.addMetaDataLine(new VCFFilterHeaderLine(javascriptFilter.filter, "Set by User in JfxNgs:"+
        				this.javascriptArea.getText().replaceAll("[\n\t\r ]+"," ")
        				));
        		}
			w.writeHeader(header2);
			iter= this.getVcfFile().iterator();
			while(iter.hasNext())
				{
				VariantContext ctx=iter.next();
				if(javascriptFilter!=null)
        			{
        			ctx = javascriptFilter.eval(ctx);
        			if(ctx==null) continue;
        			}
				w.add(ctx);
				}
			w.close();
			}
		catch(final Exception err)
			{
			JfxNgs.showExceptionDialog(this, err);
			return;
			}
		finally
			{
			CloserUtil.close(w);
			}    		
		}
	
	private Tab buildJavascriptPane()
		{
		final ScrollPane scroll=new ScrollPane(super.javascriptArea);
		scroll.setFitToWidth(true);
		scroll.setFitToHeight(true);
		final BorderPane pane=new BorderPane(scroll);
		if(this.owner.javascriptCompiler.isPresent()) {
    		final HBox flowPane = new HBox(
    				new Label("set following FILTER on rejection:"),
    				javascriptFILTERfield
    				);
    		this.javascriptFILTERfield.setPrefColumnCount(30);
    		this.javascriptFILTERfield.setPromptText("If not empty , don't discard the variant but set the FILTER");
    		flowPane.getChildren().addAll(super.makeJavascriptButtons());
    		pane.setTop(flowPane);
    		}
		final Label helpLabel=new Label("The script injects:\n"+
				"* header ( htsjdk.variant.vcf.VCFHeader )\n"+
				"* variant ( htsjdk.variant.variantcontext.VariantContext )\n"+
				"The script should return a boolean: true (accept variant) or false (reject variant)"
				);
		helpLabel.setWrapText(true);
		pane.setBottom(helpLabel);
		
		final Tab tab=new Tab(JAVASCRIPT_TAB_KEY,pane);
		tab.setClosable(false);
		return tab;
		}
	
	/** build FILTER table */
	private TableView<String> buildFilterTable()
		{
		final TableView<String> table=new TableView<>();
		final TableColumn<String,String>  scol = new TableColumn<>("Filter");
		scol.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<String,String>, ObservableValue<String>>() {				
			@Override
			public ObservableValue<String> call(CellDataFeatures<String, String> param) {
				return new ReadOnlyObjectWrapper<String>(param.getValue());
				}
			});
        table.getColumns().add(scol);
		
        table.setPlaceholder(new Label("No Variant or Variant contains no Filter"));
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
        table.getColumns().add(formatIntegerColumn(makeColumn("POS", V->V.getStart())));
        table.getColumns().add(makeColumn("ID", V->V.hasID()?V.getID():null));
        table.getColumns().add(makeColumn("REF", V->V.getReference().getDisplayString()));
        table.getColumns().add(makeColumn("ALT", V->V.getAlternateAlleles().stream().map(A->A.getDisplayString()).collect(Collectors.joining(","))));
        table.getColumns().add(makeColumn("FILTER", V->V.getFilters().stream().collect(Collectors.joining(","))));
        table.getColumns().add(makeColumn("QUAL", V->V.hasLog10PError()?V.getPhredScaledQual():null));
        
        table.setPlaceholder(new Label("No Variant."));
        
        
        final ContextMenu ctxMenu=new ContextMenu();
        
        MenuItem menuItem=new MenuItem("dbSNP...");
        menuItem.setOnAction(AE->{
        	final VariantContext ctx=table.getSelectionModel().getSelectedItem();
        	if(ctx==null || !ctx.hasID() || !ctx.getID().matches("rs[0-9]+")) return ;
        	//http://stackoverflow.com/questions/16604341
        	VcfStage.this.owner.getHostServices().showDocument("https://www.ncbi.nlm.nih.gov/snp/"+ctx.getID().substring(2));
        });
        ctxMenu.getItems().add(menuItem);
        ctxMenu.getItems().addAll(super.buildItemsForContextMenu());
        
        for(final String build:new String[]{"human"}) {
	        menuItem=new MenuItem("Prediction Ensembl REST ["+build+"]");
	        menuItem.setOnAction(AE->{
	        	final VariantContext ctx=table.getSelectionModel().getSelectedItem();
	        	if(ctx==null) return;
	        	
	        	for(final Allele a: ctx.getAlternateAlleles())
	        		{
	        		if(a.isReference()) continue;
	        		if(a.isSymbolic()) continue;
	        		if(a.isNoCall()) continue;
	        	
		        	VcfStage.this.owner.getHostServices().showDocument(
		        		"http://rest.ensembl.org/vep/"+build+"/region/"
		        				+ JfxNgs.ContigToEnseml.apply(ctx.getContig())
		        				+"%3A"+ctx.getStart()+"-"+ctx.getEnd()+"%3A1%2F"+ a.getDisplayString() +"?content-type=text%2Fxml");
	        		}
	        	});
	        ctxMenu.getItems().add(menuItem);
	        }
    	for(final String database : new String[]{"Exac","gnomAD"}) {
        menuItem=new MenuItem("Open Variant (ALT) in "+database+" ... ");
		menuItem.setOnAction(AE->{
        	final VariantContext ctx=table.getSelectionModel().getSelectedItem();
        	if(ctx==null) return;
        	for(final Allele a: ctx.getAlternateAlleles())
	    		{
	    		if(a.isReference()) continue;
	    		if(a.isSymbolic()) continue;
	    		if(a.isNoCall()) continue;
	    		
	        	VcfStage.this.owner.getHostServices().showDocument(
	        		"http://"+ database.toLowerCase() +".broadinstitute.org/variant/"
	        				+ JfxNgs.ContigToEnseml.apply(ctx.getContig())
	        				+ "-"+ctx.getStart()+"-"
	        				+ ctx.getReference().getDisplayString()+"-"
	        				+a.getDisplayString()
	        				);
	    		}
			});
		ctxMenu.getItems().add(menuItem);
		}
        
        table.setContextMenu(ctxMenu);
        return table;
		}

	/** build INFO table */
	private TableView<Allele> buildAlleleTable()
		{
		final TableView<Allele> table=new TableView<>();
		table.getColumns().add(makeColumn("REF",A->A.isReference()?"*":null));
		table.getColumns().add(makeColumn("Sym.",A->A.isSymbolic()?"*":null));
		table.getColumns().add(makeColumn("Bases.",A->A.getBaseString()));
		table.getColumns().add(makeColumn("Length.",A->{if(A.isSymbolic()) return (Integer)null;return A.length();}));
		table.setPlaceholder(new Label("No Allele."));
		return table;
		}

	
	/** build INFO table */
	private TableView<InfoTableRow> buildInfoTableRow()
		{
		final TableView<InfoTableRow> table=new TableView<>();
		table.getColumns().add(makeColumn("Key", R->R.key));
		table.getColumns().add(makeColumn("Index", R->R.index));
		
		final TableColumn<InfoTableRow, Object> t=makeColumn("Value", R->R.value);
		t.setPrefWidth(300.0);
		table.getColumns().add(t);
		table.setPlaceholder(new Label("No INFO."));
		return table;
		}
     
	/** build table of Sample */
	private Tab buildSampleTab()
		{
		final boolean editable = this.name2sampledef.size()>1;
		final BorderPane pane=new BorderPane(this.sampleTable);
		pane.setPadding(new Insets(10));
		

		this.sampleTable.setEditable(editable);
		this.sampleTable.getColumns().add(makeColumn("Name", R->R.getName()));
		final TableColumn<SampleDef, Boolean> selectCol = new TableColumn<>("Display");
		selectCol.setCellValueFactory(new PropertyValueFactory<>("displayed"));
		selectCol.setCellFactory(col -> {
            CheckBoxTableCell<SampleDef, Boolean> cell = new CheckBoxTableCell<>(index -> {
               return sampleTable.getItems().get(index).displayedProperty();
            });
            return cell ;
        });
		selectCol.setEditable(editable);
		this.sampleTable.getColumns().add(selectCol);
		this.sampleTable.getColumns().add(makeColumn("Description", R->R.getDescription()));
		this.sampleTable.setPlaceholder(new Label("No Sample."));
		this.sampleTable.getItems().addAll(this.name2sampledef.values());
		
		final HBox top=new HBox();
		top.setPadding(new Insets(10));
		Button but= new Button("Select all");
		but.setOnAction(AE->{
			for(final SampleDef def:name2sampledef.values())
				{
				def.setDisplayed(true);
				}
			if(getCurrentSelectedItem().isPresent())
				{
				refreshGenotypeTable(getCurrentSelectedItem().get());
				}
			});
		top.getChildren().add(but);
		
		but= new Button("Unselect all");
		but.setOnAction(AE->{
			for(final SampleDef def:name2sampledef.values())
				{
				def.setDisplayed(false);
				}
			if(getCurrentSelectedItem().isPresent())
				{
				refreshGenotypeTable(getCurrentSelectedItem().get());
				}
			});
		top.getChildren().add(but);
		
		pane.setTop(top);
		final Tab tab=new Tab("Samples",pane);
		tab.setClosable(false);
		return tab;
		}
	
	private TableView<AnnPredictionParser.AnnPrediction> buildAnnTableRow(final AnnPredictionParser parser)
		{
		final TableView<AnnPredictionParser.AnnPrediction> table=new TableView<>();
		if(parser.isValid())
			{
			table.getColumns().add(makeColumn("SO", P->P.getSOTermsString()));
			table.getColumns().add(makeColumn("Allele", P->P.getAllele()));
			table.getColumns().add(makeColumn("Impact", P->P.getPutativeImpact()));
			table.getColumns().add(makeColumn("GeneId", P->P.getGeneId()));
			table.getColumns().add(makeColumn("Feature", P->P.getFeatureType()));
			table.getColumns().add(makeColumn("Biotype", P->P.getTranscriptBioType()));
			table.getColumns().add(makeColumn("HGVsc", P->P.getHGVSc()));
			table.getColumns().add(makeColumn("Rank", P->P.getRank()));
			table.getColumns().add(makeColumn("cDNA-pos", P->P.getCDNAPos()));
			table.getColumns().add(makeColumn("CDS-pos", P->P.getCDSPos()));
			table.getColumns().add(makeColumn("AA-pos", P->P.getAAPos()));
			table.getColumns().add(makeColumn("Distance", P->P.getDistance()));
			table.getColumns().add(makeColumn("Msg", P->P.getMessages()));
			}
		table.setPlaceholder(new Label("No ANN prediction available"));
		return table;
		}
	private TableView<VepPredictionParser.VepPrediction> buildVepTableRow(final VepPredictionParser parser)
		{
		final TableView<VepPredictionParser.VepPrediction> table=new TableView<>();
		if(parser.isValid())
			{
			for(final String col:parser.getCategories())
				{
				table.getColumns().add(makeColumn(col, P->P.get(col)));
				}
			}
		table.setPlaceholder(new Label("No VEP prediction available"));
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
		
		table.setPlaceholder(new Label("No Genotype."));
		return table;
		}

	
	
    private void fireSelectedVariantChanged(final VariantContext ctx)
    	{
    	refreshGenotypeTable(ctx);
    	reloadInfoTable(ctx);
    	reloadAlleleTable(ctx);
    	if(ctx!=null)
    		{
    		super.seqDictionaryCanvas.setSelectInterval(new Interval(ctx.getContig(), ctx.getStart(), ctx.getEnd()));
    		this.filterTableRow.getItems().setAll(
    				ctx.getFilters()
    				);
    		if(ctx.getNSamples()>0)
	    		{
	    		this.genotypeChartPane.setCenter(this.makeGenotypePie(ctx));
	    		}
    		}
    	else
    		{
    		super.seqDictionaryCanvas.setSelectInterval(null);
        	this.filterTableRow.getItems().clear();
        	if(getVcfFile().getHeader().getNGenotypeSamples()>0)
	        	{
	    		this.genotypeChartPane.setCenter(this.makeGenotypePie(ctx));
	    		}
    		}
    	}
	
    private PieChart makeGenotypePie(final VariantContext ctx) {
        final Counter<GenotypeType> countTypes= new Counter<>();
        if(ctx!=null) {
        	for(final Genotype g:ctx.getGenotypes())
				{
        		// ignore genotype if not displayed
        		final SampleDef sampleDef= this.name2sampledef.get(g.getSampleName());
        		if(sampleDef==null || !sampleDef.isDisplayed()) continue;

        		countTypes.incr(g.getType());
				}
        	}
    	final ObservableList<PieChart.Data> pieChartData = FXCollections.observableArrayList();
    	for(final GenotypeType t:GenotypeType.values())
    		{
    		
    		int c= (int)countTypes.count(t);
    		if(c==0) continue;
    		pieChartData.add(new PieChart.Data(
    				t.name()+" = "+c,
    				c));
    		}
        final PieChart chart = new PieChart(pieChartData);
        chart.setLegendVisible(false);
        return chart;
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
        table.setPlaceholder(new Label("No FORMAT defined."));
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
        table.setPlaceholder(new Label("No INFO defined."));
        return tab;
		}
	
	/** build a table describing the INFO column */
	private Tab buildFilterHeaderTab(final VCFHeader header)
		{
        final TableView<VCFFilterHeaderLine> table=new TableView<>(FXCollections.observableArrayList(header.getFilterLines()));
        table.getColumns().add(makeColumn("ID", F->F.getID()));
        final Tab tab=new Tab("FILTER",table);
        tab.setClosable(false);
        
        table.setPlaceholder(new Label("No FILTER defined."));
        return tab;
		}
	
	@Override
    void reloadData()
		{
    	final int max_items= super.maxItemsLimitSpinner.getValue();
    	final List<VariantContext> L= new ArrayList<>(max_items);
    	final String location = this.gotoField.getText().trim();
    	final CloseableIterator<VariantContext> iter;
    	
    	try {
	    	if(location.isEmpty())
	    		{
	    		iter = this.getVcfFile().iterator();
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
					iter= this.getVcfFile().iterator(
							interval.getContig(),
							interval.getStart(),
							interval.getEnd()
							);
	    			}
	    		}
    		}
    	catch(final IOException err)
    		{
    		JfxNgs.showExceptionDialog(this, err);
    		return;
    		}
    	VcfJavascripFilter javascripFilter=null;
    	if(this.owner.javascriptCompiler.isPresent() &&
    		!this.javascriptArea.getText().trim().isEmpty())
    		{
    		try
    			{
    			javascripFilter = new VcfJavascripFilter(
    					this.javascriptFILTERfield.getText().trim(),
    					this.getVcfFile().getHeader(),
    					Optional.of(this.owner.javascriptCompiler.get().compile(this.javascriptArea.getText()))
    					);
    			}
    		catch(final Exception err)
    			{
    			updateStatusBar(AlertType.ERROR, err);
    			LOG.warning(err.getMessage());
    			}
    		}
    	int count_items=0;
    	while(iter!=null && iter.hasNext() && count_items<max_items)
    		{
    		VariantContext rec = iter.next();
    		++count_items;
    		if(javascripFilter!=null)
    			{
    			if(javascripFilter.encounteredException.isPresent()) break;
    			rec= javascripFilter.eval(rec);
    			if(rec==null) continue;
    			}
    		L.add(rec);
    		}
    	if(iter!=null) iter.close();
    	this.variantTable.getItems().setAll(L);
    	if(javascripFilter!=null && javascripFilter.encounteredException.isPresent())
    		{
    		JfxNgs.showExceptionDialog(this, javascripFilter.encounteredException.get());
    		}
    	
    	
    	if(!this.variantTable.getItems().isEmpty())
			{
			super.seqDictionaryCanvas.setItemsInterval(
					new ContigPos(this.variantTable.getItems().get(0).getContig(), this.variantTable.getItems().get(0).getStart()),
					new ContigPos(this.variantTable.getItems().get(this.variantTable.getItems().size()-1).getContig(), this.variantTable.getItems().get(this.variantTable.getItems().size()-1).getEnd())
					);
			}
		else
			{
			super.seqDictionaryCanvas.setItemsInterval(null,null);
			}
    	
    	}
	
	private void reloadAlleleTable(final VariantContext ctx)
		{
		if(ctx==null)
			{
			this.allelesTable.getItems().clear();
			}
		else
			{
			this.allelesTable.getItems().setAll(ctx.getAlleles());
			}
		}	
	
	private void reloadInfoTable(final VariantContext ctx)
		{
		if(ctx==null)
			{
			this.infoTableRow.getItems().clear();
			this.annPredictionTable.getItems().clear();
			this.vepPredictionTable.getItems().clear();
			}
		else
			{
			final String filterStr=this.tfFilterInfo.getText();
			Pattern regex=null;
			if(!filterStr.trim().isEmpty())
				{
				try
					{
					regex=Pattern.compile(filterStr.trim(),Pattern.CASE_INSENSITIVE);
					}
				catch(PatternSyntaxException err)
					{
					LOG.warning("Invalid regexp :"+filterStr);
					regex=Pattern.compile(Pattern.quote(filterStr));
					}
				}
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
					final InfoTableRow itr=new InfoTableRow(key,(L.size()==1?null:x+1),L.get(x));
					if(regex!=null && !itr.accept(regex)) {
						continue;
						}
					infos.add(itr);
					
					
					}
				}
			this.infoTableRow.getItems().setAll(infos);
			this.annPredictionTable.getItems().setAll(this.annPredictionParser.getPredictions(ctx));
			this.vepPredictionTable.getItems().setAll(this.vepPredictionParser.getPredictions(ctx));
			}
		
		}

	private void refreshGenotypeTable(final VariantContext ctx)
		{
		if(ctx==null || ctx.getNSamples()==0)
			{
        	this.genotypeTable.getItems().clear();
			}
		else
			{
        	this.genotypeTable.getItems().setAll(
        			ctx.getGenotypes().stream().
        			filter(G->cboxShowNoCall.isSelected() || G.isCalled()).
        			filter(G->cboxShowHomRef.isSelected() || !G.isHomRef()).
        			filter(G->cboxShowFiltered.isSelected() || !G.isFiltered()).
        			filter(G->{SampleDef d=name2sampledef.get(G.getSampleName()); return d==null ||d.isDisplayed();}).
        			collect(Collectors.toList())
        			);
			}
		}

	@Override
	protected void doMenuShowWholeStats(final ChartFactory<?> factory) {
			Optional<CompiledScript> compiledScript=Optional.empty();
        	if( VcfStage.this.owner.javascriptCompiler.isPresent() &&
        		!VcfStage.this.javascriptArea.getText().trim().isEmpty())
        		{
        		try
        			{
        			compiledScript = Optional.of(VcfStage.this.owner.javascriptCompiler.get().compile(VcfStage.this.javascriptArea.getText()));
        			}
        		catch(final Exception err)
        			{
        			LOG.warning(err.getMessage());
        			updateStatusBar(AlertType.ERROR, err);
        			return;
        			}
        		}
        	
        	VcfFile copy=null;
        	try {
        		copy=(VcfFile)this.getVcfFile().reOpen();
	        	final VcfQualityStage qcstage=new VcfQualityStage(
	        			(VariantContextChartFactory)factory,
	        			copy,
	        			compiledScript,
	        			V->true
	        			);
				qcstage.show();
        	} catch(final IOException err)
        		{
        		CloserUtil.close(copy);
        		JfxNgs.showExceptionDialog(this, err);
        		}
			}

    @Override
    protected final String getSnippetResourcePath() {
    	return "/com/github/lindenb/jvarkit/tools/vcfviewgui/vcf.snippets.xml";
    	}
	
    private VcfFile getVcfFile() {
    	return VcfFile.class.cast(super.getNgsFile());
    }
    
    @Override
    protected Optional<VariantContext> getCurrentSelectedItem()
    	{
    	return Optional.ofNullable(this.variantTable.getSelectionModel().getSelectedItem());
    	}

    
	}
