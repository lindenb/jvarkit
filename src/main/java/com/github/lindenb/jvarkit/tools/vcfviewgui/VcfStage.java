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
import com.github.lindenb.jvarkit.util.vcf.predictions.PredictionParserFactory;

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
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import javafx.beans.property.ReadOnlyObjectWrapper;
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
import javafx.scene.control.Label;
import javafx.scene.control.MenuItem;
import javafx.scene.control.ScrollPane;
import javafx.scene.control.Separator;
import javafx.scene.control.SeparatorMenuItem;
import javafx.scene.control.SpinnerValueFactory;
import javafx.scene.control.SplitPane;
import javafx.scene.control.Tab;
import javafx.scene.control.TabPane;
import javafx.scene.control.TableColumn;
import javafx.scene.control.TableView;
import javafx.scene.control.TextField;
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
    private static final int DEFAULT_VCF_RECORDS_COUNT=Integer.parseInt(System.getProperty("jxf.ngs.default.vcf", "100"));
    static final ExtensionFilter EXTENSION_FILTER=new ExtensionFilter("Variant Files", ".vcf",".vcf.gz");
    private static final String SPINNER_VALUE_KEY="vcf.spinner.value";
    
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
    
        
    
    /** Misc VCF tools that will be injected in the nashorn context */
    public static class VcfTools
    	{
    	private final AnnPredictionParser annParser;
    	VcfTools(final VCFHeader header) {
    		this.annParser  = new PredictionParserFactory().header(header).buildAnnPredictionParser();
    	}
    	
    	public List<AnnPredictionParser.AnnPrediction> getAnnPredictions(final VariantContext ctx) {
    		return this.annParser.getPredictions(ctx);
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
    
    private static class VcfJavascripFilter
	extends JavascriptFilter<VCFHeader,VariantContext>
		{
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
			boolean ok=super.accept();
			if(!ok )
				{
				if(filter==null || filter.trim().isEmpty()) return null;
				return new VariantContextBuilder(ctx).filter(this.filter).make();
				}
			return ctx;
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
	/** annotation SNPEFF table ROW */
	private final TableView<AnnPredictionParser.AnnPrediction> annTableRow;
	private final TableView<String> filterTableRow;
	private final BorderPane genotypeChartPane;
	private final CheckBox cboxShowHomRef=new CheckBox("HomRef");
	private final CheckBox cboxShowNoCall=new CheckBox("NoCall");
	private final TextField tfFilterInfo=new TextField();
	private final AnnPredictionParser annPredictionParser;
	VcfStage(final JfxNgs owner,final VcfFile vcfFile) throws IOException {
		super(owner,vcfFile);
		final VCFHeader header= vcfFile.getHeader();
		
        this.annPredictionParser=new PredictionParserFactory().
        			header(header).
        			buildAnnPredictionParser();
      
        
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
        
        super.maxItemsLimitSpinner.setValueFactory(
        		new SpinnerValueFactory.IntegerSpinnerValueFactory(0,100000,
        				owner.preferences.getInt(SPINNER_VALUE_KEY, DEFAULT_VCF_RECORDS_COUNT)
        				));;

        top1.getChildren().add(super.maxItemsLimitSpinner);
        top1.getChildren().add(new Separator(Orientation.VERTICAL));
		
        top1.getChildren().add(createIgvButton());
        
        this.gotoField.setOnAction(A->reloadData());
		final TabPane tabPane=new TabPane();
		tabPane.setPadding(new Insets(10, 10, 10, 10));
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
		FlowPane flow3 = new FlowPane(this.cboxShowHomRef,this.cboxShowNoCall);
		this.cboxShowHomRef.setSelected(true);
		this.cboxShowNoCall.setSelected(true);
		EventHandler<ActionEvent> repaintGTTable=new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent event) {
				refreshGenotypeTable(variantTable.getSelectionModel().getSelectedItem());
				}
			};
		this.cboxShowHomRef.setOnAction(repaintGTTable);
		this.cboxShowNoCall.setOnAction(repaintGTTable);

		
		if(header.getNGenotypeSamples()>0)
			{
			final BorderPane pane2=new BorderPane(this.genotypeTable);
			pane2.setTop(flow3);
			//GridPane.setConstraints(pane2,5, 0,5,14); 
			split1.getItems().add(pane2);
			}
		
		final SplitPane split2=new SplitPane();
		split2.setOrientation(Orientation.HORIZONTAL);

		/* filter table */
		this.filterTableRow = this.buildFilterTable();
		//GridPane.setConstraints( this.filterTableRow,0, 10,3,1);
		//gridPane.getChildren().add(this.filterTableRow);
		split2.getItems().add(this.filterTableRow);
		/* genotype pane */
		
		
		this.genotypeChartPane = new BorderPane(this.makeGenotypePie(null));
		if(header.getNGenotypeSamples()>0)
			{
			this.genotypeChartPane.setPadding(new Insets(5));
			//GridPane.setConstraints( this.genotypeChartPane,0, 12,3,4);
	    	//gridPane.getChildren().add(this.genotypeChartPane);
			split2.getItems().add(this.genotypeChartPane);
			}
		
		/* build info Table table */
		this.infoTableRow = this.buildInfoTableRow();
		this.annTableRow = this.buildAnnTableRow();
		flow3 = new FlowPane(new Label("Filter:"),tfFilterInfo);
		tfFilterInfo.setPrefColumnCount(10);
		tfFilterInfo.setOnAction(AE->reloadInfoTable(variantTable.getSelectionModel().getSelectedItem()));
		final TabPane tabPane2=new TabPane();
		Tab tab=new Tab("INFO",this.infoTableRow);
		tab.setClosable(false);
		tabPane2.getTabs().add(tab);
		tab=new Tab("ANN",this.annTableRow);
		tab.setClosable(false);
		tabPane2.getTabs().add(tab);
		
		BorderPane pane2=new BorderPane(tabPane2);
		pane2.setTop(flow3);
		//GridPane.setConstraints( this.infoTableRow,3, 10,8,5); // column=3 row=1
		//gridPane.getChildren().add(this.infoTableRow);
		
		
		split2.getItems().add(pane2);
				
       
		
        
        //vbox1.getChildren().add(gridPane);
		final SplitPane split3=new SplitPane();
		split3.setOrientation(Orientation.VERTICAL);
		split3.getItems().addAll(split1,split2);

		tab=new Tab("Variants", split3);
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
        		menuForSavingTable("Genotype",this.genotypeTable)
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
        table.getColumns().add(formatIntegerColumn(makeColumn("POS", V->V.getStart())));
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
		final TableView<InfoTableRow> table=new TableView<>();
		table.getColumns().add(makeColumn("Key", R->R.key));
		table.getColumns().add(makeColumn("Index", R->R.index));
		
		final TableColumn<InfoTableRow, Object> t=makeColumn("Value", R->R.value);
		t.setPrefWidth(300.0);
		table.getColumns().add(t);
		
		return table;
		}
     
	private TableView<AnnPredictionParser.AnnPrediction> buildAnnTableRow()
		{
		final TableView<AnnPredictionParser.AnnPrediction> table=new TableView<>();
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
    	refreshGenotypeTable(ctx);
    	reloadInfoTable(ctx);
    	if(ctx!=null)
    		{
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
    			javascripFilter=null;
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
    	}
	
	private void reloadInfoTable(final VariantContext ctx)
		{
		if(ctx==null)
			{
			this.infoTableRow.getItems().clear();
			this.annTableRow.getItems().clear();
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
			this.annTableRow.getItems().setAll(this.annPredictionParser.getPredictions(ctx));
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
    
	}
