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

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.function.Supplier;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;
import java.util.stream.Collectors;

import javax.script.CompiledScript;

import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.AFByPopulationChartFactory;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.AFBySexChartFactory;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.AlleleFrequencyChartFactory;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.ChartFactory;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.GenotypeTypeChartFactory;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.TiTvChartFactory;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.VariantContextChartFactory;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.VariantDepthChartFactory;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.VariantQualChartFactory;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.VariantTypeChartFactory;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree.Term;
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
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.chart.PieChart;
import javafx.scene.control.Alert;
import javafx.scene.control.Button;
import javafx.scene.control.ButtonType;
import javafx.scene.control.CheckBox;
import javafx.scene.control.ContextMenu;
import javafx.scene.control.Label;
import javafx.scene.control.Menu;
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
import javafx.scene.paint.Color;
import javafx.scene.paint.Paint;
import javafx.scene.shape.ArcType;
import javafx.scene.text.Font;
import javafx.scene.text.Text;
import javafx.scene.text.TextFlow;
import javafx.stage.FileChooser;
import javafx.stage.Screen;
import javafx.stage.WindowEvent;
import javafx.stage.FileChooser.ExtensionFilter;
import javafx.util.Callback;

public class VcfStage extends NgsStage<VCFHeader,VariantContext> {
    static final String SPINNER_VALUE_KEY="vcf.spinner.value";
    static final String VARIANT_CONTEXT_KEY="variant";
	static final int DEFAULT_VCF_RECORDS_COUNT= 100;
    static final List<ExtensionFilter> EXTENSION_FILTERS= Arrays.asList(
    		new ExtensionFilter("Indexed VCF File", "*.vcf", "*.vcf.gz"),
    		new ExtensionFilter("Tabix-Indexed VCF Files", "*.vcf.gz"),
    		new ExtensionFilter("Tribble-Indexed VCF File", "*.vcf")
    		);
    
    /** variant oriented chart-factories */
    private static final List<Supplier<ChartFactory<VCFHeader,VariantContext>>> VARIANT_CHART_FACTORIES= Arrays.asList(
		()->new VariantTypeChartFactory(),
		()->new GenotypeTypeChartFactory(),
		()->new AlleleFrequencyChartFactory(),
		()->new VariantQualChartFactory(),
		()->new TiTvChartFactory(),
		()->new VariantDepthChartFactory(),
		()->new AFByPopulationChartFactory(),
		()->new AFBySexChartFactory()
		);
    
    /** class used to display samples */
    public static class SampleDef
    	{
    	private final String name;
    	private final String description;
        private final BooleanProperty _displayedProperty = new SimpleBooleanProperty(true);
        private final Optional<PedFile.Sample> pedsample;
        public BooleanProperty displayedProperty()  { return _displayedProperty; }
        public boolean isDisplayed() { return displayedProperty().get(); }
        public void setDisplayed(boolean v) { displayedProperty().set(v); }
        public Optional<PedFile.Sample> getPedSample() { return this.pedsample;}
    	SampleDef(final String name,final String description,final PedFile.Sample pedsample)
	    	{
	    	this.name=name;
	    	this.description=description;
	    	this.pedsample = Optional.ofNullable(pedsample);
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
    	private final SequenceOntologyTree sequenceOntologyTree= SequenceOntologyTree.getInstance(); 
    	private final AnnPredictionParser annParser;
    	private final VepPredictionParser vepParser;
    	@SuppressWarnings("unused")
		private final PedFile pedfile;
    	VcfTools(final VCFHeader header,final PedFile pedfile) {
    		this.pedfile=pedfile;
    		this.annParser  = new AnnPredictionParserFactory(header).get();
    		this.vepParser  = new VepPredictionParserFactory(header).get();
    	}
    	
    	public SequenceOntologyTree getSequenceOntologyTree() {
			return sequenceOntologyTree;
		}
    	
    	public List<AnnPredictionParser.AnnPrediction> getAnnPredictions(final VariantContext ctx) {
    		return this.annParser.getPredictions(ctx);
    		}
    	
    	public List<VepPredictionParser.VepPrediction> getVepPredictions(final VariantContext ctx) {
    		return this.vepParser.getPredictions(ctx);
    		}

    	/** return true if variant has any prediction with a SO term (or its children) with this label */
    	public boolean hasSequenceOntologyLabel(final VariantContext ctx,final String lbl)
    		{
    		if(lbl==null) return false;
    		final Term t= this.getSequenceOntologyTree().getTermByLabel(lbl);
    		if(t==null) LOG.warning("don't know SO.label "+lbl);
    		return hasSequenceOntologyTerm(ctx,t);
    		}
    	/** return true if variant has any prediction with a SO term (or its children) with this accession */
    	public boolean hasSequenceOntologyAccession(final VariantContext ctx,final String acn)
			{
			if(acn==null) return false;
			final Term t= this.getSequenceOntologyTree().getTermByAcn(acn);
    		if(t==null) LOG.warning("don't know SO.acn "+acn);
			return hasSequenceOntologyTerm(ctx,t);
			}

    	/** return true if variant has any prediction with a SO term (or its children) */
    	public boolean hasSequenceOntologyTerm(final VariantContext ctx,final Term t)
			{
			if(t==null) return false;
			final Set<Term> children=t.getAllDescendants();
			for(AnnPredictionParser.AnnPrediction a: getAnnPredictions(ctx)) {
				if(!Collections.disjoint(a.getSOTerms(),children)) return true;
				}
			for(VepPredictionParser.VepPrediction a: getVepPredictions(ctx)) {
				if(!Collections.disjoint(a.getSOTerms(),children)) return true;
				}
			return false;
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
					if(child.sameGenotype(sim, true)) return false;
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
				final VcfFile vcfFile,
				final Optional<CompiledScript> compiledScript
				)
			{
			super(vcfFile.getHeader(),compiledScript);
			super.bindings.put(TOOL_CONTEXT_KEY, new VcfTools(vcfFile.getHeader(),vcfFile.getPedigree()));
			super.bindings.put(PEDIGREE_CONTEXT_KEY, vcfFile.getPedigree());
			
			this.filter=(filter==null?"":filter);
			}
		@Override public VariantContext eval(final VariantContext ctx)
			{
			if(!super.compiledScript.isPresent()) return ctx;
			super.bindings.put(VARIANT_CONTEXT_KEY, ctx);
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
    			ScanVcfThread(final ChartFactory<VCFHeader,VariantContext> factory,
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
    									new VcfJavascripFilter("",
    											getVcfFile(),
    											compiledScript
    											));
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
    		
    		VcfQualityStage(final ChartFactory<VCFHeader,VariantContext> factory,
    				final VcfFile vcfinput,
    				final Optional<CompiledScript> compiledScript,
    				final Predicate<VariantContext> otherFilters
    				)
	    		{
    			super(factory,vcfinput,compiledScript,otherFilters);
	    		}

			@Override
			protected ScanVcfThread createThread(
					final ChartFactory<VCFHeader,VariantContext> factory,
					final NgsFile<VCFHeader,VariantContext> ngsfile,
					final Optional<CompiledScript> compiledScript,
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
	private final TableView<PedFile.TrioGenotype> triosTable;
	private final BorderPane genotypeChartPane;
	private final CheckBox canvasEvenlySpaced=new CheckBox("Evenly Spaced");
	private final CheckBox cboxShowHomRef=new CheckBox("HomRef");
	private final CheckBox cboxShowNoCall=new CheckBox("NoCall");
	private final CheckBox cboxShowFiltered=new CheckBox("Filtered");
	private final TextField tfFilterInfo=new TextField();
	private final AnnPredictionParser annPredictionParser;
	private final VepPredictionParser vepPredictionParser;
	/* don't display allele if it's too big */
	private final Function<Allele,String> allele2stringConverter;
	
	/** canvas displaying variants */
	private final ResizableCanvas drawingArea;
	
    /** bioalcidae instance */
    private final AbstractAwkLike bioalcidae=new  AbstractAwkLike()
    		{
    		@Override
			protected TextFlow getHelpString() {
    			final TextFlow tf= super.getHelpString();
    			
				 tf.getChildren().addAll(new Text("\n"+
						"* '"+ HEADER_CONTEXT_KEY+"' an instance of java class  "),javadocFor(VCFHeader.class),new Text("\n"+
						"* '"+ ITER_CONTEXT_KEY+"' an Iterator over instances of java class "),javadocFor(VariantContext.class),new Text("\n"+
						"* '"+ TOOL_CONTEXT_KEY +"' an instance of "),javadocFor(VcfTools.class),new Text("\n"+
						"* '"+ PEDIGREE_CONTEXT_KEY +"' an instance of "),javadocFor(PedFile.class),new Text("\n"
						));
				 return tf;
				}
    	
    		@Override
    		protected javax.script.SimpleBindings completeBindings(javax.script.SimpleBindings sb, final VCFHeader h)
    			{
    			sb.put(TOOL_CONTEXT_KEY, new VcfTools(h, getPedigree()));
    			sb.put(PEDIGREE_CONTEXT_KEY, getPedigree());
    			return sb;
    			}
    		};

	
	VcfStage(final JfxNgs owner,final VcfFile vcfFile) throws IOException {
		super(owner,vcfFile);
		final VCFHeader header= vcfFile.getHeader();
		
        this.annPredictionParser=new AnnPredictionParserFactory(header).get();
        this.vepPredictionParser=new VepPredictionParserFactory(header).get();
        
        /* create allele2stringConverter , it won't display read length having more than xxx bases*/
        {
        	int prefnum=0;
        	try {
        		prefnum=Integer.parseInt(this.owner.preferences.get(this.owner.pref_vcf_max_allele_length_displayed.key,"100"));
        	} catch(final NumberFormatException err) {
        		prefnum=100;
        	}
        	final int vcf_max_allele_length = prefnum;
        	this.allele2stringConverter = new Function<Allele, String>() {
				@Override
				public String apply(final Allele t) {
					if(t==null) return null;
					if(t.isNoCall()) return Allele.NO_CALL_STRING; 
					if(t.isSymbolic()) return t.getDisplayString();
					if(vcf_max_allele_length>=0 && t.length() > vcf_max_allele_length) {
						return t.getDisplayString().substring(0, vcf_max_allele_length)+"... (len="+t.length()+")";
					}
					return t.getDisplayString();
				}
			};
        }
        
        
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
				final SampleDef sampledef = new SampleDef(
						sampleName,
						desc,
						vcfFile.getPedigree().get(sampleName)
						);
				sampledef.displayedProperty().addListener(CL->{paintDrawingArea();});
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
		
		this.triosTable = buildTrioTable();
		if(vcfFile.getHeader().getNGenotypeSamples()>0 && !vcfFile.getPedigree().isEmpty()) {
			tabPane2.getTabs().add(new Tab("TRIOS",this.triosTable));
			}
		
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
        
        /* build canvas pane */
        this.drawingArea = new ResizableCanvas(900, 400)
        		{
        		public void repaintCanvas() {paintDrawingArea();};
        		};
        BorderPane canvasPane=new BorderPane(this.drawingArea);
        this.canvasEvenlySpaced.setSelected(false);
        this.canvasEvenlySpaced.setOnAction(AE->{paintDrawingArea();});
        canvasPane.setTop(new FlowPane(this.canvasEvenlySpaced));
        tab=new Tab("Canvas",canvasPane);
        tab.setClosable(false);
		tabPane.getTabs().add(tab);
        
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
        
        for(final Supplier<ChartFactory<VCFHeader,VariantContext>> supplier: VARIANT_CHART_FACTORIES)
	        {
        	final ChartFactory<VCFHeader,VariantContext> factory = supplier.get();
        	final MenuItem menuItem=new MenuItem("Local "+factory.getName());
        	menuItem.setOnAction(new EventHandler<ActionEvent>() {
				@Override
				public void handle(ActionEvent event) {
					doMenuShowLocalStats(supplier.get(), variantsProvider);
				}
			});
        	statsMenu.getItems().add(menuItem);
	        }
        super.statsMenu.getItems().add(new SeparatorMenuItem());
        for(final Supplier<ChartFactory<VCFHeader,VariantContext>> supplier: VARIANT_CHART_FACTORIES)
	        {
	    	final ChartFactory<VCFHeader,VariantContext> factory = supplier.get();
	    	
	    	final MenuItem menuItem=new MenuItem("Whole"+factory.getName());
	    	menuItem.setOnAction(new EventHandler<ActionEvent>() {
				@Override
				public void handle(ActionEvent event) {
					doMenuShowWholeStats(supplier.get());
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
    	fc.getExtensionFilters().addAll(EXTENSION_FILTERS);
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
    				this.getVcfFile(),
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
			vcwb.setOutputFile(saveAs);
			vcwb.setOutputFileType(OutputType.BLOCK_COMPRESSED_VCF);
			vcwb.setReferenceDictionary(this.getNgsFile().getSequenceDictionary());
			w= vcwb.build();
			final VCFHeader header2= new VCFHeader(this.getVcfFile().getHeader());
        	if(javascriptFilter!=null && !(javascriptFilter.filter==null || javascriptFilter.filter.isEmpty()))
        		{
        		header2.addMetaDataLine(new VCFFilterHeaderLine(javascriptFilter.filter, "Set by User in JfxNgs:"+
        				this.javascriptArea.getText().replaceAll("[\n\t\r ]+"," ")
        				));
        		}
			w.writeHeader(header2);
			iter= new LogCloseableIterator(this.getVcfFile().iterator(),null);
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
			w=null;
			iter.close();iter=null;
    		final Alert alert = new Alert(AlertType.CONFIRMATION, "Done", ButtonType.OK);
			alert.showAndWait();
			}
		catch(final Exception err)
			{
			err.printStackTrace();
			JfxNgs.showExceptionDialog(this, err);
			return;
			}
		finally
			{
			CloserUtil.close(iter);
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
		final FlowPane bottom=new FlowPane(new TextFlow(new Text("The script injects:\n"+
				"* '"+HEADER_CONTEXT_KEY+"' an instance of  "),javadocFor(VCFHeader.class),new Text("\n"+
				"* '"+VARIANT_CONTEXT_KEY+"' an instance of "),javadocFor(VariantContext.class),new Text(" )\n"+
				"* '"+PEDIGREE_CONTEXT_KEY+"' an instance of "),javadocFor(PedFile.class),new Text("\n"+
				"* '"+TOOL_CONTEXT_KEY+"' an instance of "),javadocFor(VcfTools.class),new Text("\n"+
				"The script should return a boolean: true (accept variant) or false (reject variant)"
				)));
		pane.setBottom(bottom);
		
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
		updateStatusBar(AlertType.NONE,"");
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
        table.getColumns().add(makeColumn("REF", V->allele2stringConverter.apply(V.getReference())));
        table.getColumns().add(makeColumn("ALT", V->V.getAlternateAlleles().stream().map(A->allele2stringConverter.apply(A)).collect(Collectors.joining(","))));
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
		table.getColumns().add(makeColumn("Bases.",A->allele2stringConverter.apply(A)));
		table.getColumns().add(makeColumn("Length.",A->{if(A.isSymbolic()) return (Integer)null;return A.length();}));
		table.setPlaceholder(new Label("No Allele."));
		return table;
		}

	/** build Trio table */
	private TableView<PedFile.TrioGenotype> buildTrioTable()
		{

		final TableView<PedFile.TrioGenotype> table=new TableView<>();
		if(getVcfFile().getHeader().getNGenotypeSamples()>0 && !getPedigree().isEmpty()) {
			final Function<Genotype,String> gt2str=new Function<Genotype,String>()
				{
				@Override
				public String apply(final Genotype gt)
					{
					if(gt==null || !gt.isCalled()) return null;
					return gt.getAlleles().stream().map(S->allele2stringConverter.apply(S)).collect(Collectors.joining(gt.isPhased()?"|":"/"));
					}
				};
			
			table.getColumns().add(makeColumn("Child",T->T.getChildren()==null?null:T.getChildren().getSampleName()));
			table.getColumns().add(makeColumn("Child-GT",T->gt2str.apply(T.getChildren())));
			table.getColumns().add(makeColumn("Father",T->T.getFather()==null?null:T.getFather().getSampleName()));
			table.getColumns().add(makeColumn("Father-GT",T->gt2str.apply(T.getFather())));
			table.getColumns().add(makeColumn("Mother",T->T.getMother()==null?null:T.getMother().getSampleName()));
			table.getColumns().add(makeColumn("Mother-GT",T->gt2str.apply(T.getMother())));
			table.getColumns().add(makeColumn("Violation",T->T.isMendelianIncompatibility()));
			}
		table.setPlaceholder(new Label("No Trio."));
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
		
		if(!getVcfFile().getPedigree().isEmpty())
			{
			this.sampleTable.getColumns().add(makeColumn("Family", R->{
				if(R.getPedSample().isPresent())
					{
					return R.getPedSample().get().getFamily();
					}
				return null;
				}));
			}
		
		
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
		
		
		if(!getVcfFile().getPedigree().isEmpty())
			{
			this.sampleTable.getColumns().add(makeColumn("Father", R->{
				if(R.getPedSample().isPresent())
					{
					final PedFile.Sample parent= R.getPedSample().get().getFather();
					return parent==null?null:parent.getName();
					}
				return (String)null;
				}));
			
			this.sampleTable.getColumns().add(makeColumn("Mother", R->{
				if(R.getPedSample().isPresent())
					{
					final PedFile.Sample parent= R.getPedSample().get().getMother();
					return parent==null?null:parent.getName();
					}
				return (String)null;
				}));

			this.sampleTable.getColumns().add(makeColumn("Sex", R->{
				if(R.getPedSample().isPresent())
					{
					return R.getPedSample().get().getSex().name();
					}
				return (String)null;
				}));
			
			this.sampleTable.getColumns().add(makeColumn("Status", R->{
				if(R.getPedSample().isPresent())
					{
					return R.getPedSample().get().getStatus().name();
					}
				return (String)null;
				}));
			}
		
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
			table.getColumns().add(makeColumn("GeneName", P->P.getGeneName()));
			table.getColumns().add(makeColumn("GeneId", P->P.getGeneId()));
			table.getColumns().add(makeColumn("Feature", P->P.getFeatureType()));
			table.getColumns().add(makeColumn("FeatureId", P->P.getFeatureId()));
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
					final String delim;
					Object o;
					if(param.getTableColumn().getText().equals(VCFConstants.GENOTYPE_KEY))
						{
						delim =  param.getValue().isPhased()?"|":"/";
						o = param.getValue().
								getAlleles().
								stream().map(A->allele2stringConverter.apply(A)).
								collect(Collectors.toList());
						}
					else
						{
						delim=",";
						o = param.getValue().getAnyAttribute(param.getTableColumn().getText());
						}
					
					if(o==null)
						{
						return new ReadOnlyObjectWrapper<String>(null);
						}
					if(o instanceof List)
						{
						List<?> L=(List<?>)o;
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
    	refreshTrioTable(ctx);
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
        table.getColumns().add(makeColumn("Description", F->F.getDescription()));
        final Tab tab=new Tab("FILTER",table);
        tab.setClosable(false);
        
        table.setPlaceholder(new Label("No FILTER defined."));
        return tab;
		}
	
	@Override
    void reloadData()
		{
		updateStatusBar(AlertType.NONE,"");
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
    					this.getVcfFile(),
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
    		int last_index=this.variantTable.getItems().size()-1;
			super.seqDictionaryCanvas.setItemsInterval(
					new ContigPos(this.variantTable.getItems().get(0).getContig(), this.variantTable.getItems().get(0).getStart()),
					new ContigPos(this.variantTable.getItems().get(last_index).getContig(), this.variantTable.getItems().get(last_index).getEnd())
					);
			if(this.variantTable.getItems().get(0).getContig().equals(
				this.variantTable.getItems().get(last_index).getContig()))
				{
				this.gotoField.setText(this.variantTable.getItems().get(0).getContig()+":"+this.variantTable.getItems().get(0).getStart()+"-"+this.variantTable.getItems().get(last_index).getEnd());
				}
			
			}
		else
			{
			super.seqDictionaryCanvas.setItemsInterval(null,null);
			}
    	paintDrawingArea();
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
			for(final String key: new TreeSet<String>(atts.keySet()))
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

	
	private void refreshTrioTable(final VariantContext ctx)
		{
		this.triosTable.getItems().clear();
		for(PedFile.TrioGenotype gt: getPedigree().getTrios(ctx))
			{
			this.triosTable.getItems().add(gt);
			}
		}
	
	
	
	@Override
	protected void doMenuShowWholeStats(
			final ChartFactory<VCFHeader,VariantContext> factory) {
			updateStatusBar(AlertType.NONE,"");
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
        		factory.setHeader(copy.getHeader());
        		factory.setPedigree(getPedigree());
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
    public PedFile getPedigree() {
    	return getVcfFile().getPedigree();
    	}
    
    @Override
    protected Optional<VariantContext> getCurrentSelectedItem()
    	{
    	return Optional.ofNullable(this.variantTable.getSelectionModel().getSelectedItem());
    	}

    @Override
    protected Menu createJavascriptSnippetMenu() {
    	final Menu menu = super.createJavascriptSnippetMenu();
    	menu.getItems().add(new SeparatorMenuItem());
    	
    	MenuItem item= new MenuItem("Insert Samples' name as Array");
    	menu.getItems().add(item);
    	item.setOnAction(AE->{
    		final int caret = super.javascriptArea.getCaretPosition();
        	final StringBuilder sb=new StringBuilder("var sampleNames=[").
        			append(getVcfFile().getHeader().getSampleNamesInOrder().stream().map(S->"\""+S+"\"").collect(Collectors.joining(","))).
        			append("];\n")
        			;
        	super.javascriptArea.insertText(caret, sb.toString());
    		});
    	
    	return menu;
    	}
    /** invoke bioalcidae */
   @Override
   void invokeBioalcidae()
    	{
    	this.bioalcidae.show();
    	}
   
   
   private void paintDrawingArea()
	   	{
    	final double canvaswidth= this.drawingArea.getCanvas().getWidth();
    	if(canvaswidth<1) return;
    	final double canvasheight= this.drawingArea.getCanvas().getHeight();
    	if(canvasheight<1) return;
    	final GraphicsContext gc=this.drawingArea.getCanvas().getGraphicsContext2D();
    	gc.setGlobalAlpha(1.0);
    	gc.setFill(Color.WHITE);
    	gc.fillRect(0, 0, canvaswidth, canvasheight);
		if(this.variantTable.getItems().isEmpty()) return;
		
		gc.setFill(Color.BLACK);
		final int MAX_VARIANTS=(int)(canvaswidth/4.0);
		if(this.variantTable.getItems().size() >MAX_VARIANTS)
			{
			gc.setFont(Font.font ("Verdana", 24));
			gc.fillText("Sorry. Too Many Variants (N="+this.variantTable.getItems().size()+">"+MAX_VARIANTS+")",2,50);
			return;
			}
		
		int nrows=1;//genome track
		for(final SampleDef sampleDef : this.name2sampledef.values())
			{
			if(!sampleDef.isDisplayed()) continue;
			nrows++;
			}
    	double rowheight = canvasheight/nrows;
    	if(rowheight<5) rowheight=5;
		
		final double marginleft = (nrows<=1?0.0:canvaswidth/10.0);//no margin if no sample
		
		final Function<Long,Double> convertGenomicIndexToPixel ;
		final Function<Integer,Double> convertListIndexToPixel ;
		
		if( !this.canvasEvenlySpaced.isSelected() )
			{
			Long minGenomicIndex;
			Long maxGenomicIndex;
		
				{
				final VariantContext ctx = this.variantTable.getItems().get(0);
				minGenomicIndex = convertContigPosToGenomicIndex(ctx.getContig(), ctx.getStart());
				maxGenomicIndex = minGenomicIndex;
				}
		
				{
				final VariantContext ctx = this.variantTable.getItems().get(this.variantTable.getItems().size()-1);
				maxGenomicIndex = convertContigPosToGenomicIndex(ctx.getContig(), ctx.getEnd());
				}
			if(minGenomicIndex==null || maxGenomicIndex==null) return;
		
		
			long extend = (long)((maxGenomicIndex-minGenomicIndex)*0.1);
			minGenomicIndex = minGenomicIndex-extend;
			maxGenomicIndex = maxGenomicIndex+extend;
		
			final long genomicIndexStart = minGenomicIndex;
			final long genomicIndexLength = (maxGenomicIndex-minGenomicIndex);
			
			convertGenomicIndexToPixel = (GI)->
				{
					double x=	marginleft+(canvaswidth-marginleft)*((double)(GI-genomicIndexStart))/((double)(genomicIndexLength));
				
					return x;
				};
			convertListIndexToPixel = null;
			}
		else
			{
			final double nItems = 2.0 + this.variantTable.getItems().size();
			convertGenomicIndexToPixel = null;
			convertListIndexToPixel = (IDX)-> {
					return marginleft + ((1+IDX)/nItems)*(canvaswidth-marginleft);
				};
			}
		
		//draw vertical lines
		String prev_chr="";
		Paint fill_chr1=Color.BLACK;
		Paint fill_chr2=Color.DARKGRAY;
		Paint fill_current = fill_chr1;
		gc.setFont(Font.font ("Verdana", 7));
		for(int idx= 0;idx < this.variantTable.getItems().size(); ++idx )
			{
			final VariantContext ctx = this.variantTable.getItems().get(idx);
			
			if(!ctx.getContig().equals(prev_chr)) {
				prev_chr=ctx.getContig();
				fill_current =(fill_current==fill_chr1?fill_chr2:fill_chr1); 
				}
			gc.setFill(fill_current);
			
			double x0;
			double x1;
			
			if( convertListIndexToPixel==null)
				{
				x0 = convertGenomicIndexToPixel.apply(convertContigPosToGenomicIndex(ctx.getContig(), ctx.getStart()));
			    x1 = convertGenomicIndexToPixel.apply(convertContigPosToGenomicIndex(ctx.getContig(), ctx.getEnd()));
			    if(x0<=x1) x1=x0+1;
				}
			else
				{
				x0 = convertListIndexToPixel.apply(idx);
				x1=x0+1;
				}
			
			gc.setGlobalAlpha(0.4);
			gc.fillRect(x0, 0, (x1-x0), canvasheight);
			
			// http://stackoverflow.com/questions/39524792/
			gc.setGlobalAlpha(1.0);
			gc.setFill(Color.BLACK);
			gc.save();
		    gc.translate(x0, 0);
		    gc.rotate(90);
		    gc.fillText(ctx.getContig()+":"+ctx.getStart(), 5, 0,rowheight);
		    gc.restore();
			}
    	final double radius=3.0;


		
		gc.setFont(Font.font ("Verdana", 12));
    	double y = rowheight;
    	for(final SampleDef sampleDef : this.name2sampledef.values())
			{
			if(!sampleDef.isDisplayed()) continue;
			gc.setLineWidth(1.0);
			gc.setGlobalAlpha(1.0);
    		gc.setFill(Color.BLACK);
    		gc.setStroke(Color.BLACK);
    		gc.strokeLine(0, y, canvaswidth, y);
    		gc.fillText(sampleDef.getName(),1,y+gc.getFont().getSize(),marginleft);
    		gc.setLineWidth(0.5);
    		final double midy= y+ rowheight/2.0;
    		
    		gc.setGlobalAlpha(0.8);
    		for(int idx=0;idx < this.variantTable.getItems().size();++idx)
    			{
    			final VariantContext ctx= this.variantTable.getItems().get(idx);
    			
    			double x0;
    			if(convertListIndexToPixel==null)
    				{
    				x0 = convertGenomicIndexToPixel.apply(convertContigPosToGenomicIndex(ctx.getContig(), ctx.getStart()));
    				}
    			else
    				{
    				x0 = convertListIndexToPixel.apply(idx);
    				}
    			final Genotype g= ctx.getGenotype(sampleDef.getName());
    			if(g==null || g.isNoCall()) continue;
    			if(g.isHomRef()) {
    				gc.setStroke(Color.DARKGRAY);
    				gc.strokeOval(x0-radius, midy-radius, radius*2, radius*2);
    				}
    			else if(g.isHet() && !g.isHetNonRef())
    				{
    				gc.setStroke(Color.DARKGRAY);
    				gc.strokeOval(x0-radius, midy-radius, radius*2, radius*2);
    				gc.setFill(Color.RED);
    				gc.fillArc(x0-radius, midy-radius, radius*2, radius*2,0,180,ArcType.CHORD);
    				}
    			else if(g.isHetNonRef())
					{
					gc.setFill(Color.PINK);
					gc.fillOval(x0-radius, midy-radius, radius*2, radius*2);
					gc.setFill(Color.RED);
					gc.fillArc(x0-radius, midy-radius, radius*2, radius*2,0,180,ArcType.CHORD);
					}
    			else
    				{
    				gc.setFill(Color.RED);
    				gc.fillOval(x0-radius, midy-radius, radius*2, radius*2);
    				}
    			}
    		
    		y += rowheight;
    		if(y> canvasheight) break;
			}
		
		
	   	}
   
	}
