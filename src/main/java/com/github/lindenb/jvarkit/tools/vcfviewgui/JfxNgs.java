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
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.net.URL;
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
import java.util.Optional;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.prefs.BackingStoreException;
import java.util.prefs.Preferences;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import javax.script.Bindings;
import javax.script.Compilable;
import javax.script.CompiledScript;
import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;
import javax.script.ScriptException;
import javax.script.SimpleBindings;
import javax.xml.namespace.QName;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;

import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.AlleleFrequencyChartFactory;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.BasesPerPositionChartFactory;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.GCPercentChartFactory;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.GenotypeTypeChartFactory;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.MapqChartFactory;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.QualityPerPositionChartFactory;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.ReadChartFactory;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.ReadLengthChartFactory;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.ReadQualityChartFactory;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.SamFlagsChartFactory;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.TiTvChartFactory;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.VariantContextChartFactory;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.VariantDepthChartFactory;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.VariantQualChartFactory;
import com.github.lindenb.jvarkit.tools.vcfviewgui.chart.VariantTypeChartFactory;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.hershey.JfxHershey;
import com.github.lindenb.jvarkit.util.igv.IgvSocket;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.BAMIndex;
import htsjdk.samtools.BAMIndexer;
import htsjdk.samtools.BamFileIoUtils;
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
import htsjdk.samtools.SamReaderFactory.Option;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.index.IndexFactory.IndexType;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndex;
import htsjdk.tribble.index.tabix.TabixIndexCreator;
import htsjdk.tribble.util.TabixUtils;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFCodec;
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
import javafx.collections.ObservableList;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.geometry.Insets;
import javafx.geometry.Orientation;
import javafx.geometry.Rectangle2D;
import javafx.scene.Scene;
import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.chart.Chart;
import javafx.scene.chart.PieChart;
import javafx.scene.control.Alert;
import javafx.scene.control.Alert.AlertType;
import javafx.scene.control.Button;
import javafx.scene.control.ButtonType;
import javafx.scene.control.CheckBox;
import javafx.scene.control.CheckMenuItem;
import javafx.scene.control.Dialog;
import javafx.scene.control.Label;
import javafx.scene.control.Menu;
import javafx.scene.control.MenuBar;
import javafx.scene.control.MenuItem;
import javafx.scene.control.ScrollBar;
import javafx.scene.control.ScrollPane;
import javafx.scene.control.Separator;
import javafx.scene.control.SeparatorMenuItem;
import javafx.scene.control.Spinner;
import javafx.scene.control.Tab;
import javafx.scene.control.TabPane;
import javafx.scene.control.TableColumn;
import javafx.scene.control.TableColumn.CellDataFeatures;
import javafx.scene.control.TableView;
import javafx.scene.control.TextArea;
import javafx.scene.control.TextField;
import javafx.scene.control.TextInputDialog;
import javafx.scene.control.Tooltip;
import javafx.scene.layout.Border;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.FlowPane;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.scene.layout.VBox;
import javafx.scene.paint.Color;
import javafx.stage.FileChooser;
import javafx.stage.FileChooser.ExtensionFilter;
import javafx.stage.Modality;
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
@SuppressWarnings("unused")
public class JfxNgs extends Application {
    private static final Logger LOG= Logger.build(JfxNgs.class).make();
    final Preferences preferences ;
    final Optional<Compilable> javascriptCompiler;
    private static final String LAST_USED_DIR_KEY="last.used.dir";
    private final List<NgsStage<?,?>> all_opened_stages=new ArrayList<>();

    /** class used to store key/value for preferences */
    static class PrefItem
    	{
    	final String key;
    	final String label;
    	final TextField textField;
    	String defaultValue="";
    	PrefItem(final String key, final String label,String description) {
    		this.key = key;
    		this.label= label;
    		if(description==null) description=label;
    		this.textField =new TextField();
    		this.textField.setPrefColumnCount(25);
    		this.textField.setTooltip(new Tooltip(description));
    		this.textField.setPromptText(description);
    		}
    	PrefItem setDefault(final String s) { this.defaultValue=s;this.textField.setText(s); return this;}
    	String getValue() { return this.textField.getText().trim();}
    	}

	private final PrefItem pref_httpHost = new PrefItem("http.proxyHost", "Http Proxy Host",null);
	private final PrefItem pref_httpPort = new PrefItem("http.proxyPort", "Http Proxy Port",null);
	private final PrefItem pref_httpsHost = new PrefItem("https.proxyHost", "Https Proxy Host",null);
	private final PrefItem pref_httpsPort = new PrefItem("https.proxyPort", "Https Proxy Port",null);
	private final PrefItem pref_max_sam_items = new PrefItem(BamStage.SPINNER_VALUE_KEY, "Default number of Reads",null);
	private final PrefItem pref_max_vcf_items = new PrefItem(VcfStage.SPINNER_VALUE_KEY, "Default number of Variants",null);
	
	final PrefItem pref_bam_max_seq_length_displayed = new PrefItem("bam.max.seq.length.displayed", "Max sequence length to be displayed",null).setDefault("1000");
	final PrefItem pref_bam_max_cigar_items_displayed = new PrefItem("bam.cigar.max.items", "Max number of cigar elements to be displayed",null).setDefault("50");
	final PrefItem pref_vcf_max_allele_length_displayed = new PrefItem("allele.max.length", "Max Allele size to be displayed",null).setDefault("20");
	
	final PrefItem pref_bioalcidae_max_stream = new PrefItem("bioalcidae.max.stream", "Max number of bytes to be written by bioalcidae using stream. -1=no limit",null).setDefault("-1");
	final PrefItem pref_bioalcidae_max_string = new PrefItem("bioalcidae.max.string", "Max number of bytes to be written by bioalcidae using textpane",null).setDefault("10000");
	
	
	
	private final PrefItem all_preferences[]=new PrefItem[]{
			pref_httpHost,pref_httpPort,
			pref_httpsHost,pref_httpsPort,
			pref_max_sam_items,
			pref_max_vcf_items,
			pref_bam_max_seq_length_displayed,
			pref_bam_max_cigar_items_displayed,
			pref_vcf_max_allele_length_displayed,
			pref_bioalcidae_max_stream,
			pref_bioalcidae_max_string
		};

	

    /** utility Function to convert base to Color */
    public static final Function<Character, Color> BASE2COLOR= new Function<Character, Color>() {
		@Override
		public Color apply(final Character c) {
			switch(c)
			{
			case 'A':case 'a': return Color.BLUE;
			case 'T': case 't' : return Color.GREEN;
			case 'C': case 'c': return Color.ORANGE;
			case 'G': case 'g': return Color.RED;
			default: break;
			}
		return Color.BLACK;
		}
	};


	/** utility to convert UCSC chrom to Ensembl , used web opening web browser*/
    public static final Function<String, String> ContigToEnseml= new Function<String, String>() {
		@Override
		public String apply(String c) {
			if(c==null) return null;

			if(c.toLowerCase().endsWith("_random") &&
				c.toLowerCase().startsWith("chr") &&
				c.toLowerCase().contains("_gl")
				)
				{
				c=c.split("[_]")[1];
				if(c.startsWith("gl")) return "GL"+c.substring(2)+".1";
				}
			else if(c.toLowerCase().startsWith("chrUn_gl"))
				{
				return "GL"+c.substring(8)+".1";
				}
			else if(c.toLowerCase().startsWith("chr"))
				{
				c=c.substring(3);
				}
			if(c.equals("M")) c="MT";

			return c;
			}
		};

	/** utility to convert Ensembl chrom to UCSC  , used web opening web browser*/
    public static final Function<String, String> ContigToUCSC = new Function<String, String>() {
		@Override
		public String apply(String c) {
			if(c==null) return null;
			if(!c.toLowerCase().startsWith("chr")) {
				c="chr"+c;
				}
			if(c.equals("chrMT")) c="chrM";

			return c;
			}
		};


    public JfxNgs()
		{
		this.preferences = Preferences.userNodeForPackage(JfxNgs.class);
		Compilable engine=null;
		try {
			final ScriptEngineManager manager = new ScriptEngineManager();
			final ScriptEngine scriptEngine = manager.getEngineByName("js");

			if(scriptEngine!=null)
				{
				if(!(scriptEngine instanceof Compilable)) {
					LOG.info("cannot find compilable nashorn");;
					}
				else
					{
					engine = (Compilable)scriptEngine;
					}
				}
			else
				{
				LOG.info("Cannot get instance of nashorn");
				}
			}
		catch(final Exception err)
			{
			err.printStackTrace();
			engine=null;
			LOG.warning("Cannot get Compilable JS engine "+err.getMessage());
			}
		this.javascriptCompiler = Optional.ofNullable(engine);
		/* init some prefs */
		if( this.preferences.get(pref_max_sam_items.key,null)==null)
			{
			this.preferences.put(pref_max_sam_items.key,String.valueOf(BamStage.DEFAULT_BAM_RECORDS_COUNT));
			}
		if( this.preferences.get(pref_max_vcf_items.key,null)==null)
			{
			this.preferences.put(pref_max_vcf_items.key,String.valueOf(VcfStage.DEFAULT_VCF_RECORDS_COUNT));
			}
		applyPreferences();
		}

    /** obtain a new file chooser, update the preferences if needed */
    FileChooser newFileChooser()
    	{
    	final FileChooser fc=new FileChooser();
    	final String lastDirStr= preferences.get(LAST_USED_DIR_KEY, null);
    	if(lastDirStr!=null && !lastDirStr.isEmpty())
    		{
    		fc.setInitialDirectory(new File(lastDirStr));
    		}
    	return fc;
    	}
    /** update the last user directory used by newFileChooser() */
    File updateLastDir(File f)
    	{
    	if(f==null || !f.exists()) return f;
    	File dir=f;
    	if(f.isFile()){
    		dir=f.getParentFile();
    		if(dir==null) return f;
    		}
    	this.preferences.put(LAST_USED_DIR_KEY, dir.getPath());
    	return f;
    	}

    /** at stop, save the preferences */
    @Override
    public void stop() throws Exception
    	{
    	try {
    		this.preferences.flush();
    		}
    	catch(final BackingStoreException err)
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

        menu.getItems().addAll(createCommonMenuItems(primaryStage));

        menu.getItems().add(new SeparatorMenuItem());
        MenuItem menuItem=new MenuItem("Quit...");
        menuItem.setOnAction(AE->doMenuQuit());
        menu.getItems().add(menuItem);

        MenuBar bar=new MenuBar(menu);
        FlowPane flow=new FlowPane(5,5);
        flow.setPadding(new Insets(10));
        flow.getChildren().add(new Label("Set Location of all frames to:"));
        final TextField textField=new TextField();
        textField.setPrefColumnCount(25);
        textField.setPromptText("Location. e:g '2:1234-5678'");
        flow.getChildren().add(textField);
        Button button=new Button("Go");
        flow.getChildren().add(button);

        textField.setTooltip(new Tooltip(
        		"set genomic location can be: empty, 'contig', 'contig:pos', 'contig:start-end' and (\"unmapped\" for bam)"
        		));

        final EventHandler<ActionEvent> handler=new EventHandler<ActionEvent>() {
			@Override
			public void handle(final ActionEvent event) {
				final String loc=textField.getText().trim();
				LOG.info("moveTo all to "+loc);
				for(final NgsStage<?,?> sc:all_opened_stages )
					{
					LOG.info("moveTo "+sc.getTitle()+" to "+loc);
					sc.moveTo(loc);
					}
				}
		};
        button.setOnAction(handler);
        button.setTooltip(new Tooltip("Go the specified genomic location."));
        textField.setOnAction(handler);


        BorderPane pane=new BorderPane();
        pane.setPadding(new Insets(5));

        pane.setBottom(new Label("Author: Pierre Lindenbaum PhD."));


        VBox vbox1= new VBox(bar,flow,pane);

        final Scene scene = new Scene(vbox1,500,300);

        primaryStage.setScene(scene);

        Exception lastException=null;

        primaryStage.addEventHandler(
    			WindowEvent.WINDOW_SHOWING ,new EventHandler<WindowEvent>() {
                    @Override
                    public void handle(final WindowEvent event) {
                    	final List<String> unnamedParams = new ArrayList<>(params.getUnnamed());
                    	String startPos="";
                    	int optind=0;
                    	while(optind+1< unnamedParams.size())
                    		{
                    		if(unnamedParams.get(optind).equals("-h") || unnamedParams.get(optind).equals("--help"))
	                			{
                    			unnamedParams.remove(optind);
	                			System.out.println("JfxNgs : Pierre Lindenbaum PhD 2017");
	                			System.out.println("Options:");
	                			System.out.println(" -h|--help this screen.");
	                			System.out.println(" -p|--position (string) the starting position");
	                			Platform.exit();
	                			}
                    		else if(unnamedParams.get(optind).equals("-p") || unnamedParams.get(optind).equals("--position"))
                    			{
                    			startPos=unnamedParams.get(optind+1);
                    			unnamedParams.remove(optind+1);
                    			unnamedParams.remove(optind);
                    			}
                    		else
                    			{
                    			optind++;
                    			}
                    		}
                    	
                    	for(final String arg: unnamedParams )
            	        	{
                    		VcfFile vcfin=null;
                    		BamFile bamin=null;

                    		try {
                    			if(IOUtil.isUrl(arg))
	                    			{
	                    			if(arg.endsWith(".bam")) {
	                    				bamin=BamFile.newInstance(arg);
	                    				}
	                    			else if(arg.endsWith(".vcf.gz")) {
	                    				vcfin=VcfFile.newInstance(arg);
	                    				}
	                    			}
                    			else
	                    			{
	                    			final File f=new File(arg);

									if(fileMatchExtensionFilter(f.getName(), BamStage.EXTENSION_FILTERS))
										{
										bamin=BamFile.newInstance(f);
										}
									else if(fileMatchExtensionFilter(f.getName(), VcfStage.EXTENSION_FILTERS))
										{
										vcfin=VcfFile.newInstance(f);
										}
									else
										{
										JfxNgs.showExceptionDialog(primaryStage,"Cannot open "+f);
										}
	                    			}
                    			if(vcfin!=null) {
                    				new VcfStage(JfxNgs.this, vcfin).setLocationOnOpen(startPos).show();
                    				}
                    			else if(bamin!=null) {
                    				new BamStage(JfxNgs.this, bamin).show();
                    				}
							} catch (Exception e) {
								CloserUtil.close(vcfin);
								CloserUtil.close(bamin);
								showExceptionDialog(primaryStage, e);
								}
            	        	}
                        }
                    });
        primaryStage.show();
        }

    void unregisterStage(final NgsStage<?,?> s) {
    	this.all_opened_stages.remove(s);
    }
    void registerStage(final NgsStage<?,?> s) {
    	this.all_opened_stages.add(s);
    }


    static void showExceptionDialog(final Window owner,Object err)
    	{
    	final Alert alert = new Alert(AlertType.WARNING);
    	alert.setTitle("Error");


		if(err!=null && err instanceof Throwable)
			{

			alert.setHeaderText("Error");

			final Throwable error=(Throwable) err;
			alert.setContentText(String.valueOf(error.getMessage()));

			error.printStackTrace();

			StringWriter sw = new StringWriter();
			PrintWriter pw = new PrintWriter(sw);
			error.printStackTrace(pw);
			String exceptionText = sw.toString();

			final Label label = new Label("The exception stacktrace was:");

			final TextArea textArea = new TextArea(exceptionText);
			textArea.setEditable(false);
			textArea.setWrapText(false);


			textArea.setMaxWidth(Double.MAX_VALUE);
			textArea.setMaxHeight(Double.MAX_VALUE);
			GridPane.setVgrow(textArea, Priority.ALWAYS);
			GridPane.setHgrow(textArea, Priority.ALWAYS);

			final BorderPane expContent = new BorderPane(new ScrollPane(textArea));
			expContent.setTop(label);
			expContent.setMinHeight(500);
			expContent.setMinWidth(800);
			expContent.setPrefWidth(1000);
			expContent.setPrefHeight(900);


			// Set expandable Exception into the dialog pane.
			alert.getDialogPane().setExpandableContent(expContent);
			}
		else
			{
			final String msg=String.valueOf(err);
			alert.setHeaderText(msg);
			alert.setContentText(msg);
			alert.getDialogPane().setExpandableContent(new Label(msg));
			}
		alert.showAndWait();
    	}

    private void doMenuQuit()
    	{
    	Platform.exit();
    	}


    void openBamUrl(final Window owner)
		{
    	final String lastkey="last.bam.url";
    	final TextInputDialog dialog=new TextInputDialog(this.preferences.get(lastkey, ""));
    	dialog.setContentText("URL:");
    	dialog.setTitle("Get BAM URL");
    	dialog.setHeaderText("Input BAM URL");
    	final Optional<String> choice=dialog.showAndWait();
    	if(!choice.isPresent()) return;
    	BamFile input=null;
    	try {
			input = BamFile.newInstance(choice.get());
			this.preferences.put(lastkey, choice.get());
			final BamStage stage=new BamStage(JfxNgs.this, input);
			stage.show();
		} catch (final Exception err) {
			CloserUtil.close(input);
			showExceptionDialog(owner, err);
			}

		}

    void openVcfUrl(final Window owner)
		{
    	final String lastkey="last.vcf.url";
		final TextInputDialog dialog=new TextInputDialog(this.preferences.get(lastkey, ""));
		dialog.setContentText("URL:");
		dialog.setTitle("Get VCF URL");
		dialog.setHeaderText("Input VCF URL");
		final Optional<String> choice=dialog.showAndWait();
		if(!choice.isPresent()) return;
    	VcfFile input=null;
    	try {
			input = VcfFile.newInstance(choice.get());
			this.preferences.put(lastkey, choice.get());
			VcfStage stage=new VcfStage(JfxNgs.this, input);
			stage.show();
		} catch (final Exception err) {
			CloserUtil.close(input);
			showExceptionDialog(owner, err);
			}
		}


    private void _openVcfAnySource(final Window owner,final String src)
    	{
    	VcfFile input=null;
    	try {
			input = VcfFile.newInstance(src);
			final VcfStage stage=new VcfStage(JfxNgs.this, input);
			stage.show();
		} catch (final Exception err) {
			CloserUtil.close(input);
			showExceptionDialog(owner, err);
			}
    	}

    private void doMenuOpenInExcel(final Window owner)
    	{
    	showExceptionDialog(owner,"DO NOT USE EXCEL");
    	}

    /** open index a BAM file */
    private void doMenuIndexBam(final Window owner)
		{
    	final FileChooser fc = newFileChooser();
    	fc.getExtensionFilters().addAll(BamStage.EXTENSION_FILTERS);
    	final List<File> files = fc.showOpenMultipleDialog(owner);
    	if(files==null) return;
    	for(final File file : files)
	    	{
    		LOG.info("indexing "+file);
    		updateLastDir(file);
	    	SamReader bam=null;
	    	try
	    		{
	    		final File output=  new File(file.getAbsolutePath() + BAMIndex.BAMIndexSuffix);
	    		if(output.exists())
	    			{
	    			throw new IOException("Bam index "+output+" already exists.");
	    			}
	    		bam = SamReaderFactory.makeDefault()
	    				.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS)
	    				.validationStringency(ValidationStringency.LENIENT)
	                    .open(file);
	    				;
	    		BAMIndexer.createIndex(bam, output);
	    		final Alert alert = new Alert(AlertType.CONFIRMATION, "Done. ?", ButtonType.OK);
				alert.showAndWait();
	    		}
	    	catch(final Exception err)
	    		{
				showExceptionDialog(owner, err);
				break;
	    		}
	    	finally
	    		{
	    		CloserUtil.close(bam);
	    		}
	    	}
		}



    List<MenuItem> createCommonMenuItems(final Stage stage) {
    	final List<MenuItem> L=new ArrayList<>();
    	L.add(createMenuItem("About...",()->doMenuAbout(stage)));
    	L.add(createMenuItem("Preferences...",()->showPreferencesDialog(stage)));
    	L.add(createMenuItem("Open in Microsoft Excel...",()->doMenuOpenInExcel(stage)));
        L.add(new SeparatorMenuItem());
    	L.add(createMenuItem("Open VCF/BAM File...",()->openNgsFile(stage)));
        L.add(createMenuItem("Open Remote BAM...",()->openBamUrl(stage)));
        L.add(createMenuItem("Open Remote VCF...",()->openVcfUrl(stage)));
        L.add(new SeparatorMenuItem());
        L.add(createMenuItem("Tool: index BAM file...",()->doMenuIndexBam(stage)));
        L.add(createMenuItem("Tool: index VCF file...",()->doMenuIndexVcf(stage)));
        L.add(new SeparatorMenuItem());
        L.add(createMenuItem("Close",()->stage.hide()));
        return L;
    	}

    private File tryDownloadBamIndex(final String url)
    	throws IOException
    	{
    	if(!url.endsWith(BamFileIoUtils.BAM_FILE_EXTENSION)) throw new IllegalStateException("bam suffix not checked");
    	final File baiFile=File.createTempFile("tmp.", ".bai");
    	baiFile.deleteOnExit();
    	for(int i=0;i< 2;++i)
    		{
    		final String baiurl=(i==0?
    				url:url.substring(0, url.length()-BamFileIoUtils.BAM_FILE_EXTENSION.length())
    				)+".bai";
    		InputStream in=null;
    		FileOutputStream out=null;
    		try {
    			LOG.info("trying "+baiurl);
				in = new URL(baiurl).openStream();
				out = new FileOutputStream(baiFile);
				IOUtil.copyStream(in, out);
				out.flush();
				out.close();
				in.close();
				return baiFile;
			} catch (final IOException err) {
				baiFile.delete();
				LOG.warning("Cannot fetch "+baiurl+" : "+err.getMessage());
				}
    		finally
    			{
    			CloserUtil.close(out);
    			CloserUtil.close(in);
    			}
    		}
    	throw new IOException("cannot get bam index for "+url);
    	}


    void openNgsFile(final Window owner)
		{
		final FileChooser fc = newFileChooser();
		
		final List<String> suffixes=new ArrayList<>();
		for(FileChooser.ExtensionFilter e:VcfStage.EXTENSION_FILTERS)
			suffixes.addAll(e.getExtensions());
		for(FileChooser.ExtensionFilter e:BamStage.EXTENSION_FILTERS)
			suffixes.addAll(e.getExtensions());
		
		
		fc.getExtensionFilters().add(new FileChooser.ExtensionFilter("BAM/VCF",suffixes));
		final List<File> selFiles = fc.showOpenMultipleDialog(owner);
		if(selFiles==null) return ;
		for(final File file: selFiles) {
			updateLastDir(file);
			if(fileMatchExtensionFilter(file.getName(), VcfStage.EXTENSION_FILTERS))
				{
				VcfFile input=null;
		    	try {
					input = VcfFile.newInstance(file);
					VcfStage stage=new VcfStage(JfxNgs.this, input);
					stage.show();
				} catch (final Exception err) {
					CloserUtil.close(input);
					showExceptionDialog(owner, err);
					}
				}
			else if(fileMatchExtensionFilter(file.getName(), BamStage.EXTENSION_FILTERS))
				{
				BamFile input=null;
		    	try {
					input = BamFile.newInstance(file);
					final BamStage stage=new BamStage(JfxNgs.this, input);
					stage.show();
				} catch (final Exception err) {
					CloserUtil.close(input);
					showExceptionDialog(owner, err);
					}
		    	}
			else
				{
				LOG.info("Cannot open "+file);;
				}
			}
		}

    /** open index a VCF file */
    private void doMenuIndexVcf(final Window owner)
		{
    	final FileChooser fc = newFileChooser();
    	fc.getExtensionFilters().addAll(VcfStage.EXTENSION_FILTERS);
    	final List<File> files = fc.showOpenMultipleDialog(owner);
    	if(files==null) return;
    	for(final File file : files)
	    	{
    		updateLastDir(file);
    		if(file.getName().endsWith(".vcf.gz"))
    			{
    			LOG.info("writing tabix index for "+file);
    			final File output=new File(file.getAbsolutePath()+TabixUtils.STANDARD_INDEX_EXTENSION);
    			try
    				{
    				if(output.exists())
		    			{
		    			throw new IOException("Tabix index "+output+" already exists.");
		    			}
    				final TabixIndex index=IndexFactory.createTabixIndex(file,new VCFCodec(),(SAMSequenceDictionary)null);
    				index.write(output);
					final Alert alert = new Alert(AlertType.CONFIRMATION, "Done. ?", ButtonType.OK);
					alert.showAndWait();
    				}
    			catch(final Exception err)
    				{
    				showExceptionDialog(owner, err);
    				break;
    				}
    			}
    		else if(file.getName().endsWith(".vcf"))
    			{
    			LOG.info("writing tribble index for "+file);
    			final File output=new File(file.getAbsolutePath()+Tribble.STANDARD_INDEX_EXTENSION);
    			try
					{
					if(output.exists())
		    			{
		    			throw new IOException("Tribble index "+output+" already exists.");
		    			}
					final Index index=IndexFactory.createIndex(file,new VCFCodec(),IndexType.LINEAR);
					index.writeBasedOnFeatureFile(file);
					final Alert alert = new Alert(AlertType.CONFIRMATION, "Done. ?", ButtonType.OK);
					alert.showAndWait();
					}
    			catch(final Exception err)
					{
					showExceptionDialog(owner, err);
					break;
					}
				}
    		else
    			{
    			showExceptionDialog(owner,"Cannot index file "+file);
    			break;
    			}
	    	}
		}




    static void doMenuAbout(final Stage stage)
    	{
    	final Alert alert=new Alert(AlertType.INFORMATION);
		alert.setHeaderText("JFXNGS");
		alert.setContentText("Pierre Lindenbaum PhD. 2017.\n"
				+ "@yokofakun\n"
				+ "Institut du Thorax - Nantes - France\n"
				+ "https://github.com/lindenb/jvarkit"
				);
		alert.showAndWait();
    	}

    /** utiliti for create a menuItem */
    static MenuItem createMenuItem(final String label,final Runnable runner)
    	{
    	final MenuItem menu=new MenuItem(label);
    	menu.setOnAction(AE->runner.run());
    	return menu;
    	}

    private  void showPreferencesDialog(final Stage parentStage)
    	{
    	Stage dialog = new Stage();
    	dialog.setTitle("Preferences...");
    	// populate dialog with controls.
    	final VBox vbox=new VBox();
    	for(final PrefItem pref: all_preferences)
    		{
    		final Label label=new Label(pref.label+" :");
        	final HBox hbox=new HBox(label,pref.textField);
        	hbox.setPadding(new Insets(10));
        	pref.textField.setText(this.preferences.get(pref.key, pref.defaultValue));
        	vbox.getChildren().add(hbox);
    		}
    	vbox.setPadding(new Insets(5));
    	dialog.setScene(new Scene(vbox));
    	dialog.initOwner(parentStage);
    	dialog.initModality(Modality.APPLICATION_MODAL);
    	dialog.showAndWait();

    	for(final PrefItem pref: this.all_preferences)
    		{
    		this.preferences.put(pref.key, pref.getValue());
    		}
    	applyPreferences();
    	}

    /** apply preferences after application is launched of after the pref dialog was opened*/
    private void applyPreferences()
    	{
    	for(final PrefItem pref:new PrefItem[]{pref_httpHost,pref_httpPort,pref_httpsHost,pref_httpsPort})
			{
			System.setProperty(pref.key,this.preferences.get(pref.key, ""));
			}
    	}
    
    /** utility function to test if a file match FileChooser.ExtFilter*/
    static boolean fileMatchExtensionFilter(final String f,final Collection<FileChooser.ExtensionFilter> extensionFilters)
    	{
    	if(f==null) return false;
    	for(final FileChooser.ExtensionFilter extensionFilter:extensionFilters)
    		{
    		for(String ext:extensionFilter.getExtensions()) {
    			if(ext.equals("*.*")) continue;
    			if(ext.startsWith("*.")) ext=ext.substring(2);
    			if(f.endsWith(ext)) return true;
    			}
    		}
    	return false;
    	}
    
    public static void main(String[] args) {
    	launch(args);
    }

}
