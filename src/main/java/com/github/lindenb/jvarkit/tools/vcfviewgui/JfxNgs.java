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
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
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
@SuppressWarnings("unused")
public class JfxNgs extends Application {
    private static final Logger LOG= Logger.getLogger("JfxNgs");
    private final Preferences preferences ;
    final Compilable javascriptEngine;
    private static final String LAST_USED_DIR_KEY="last.used.dir";
    private final List<NgsStage> all_opened_stages=new ArrayList<>();
    
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
    	preferences.get(LAST_USED_DIR_KEY, dir.getPath());
    	return f;
    	}
    
    /** at stop, save the preferences */
    @Override
    public void stop() throws Exception
    	{
    	try {
    		LOG.info("flush preferences");
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
				for(final NgsStage sc:all_opened_stages )
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

        final Scene scene = new Scene(vbox1,500,300);
       
        primaryStage.setScene(scene);
        primaryStage.setOnHidden(e -> doMenuQuit());
       
        Collection<NgsStage> newStages=new ArrayList<>();
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
        for(NgsStage sc:newStages)
	    	{
	    	sc.show();
	    	}
        }

    void unregisterStage(NgsStage s) {
    	this.all_opened_stages.remove(s);
    	LOG.info("unregister nbr.win"+this.all_opened_stages.size());

    }
    void registerStage(NgsStage s) {
    	this.all_opened_stages.add(s);  
    	LOG.info("register nbr.win"+this.all_opened_stages.size());
    }
    
    
    static void showExceptionDialog(final Window owner,Throwable error)
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
    
    void openNgsFiles(final Window owner)
    	{
    	final FileChooser fc=newFileChooser();
    	
    	fc.setSelectedExtensionFilter(new ExtensionFilter("NGS Files", "*.bam","*.vcf","*.vcf.gz","*.list"));
    	final List<File> selFiles = fc.showOpenMultipleDialog(owner);
    	if(selFiles==null || selFiles.isEmpty()) return ;
    	List<NgsStage> stages =new ArrayList<>();
    	try 
    		{
	    	for(final File f:selFiles)
	    		{
	    		stages.addAll( openNgsFiles(f.getPath()) );
	    		updateLastDir(f.getParentFile());
	    		}
	    	for(NgsStage sc:stages)
	    		{	
	    		sc.show();
	    		}
    		}
    	catch(final Exception err)
    		{
    		for(NgsStage sc:stages) sc.closeNgsResource();
    		showExceptionDialog(owner, err);
    		}		
    	}

    Collection<NgsStage> openNgsFiles(final String path0) throws Exception
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
		final List<NgsStage> stages =new ArrayList<>();

    	for(final String uri:pathList)
	    	{
	    	//try as BAM
	    	if(uri.endsWith(".bam") )
		    	{
	    		final NgsStage sc=openBam(uri);
		    	if(sc!=null) stages.add(sc);
		    	}
	    	else if(uri.endsWith(".vcf") || uri.endsWith(".bcf")  || uri.endsWith(".vcf.gz") )
		    	{
	    		final NgsStage sc=openVcf(uri);
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
    
    private NgsStage openBam(final String uri)
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
			return new BamStage(this,uri);
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
    		File file=new File(uri);
    		vcfIn = new VCFFileReader(file, true);
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
				return new VcfStage(this,file);
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
    
    static MenuItem createMenuItem(final String label,final Runnable runner)
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

