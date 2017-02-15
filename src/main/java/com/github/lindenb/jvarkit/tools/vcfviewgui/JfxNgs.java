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
import java.util.logging.Level;
import java.util.logging.Logger;
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

import com.github.lindenb.jvarkit.tools.vcfviewgui.VcfStage.VariantFileReader;
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
    final Preferences preferences ;
    final Optional<Compilable> javascriptCompiler;
    private static final String LAST_USED_DIR_KEY="last.used.dir";
    private final List<NgsStage<?,?>> all_opened_stages=new ArrayList<>();
    
    /** utility Function to convert base to Color */
    public static final Function<Character, Color> BASE2COLOR= new Function<Character, Color>() {
		@Override
		public Color apply(final Character c) {
			switch(c)
			{
			case 'A':case 'a': return Color.BLUE;
			case 'T': case 't' : return Color.GREEN;
			case 'C': case 'c': return Color.YELLOW;
			case 'G': case 'g': return Color.RED;
			default: break;
			}    	
		return Color.BLACK;
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
        
        menu.getItems().addAll(createCommonMenuItems(primaryStage));
        
        menu.getItems().add(new SeparatorMenuItem());
        MenuItem menuItem=new MenuItem("Quit...");
        menuItem.setOnAction(AE->doMenuQuit());
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
        button.setTooltip(new Tooltip("Go the the specified genomic location."));
        textField.setOnAction(handler);
        
        
        BorderPane pane=new BorderPane();
        
        
        pane.setBottom(new Label("Author: Pierre Lindenbaum PhD."));
        
       
        VBox vbox1= new VBox(bar,flow,pane);

        final Scene scene = new Scene(vbox1,500,300);
       
        primaryStage.setScene(scene);
        primaryStage.setOnHidden(e -> doMenuQuit());
       
        Exception lastException=null;
        
        primaryStage.addEventHandler(
    			WindowEvent.WINDOW_SHOWING ,new EventHandler<WindowEvent>() {
                    @Override
                    public void handle(final WindowEvent event) {
                    	for(final String arg: params.getUnnamed())
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
	                        		
									if(BamStage.EXTENSION_FILTER.getExtensions().stream().filter(S->f.getName().endsWith(S)).findAny().isPresent())
										{
										bamin=BamFile.newInstance(f); 
										}
									else if(VcfStage.EXTENSION_FILTER.getExtensions().stream().filter(S->f.getName().endsWith(S)).findAny().isPresent())
										{
										vcfin=VcfFile.newInstance(f);
										}
									else
										{
										JfxNgs.showExceptionDialog(primaryStage,"Cannot open "+f);
										}
	                    			}
                    			if(vcfin!=null) {
                    				new VcfStage(JfxNgs.this, vcfin).show();
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
			error.printStackTrace();
			
			StringWriter sw = new StringWriter();
			PrintWriter pw = new PrintWriter(sw);
			error.printStackTrace(pw);
			String exceptionText = sw.toString();
	
			final Label label = new Label("The exception stacktrace was:");
	
			final TextArea textArea = new TextArea(exceptionText);
			textArea.setEditable(true);
			textArea.setWrapText(true);
	
			textArea.setMaxWidth(Double.MAX_VALUE);
			textArea.setMaxHeight(Double.MAX_VALUE);
			GridPane.setVgrow(textArea, Priority.ALWAYS);
			GridPane.setHgrow(textArea, Priority.ALWAYS);
	
			final GridPane expContent = new GridPane();
			expContent.setMaxWidth(Double.MAX_VALUE);
			expContent.add(label, 0, 0);
			expContent.add(textArea, 0, 1);
	
			// Set expandable Exception into the dialog pane.
			alert.getDialogPane().setExpandableContent(expContent);
			}
		else
			{
			String msg=String.valueOf(err);
			alert.setHeaderText(msg);
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
    	final TextInputDialog dialog=new TextInputDialog();
    	dialog.setTitle("Get URL");
    	dialog.setHeaderText("Input BAM URL");
    	final Optional<String> choice=dialog.showAndWait();
    	if(!choice.isPresent()) return;
    	BamFile input=null;
    	try {
			input = BamFile.newInstance(choice.get());
			BamStage stage=new BamStage(JfxNgs.this, input);
			stage.show();
		} catch (final Exception err) {
			CloserUtil.close(input);
			showExceptionDialog(owner, err);
			}

		}
    
    void openVcfUrl(final Window owner)
		{
		final TextInputDialog dialog=new TextInputDialog();
		dialog.setTitle("Get URL");
		dialog.setHeaderText("Input VCF URL");
		final Optional<String> choice=dialog.showAndWait();
		if(!choice.isPresent()) return;
    	VcfFile input=null;
    	try {
			input = VcfFile.newInstance(choice.get());
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
			VcfStage stage=new VcfStage(JfxNgs.this, input);
			stage.show();
		} catch (final Exception err) {
			CloserUtil.close(input);
			showExceptionDialog(owner, err);
			}
    	}
    
    
    
    /** open FileChoser to open a BAM file */
    void openBamFile(final Window owner)
    	{
    	final FileChooser fc = newFileChooser();
    	fc.setSelectedExtensionFilter(BamStage.EXTENSION_FILTER);
    	final File file = updateLastDir(fc.showOpenDialog(owner));
    	if(file==null) return;
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

    List<MenuItem> createCommonMenuItems(final Stage stage) {
    	final List<MenuItem> L=new ArrayList<>();
    	L.add(createMenuItem("About...",()->doMenuAbout(stage)));
        L.add(new SeparatorMenuItem());
    	L.add(createMenuItem("Open VCF File...",()->openVcfFile(stage)));
        L.add(createMenuItem("Open BAM File...",()->openBamFile(stage)));
        L.add(createMenuItem("Open Remote BAM...",()->openBamUrl(stage)));
        L.add(createMenuItem("Open Remote VCF...",()->openVcfUrl(stage)));
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
    
    
    void openVcfFile(final Window owner)
		{
		final FileChooser fc = newFileChooser();
		fc.setSelectedExtensionFilter(VcfStage.EXTENSION_FILTER);
		final File file = updateLastDir(fc.showOpenDialog(owner));
		if(file==null) return;
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

    public static void main(String[] args) {
    	launch(args);
    }

}

