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
package com.github.lindenb.jvarkit.tools.igvreview;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.Socket;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Vector;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.prefs.BackingStoreException;
import java.util.prefs.Preferences;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.jfx.components.JfxUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.igv.IgvConstants;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.JexlVariantPredicate;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import javafx.application.Application;
import javafx.application.Platform;
import javafx.beans.property.SimpleStringProperty;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.geometry.Insets;
import javafx.scene.Scene;
import javafx.scene.control.ContextMenu;
import javafx.scene.control.Menu;
import javafx.scene.control.MenuBar;
import javafx.scene.control.MenuItem;
import javafx.scene.control.SeparatorMenuItem;
import javafx.scene.control.TableColumn;
import javafx.scene.control.TableView;
import javafx.scene.control.cell.TextFieldTableCell;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.VBox;
import javafx.stage.FileChooser;
import javafx.stage.Stage;
import javafx.stage.Window;

public class IgvReview extends Application
	{
	private static final Logger LOG = Logger.build(IgvReview.class).make();
	private static final String LAST_USED_DIR_KEY="last.used.dir";
    final Preferences preferences ;
    private List<VariantContext> variants = new Vector<>();
	private TableView<VariantAndGenotype> genotypeTable = null;
	private File saveAsFile=null;
    private VCFHeader vcfHeader = null;
    private Map<String,File> sample2bamFile=new HashMap<>();
    private final VCFFormatHeaderLine reviewFormat;
    private static class IgvConnection 
    	implements Closeable
    	{
    	private final Socket socket;
    	private final PrintWriter out;
    	private final  BufferedReader in;
    	
    	IgvConnection() throws IOException {
    		this.socket = new Socket(
    				IgvConstants.DEFAULT_HOST,
    				IgvConstants.DEFAULT_PORT
    				);
            this.out = new PrintWriter(socket.getOutputStream(), true);
            this.in = new BufferedReader(new InputStreamReader(socket.getInputStream()));
    		}
    	
    	void send(final String cmd) {
    		try {
    			LOG.info("sending "+cmd);
    			this.out.println(cmd);
    			final String msg=in.readLine();
    			LOG.info("IGV replied :"+msg);
    		} catch(IOException err)
    			{
    			LOG.error(err);
    			}
    		}
    	
    	@Override
    	public void close() throws IOException {
    		CloserUtil.close(this.out);
    		CloserUtil.close(this.in);
    		CloserUtil.close(this.socket);
    		}
    	}
    
 
    
    public class VariantAndGenotype
    	{
    	final VariantContext ctx;
    	final Genotype g;
    	private SimpleStringProperty reviewProperty;
    	VariantAndGenotype(final VariantContext ctx,final Genotype g) {
    		this.ctx = ctx;
    		this.g = g;
    		
    		String att=this.g.getAttributeAsString(IgvReview.this.reviewFormat.getID(),"");
    		this.reviewProperty = new SimpleStringProperty(att);
    		}
    	
    	public SimpleStringProperty getReviewProperty() {
			return reviewProperty;
		}
    	
    	public String getReview() {
			return getReviewProperty().get();
			}
    	public void setReview(String review) {
    		getReviewProperty().set(review);
			}
    	public Genotype makeGenotype() {
    		GenotypeBuilder gb= new GenotypeBuilder(this.g);
    		String fmt = getReview(); if(fmt==null) fmt="";
    		fmt= fmt.trim().replaceAll("[ \n\t\\:_]+","_");
        	gb.attribute(IgvReview.this.reviewFormat.getID(), fmt);
        	return gb.make();
    		}
    	
    	}
    
	public IgvReview()
		{
		this.preferences = Preferences.userNodeForPackage(IgvReview.class);
		
		final String userName=System.getProperty("user.name", "user").replace(" ", "_");
		this.reviewFormat = new VCFFormatHeaderLine(userName+"_REVIEW",1,VCFHeaderLineType.String,"Review genotypes by "+userName);
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
    	CloserUtil.close(this.igvConnection);
    	super.stop();
    	}

	
	@Override
	public void start(final Stage stage) throws Exception {
		stage.setTitle(getClass().getSimpleName());
		Predicate<VariantContext> ctxFilter;
		
		Map<String,String> params=super.getParameters().getNamed();
		if(params.containsKey("--filter")) {
			ctxFilter = JexlVariantPredicate.create(params.get("--filter"));
		} else
			{
			ctxFilter = V->true;
			}
		
		
		final List<String> args=super.getParameters().getUnnamed();
		final File configFile;
		if(args.isEmpty())
			{
	    	final FileChooser fc=new FileChooser();
	    	final String lastDirStr= preferences.get(LAST_USED_DIR_KEY, null);
	    	if(lastDirStr!=null && !lastDirStr.isEmpty())
	    		{
	    		fc.setInitialDirectory(new File(lastDirStr));
	    		}
	    	fc.getExtensionFilters().addAll(Collections.singletonList(new FileChooser.ExtensionFilter("Config file", "*.config","*.cfg","*.list")));
	    	configFile = fc.showOpenDialog(stage);
			}
		else if(args.size()==1)
			{
			configFile=new File(args.get(0));
			}
		else
			{
			configFile=null;
			}
		if(configFile==null || !configFile.exists()) {
			JfxUtils.dialog().cause("Illegal number of arguments or file doesn't exists.").show(stage);
			Platform.exit();
			return;
			}
		
    	if(configFile.isFile() && configFile.getParentFile()!=null){
    		this.preferences.put(LAST_USED_DIR_KEY, configFile.getParentFile().getPath());
    		}
    	
		
		final List<String> configLines = Files.readAllLines(configFile.toPath());
		final Predicate<String> ignoreComment=(S)->!S.startsWith("#");
		final Predicate<String> predVcf=S->S.endsWith(".vcf") || S.endsWith(".vcf.gz");
		
		
		if(configLines.stream().filter(ignoreComment).filter(predVcf).count()!=1)
			{
			JfxUtils.dialog().cause("Found more than one vcf file in "+configFile).show(stage);
			Platform.exit();
			return;
			}
		final File vcfFile = configLines.stream().filter(ignoreComment).filter(predVcf).
				map(S->new File(S)).
				findFirst().get();
		
		
		
		LOG.info("Opening vcf file and loading in memory");
		VCFFileReader vfr=null;
		CloseableIterator<VariantContext> iter=null;
		final Set<String> sampleNames;
		try {
			this.variants.clear();
			vfr=new VCFFileReader(vcfFile, false);
			
			this.vcfHeader = vfr.getFileHeader();
			sampleNames = new HashSet<>(this.vcfHeader.getSampleNamesInOrder());
			if(sampleNames.isEmpty())
				{
				JfxUtils.dialog().cause("No Genotypes in "+vcfFile).show(stage);
				Platform.exit();
				return;
				}
			
			iter = vfr.iterator();
			this.variants.addAll(
					iter.stream().
					filter(ctxFilter).
					filter(CTX->CTX.isVariant()).
					collect(Collectors.toList())
					);
			}
		catch(final Exception err)
			{
			JfxUtils.dialog().cause(err).show(stage);
			Platform.exit();
			return;
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(vfr);
			}
		if(this.variants.isEmpty())
			{
			JfxUtils.dialog().cause("No Variants").show(stage);
			Platform.exit();
			return;
			}
		
		final SAMSequenceDictionary dict = this.vcfHeader.getSequenceDictionary();
		if(dict==null || dict.isEmpty())
			{
			JfxUtils.dialog().cause(JvarkitException.VcfDictionaryMissing.getMessage(vcfFile.getPath())).show(stage);
			Platform.exit();
			return;
			}
		
		final SamReaderFactory srf =SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT);
		
		
		
		configLines.stream().filter(ignoreComment).
			filter(S->S.endsWith(".bam")).
			map(S->new File(S)).
			forEach(F->{
				final SamReader samIn = srf.open(F);
				final SAMFileHeader header = samIn.getFileHeader();
				CloserUtil.close(samIn);
				String sample=null;
				for(final SAMReadGroupRecord rg:header.getReadGroups())
					{
					String s = rg.getSample();
					if(s==null) continue;
					if(sample==null) {
						sample=s;
						}
					else if(!sample.equals(s))
						{
						JfxUtils.dialog().cause("Two samples in "+F).show(stage);
						Platform.exit();
						return;
						}
					}
				if(sample==null) {
					JfxUtils.dialog().cause("No sample in "+F+". Ignoring").show(stage);
					return;
					}
				if(!sampleNames.contains(sample))
					{
					JfxUtils.dialog().cause("Not in VCF header "+sample+" / "+F+". Ignoring").show(stage);
					return;
					}
				this.sample2bamFile.put(sample, F);
			});
		
		if(this.sample2bamFile.isEmpty())
			{
			JfxUtils.dialog().cause("No valid bam file in "+configFile).show(stage);
			return;
			}
		sampleNames.retainAll(this.sample2bamFile.keySet());
		if(sampleNames.isEmpty())
			{
			JfxUtils.dialog().cause("No Sample associated to bam").show(stage);
			return;
			}
		
		 ObservableList<VariantAndGenotype> genotypes = FXCollections.observableArrayList(
			this.variants.stream().
				flatMap(CTX->CTX.getGenotypes().stream().
						filter(G->sampleNames.contains(G.getSampleName()))
						.map(G->new VariantAndGenotype(CTX, G))).
				collect(Collectors.toList())
				);
		if(genotypes.isEmpty())
			{
			JfxUtils.dialog().cause("No Genotype to show").show(stage);
			return;
			}
		 Menu menu=new Menu("File");
		 MenuItem menuItem = new MenuItem("Save as...");
		 menuItem.setOnAction(AE->{
			 saveVariantsAs(stage);
			 
		 });
		 menu.getItems().add(menuItem);
		 menuItem = new MenuItem("Save");
		 menuItem.setOnAction(AE-> {
			 if(this.saveAsFile!=null)
				 {
				 saveVariants(stage,this.saveAsFile);
				 }
			 else
			 	{
				saveVariantsAs(stage); 
			 	}
		 	});
		 menu.getItems().add(menuItem);

		 menu.getItems().add(new SeparatorMenuItem());
		 menuItem = new MenuItem("Quit");
		 menuItem.setOnAction(AE->{Platform.exit();});
		 menu.getItems().add(menuItem);
		 MenuBar bar=new MenuBar(menu);
		 
		 

		 this.genotypeTable = new TableView<>(genotypes);
		 this.genotypeTable.getColumns().add(makeColumn("CHROM", G->G.ctx.getContig()));
		 this.genotypeTable.getColumns().add(makeColumn("POS", G->G.ctx.getStart()));
		 this.genotypeTable.getColumns().add(makeColumn("ID", G->G.ctx.getID()));
		 this.genotypeTable.getColumns().add(makeColumn("REF", G->G.ctx.getReference().getDisplayString()));
		 this.genotypeTable.getColumns().add(makeColumn("ALT", G->G.ctx.getAlternateAlleles().stream().map(A->A.getDisplayString()).collect(Collectors.joining(","))));
		 this.genotypeTable.getColumns().add(makeColumn("Sample",G->G.g.getSampleName()));
		 this.genotypeTable.getColumns().add(makeColumn("Type",G->G.g.getType().name()));
		 this.genotypeTable.getColumns().add(makeColumn("Alleles", G->G.g.getAlleles().stream().map(A->A.getDisplayString()).collect(Collectors.joining(","))));
		 
		 TableColumn<VariantAndGenotype,String> reviewCol = new TableColumn<>("Review");
		 reviewCol.setCellValueFactory(C->C.getValue().getReviewProperty());
		 reviewCol.setCellFactory(TextFieldTableCell.forTableColumn());
		 reviewCol.setOnEditCommit((E)->{
			 int y=E.getTablePosition().getRow();
			 this.genotypeTable.getItems().get(y).setReview(E.getNewValue());
		 	});
		 reviewCol.setEditable(true);
		 
		 
		 this.genotypeTable.getColumns().add(reviewCol);
		 this.genotypeTable.getSelectionModel().cellSelectionEnabledProperty().set(true);
		 this.genotypeTable.setEditable(true);
		 
		 final ContextMenu cm = new ContextMenu();
		 MenuItem mi1 = new MenuItem("Menu 1");
	     cm.getItems().add(mi1);
	     MenuItem mi2 = new MenuItem("Menu 2");
	     cm.getItems().add(mi2);

		 
		 this.genotypeTable.setOnMousePressed(event->{
	        if (event.isPrimaryButtonDown() && (event.getClickCount() == 3 || event.isShiftDown())) {
	        	
	            moveIgvTo(stage,genotypeTable.getSelectionModel().getSelectedItem());                   
	        	}
	        else if(event.isSecondaryButtonDown())
	            {
	        	cm.show(genotypeTable ,event.getScreenX() ,event.getScreenY());
	            }
		 	});
		 
		 final BorderPane pane2=new BorderPane(this.genotypeTable);
		 pane2.setPadding(new Insets(10, 10, 10, 10));
		 
		VBox vbox1= new VBox(bar,pane2);
		final Scene scene = new Scene(vbox1,500,300);
		stage.setScene(scene);
		stage.show();
		}
    protected <T,R> TableColumn<T,R> makeColumn(final String tag,final Function<T,R> supplier)
		{
	    return JfxUtils.makeTableColumn(tag, supplier);
		}
    
    
    private File previousBamDisplayedInIgv=null;
    private IgvConnection igvConnection=null;
    private void moveIgvTo(Window win,final VariantAndGenotype vg) {
    	if(vg==null) return;
    	if(igvConnection==null) {
    		try {
    			LOG.info("trying to connect to IGV...");
    			this.igvConnection = new IgvConnection();
    			}
    		catch(Exception err)
    			{
    			this.igvConnection = null;
    			JfxUtils.dialog().cause(err).show(win);
    			return;
    			}
    		}
    	final File bam = this.sample2bamFile.get(vg.g.getSampleName());
    	
    	
    	if(previousBamDisplayedInIgv==null || !previousBamDisplayedInIgv.equals(bam))
    		{
    		if(previousBamDisplayedInIgv!=null) 
    			{
    			this.igvConnection.send("remove "+previousBamDisplayedInIgv.getName());
        		sleep(2000);
    			this.igvConnection.send("remove \""+previousBamDisplayedInIgv.getName()+" Coverage\"");
        		sleep(1000);
    			}
    		this.igvConnection.send("load "+bam.getPath()+"");
    		sleep(3000);
    		}
    	final int extend=5;
    	this.igvConnection.send("goto "+vg.ctx.getContig()+":"+
    			Math.max(1, vg.ctx.getStart()-extend)+
    			"-"+
    			(vg.ctx.getEnd()+extend)
    			);
    	sleep(1000);
    	
    	previousBamDisplayedInIgv = bam;
    	}
    
    private void sleep(final long tms) {
    	try {
    		Thread.sleep(tms);
    	} catch(Exception err) {LOG.error(err);}
    }
    private File saveVariantsAs(Window stage) {
		 FileChooser fc = new FileChooser();
		 fc.getExtensionFilters().addAll(JfxUtils.getIndexedVcfExtensionFilters());
		 File vcfout=fc.showSaveDialog(stage);
		 if(vcfout==null) return null;
		 this.saveAsFile=saveVariants(stage,vcfout);
		 return this.saveAsFile;
    	}
    
    private File saveVariants(Window owner,final File f) {
    	if(f==null) return null;
    	VariantContextWriterBuilder vcb = new VariantContextWriterBuilder();
    	VariantContextWriter w=null;
    	try {
    		SAMSequenceDictionary dict=this.vcfHeader.getSequenceDictionary();
    		if(dict!=null)
    			{
    			vcb.setReferenceDictionary(dict);
    			}
    		vcb.setOutputFile(f);
    		final VCFHeader header2=new VCFHeader(this.vcfHeader);
    		if(header2.getFormatHeaderLine(this.reviewFormat.getID())==null)
				{
    			header2.addMetaDataLine(this.reviewFormat);
				}
    		w = vcb.build();
    		w.writeHeader(header2);
    		int x=0;
    		while(x<this.genotypeTable.getItems().size())
    			{
    			VariantContext v1= this.genotypeTable.getItems().get(x).ctx;
    			List<Genotype> genotypes = new ArrayList<>();
    			genotypes.add(this.genotypeTable.getItems().get(x).makeGenotype());

    			int y=x+1;
    			while(y<this.genotypeTable.getItems().size())
    				{
        			VariantContext v2= this.genotypeTable.getItems().get(y).ctx;
        			if(v2!=v1) break;//yes '!=' is enough
        			genotypes.add(this.genotypeTable.getItems().get(y).makeGenotype());
    				y++;
    				}
    			VariantContextBuilder vb=new VariantContextBuilder(v1);
    			vb.genotypes(genotypes);
    			w.add(vb.make());
    			x=y;
    			}
    		w.close();
    		return f;
    		}
    	catch(final Exception err) {
    		JfxUtils.dialog().cause(err).show(owner);
    		return null;
    		}
    	finally
    		{	
    		CloserUtil.close(w);
    		}
    }
    
    public static class IgvReviewLauncher extends Launcher {

    	@Override
    	public int doWork(final List<String> args) {
    		Application.launch(IgvReview.class,args.toArray(new String[args.size()]));
    		return 0;
    		}
    	}
    
	public static void main(final String[] args) {
		new IgvReviewLauncher().instanceMain(args);//not instanceMainWithExit
		}
	}

