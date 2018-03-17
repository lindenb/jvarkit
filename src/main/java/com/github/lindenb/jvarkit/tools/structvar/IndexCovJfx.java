/*
The MIT License (MIT)

Copyright (c) 2018 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.structvar;
import java.io.BufferedReader;
import java.io.File;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Optional;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.Hershey;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.StringUtil;
import javafx.application.Application;
import javafx.application.Platform;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.geometry.Insets;
import javafx.geometry.Rectangle2D;
import javafx.scene.Scene;
import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.control.Alert;
import javafx.scene.control.Alert.AlertType;
import javafx.scene.control.ButtonType;
import javafx.scene.control.Label;
import javafx.scene.control.ListView;
import javafx.scene.control.Menu;
import javafx.scene.control.MenuBar;
import javafx.scene.control.MenuItem;
import javafx.scene.control.ScrollPane;
import javafx.scene.control.ScrollPane.ScrollBarPolicy;
import javafx.scene.control.SelectionMode;
import javafx.scene.control.SeparatorMenuItem;
import javafx.scene.control.Spinner;
import javafx.scene.control.TextInputDialog;
import javafx.scene.input.KeyCode;
import javafx.scene.input.KeyCodeCombination;
import javafx.scene.input.KeyCombination;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.VBox;
import javafx.scene.paint.Color;
import javafx.stage.Screen;
import javafx.stage.Stage;
/**
BEGIN_DOC

## Input

input is a tab-delimited file created by e.g: indexcov (https://github.com/brentp/goleft/tree/master/indexcov)

```
#chrom  start  end     SampleBB  SampleBC  SampleBD  SampleBE  SampleBF  SampleBG  SampleBH
chr1    23778  40778   1.59      1.31      1.67      1.61      1.83      1.52      1.48
chr1    29106  46106   1.9       1.54      1.72      1.97      1.88      1.53      1.95
chr1    84581  101581  0.764     0.841     1.2       1.16      1.18      1.13      1.23
chr1    15220  32220   0.355     0.704     1.09      0.784     0.81      1.37      0.954
chr1    58553  75553   0.353     0.436     0.912     0.836     1.16      1.09      0.611
chr1    19347  36347   0.381     0.411     0.811     0.795     1.16      1.22      0.495
chr1    81062  98062   1.09      0.972     1.35      1.22      1.66      1.76      1.1
chr1    17353  34353   1.06      1.06      1.23      1.26      1.44      1.43      1.03
chr1    48498  65498   1.08      0.996     1.28      1.44      1.52      1.57      1.05
```

## Screenshot

![https://pbs.twimg.com/media/DYbK3f1X0AECEfw.jpg:large](https://pbs.twimg.com/media/DYbK3f1X0AECEfw.jpg:large)

## See also

* https://pbs.twimg.com/media/DYbK3f1X0AECEfw.jpg
* https://github.com/brentp/goleft/tree/master/indexcov

END_DOC
 */
@Program(
		name="indexcovjfx",
		description="display indexcov data in a jfx client",
		keywords={"cnv","jfx","duplication","deletion","sv"}
		)
public class IndexCovJfx extends Application{
	private static final Logger LOG = Logger.build(IndexCovJfx.class).make();

	
	private final ObservableList<Sample> sampleNames = FXCollections.observableArrayList();
	private final ObservableList<IndexCovRow> visibleIndexCovRows = FXCollections.observableArrayList();
	private Canvas canvas = null;
	private ScrollPane canvasSrollPane = null;
	private ListView<Sample> sampleListView = null;
	private final Map<String,Color> contig2color = new HashMap<>();
	private final int CHUNK_WIDTH=100;
	private final Hershey hershey = new Hershey();
	private static final float DEFAULT_deletionTreshold = 0.4f;
	private static final float DEFAULT_duplicationTreshold = 1.9f;
	private Spinner<Double> deletionSpinner; 
	private Spinner<Double> duplicationSpinner; 
	
	
	@ParametersDelegate
	private Launcher.UsageBuider usageBuider=null;
	@Parameter(description = "Files")
	private List<String> args = new ArrayList<>();
	@Parameter(names="--testng",description = "testng",hidden=true)
	private boolean this_is_a_unit_test = false;
	
	private class Sample
		{
		final String name;
		Sample(final String name) {
			this.name = name;
			}
		
		@Override
		public String toString() {
			return name;
			}
		}
	
	private class IndexCovRow
		implements Locatable
		{
		final String contig;
		final int start;
		final int end;
		final float folds[];
		IndexCovRow(final List<String> tokens) {
			this.contig = tokens.get(0);
			this.start = Integer.parseInt(tokens.get(1));
			this.end = Integer.parseInt(tokens.get(2));
			this.folds = new float[sampleNames.size()];
			for(int i=3;i<tokens.size();i++) {
				this.folds[i-3] = Float.parseFloat(tokens.get(i));
				}
			}
		@Override
		public String getContig() { return contig; }
		@Override
		public int getStart() { return start;}
		@Override
		public int getEnd() { return end;}
		
		private String normContig( String s) {
			s=s.toLowerCase();
			if(s.startsWith("chr")) s=s.substring(3);
			return s;
 			}
		boolean hasContig(final String s) {
			return getContig().equals(s) &&
				normContig(this.getContig()).equals(normContig(s));
		}
		
		boolean overlaps(final String c,int p)
			{
			if(p<getStart() || p>getEnd()) return false;
			if(hasContig(c)) return true;
			return false;
			}
		
		}
	
	public IndexCovJfx()
	{
		this.usageBuider = new Launcher.UsageBuider(this.getClass());
	}
	
	private Color ligther(final Color c) {
		final double r=0.95;
		return new Color(c.getRed()*r,c.getGreen()*r,c.getBlue()*r,1);
	}
	
	private Color gray(double g) {
		return new Color(g,g,g,1);
	}
	@Override
	public void start(final Stage primaryStage) throws Exception {
			final Rectangle2D screen=Screen.getPrimary().getVisualBounds();

			final Pattern tab = Pattern.compile("[\t]");
			final JCommander jCommander = new JCommander(this);
			
			

			
			final List<String> unammed=this.getParameters().getUnnamed();
			jCommander.parse(unammed.toArray(new String[unammed.size()]));
			
			if(this.usageBuider.shouldPrintUsage())
				{
				this.usageBuider.usage(jCommander);
				Platform.exit();
				return;
				}
			
			BufferedReader r = null;
			try {
				if(this.args.isEmpty())
					{
					r = IOUtils.openStdinForBufferedReader();
					}
				else if(this.args.size()==1)
					{
					r = IOUtils.openFileForBufferedReading(new File(this.args.get(0)));
					}
				else
					{
					LOG.error("Illegal Number of arguments: " + this.args);
					Platform.exit();
					return;
					}
				String line = r.readLine();
				if(line==null) {
					LOG.error("Cannot read first line.");
					Platform.exit();
					return;
					}
				String tokens[] = tab.split(line);
				if(tokens.length<4 ||
					!tokens[0].equals("#chrom") ||
					!tokens[1].equals("start") ||
					!tokens[2].equals("end")) {
					LOG.error("bad first line "+line);
					Platform.exit();
					return;
					}
				this.sampleNames.addAll(Arrays.asList(tokens).
						subList(3,tokens.length).
						stream().
						map(S->new Sample(S)).
						collect(Collectors.toList())
						);
				
				this.sampleListView = new ListView<>(this.sampleNames);
				this.sampleListView.getSelectionModel().setSelectionMode(SelectionMode.MULTIPLE);
				//this.sampleListView.setPrefWidth(200);
				
				
				this.visibleIndexCovRows.addAll(r.lines().
					filter(L->!StringUtil.isBlank(L)).
					map(L->Arrays.asList(tab.split(L))).
					map(T->new IndexCovRow(T)).
					sorted((A,B)->{
						int i=  A.getContig().compareTo(B.getContig());
						if(i!=0) return i;
						return A.getStart() - B.getStart();
					}).
					collect(Collectors.toList())
					);
				
				String lastContig=null;
				for(final IndexCovRow row: this.visibleIndexCovRows) 
					{
					if(lastContig==null || !lastContig.equals(row.getContig()))
						{
						this.contig2color.put(row.getContig(),
								this.contig2color.size()%2==0?
								gray(0.96):
								gray(0.98)
								);
						lastContig=row.getContig();
						}
					}
				
		        
				
				this.canvas = new Canvas(
						screen.getWidth()-200, 
						screen.getHeight()-200
						);
				
				this.canvasSrollPane =new ScrollPane(canvas);
				
				this.canvasSrollPane.setFitToHeight(true);
				this.canvasSrollPane.setFitToWidth(true);
				this.canvasSrollPane.setHbarPolicy(ScrollBarPolicy.ALWAYS);
				this.canvasSrollPane.setHmin(0);
				// NOT HERE: see adjustScollPane();
				//this.canvasSrollPane.setHmax(this.visibleIndexCovRows.size()*CHUNK_WIDTH);
				//this.canvasSrollPane.setHvalue(0);
				
				this.canvasSrollPane.hvalueProperty().addListener(E->repaintCanvas());
				
				
				
				final VBox root = new VBox();
				final MenuBar menuBar = new MenuBar();
			    Menu menu = new Menu("Tools");
			    MenuItem item=new MenuItem("Goto");
			    item.setOnAction(AE->askGoto());
			    menu.getItems().add(item);
			    menu.getItems().add(new SeparatorMenuItem());
			    item=new MenuItem("Cleanup: remove data > DEL && data < DUP");
			    item.setOnAction(AE->askCleanup());
			    menu.getItems().add(item);
			    
			    item=new MenuItem("Cleanup: remove ALL samples <= DEL or ALL samples >= DUP");
			    item.setOnAction(AE->askEveryWhere());
			    menu.getItems().add(item);
			    menu.getItems().add(new SeparatorMenuItem());
			    
			    item=new MenuItem("Next Interesting");
			    item.setOnAction(AE->goToNextInteresting(1));
			    item.setAccelerator(new KeyCodeCombination(KeyCode.N, KeyCombination.CONTROL_DOWN));
			    menu.getItems().add(item);
			    item=new MenuItem("Previous Interesting");
			    item.setOnAction(AE->goToNextInteresting(-1));
			    item.setAccelerator(new KeyCodeCombination(KeyCode.P, KeyCombination.CONTROL_DOWN));
			    menu.getItems().add(item);

			    
			    menu.getItems().add(new SeparatorMenuItem());
			    item=new MenuItem("Quit");
			    menu.getItems().add(item);
			    item.setOnAction(AE->{Platform.exit();});
			    menuBar.getMenus().add(menu);
			    root.getChildren().add(menuBar);
			    
			    final HBox toolboxPane = new HBox();
			    Label label = new Label("DEL when \u2264 :");
			    toolboxPane.getChildren().add(label);
			    this.deletionSpinner = new Spinner<>(0.0,0.9,DEFAULT_deletionTreshold,0.05);
			    label.setLabelFor(this.deletionSpinner);
			    this.deletionSpinner.valueProperty().addListener((a,b,c)->repaintCanvas());
			    toolboxPane.getChildren().add(this.deletionSpinner );
			    
			    label = new Label("DUP when \u2265 :");
			    toolboxPane.getChildren().add(label);
			    this.duplicationSpinner = new Spinner<>(1.1,10.0,DEFAULT_duplicationTreshold,0.05);
			    this.duplicationSpinner.valueProperty().addListener((a,b,c)->repaintCanvas());
			    label.setLabelFor(this.duplicationSpinner);
			    toolboxPane.getChildren().add(this.duplicationSpinner );
			    root.getChildren().add(toolboxPane);
			    
				
				//HBox hbox = new HBox(sampleListView,this.canvasSrollPane);
				GridPane grid = new GridPane();
				grid.setHgap(10);
				grid.setVgap(10);
				grid.setPadding(new Insets(5, 5, 5, 5));
				grid.add(this.sampleListView, 0,0, 1,1);
				grid.add(this.canvasSrollPane, 1,0, 9,1);
						        
				
				//final StackPane root = new StackPane();
				root.getChildren().add(grid);
				
				final Scene scene = new Scene(root);
				primaryStage.setTitle(IndexCovJfx.class.getSimpleName()+
						" " + this.sampleNames.size()+" Sample(s) " +
						this.visibleIndexCovRows.size()+" Point(s)."
						);
				primaryStage.setOnShown(E->{
					adjustScollPane();
					canvasSrollPane.requestFocus();
					repaintCanvas();
					if(this.this_is_a_unit_test) {
						Platform.exit();
					}
				});

				canvasSrollPane.setOnKeyPressed(e -> {
				    if (e.getCode() == KeyCode.LEFT) {
				    	canvasSrollPane.setHvalue(Math.max(canvasSrollPane.getHmin(), canvasSrollPane.getHvalue()-CHUNK_WIDTH));
				    }
				    if (e.getCode() == KeyCode.RIGHT) {
				    	canvasSrollPane.setHvalue(Math.min(canvasSrollPane.getHmax(), canvasSrollPane.getHvalue()+CHUNK_WIDTH));
				    }
				});
				
				
		        primaryStage.setScene(scene);
		        primaryStage.show();
		        
		        this.sampleListView.getSelectionModel().selectedItemProperty().addListener(E->{
		        	repaintCanvas();
					});
		        this.sampleListView.getSelectionModel().selectedItemProperty().addListener(E->repaintCanvas());
				}
			catch(final Exception err) {
				LOG.error(err);
				Platform.exit();
				return;
				}
			finally {
				CloserUtil.close(r);
				}
		
			}
		
		private float getDeletionTreshold() {
			Double d=this.deletionSpinner.getValue();
			return d==null?DEFAULT_deletionTreshold:d.floatValue();
		}
		
		private float getDuplicationTreshold() {
			Double d=this.duplicationSpinner.getValue();
			return d==null?DEFAULT_duplicationTreshold:d.floatValue();
		}
	
		private  void repaintCanvas() {
			
			final int samplesIndices[];
			if(this.sampleListView.getSelectionModel().isEmpty())
				{
				samplesIndices = IntStream.range(0, this.sampleNames.size()).
						toArray();
				}
			else
				{
				samplesIndices = this.sampleListView.
						getSelectionModel().
						getSelectedIndices().stream().
						mapToInt(I->I.intValue()).
						toArray();
				}
			
			
			final GraphicsContext gc = canvas.getGraphicsContext2D();
			gc.setFill(Color.WHITESMOKE);
			gc.fillRect(0, 0, canvas.getWidth(), canvas.getHeight());
			
			double x= - this.canvasSrollPane.getHvalue();
			for(int row_idx=0;row_idx<this.visibleIndexCovRows.size();++row_idx)
				{
				if(x+CHUNK_WIDTH < 0) {
					x+=CHUNK_WIDTH;
					continue;
					}
				final IndexCovRow row = this.visibleIndexCovRows.get(row_idx);

				if(x > this.canvasSrollPane.getWidth()) break;
				
				Color bckg = contig2color.getOrDefault(row.getContig(), Color.ALICEBLUE);
				if(row_idx%2==0) bckg = ligther(bckg);
				gc.setFill(bckg);
				gc.fillRect(x, 0, CHUNK_WIDTH, canvas.getHeight());
				
				float minV = 0;
				float maxV = 2.1f;
				for(int sampleIdx:samplesIndices)
					{
					if(sampleIdx<0 || sampleIdx>=row.folds.length) continue;

					float v = row.folds[sampleIdx];
					maxV = Math.max(v,maxV);
					}
				
				
				double x_v1 = x + ((1.0-minV)/(maxV-minV)) * CHUNK_WIDTH;
				double x_v0_5 = x + ((0.5-minV)/(maxV-minV)) * CHUNK_WIDTH;
				double x_v2= x+ ((2-minV)/(maxV-minV)) * CHUNK_WIDTH;
				
				double y =0;
				gc.setStroke(Color.DARKGRAY);
				
				
				this.hershey.paint(gc, String.valueOf(row.getContig()),
						x+1, y,
						CHUNK_WIDTH-2,
						11
						);
				y+=12;
				gc.setStroke(Color.BLACK);
				this.hershey.paint(gc,
						NumberFormat.getNumberInstance(Locale.US).format(row.getStart()),
						x+1, y,
						CHUNK_WIDTH-2,
						11
						);
				y+=12;
				

				
				final double topY= y;
				final double sampleHeight = (canvas.getHeight()-topY)/samplesIndices.length;
				final float delLimit = this.getDeletionTreshold();
				final float dupLimit = this.getDuplicationTreshold();
				
				for(int sampleIdx:samplesIndices)
					{
					if(sampleIdx<0 || sampleIdx>=row.folds.length) continue;
					float v = row.folds[sampleIdx];
					double sample_x = x + ((v-minV)/(maxV-minV)) * CHUNK_WIDTH;
					Rectangle2D rect; 
					if(sample_x< x_v1) {
						rect = new Rectangle2D(
							sample_x, y,
							x_v1-sample_x,
							sampleHeight
							);
						if(v <= delLimit ) {
							gc.setStroke(Color.DARKGRAY);
							this.hershey.paint(gc,
									this.sampleNames.get(sampleIdx).name , 
									x_v1,
									y,
									CHUNK_WIDTH-(x_v1-x),
									Math.min(10,sampleHeight)
									);
							gc.setFill(Color.BLUE);
							}
						else
							{
							gc.setFill(Color.LIGHTGREY);
							}
						}
					else
						{
						rect = new Rectangle2D(
								x_v1, y,
								sample_x-x_v1,
								sampleHeight
								);
						if(v >= dupLimit)
							{
							gc.setStroke(Color.DARKGRAY);
							this.hershey.paint(gc,
									this.sampleNames.get(sampleIdx).name , 
									x,
									y,
									(x_v1-x),
									Math.min(10,sampleHeight)
									);
							gc.setFill(Color.RED);
							}
						else
							{
							gc.setFill(Color.LIGHTGREY);
							}
						
						}
					gc.fillRect(rect.getMinX(),rect.getMinY(),rect.getWidth(),rect.getHeight());
					gc.setStroke(Color.DARKGRAY);
					gc.strokeRect(rect.getMinX(),rect.getMinY(),rect.getWidth(),rect.getHeight());
					y+= sampleHeight;
					}
				
				gc.setStroke(Color.ORANGE);
				gc.strokeLine(x_v0_5, topY, x_v0_5, canvas.getHeight());
				gc.setStroke(Color.MAGENTA);
				gc.strokeLine(x_v2, topY, x_v2, canvas.getHeight());
				
				gc.setStroke(Color.DARKGRAY);
				gc.strokeRect(x, 0, CHUNK_WIDTH, canvas.getHeight());
				x+=CHUNK_WIDTH;
				}
			
			}
		/** remove data > DEL && data < DUP */
		private void askCleanup() {
			final float delLimit = this.getDeletionTreshold();
			final float dupLimit = this.getDuplicationTreshold();
			final Alert alert = new Alert(AlertType.CONFIRMATION);
			alert.setTitle("Confirmation Dialog");
			alert.setHeaderText("Keep with tresholds DEL:" + 
					delLimit +"  DUP: "+ dupLimit
					);
			alert.setContentText("This will remove some data. Are you ok with this?");

			final Optional<ButtonType> result = alert.showAndWait();
			if (result.get() != ButtonType.OK) return;
			this.visibleIndexCovRows.removeIf(R->{
				for(float v: R.folds)
					{
					if(v<=delLimit) return false;
					if(v>=dupLimit) return false;
					}
				return true;
				});
			adjustScollPane();
			}
		
		/** remove data > DEL && data < DUP for all samples*/
		private void askEveryWhere() {
			final float delLimit = this.getDeletionTreshold();
			final float dupLimit = this.getDuplicationTreshold();

			final Alert alert = new Alert(AlertType.CONFIRMATION);
			alert.setTitle("Confirmation Dialog");
			alert.setHeaderText(
					"Keep with tresholds DEL:" +
					delLimit+"  DUP: " + 
					dupLimit
					);
			alert.setContentText("This will remove some data. Are you ok with this?");

			final Optional<ButtonType> result = alert.showAndWait();
			if (result.get() != ButtonType.OK) return;
			this.visibleIndexCovRows.removeIf(R->{
				
				int count=0;
				for(final float v: R.folds)
					{
					if(v <= delLimit) 
						{
						count++;
						}
					}
				if(count==R.folds.length) return true;
				count=0;
				for(final float v: R.folds)
					{
					if(v >= dupLimit)
						{
						count++;
						}
					}
				if(count==R.folds.length) return true;
				return false;
				});
			adjustScollPane();
			}
		
		private void adjustScollPane() {
			this.canvasSrollPane.setHmax(this.visibleIndexCovRows.size()*CHUNK_WIDTH);
			this.canvasSrollPane.setHvalue(0);
			repaintCanvas();
		}
		
		private void askGoto() {
			TextInputDialog dialog = new TextInputDialog("Enter position");
			dialog.setTitle("Go To...");
			dialog.setHeaderText("Enter a position");
			dialog.setContentText("chrom:pos");
			// Traditional way to get the response value.
			Optional<String> result = dialog.showAndWait();
			if (!result.isPresent()){
				LOG.error("No input");
				return;
				}
			int colon = result.get().indexOf(':');
			if(colon<=0) return;
			final String contig =  result.get().substring(0, colon).trim();
			final int pos;
			try {
				pos = Integer.parseInt( result.get().
						substring(colon+1).
						replaceAll("[, ]", "").
						trim());
				}
			catch(Exception err) {
				LOG.error(err);
				return;
			 	}
			int x=0;
			while(x<this.visibleIndexCovRows.size())
				{
				final IndexCovRow row = this.visibleIndexCovRows.get(x);
				if(row.overlaps(contig,pos) ||
						row.hasContig(contig) && row.getStart()> pos) {
					this.canvasSrollPane.setHvalue(x*CHUNK_WIDTH);
					repaintCanvas();
					return;
					}
				++x;
				}
			LOG.error("Not found "+result.get());
			}
		
		
		private void goToNextInteresting(int direction) {
			final float delLimit = this.getDeletionTreshold();
			final float dupLimit = this.getDuplicationTreshold();

			int x = (int)(this.canvasSrollPane.getHvalue()/(double)CHUNK_WIDTH);
			for(;;)
				{
				x+=direction;
				if(x<0 || x>=this.visibleIndexCovRows.size()) break;
				final IndexCovRow row = this.visibleIndexCovRows.get(x);
				
				for(float v:row.folds)
					{
					if(v <= delLimit || v >= dupLimit)
						{
						this.canvasSrollPane.setHvalue(x*CHUNK_WIDTH);
						repaintCanvas();
						return;
						}
					}
				
				}
		}
		
		public static void main(String[] args) {
			Application.launch(args);
			}

}
