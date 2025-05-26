/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.regenie;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Desktop;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Graphics2D;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.IOException;
import java.net.URI;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Vector;
import java.util.function.Predicate;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import javax.swing.AbstractAction;
import javax.swing.BoxLayout;
import javax.swing.DefaultComboBoxModel;
import javax.swing.ImageIcon;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.ListSelectionModel;
import javax.swing.SwingUtilities;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.dict.SequenceDictionaryExtractor;
import com.github.lindenb.jvarkit.io.FileHeader;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.iterator.AbstractCloseableIterator;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.primitive.FloatArray;
import com.github.lindenb.jvarkit.locatable.SimpleInterval;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.math.MathUtils;
import com.github.lindenb.jvarkit.net.UrlSupplier;
import com.github.lindenb.jvarkit.swing.AbstractGenericTableModel;
import com.github.lindenb.jvarkit.util.Algorithm;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;

/**

BEGIN_DOC




END_DOC
**/
@Program(
		name="swingregenie",
		description="view regenie output",
		creationDate="20250324",
		modificationDate="20250324",
		keywords= {"regenie","gwas","burden"},
		jvarkit_amalgamion = true
		)
public class RegenieSwing extends Launcher {
	private static final Logger LOG = Logger.of(RegenieSwing.class);
	private final static String[] COLUMN_HEADER=new String[] {"CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","N","TEST","BETA","SE","CHISQ","LOG10P","EXTRA"};	
	private static final int MANHATTAN_WIDTH = 900;
	private static final int MANHATTAN_HEIGHT = 200;
	private static final int QQPLOT_SIZE = MANHATTAN_HEIGHT;
	@Parameter(names={"-R","--reference"},description= DICTIONARY_SOURCE,required = true)
	private Path dictPath = null;

	private Path regeniePath;
	private SAMSequenceDictionary main_dictionary = null;
	
	private  static class ChromInfo {
		SAMSequenceRecord ssr=null;
		int tid;
		long start = 0;
		}
	
	private static class BrowseURL {
		final UrlSupplier.LabelledUrl url;
		BrowseURL(UrlSupplier.LabelledUrl url) {
			this.url=url;
			}
		@Override
		public String toString() {
			String s= url.getLabel()+" "+url.getUrl();
			if(s.length()>50) s=s.substring(0,47)+"...";
			return s;
			}
		}
	
	@SuppressWarnings("serial")
	private class RecordTableModel extends AbstractGenericTableModel<Record> {
		@Override
		public int getColumnCount() {
			return COLUMN_HEADER.length;
			}
		@Override
		public String getColumnName(int column) {
			return COLUMN_HEADER[column];
			}
		@Override
		public Object getValueOf(Record o, int columnIndex) {
			return o.row.get(COLUMN_HEADER[columnIndex]);
			}
		}
	
	private class Record implements  Locatable {
		FileHeader.RowMap row;
		int pos;
		SAMSequenceRecord ssr;
		
		Record(final FileHeader.RowMap row) {
			this.row = row;
		}
		public String getFreq() {
			return split()[2];
		}
		
		public String getMask() {
			return split()[1];
		}
		public String getGene() {
			return split()[0];
		}
		public String getTest() {
			return row.get("TEST");
		}
		public double getLog10P() {
			return Double.parseDouble(row.get("LOG10P"));
			}
		private String[] split() {
			final String[] array = new String[3];
			String id = row.get("ID");
			if(id.endsWith(".singleton")) {
				array[2] = "singleton";
				id = id.substring(0,id.length()-10);
				}
			else
				{
				int zero= id.lastIndexOf(".0.");
				if(zero==-1) throw new IllegalArgumentException();
				array[2] = id.substring(zero+1);
				id = id.substring(0,zero);
				}
			int dot= id.lastIndexOf(".");
			array[1] = id.substring(dot+1);
			id=id.substring(0,dot);
			array[0]=id;
			return array;
			}
		
		
		@Override
		public String getContig() {
			return this.ssr.getSequenceName();
			}
		public int getPosition() {
			return pos;
			}
		@Override
		public int getStart() {
			return getPosition();
			}
		@Override
		public int getEnd() {
			return getPosition();
			}
		}
	private static class NamedFilter implements Predicate<Record> {
		final String name;
		final Predicate<Record> pred;
		NamedFilter(final String name,final Predicate<Record> pred) {
			this.name = name;
			this.pred = pred;
			}
		@Override
		public boolean test(Record t) {
			return pred.test(t);
			}
		@Override
		public String toString() {
			return name;
			}
		}
	
	private static class NamedComparator implements Comparator<Record> {
		final String name;
		final Comparator<Record> cmp;
		NamedComparator(final String name,final Comparator<Record> cmp) {
			this.name = name;
			this.cmp = cmp;
			}
		@Override
		public int compare(Record a,Record b) {
			return cmp.compare(a,b);
			}
		@Override
		public String toString() {
			return name;
			}
		}
	
	private class RecordIterator extends AbstractCloseableIterator<Record> {
		final BufferedReader br;
		private final FileHeader fileHeader;
		final int column_chrom;
		final int column_pos;
		final int column_extra;
		final int column_log10P;
		private final Map<String,SAMSequenceRecord> aliases=new HashMap<>();
		RecordIterator() {
			try {
				for(int i=0;i< RegenieSwing.this.main_dictionary.size();i++) {
					final SAMSequenceRecord ssr = RegenieSwing.this.main_dictionary.getSequence(i);
					final String chrom = ssr.getSequenceName();
					if(!chrom.matches("(chr)?[0-9X]+")) continue;//Y not in regenie
					aliases.put(chrom, ssr);
					if(chrom.startsWith("chr")) {
						aliases.put(chrom.substring(3), ssr);
						}
					else {
						aliases.put("chr"+chrom, ssr);
						}
					if((chrom.equals("chrX") || chrom.equals("X")) && RegenieSwing.this.main_dictionary.getSequence("23")==null) {
						aliases.put("23", ssr);
						}
					}
				if(aliases.isEmpty()) throw new IllegalStateException("no aliases");
				
				this.br = IOUtils.openPathForBufferedReading(RegenieSwing.this.regeniePath);
				for(;;) {
					final String line = br.readLine();
					if(line!=null && line.startsWith("#")) continue;
					if(line==null) throw new IOException("cannot read first line of "+regeniePath);
					final Pattern ws  = Pattern.compile("[ \t]");
					this.fileHeader = new FileHeader(line, S->Arrays.asList(ws.split(S)));
					this.column_chrom = this.fileHeader.getColumnIndex("CHROM");
					this.column_pos = this.fileHeader.getColumnIndex("GENPOS");
					this.column_extra = this.fileHeader.getColumnIndex("EXTRA");
					this.column_log10P =  this.fileHeader.getColumnIndex("LOG10P");
					break;
					}
				}
			catch(IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		@Override
		protected Record advance() {
			try {
				for(;;) {
					final String line = br.readLine();
					if(line==null) return null;
					if(line.startsWith("#")|| line.startsWith(COLUMN_HEADER[0])) continue;
					final Record rec= new Record(this.fileHeader.toMap(line));
					if(rec.row.at(this.column_extra).equals("TEST_FAIL")) continue;
					if(rec.row.at(this.column_log10P).equals("NA")) continue;
					rec.pos = Integer.parseInt(rec.row.at(this.column_pos));
					rec.ssr = this.aliases.get(rec.row.at(this.column_chrom));
					if(rec.ssr==null) {
						continue;
						}
					return rec;
					}
				}
			catch(final IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		@Override
		public void close() {
			try { this.br.close();}
			catch(IOException err) {
				LOG.error(err);
				}
			}
		}
	
	@SuppressWarnings("serial")
	private class SwingFrame extends JFrame{
		
	
		private class ScanThread extends Thread {
			volatile boolean stop_flag = false;
			final FloatArray x_array = new FloatArray(SwingFrame.this.number_of_records);
			final FloatArray y_array = new FloatArray(SwingFrame.this.number_of_records);
			boolean hide_xy=false;
			boolean show_qqplot=false;
			final List<Record> best_records = new ArrayList<>();
			List<Predicate<Record>> predicates = new ArrayList<>();
			final Map<String,ChromInfo> dictionary = new LinkedHashMap<>(RegenieSwing.this.main_dictionary.size());
			long genome_size=0L;
			ScanThread() {
				}
			@Override
			public void run() {
				int count=0;
				for(int i=0;i< RegenieSwing.this.main_dictionary.size();++i) {
					final SAMSequenceRecord ssr = RegenieSwing.this.main_dictionary.getSequence(i);
					final String chrom = ssr.getSequenceName();
					if(!chrom.matches("(chr)?[0-9X]+")) continue;//Y not in regenie
					if(hide_xy && this.isXY(chrom)) continue;
					final ChromInfo ci = new ChromInfo();
					ci.ssr=ssr;
					ci.tid = this.dictionary.size();
					this.dictionary.put(chrom, ci);
					}
				
				long x1=0;
				for(ChromInfo ci:this.dictionary.values()) {
					ci.start =  x1;
					x1+=ci.ssr.getLengthOnReference();
					}
				
				this.genome_size = dictionary.values().stream().mapToLong(S->S.ssr.getLengthOnReference()).sum();
				if(this.genome_size<=0L || dictionary.isEmpty()) {
					System.err.println("NO GENOME SIZE ??");
					return;
					}
				final Algorithm<Record, Double> sort_on_log10 = new Algorithm<>(R->R.getLog10P());
					try(RecordIterator iter=new RecordIterator()) {
						while(iter.hasNext()) {
							if(stop_flag) break;
							final Record rec = iter.next();
							final ChromInfo chromInfo = dictionary.get(rec.ssr.getSequenceName());
							if(chromInfo==null) {
								continue;
								}
							if(predicates.stream().anyMatch(P->!P.test(rec))) continue;
							count++;
							final double log10p = rec.getLog10P();
							if(best_records.isEmpty() || log10p >= best_records.get(0).getLog10P())
								{
								int i = sort_on_log10.lower_bound(best_records, log10p);
								if(i<=best_records.size()) {
									best_records.add(i,rec);
									while(best_records.size()>1000) {
										best_records.remove(0);
										}
									}
								}
						
							final double x  = toPixelX(chromInfo.start+rec.getPosition());
							final double y = log10p;
							x_array.add((float)x);
							y_array.add((float)y);
							if(count%5_000_000==0) {
								updateGUI();
								}
							}
						updateGUI();
						SwingUtilities.invokeLater(()->{
							if(!stop_flag && ScanThread.this.equals(SwingFrame.this.thread)) {
								progressBar.setIndeterminate(false);
								}
						});
						}
					catch(final Throwable err) {
						err.printStackTrace();
						return ;
						}
					
					}
			
				private boolean isXY(String chrom) {
					if(chrom.equals("X") || chrom.equals("chrX") || chrom.equals("23")) return true;
					if(chrom.equals("Y") || chrom.equals("chrY") || chrom.equals("24")) return true;
					return false;
					}
				
				private double toPixelX(long genome_pos) {
					return ((genome_pos/(double)genome_size)*MANHATTAN_WIDTH);
					}
				
				
				private void updateGUI() {
					if(stop_flag)return;
					try {
					SwingUtilities.invokeAndWait(()->{
						System.err.println("Update GUI "+x_array.size());
						if(!stop_flag && ScanThread.this.equals(SwingFrame.this.thread) && !y_array.isEmpty() ) {
						BufferedImage img = new BufferedImage(MANHATTAN_WIDTH, MANHATTAN_HEIGHT, BufferedImage.TYPE_INT_ARGB);
						
						Graphics2D g=img.createGraphics();
						
						g.setColor(Color.WHITE);
						g.fillRect(0, 0,MANHATTAN_WIDTH,MANHATTAN_HEIGHT);
						for(ChromInfo ci: this.dictionary.values()) {
							g.setColor(ci.tid%2==0?Color.WHITE:Color.LIGHT_GRAY);
							g.fill(new Rectangle2D.Double(toPixelX(ci.start), 0, toPixelX(ci.ssr.getLengthOnReference()),MANHATTAN_HEIGHT));
							}
						for(ChromInfo ci: this.dictionary.values()) {
							g.setColor(Color.DARK_GRAY);
							g.fillRect((int)toPixelX(ci.start), 0,1,MANHATTAN_HEIGHT);
							}
						double max_y = 1+y_array.stream().mapToDouble(f->f).max().orElse(1.0);
						
						for(int i=6;i<=(int)(max_y+1);++i) {
							g.setColor(Color.GREEN);
							final double y = MANHATTAN_HEIGHT - (i/max_y)*MANHATTAN_HEIGHT;
							g.fillRect(0,(int) y,MANHATTAN_WIDTH,1);
							}
						
						for(int i=0;i< x_array.size();++i) {
							if(stop_flag) return;
							double x = x_array.get(i);
							double y = MANHATTAN_HEIGHT - (y_array.get(i)/max_y)*MANHATTAN_HEIGHT;
							int r=3;
							g.setColor(Color.BLUE);
							g.fillOval((int)(x-r/2.0),(int)(y-r/2.0),r, r);
							}
						
						for(ChromInfo ci: dictionary.values()) {
							g.setColor(Color.DARK_GRAY);
							g.fillRect((int)toPixelX(ci.start), 0,1,MANHATTAN_HEIGHT);
							}
						
						g.setColor(Color.DARK_GRAY);
						g.drawRect(0,0,MANHATTAN_WIDTH-1,MANHATTAN_HEIGHT-1);
						g.dispose();
						manhattan_label.setIcon(new ImageIcon(img));

						
						img = new BufferedImage(QQPLOT_SIZE, QQPLOT_SIZE, BufferedImage.TYPE_INT_ARGB);
						g=img.createGraphics();
						g.setColor(Color.WHITE);
						g.fillRect(0, 0,QQPLOT_SIZE,QQPLOT_SIZE);
						
						if(this.show_qqplot) {
							max_y = y_array.stream().mapToDouble(f->f).max().orElse(1.0);
							g.setColor(Color.GREEN);
							g.drawLine(0,QQPLOT_SIZE,  QQPLOT_SIZE,0);
							final float[] x2_array = this.y_array.toArray();
							if(x2_array.length>0 && !stop_flag) {
								Arrays.sort(x2_array);
								double[] y2_array = MathUtils.ppoints(x2_array.length);
								
								for(int i=0;i<y2_array.length && !stop_flag;++i ) {
									y2_array[i] = -Math.log10(y2_array[i]);
									}
								Arrays.sort(y2_array);
								final double min_y2=y2_array[0];
								final double max_y2=y2_array[y2_array.length-1];
								
								
								for(int i=0;i<x2_array.length && !stop_flag;++i ) {
									double cx = (x2_array[i]/max_y)*QQPLOT_SIZE;
									double cy = QQPLOT_SIZE - ((y2_array[i]-min_y2)/(max_y2-min_y2))*QQPLOT_SIZE;
									g.setColor(Color.BLUE);
									int r=3;
									g.fillOval((int)(cx-r/2.0),(int)(cy-r/2.0),r, r);
									}
								}
							}
						g.dispose();
						qqplot_label.setIcon(new ImageIcon(img));
						//System.err.println("Update GUI Done "+x_array.size());
						}
					});
					} catch(Throwable err ) {
						err.printStackTrace();
						}
					SwingUtilities.invokeLater(()->{
						if(!stop_flag && ScanThread.this.equals(SwingFrame.this.thread) ) {
							List<Record> copy = new ArrayList<>(best_records);
							Collections.reverse(copy);
							recordTableModel.setRows(copy);
							}
						});
					}
				}
		private final RecordTableModel recordTableModel;
		private final JTable bestRecordTable;
		private final JLabel manhattan_label;
		private final JLabel qqplot_label;
		private ScanThread thread=null;
		private JProgressBar progressBar;
		private JComboBox<NamedFilter> jcomboFreq;
		private JComboBox<NamedFilter> jcomboTest;
		private JComboBox<NamedFilter> jcomboMask;
		private JComboBox<NamedComparator> jcomboSortRecord;
		private final JCheckBox jCheckhideSexualChrom;
		private final JCheckBox jCheckShowQQPlot;
		private final JComboBox<BrowseURL> jcomboURLS;
		private int number_of_records = 0;
		public SwingFrame() {
			this.setTitle(regeniePath.getFileName().toString());
			this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
			final Set<String> frequency_set = new HashSet<>();
			final Set<String> test_set = new HashSet<>();
			final Set<String> mask_set = new HashSet<>();
			try(RecordIterator iter=new RecordIterator()) {
				while(iter.hasNext()) {
					final Record rec=iter.next();
					number_of_records++;
					frequency_set.add(rec.getFreq());
					test_set.add(rec.getTest());
					mask_set.add(rec.getMask());
				}
			}
			final JMenuBar menuBar = new JMenuBar();
			setJMenuBar(menuBar);
			JMenu menu = new JMenu("File");
			menuBar.add(menu);
			menu.add(new JMenuItem(new AbstractAction("Quit") {
				@Override
				public void actionPerformed(ActionEvent e) {
					SwingFrame.this.setVisible(false);
					SwingFrame.this.dispose();
					if(thread!=null) {
						thread.stop_flag=true;
						thread.interrupt();
					}
				}
			}));
			
			final JPanel contentPane = new JPanel(new BorderLayout(5,5));
			this.setContentPane(contentPane);
			final JPanel mainPane = new JPanel();
			mainPane.setLayout(new BoxLayout(mainPane,BoxLayout.Y_AXIS));
			contentPane.add(mainPane,BorderLayout.CENTER);
			final JPanel top=new JPanel(new FlowLayout());
			contentPane.add(top,BorderLayout.NORTH);
			
			final JPanel bottom=new JPanel(new FlowLayout());
			contentPane.add(bottom,BorderLayout.SOUTH);
			this.progressBar = new JProgressBar();
			bottom.add(progressBar);
			
			
			top.add(this.jCheckhideSexualChrom=new JCheckBox("Hide X/Y"));
			this.jCheckhideSexualChrom.setSelected(true);
			this.jCheckhideSexualChrom.addActionListener(AE->reloadData());
			
			top.add(this.jCheckShowQQPlot =new JCheckBox("Show QQPlot"));
			this.jCheckShowQQPlot.setSelected(false);
			this.jCheckShowQQPlot.addActionListener(AE->reloadData());
			
			Vector<NamedFilter> filters = new Vector<>();
			
			// frequency
			filters.add(new NamedFilter("*",R->true));
			for(String f: frequency_set) {
				filters.add(new NamedFilter(f,R->R.getFreq().equals(f)));
				}
			this.jcomboFreq = new JComboBox<>(filters);
			top.add(new JLabel("Freq:"));
			top.add(this.jcomboFreq);
			this.jcomboFreq.addActionListener(AE->reloadData());
			
			// test
			 filters = new Vector<>();
			filters.add(new NamedFilter("*",R->true));
			for(String f: test_set) {
				filters.add(new NamedFilter(f,R->R.getTest().equals(f)));
				}
			this.jcomboTest = new JComboBox<>(filters);
			top.add(new JLabel("Test:"));
			top.add(this.jcomboTest);
			this.jcomboTest.addActionListener(AE->reloadData());
			
			// mask
			filters = new Vector<>();
			filters.add(new NamedFilter("*",R->true));
			for(String f: mask_set) {
				filters.add(new NamedFilter(f,R->R.getMask().equals(f)));
				}
			this.jcomboMask = new JComboBox<>(filters);
			top.add(new JLabel("Mask:"));
			top.add(this.jcomboMask);
			this.jcomboMask.addActionListener(AE->reloadData());
			
			JPanel paneChart = new JPanel(new BorderLayout());
			ImageIcon icn=new ImageIcon();
		
			this.manhattan_label= new JLabel(icn);
			this.manhattan_label.setPreferredSize(new Dimension(MANHATTAN_WIDTH,MANHATTAN_HEIGHT));
			paneChart.add(this.manhattan_label,BorderLayout.CENTER);
			icn=new ImageIcon();
			this.qqplot_label= new JLabel(icn);
			this.qqplot_label.setPreferredSize(new Dimension(QQPLOT_SIZE,QQPLOT_SIZE));
			paneChart.add(this.qqplot_label,BorderLayout.EAST);
			
			mainPane.add(paneChart);
			
			final JPanel tablePane = new JPanel(new BorderLayout(5, 5));
			final JPanel top2= new JPanel(new FlowLayout());
			tablePane.add(top2,BorderLayout.NORTH);
			
			top2.add(new JLabel("Sort:"));
			Vector<NamedComparator> recordSorters=new Vector<>(Arrays.asList(
					new NamedComparator("Log10", (A,B)->Double.compare(B.getLog10P(),A.getLog10P())),
					new NamedComparator("Gene", (A,B)->A.getGene().compareTo(B.getGene())),
					new NamedComparator("Position", (A,B)->{
						int i = Integer.compare(A.ssr.getSequenceIndex(),B.ssr.getSequenceIndex());
						if(i!=0) return i;
						return Integer.compare(A.pos,B.pos);
						})));
			top2.add(jcomboSortRecord=new  JComboBox<>(recordSorters));
			
			
			top2.add(new JLabel("URL:"));
			top2.add(jcomboURLS = new JComboBox<>(new DefaultComboBoxModel<BrowseURL>()));
			jcomboURLS.addActionListener(AE->{
				BrowseURL url = (BrowseURL) jcomboURLS.getSelectedItem();
				if(url==null) return;
				 try {
                     Desktop.getDesktop().browse(new URI(url.url.getUrl()));
                     }
                 catch(final Throwable err) {
                     JOptionPane.showInputDialog(SwingFrame.this,"URL",url.url.getUrl());
                     }
				});
			
			mainPane.add(tablePane);
			this.recordTableModel = new RecordTableModel();
			this.bestRecordTable = new JTable(this.recordTableModel);
			this.bestRecordTable.getSelectionModel().setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
			tablePane.add(new JScrollPane(this.bestRecordTable),BorderLayout.CENTER);
			this.bestRecordTable.getSelectionModel().addListSelectionListener(EVT->{
				DefaultComboBoxModel<BrowseURL> model = (DefaultComboBoxModel<BrowseURL>)(jcomboURLS.getModel());
				model.removeAllElements();
				int r =bestRecordTable.getSelectedRow();
				if(r<0) return;
				r = bestRecordTable.convertRowIndexToModel(r);
				if(r<0) return;
				
				Record rec=this.recordTableModel.getElementAt(r);
				
				if(rec==null) return;
				int extend=500;
				final Locatable loc = new SimpleInterval(rec.getContig(),Math.max(1, rec.getStart()-extend),rec.getEnd()+extend);
				final UrlSupplier urlSupplier = new UrlSupplier(main_dictionary);
				final Set<UrlSupplier.LabelledUrl> unique_urls = new HashSet<>();
				unique_urls.addAll(urlSupplier.of(loc));
				if(rec.getGene().matches("ENS[TG][0-9\\.]+") || rec.getGene().matches("[A-Z][A-Z0-9\\.]*")) {
					unique_urls.addAll(urlSupplier.of(rec.getGene()));
					}
				
				for(BrowseURL bu:unique_urls.stream().map(X->new BrowseURL(X)).collect(Collectors.toList())) {
					model.addElement(bu);
					}
				});
			
			
			jcomboSortRecord.addActionListener(AE->{
				NamedComparator cmp = (NamedComparator)jcomboSortRecord.getSelectedItem();
				if(cmp==null) return;
				final List<Record> rows = new ArrayList<>(this.recordTableModel.getRows());
				Collections.sort(rows,cmp);
				this.recordTableModel.setRows(rows);
				});
			
			this.addWindowListener(new WindowAdapter() {
				@Override
				public void windowOpened(WindowEvent e) {
					reloadData();
					}
				});
			}
		
		
		private synchronized void reloadData() {
			if(thread!=null) {
				thread.stop_flag=true;
				this.progressBar.setIndeterminate(false);
				}
			this.thread=new ScanThread();
			NamedFilter flt = (NamedFilter)this.jcomboFreq.getSelectedItem();
			if(flt!=null) this.thread.predicates.add(flt);
			flt = (NamedFilter)this.jcomboMask.getSelectedItem();
			if(flt!=null) this.thread.predicates.add(flt);
			flt = (NamedFilter)this.jcomboTest.getSelectedItem();
			if(flt!=null) this.thread.predicates.add(flt);
			this.thread.hide_xy = this.jCheckhideSexualChrom.isSelected();
			this.thread.show_qqplot = this.jCheckShowQQPlot.isSelected();
			this.progressBar.setIndeterminate(true);
			this.thread.start();
			}
		}
	
	@Override
	public int doWork(final List<String> args) {
		
		this.main_dictionary = new SequenceDictionaryExtractor().extractRequiredDictionary(this.dictPath);
		this.regeniePath = Paths.get(oneAndOnlyOneFile(args));
		IOUtil.assertFileIsReadable(this.regeniePath);
		SwingUtilities.invokeLater(()->{
			JFrame.setDefaultLookAndFeelDecorated(true);
			final SwingFrame f = new SwingFrame();
			final Dimension d = Toolkit.getDefaultToolkit().getScreenSize();
			f.setBounds(50, 50, d.width-100	, d.height-100);
			f.setVisible(true);
			});
		return 0;
		}
	public static void main(final String[] args) {
		new RegenieSwing().instanceMain(args);
		}

}
