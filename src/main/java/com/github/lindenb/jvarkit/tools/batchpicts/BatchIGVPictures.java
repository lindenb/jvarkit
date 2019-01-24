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
package com.github.lindenb.jvarkit.tools.batchpicts;

import htsjdk.samtools.util.CloserUtil;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.GridLayout;
import java.awt.RenderingHints;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.io.Reader;
import java.io.StringReader;
import java.io.StringWriter;
import java.net.URL;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.prefs.Preferences;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import javax.imageio.ImageIO;
import javax.swing.AbstractAction;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;
import javax.swing.border.TitledBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.text.BadLocationException;
import javax.xml.namespace.QName;
import javax.xml.stream.XMLEventFactory;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLEventWriter;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamWriter;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.Namespace;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.util.igv.IgvSocket;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

/*
BEGIN_DOC
## Motivation

Takes IGV pictures in batch. Save as HTML+png images

## Screenshot

![screenshot](http://i.imgur.com/pasROkt.jpg)

## History

2019: removed open office support

END_DOC
*/

@SuppressWarnings("serial")
@Program(name="batchigvpictures",
		description="Takes IGV pictures in batch. Save as HTML+png image",
		keywords={"gui","igv","visualization"}
		)
class BatchIGVPicturesFrame extends JFrame
	{
	private static final String PREF_EXTEND="snp.extend";
	private static final String PREF_IGV_PORT="igv.port";
	private static final String PREF_IGV_HOST="igv.host";
	private static final String PREF_SAVE_AS_FILE="save.as.file";
	private static final String PREF_IMG_WIDTH="img.width";
	private static final String PREF_LOAD_VCF_FILE="load.vcf.file";
	private static final String PREF_IGV_SLEEP="igv.sleep";
	private static final String PREF_FORMAT="out.fmt";
	private static final int FORMAT_HTML=0;
	
	String about="";
	private Preferences preferences;
	private AbstractSlave currentSlave=null;
	private JTextArea snpArea;
	
	private SpinnerNumberModel spinnerExtendSize;
	private SpinnerNumberModel spinnerImageWidth;
	private SpinnerNumberModel spinnerIgvWait;
	private JTextField tfIgvHost;
	private SpinnerNumberModel spinnerIgvPort;
	private JLabel progressLbl=null;
	private JComboBox<String> cboxFormat;
	
	
	private static class Mutation
		implements Comparable<Mutation>
		{
		String chrom;
		int position;
		@Override
		public int compareTo(Mutation o) {
			int i=chrom.compareTo(o.chrom);
			if(i!=0) return i;
			return position-o.position;
			}
		@Override
		public boolean equals(Object obj) {
			if(this==obj) return true;
			return compareTo(Mutation.class.cast(obj))==0;
			}
		@Override
		public int hashCode() {
			return chrom.hashCode()*31+position;
			}
		@Override
		public String toString() {
			return chrom+":"+position;
			}
		}
	
	
	
	private abstract class AbstractSlave
		extends Thread
		{
		String pathOut;
		int igvPort;
		String igvHost;
		int extend;
		int imageWidth;
		int waitMiliSecs;
		Set<Mutation> mutations=new TreeSet<Mutation>();
		
		protected void log(final String labelStr)
			{
			try {
				SwingUtilities.invokeAndWait(new Runnable()
					{
					@Override
					public void run() {
						progressLbl.setText(labelStr);
					}
				});
				} 
			catch (Exception e) {
				
				}
			}
		protected IgvSocket createIgvSocket()
			{
			IgvSocket igvSocket=new IgvSocket();
			igvSocket.setHost(this.igvHost);
			igvSocket.setPort(this.igvPort);
			return igvSocket;
			}
	
		protected BufferedImage scale(File pngFile)
			throws IOException
			{
			BufferedImage src=ImageIO.read(pngFile);
			if(src==null)
				{
				throw new IOException("Cannot read Image "+pngFile.getName());
				}
			
			
			BufferedImage dest=new BufferedImage(
					this.imageWidth,
					(int)((this.imageWidth/(double)src.getWidth())*src.getHeight()),
					src.getType()
					);
			Graphics2D g=(Graphics2D)dest.getGraphics();
			g.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
			g.drawImage(src, 0, 0, dest.getWidth(), dest.getHeight(), null);
			g.dispose();
			return dest;
			}

		}
	
	
	private class HtmlSlave
		extends AbstractSlave
		{
		@Override
		public void run()
			{
			File baseFile=new File(pathOut);
			File resultFileOut=new File(pathOut+".html");
				
			File pngFile=null;
			PrintWriter out =null;
			BufferedReader in =null;
			XMLStreamWriter w=null;
			FileOutputStream fout=null;
			IgvSocket igvSocket=null;
			PrintWriter igvLog=null;
			try {
				XMLOutputFactory xof=XMLOutputFactory.newFactory();
				igvLog=new PrintWriter(new NullOuputStream());
				fout=new FileOutputStream(resultFileOut);
				w=xof.createXMLStreamWriter(fout, "UTF-8");
				w.writeStartElement("html");
				w.writeStartElement("body");
				pngFile=File.createTempFile("igv_", ".png");
				igvSocket=createIgvSocket();
				out=igvSocket.getWriter();
				in=igvSocket.getReader();
				
				out.println("snapshotDirectory "+pngFile.getParentFile());				 
				igvLog.println(in.readLine());
				out.println("setSleepInterval "+(waitMiliSecs+1));				 
				igvLog.println(in.readLine());

				int index=0;
				for(Mutation m:this.mutations)
					{
					if(this!=BatchIGVPicturesFrame.this.currentSlave) break;
					log("Drawing "+(++index)+"/"+this.mutations.size()+" "+m);
				
					int chromStart= Math.max(1,m.position - extend);
					int chromEnd =  m.position + extend;
					out.println("goto "+m.chrom+":"+chromStart+"-"+chromEnd);
					out.flush();
					
					igvLog.println(in.readLine()); 
					out.println("snapshot "+pngFile.getName());
					out.flush();
					igvLog.println(in.readLine());
					
					BufferedImage dest=scale(pngFile);
					String imgName=baseFile.getName()+"_"+m.chrom+"_"+m.position+".png";
					
					w.writeStartElement("h3");
					w.writeCharacters(m.toString());
					w.writeEndElement();
					w.writeEmptyElement("img");
					w.writeAttribute("src",imgName);
					w.writeAttribute("alt",m.toString());
					w.writeAttribute("width",String.valueOf(dest.getWidth()));
					w.writeAttribute("height",String.valueOf(dest.getHeight()));
					
					
					
					File destScaledImage=new File(baseFile.getParentFile(),imgName);
					ImageIO.write(dest, "png",destScaledImage);
					
					dest=null;
					log("Sleep "+waitMiliSecs+" millisecs...");
					Thread.sleep(waitMiliSecs);
					}
				log("Done.");
				
				w.writeEmptyElement("hr");
				w.writeStartElement("div");
				w.writeCharacters(about);
				w.writeEndElement();//div
				
				w.writeEndElement();//body
				w.writeEndElement();//html
				w.flush();
				fout.flush();
				
				} 
			catch (Exception e)
				{
				e.printStackTrace();
				log("Error "+e.getMessage());
				}
			finally
				{
				CloserUtil.close(igvSocket);
				CloserUtil.close(out);
				CloserUtil.close(in);
				CloserUtil.close(w);
				CloserUtil.close(fout);
				CloserUtil.close(igvLog);
				pngFile.delete();
				}
			}
		}


	
	
	BatchIGVPicturesFrame()
		{
		setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);
		try {
			this.preferences=Preferences.userNodeForPackage(BatchIGVPictures.class);
			}	 
			catch (Exception e) {
			}
		addWindowListener(new WindowAdapter()
			{
			@Override
			public void windowClosing(WindowEvent e) {
				doMenuClose();
				}
			});
		JMenuBar menubar=new JMenuBar();
		setJMenuBar(menubar);
		JMenu menu=new JMenu("File");
		menubar.add(menu);
		menu.add(new AbstractAction("Quit")
			{
			@Override
			public void actionPerformed(ActionEvent e) {
				doMenuClose();
				}
			});
		menu=new JMenu("Help");
		menubar.add(menu);
		menu.add(new AbstractAction("About")
			{
			@Override
			public void actionPerformed(ActionEvent e) {
				JOptionPane.showMessageDialog(BatchIGVPicturesFrame.this, about);
				}
			});
		JPanel contentPane=new JPanel(new BorderLayout(5, 5));
		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));

		setContentPane(contentPane);
		
		JPanel top0=new JPanel(new GridLayout(0,1,5,5));
		contentPane.add(top0,BorderLayout.NORTH);
		
		
		JPanel top=new JPanel(new FlowLayout(FlowLayout.LEADING));
		top0.add(top);

		/* extend */
		JLabel lbl=new JLabel("Position Extend:",JLabel.RIGHT);
		top.add(lbl);
		JSpinner spin=new JSpinner(this.spinnerExtendSize=new SpinnerNumberModel(
				preferences.getInt(PREF_EXTEND, 10), 1, 100, 1));
		lbl.setLabelFor(spin);
		top.add(spin);
		this.spinnerExtendSize.addChangeListener(new ChangeListener()
			{
			@Override
			public void stateChanged(ChangeEvent e) {
				preferences.putInt(PREF_EXTEND,SpinnerNumberModel.class.cast(e.getSource()).getNumber().intValue());
				}
			});
		/* image width  */
		lbl=new JLabel("Image Width:",JLabel.RIGHT);
		top.add(lbl);
		spin=new JSpinner(this.spinnerImageWidth=new SpinnerNumberModel(
				preferences.getInt(PREF_IMG_WIDTH, 640), 1, 2000, 1));
		lbl.setLabelFor(spin);
		top.add(spin);
		this.spinnerImageWidth.addChangeListener(new ChangeListener()
			{
			@Override
			public void stateChanged(ChangeEvent e) {
				preferences.putInt(PREF_IMG_WIDTH,SpinnerNumberModel.class.cast(e.getSource()).getNumber().intValue());
				}
			});
		
		/* image width  */
		lbl=new JLabel("Wait ms:",JLabel.RIGHT);
		top.add(lbl);
		spin=new JSpinner(this.spinnerIgvWait=new SpinnerNumberModel(
				preferences.getInt(PREF_IGV_SLEEP, 100), 1, 10000, 1));
		lbl.setLabelFor(spin);
		top.add(spin);
		this.spinnerIgvWait.addChangeListener(new ChangeListener()
			{
			@Override
			public void stateChanged(ChangeEvent e) {
				preferences.putInt(PREF_IMG_WIDTH,SpinnerNumberModel.class.cast(e.getSource()).getNumber().intValue());
				}
			});

		
		
		/* IGV Host */
		lbl=new JLabel("IGV Host:",JLabel.RIGHT);
		top.add(lbl);
		this.tfIgvHost = new JTextField(
				preferences.get(PREF_IGV_HOST,IgvSocket.DEFAULT_HOST),8);
		top.add(this.tfIgvHost);
		this.tfIgvHost.getDocument().addDocumentListener(new DocumentListener()
			{
			@Override
			public void removeUpdate(DocumentEvent e) {
				changedUpdate(e);
				}
			
			@Override
			public void insertUpdate(DocumentEvent e) {
				changedUpdate(e);
				}
			
			@Override
			public void changedUpdate(DocumentEvent e) {
				String s;
				try {
					s = e.getDocument().getText(0, e.getDocument().getLength()).trim();
					preferences.put(PREF_IGV_HOST,s);
				} catch (BadLocationException e1) {
					e1.printStackTrace();
					}
				}
			});
		/* IGV port */
		lbl=new JLabel("IGV PORT:",JLabel.RIGHT);
		top.add(lbl);
		spin=new JSpinner(this.spinnerIgvPort=new SpinnerNumberModel(
				preferences.getInt(PREF_IGV_PORT, IgvSocket.DEFAULT_PORT), 1, 1000000, 1));
		lbl.setLabelFor(spin);
		top.add(spin);
		this.spinnerIgvPort.addChangeListener(new ChangeListener()
			{
			@Override
			public void stateChanged(ChangeEvent e) {
				preferences.putInt(PREF_IGV_PORT,SpinnerNumberModel.class.cast(e.getSource()).getNumber().intValue());
				}
			});
		
		
		top=new JPanel(new FlowLayout(FlowLayout.LEADING));
		top0.add(top);
		/* FORMAT */
		lbl=new JLabel("Format:",JLabel.RIGHT);
		top.add(lbl);
		this.cboxFormat =new JComboBox<String>(new String[]{"HTML"});
		this.cboxFormat.setSelectedIndex(preferences.getInt(PREF_FORMAT, FORMAT_HTML));
		lbl.setLabelFor(this.cboxFormat);
		top.add(this.cboxFormat);
		this.cboxFormat.addActionListener(new ActionListener()
			{
			@Override
			public void actionPerformed(ActionEvent e) {
				preferences.putInt(PREF_FORMAT,Math.max(cboxFormat.getSelectedIndex(),0));
				}
			});
		
		
		/* paint */
		JButton but=new JButton(new AbstractAction("Load Mutation File")
			{
			@Override
			public void actionPerformed(ActionEvent e) {
				loadVCF();
				}
			});
	
		top.add(but);
	

		
		
		but=new JButton(new AbstractAction("Stop")
			{
			@Override
			public void actionPerformed(ActionEvent e) {
				stopPaint();
				}
			});
		
		but.setFont(new Font("Dialog",Font.BOLD,18));
		but.setBackground(Color.RED);
		but.setForeground(Color.WHITE);
		top.add(but);
		
		but=new JButton(new AbstractAction("Draw")
			{
			@Override
			public void actionPerformed(ActionEvent e) {
				doMenuPaint();
				}
			});
		
		but.setFont(new Font("Dialog",Font.BOLD,18));
		but.setBackground(Color.GREEN);
		but.setForeground(Color.WHITE);
		top.add(but);

		
		
		
		JPanel pane2=new JPanel(new BorderLayout(5,5));
		contentPane.add(pane2,BorderLayout.CENTER);
		
		
		
		pane2.add(new JScrollPane(this.snpArea=new JTextArea(10,80)),BorderLayout.CENTER);
		pane2.setBorder(new TitledBorder("SNP: CHROM(tab)POSITION"));
		/*
		this.snpArea.setText(
				"chr6\t85356177\n"+
				"chr6\t85356228\n"+
				"chr6\t85356328\n"	
				);
		*/
		
		this.progressLbl=new JLabel("");
		contentPane.add(this.progressLbl,BorderLayout.SOUTH);
		}
	
	
	
	private void stopPaint()
		{
		if(this.currentSlave!=null)
			{
			try {
				this.currentSlave.interrupt();
				} 
			catch (Exception e) {
				
				}
			}
		this.currentSlave=null;
		this.progressLbl.setText("");
		}
	
	private void doMenuClose()
		{
		stopPaint();
		try
			{
			this.preferences.flush();
			}
		catch(Exception err)
			{
			
			}
		this.setVisible(false);
		this.dispose();
		}
	
	private void loadVCF()
		{
		String prevFile=this.preferences.get(PREF_LOAD_VCF_FILE, null);
		File selFile=prevFile==null?null:new File(prevFile);
		JFileChooser fc=new JFileChooser();
		if(selFile!=null)
			{
			fc.setSelectedFile(selFile);
			}
		if(fc.showOpenDialog(this)!=JFileChooser.APPROVE_OPTION)
			{
			return ;
			}
		BufferedReader r=null;
		try
			{
			r=IOUtils.openFileForBufferedReading(fc.getSelectedFile());
			Set<Mutation> mutations=this.loadMutations(r);
			StringWriter sw=new StringWriter();
			for(Mutation m:mutations)
				{
				sw.write(m.chrom+"\t"+m.position+"\n");
				}
			this.snpArea.setText(sw.toString());
			this.snpArea.setCaretPosition(0);
			this.preferences.put(PREF_LOAD_VCF_FILE, fc.getSelectedFile().getPath());
			}
		catch(Exception err)
			{
			JOptionPane.showMessageDialog(this,
				"Cannot pars data in "+fc.getSelectedFile()+" "+err.getMessage());
			}
		finally
			{
			CloserUtil.close(r);
			}
		}
	
	private Set<Mutation> loadMutations(Reader reader)
		throws IOException
		{
		Set<Mutation> mutations=new TreeSet<Mutation>();
		BufferedReader in=new BufferedReader(reader);
		String line;
		while(((line=in.readLine()))!=null)
			{
			if(line.isEmpty() || line.startsWith("#")) continue;
			String tokens[]=line.split("[\t]",3);
			if(tokens.length<2) throw new IOException("Bad line "+line);
			Mutation m=new Mutation();
			m.chrom=tokens[0].trim();
			if(!m.chrom.matches("[A-Z0-9a-z_\\.]+"))
				{
				throw new IOException("Bad chrom in "+line);
				}
			try {
				m.position=Integer.parseInt(tokens[1].replace(",", "").replace(".", ""));
			} catch (Exception e) {
				throw new IOException("Bad Position in in "+line);
				}
			if(m.position<0) throw new IOException("Bad Position in in "+line);
			mutations.add(m);
			}
		return mutations;
		}

	
	private void doMenuPaint()
		{
		stopPaint();
		int fotmat=Math.max(this.cboxFormat.getSelectedIndex(),0);
		
		this.currentSlave= new HtmlSlave();
		

		this.currentSlave.igvHost= this.tfIgvHost.getText().trim();
		this.currentSlave.igvPort= this.spinnerIgvPort.getNumber().intValue();
		this.currentSlave.extend= this.spinnerExtendSize.getNumber().intValue();
		this.currentSlave.imageWidth=this.spinnerImageWidth.getNumber().intValue();
		this.currentSlave.waitMiliSecs=this.spinnerIgvWait.getNumber().intValue();
		try {
			
			this.currentSlave.mutations=loadMutations(new StringReader(this.snpArea.getText()));
			}
		catch (final Exception e) {
			JOptionPane.showMessageDialog(this,"Cannot parse SNPs area:"+e.getMessage());
			stopPaint();
			return;
			}
		
		if(this.currentSlave.mutations.isEmpty())
			{
			stopPaint();
			return;
			}
		
		File file=null;
		String fileStr=this.preferences.get(PREF_SAVE_AS_FILE, null);
		if(fileStr!=null)
			{
			file=new File(fileStr);
			}
		
			
		JFileChooser chooser=new JFileChooser(file!=null?file.getParentFile():null);
		if(file!=null )
			{
			chooser.setSelectedFile(file);
			}
		if(chooser.showSaveDialog(this)!=JFileChooser.APPROVE_OPTION)
			{
			return;
			}
		
		File selectedFile=chooser.getSelectedFile();
		
		
	
		this.currentSlave.pathOut=selectedFile.getPath();
		this.preferences.put(PREF_SAVE_AS_FILE, selectedFile.getPath());
		if(fotmat==FORMAT_HTML)
			{
			if(this.currentSlave.pathOut.endsWith(".html"))
				{
				this.currentSlave.pathOut=this.currentSlave.pathOut.substring(0,this.currentSlave.pathOut.length()-5);
				}
			else if(this.currentSlave.pathOut.endsWith(".htm"))
				{
				this.currentSlave.pathOut=this.currentSlave.pathOut.substring(0,this.currentSlave.pathOut.length()-4);
				}
			selectedFile=new File(this.currentSlave.pathOut+".html");

			}
		
		else
			{
			throw new IllegalStateException();
			}
		if(selectedFile!=null && selectedFile.exists() &&
				JOptionPane.showConfirmDialog(this, selectedFile.toString()+" already exists. Overwrite ?", "Overwrite", JOptionPane.YES_NO_OPTION, JOptionPane.WARNING_MESSAGE, null)!=JOptionPane.OK_OPTION
				)
			{
			stopPaint();
			return;
			}
		this.currentSlave.start();
		}

	
	}

public class BatchIGVPictures extends Launcher
	{
	private static final Logger LOG=Logger.build(BatchIGVPictures.class).make();
	
		@Override
		public int doWork(final List<String> args) {

			final BatchIGVPicturesFrame app=new BatchIGVPicturesFrame();
			
			app.about=" Author: Pierre Lindenbaum "+
						" Version:"+getVersion()+
						" Date:"+getCompileDate()
						;
			try
				{
				SwingUtilities.invokeAndWait(new Runnable() {
					@Override
					public void run()
						{
						Dimension win=Toolkit.getDefaultToolkit().getScreenSize();
						app.setTitle(getProgramName());
						app.pack();
						Dimension ps=app.getPreferredSize();
						app.setLocation(
								(win.width-ps.width)/2,
								(win.height-ps.height)/2
								);
						app.setVisible(true);
						}
					});
				return RETURN_OK;
				}
			catch(Exception err)
				{
				LOG.error(err);
				return -1;
				}
			finally
				{
				}
			}
		
	
	public static void main(String[] args) {
		new BatchIGVPictures().instanceMain(args);
		}
	}
