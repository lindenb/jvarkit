package com.github.lindenb.jvarkit.util.jcommander;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.function.Function;
import java.util.jar.Manifest;
import java.util.prefs.Preferences;
import java.util.stream.Collectors;
import java.util.zip.Deflater;

import javax.swing.AbstractAction;
import javax.swing.ActionMap;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JToolBar;
import javax.swing.SwingUtilities;

import com.beust.jcommander.IStringConverter;
import com.beust.jcommander.IStringConverterFactory;
import com.beust.jcommander.IValueValidator;
import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParameterException;
import com.beust.jcommander.ParametersDelegate;
import com.beust.jcommander.converters.IntegerConverter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;



public class Launcher {
private static final Logger LOG=Logger.build().
			prefix("Launcher").
			make();
public static final String[]OUTPUT_OPTIONS={"-o","--out"};
protected static final int RETURN_OK=0;
public enum Status { OK, PRINT_HELP,PRINT_VERSION,EXIT_SUCCESS,EXIT_FAILURE};


/** custom instance of jcommander, don't add same command twice. */
private class MyJCommander extends JCommander
	{
	/** when registering the option for jcommander, we take care of not visiting object twice */
	private Collection<Object> ojectsVisitedByJCommander = new ArrayList<>();

	
	@Override
	public void usage(final StringBuilder sb) {
		final Class<?> clazz=Launcher.this.getClass();
		
		final Program programdesc=clazz.getAnnotation(Program.class);
		
		if(programdesc!=null){
			this.setProgramName(programdesc.name());
		} else
			{
			this.setProgramName(Launcher.this.getProgramName());
			}
		
		if(print_markdown_help) sb.append("\n```\n");
		super.usage(sb);
		if(print_markdown_help) sb.append("\n```\n\n");

		if(programdesc!=null){
			if(!programdesc.deprecatedMsg().isEmpty())
				{
				sb.append("\n##DEPRECATED\n\n").
					append(programdesc.deprecatedMsg()).
					append("\n");
				}
			
			sb.append("\n##Description\n\n").
				append(programdesc.description()).
				append("\n");
			
			if(programdesc.keywords()!=null && programdesc.keywords().length>0) {
				sb.append("\n##Keywords\n\n");
				for(String sk:programdesc.keywords()) sb.append(" * "+sk+"\n");
				sb.append("\n");
			}
			if(programdesc.biostars()!=null && programdesc.biostars().length>0) {
				sb.append("\n## See also in Biostars\n\n");
				for(int postid:programdesc.biostars()) sb.append(" * https://www.biostars.org/p/"+postid+"\n");
				sb.append("\n");
			}
			
		}
		
		if(print_markdown_help)
			{
			final String progName=(programdesc==null?"software":programdesc.name());
			sb.append("##Compilation\n");
			sb.append("\n");
			sb.append("### Requirements / Dependencies\n");
			sb.append("\n");
			sb.append("* java compiler SDK 1.8 http://www.oracle.com/technetwork/java/index.html (**NOT the old java 1.7 or 1.6**) . Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )\n");
			sb.append("* GNU Make >= 3.81\n");
			sb.append("* curl/wget\n");
			sb.append("* git\n");
			sb.append("* xsltproc http://xmlsoft.org/XSLT/xsltproc2.html (tested with \"libxml 20706, libxslt 10126 and libexslt 815\")\n");
			sb.append("\n");
			sb.append("\n");
			sb.append("### Download and Compile\n");
			sb.append("\n");
			sb.append("```bash\n");
			sb.append("$ git clone \"https://github.com/lindenb/jvarkit.git\"\n");
			sb.append("$ cd jvarkit\n");
			sb.append("$ make "+progName+"\n");
			sb.append("```\n");
			sb.append("\n");
			sb.append("The *.jar libraries are not included in the main jar file, so you shouldn\'t move them (https://github.com/lindenb/jvarkit/issues/15#issuecomment-140099011 ).");
			sb.append("\n");
			sb.append("The required libraries will be downloaded and installed in the `dist` directory.\n");
			sb.append("\n");
			sb.append("### edit \'local.mk\' (optional)\n");
			sb.append("\n");
			sb.append("The a file **local.mk** can be created edited to override/add some definitions.\n");
			sb.append("\n");
			sb.append("For example it can be used to set the HTTP proxy:\n");
			sb.append("\n");
			sb.append("```\n");
			sb.append("http.proxy.host=your.host.com\n");
			sb.append("http.proxy.port=124567\n");
			sb.append("```\n");
			
			sb.append("## Source code \n\n");
			sb.append("https://github.com/lindenb/jvarkit/tree/master/src/main/java/").
				append(clazz.getName().replace('.','/')).append(".java\n");
			sb.append("\n");
			
			sb.append("## Contribute\n");
			sb.append("\n");
			sb.append("- Issue Tracker: http://github.com/lindenb/jvarkit/issues\n");
			sb.append("- Source Code: http://github.com/lindenb/jvarkit\n");
			sb.append("\n");
			sb.append("## License\n");
			sb.append("\n");
			sb.append("The project is licensed under the MIT license.\n");
			sb.append("\n");
			sb.append("## Citing\n");
			sb.append("\n");
			sb.append("Should you cite **"+progName +"** ? https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md\n");
			sb.append("\n");
			sb.append("The current reference is:\n");
			sb.append("\n");
			sb.append("http://dx.doi.org/10.6084/m9.figshare.1425030\n");
			sb.append("\n");
			sb.append("> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.\n");
			sb.append("> http://dx.doi.org/10.6084/m9.figshare.1425030\n");
			sb.append("\n");
			}
		InputStream in=null;
		try {
			String className=clazz.getName();
			int dollar=className.indexOf('$');
			if(dollar!=-1) className=className.substring(0, dollar);
			className=className.replace('.', '/')+".java";
			in=clazz.getResourceAsStream("/"+className);
			if(in!=null){
				BufferedReader r=new BufferedReader(new InputStreamReader(in));
				String line;
				boolean ok=false;
				while((line=r.readLine())!=null)
					{
					if(line.contains("BEGIN"+"_DOC"))
						{
						ok=true;
						}
					else if(line.contains("END"+"_DOC"))
						{
						ok=false;
						}
					else if(ok)
						{
						sb.append(line).append("\n");
						}
					}
				r.close();
				}
			else
				{
				LOG.debug("cannot find java code for "+className);
				}
			}
		catch(final Exception err) {
			
			}
		finally
			{
			CloserUtil.close(in);
			}
		}
	
	@Override
	public void addCommand(String name, Object object, String... aliases) 
			{
			if(this.ojectsVisitedByJCommander.stream().anyMatch(O->O==object)) return;
			this.ojectsVisitedByJCommander.add(object);
			super.addCommand(name, object, aliases);
			}
	}
	

/**
 * Special converter for Zip compression. Bound the values between 0 and 9
 * "best" is interpreted as BEST_COMPRESSION
 * "none" is no compression
 */
public static class CompressionConverter
extends IntegerConverter implements Function<String, Integer> {
	public CompressionConverter() {
		super("");
		}
	public CompressionConverter(final String arg) {
		super(arg);
		}

	@Override
	public final Integer apply(String t) {
		return convert(t);
		}
	
	@Override
	public Integer convert(final String s) {
		if(s!=null) {
			if(s.equals("best")) return Deflater.BEST_COMPRESSION;
			if(s.equals("none")) return Deflater.NO_COMPRESSION;
		}
		final Integer n = super.convert(s);
		if(n!=null) {
			if(n<0) return Deflater.NO_COMPRESSION;
			if(n>9) return Deflater.BEST_COMPRESSION;
		}
		return n;
	}
	@Override
	public String toString() {
		return "Compression converter";
		}
	}




/** original arc/argv */
private List<String> argcargv=Collections.emptyList();
private final JCommander jcommander = new MyJCommander();
/** git hash in the manifest */
private String gitHash = null;
/** compile date in the manifest */
private String compileDate = null;

@Parameter(names = {"-h","--help"},description="print help and exits", help = true)
private boolean print_help = false;
@Parameter(names = {"--markdownhelp"},description="print Markdown help and exits", help = true,hidden=true)
private boolean print_markdown_help = false;

@Parameter(names = {"--version"}, help = true,description="print version and exits")
private boolean print_version = false;
@Parameter(description = "Files")
private List<String> files = new ArrayList<>();

private String programName="";

public class CompressionArgs
	{
	@Parameter(names={"--compression"},description="Compression Level.",converter=CompressionConverter.class)
	public int compressionLevel=5;
	}
public CompressionArgs compressionArgs=new CompressionArgs();

public static class DirectoryExists implements IValueValidator<File> {
	@Override
	public void validate(String arg, final File dir) throws ParameterException {
		if(dir==null || !dir.exists() || !dir.isDirectory()) {
			throw new ParameterException("option :"+arg+": this is not a directory or it doesn't exists: "+dir);
			}
		}
	}

public static class TmpDirectoryArgs
	{
	@Parameter(names={"--tmpDir"},description="Temporary Directory.",validateValueWith=DirectoryExists.class)
	public File tmpDir=new File(System.getProperty("java.io.tmpdir","."));
	}


public static class SortingCollectionArgs
	{
	@ParametersDelegate
	private TmpDirectoryArgs tmpDirArgs;
	SortingCollectionArgs(TmpDirectoryArgs tmpDirArgs) {
		this.tmpDirArgs=tmpDirArgs;
		}
	@Parameter(names={"--xx"},description="Compression Level.",converter=CompressionConverter.class)
	public int compressionLevel=5;
	}
public class WritingSortingCollection
	{
	@Parameter(names={"--maxRecordsInRam"},description="When writing  files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file  handles needed to sort a file, and increases the amount of RAM needed")
	public int maxRecordsInRam=50000;
	
	@Parameter(names={"--tmpDir"},description= "tmp working directory. Default: java.io.tmpDir")
	private List<File> tmpDirs=new ArrayList<>();
	
	
	public WritingSortingCollection maxRecordsInRam(int n)
		{
		this.maxRecordsInRam = n;
		return this;
		}
	public int getMaxRecordsInRam() { return this.maxRecordsInRam;}
	public List<File> getTmpDirectories() {
		final List<File> L= new ArrayList<>(this.tmpDirs);
		if(L.isEmpty() )
			{
			L.add(new File(System.getProperty("java.io.tmpdir")));
			}
		return L;
		}
	
	}	

public class WritingBamArgs
	{
	private File referenceFile = null;
	
	@Parameter(names={"--bamcompression"},description="Compression Level.")
	public int compressionLevel=5;
	@Parameter(names={"--samoutputformat"},description="Sam output format.")
	public htsjdk.samtools.SamReader.Type samoutputformat = htsjdk.samtools.SamReader.Type.SAM_TYPE;
	
	/** creates a SAMFileWriterFactory */
	public htsjdk.samtools.SAMFileWriterFactory createSAMFileWriterFactory() {
		final SAMFileWriterFactory sfw =  new SAMFileWriterFactory();
		int n= this.compressionLevel;
		if(n<0) n= Deflater.NO_COMPRESSION;
		if(n>9) n= Deflater.BEST_COMPRESSION;
		sfw.setCompressionLevel(n);
		
		return sfw;
		}
	
	public WritingBamArgs setReferenceFile(final File referenceFile) {
		this.referenceFile = referenceFile;
		return this;
		}
	
	public WritingBamArgs setCompressionLevel(final int compressionLevel) {
		this.compressionLevel = compressionLevel;
		return this;
		}
	public WritingBamArgs setBestCompression() {
		return setCompressionLevel(Deflater.BEST_COMPRESSION);
		}
	/** return reference file. default implementation returns null */
	public File getReferenceFile() {
		return this.referenceFile;
		}
	
	
	
	public SAMFileWriter openSAMFileWriter(File outputFileOrNull,SAMFileHeader header,boolean presorted)
	{
		final htsjdk.samtools.SAMFileWriterFactory sfw= this.createSAMFileWriterFactory();
		
		
		if(outputFileOrNull==null)
			{
			if( this.samoutputformat!=null &&
				this.samoutputformat.equals(htsjdk.samtools.SamReader.Type.BAM_TYPE))
				{
				return sfw.makeBAMWriter(header, presorted, stdout());
				}
			else if(this.samoutputformat==null || this.samoutputformat.equals(htsjdk.samtools.SamReader.Type.SAM_TYPE))
				{
				return sfw.makeSAMWriter(header, presorted, stdout());
				}
			else if(this.samoutputformat==null || this.samoutputformat.equals(htsjdk.samtools.SamReader.Type.CRAM_TYPE))
				{
				return sfw.makeCRAMWriter(header,stdout(),getReferenceFile());
				}
			else
				{
				throw new IllegalStateException("Bad output format");
				}
			}
		else
			{
			return sfw.makeWriter(header, presorted, outputFileOrNull, getReferenceFile());
			}

		}
	
	}

public static class PrintWriterOnDemand
	extends PrintWriter
	{
	private final File fileout;
	private boolean delegate_created;
	public PrintWriterOnDemand(final File out) {
		super(new NullOuputStream());
		this.fileout=out;
		}
	public PrintWriterOnDemand() {
		this(null);
		}
	@Override
	public void write(char[] buf, int off, int len) {
		if(!this.delegate_created)
			{
			 try {
		         synchronized (lock)
		          	{
		        	if(this.delegate_created) {
		        		}
		        	else if(this.fileout==null)
					   {
					   super.out=new PrintWriter(System.out);
					   }
		        	else
						{
						super.out =IOUtils.openFileForPrintWriter(fileout);
						}
		          	}
		         this.delegate_created=true;
				}    
		    catch(final Exception err)
		    	{
		    	throw new RuntimeIOException(err);
		    	}
			}
		super.write(buf, off, len);
		}
	}

public static class VcfWriterOnDemandConverter
	implements IStringConverter<VcfWriterOnDemand> {
	@Override
	public VcfWriterOnDemand convert(String s) {
		if(s.equals("-") || s.equals("stdin") || s.isEmpty()) {
			return new VcfWriterOnDemand();
			}
		else
			{
			return new VcfWriterOnDemand(new File(s));
			}
		}
	}

public static class DimensionConverter
	implements IStringConverter<Dimension>
{
	@Override
	public Dimension convert(String v) {
		int x=v.indexOf('x');
		if(x<1)
			{
			throw new ParameterException("bad size. Expected (width)x(heigh) "+v);
			}
		int width=Integer.parseInt(v.substring(0, x));
		int height=Integer.parseInt(v.substring(x+1));
		return new Dimension(width, height);
		}
}

public static class VcfWriterOnDemand
	implements VariantContextWriter
	{
	private VariantContextWriter delegate=null;
	private final VariantContextWriterBuilder vcb=new VariantContextWriterBuilder().
				clearIndexCreator().
				clearOptions()
				;
	private List<VCFHeaderLine> extraHeaderLines=new ArrayList<>();
	public VcfWriterOnDemand() {
		vcb.setOutputVCFStream(System.out);
		}
	VcfWriterOnDemand(final File file)
		{
		vcb.setOutputFile(file);
		}
	@Override
	public void writeHeader(final VCFHeader header) {
		if(this.delegate==null) {
			this.delegate =vcb.build(); 
			}
		VCFHeader header2 = header;
		if(!this.extraHeaderLines.isEmpty()) {
			header2= new VCFHeader(header);
			for(final VCFHeaderLine hl:this.extraHeaderLines)
				{
				header2.addMetaDataLine(hl);
				}
		}
		this.delegate.writeHeader(header2);
		}
	
	
	@Override
	public void add(final VariantContext ctx) {
		if(this.delegate==null) throw new JvarkitException.ProgrammingError("null delegate");
		this.delegate.add(ctx);
		}
	@Override
	public boolean checkError() {
		return (this.delegate==null?false:delegate.checkError());
		}
	@Override
	public void close() {
		CloserUtil.close(this.delegate);
		}
	@Override
	public String toString() {
		return "Default VCF writer (stdout)";
		}
	}

public static class JConsole extends JFrame
	{
	private final Preferences prefs;
	private File lastDir=null;
	protected final ActionMap actionMap=new ActionMap();
	protected final JTextArea textArea;
	private final Class<? extends Launcher> clazz;
	private volatile Runner runner=null;
	private  class Runner
		extends Thread
		{
		Launcher instance;
		String args[];
		int returnStatus=0;
		@Override
		public void run() {
			try
				{
				int ret= instance.instanceMain(args);
				this.returnStatus=ret;
				}
			catch(Exception err)
				{
				LOG.error(err);
				try { SwingUtilities.invokeAndWait(()->{
					if(Runner.this != JConsole.this.runner) return;
					 JConsole.this.runner=null;
					JOptionPane.showMessageDialog(JConsole.this, "An error occured "+err.getMessage());
					}); } catch(Throwable err2) {
						LOG.error(err2);
					}
				returnStatus=-1;
				return;
				}
			
			try { SwingUtilities.invokeAndWait(()->{
				if(Runner.this != JConsole.this.runner) return;
				 JConsole.this.runner=null;
				JOptionPane.showMessageDialog(
						JConsole.this,
						"Program exited with status "+returnStatus,
						"End",
						returnStatus==0?JOptionPane.PLAIN_MESSAGE:JOptionPane.ERROR_MESSAGE
						);
				}); } catch(Throwable err2) {
					LOG.error(err2);
				}
			}
		
		}
	
	@SuppressWarnings("serial")
	JConsole(Class<? extends Launcher> clazz) {
		super("JConsole");
		
		this.clazz = clazz;
		final Program programAnnot = clazz.getAnnotation(Program.class);
		if(programAnnot!=null)
			{
			this.setTitle(programAnnot.name());
			}
		this.prefs = Preferences.userNodeForPackage(JConsole.class);
		String pref= this.prefs.get(this.clazz.getName()+".lastDir", null);
		if(pref!=null)
			{
			this.lastDir=new File(pref);
			}
		
		this.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		final JMenuBar menubar=new JMenuBar();
		setJMenuBar(menubar);
		this.actionMap.put("launcher.jconsole.about", new AbstractAction("About...") {
			@Override
			public void actionPerformed(ActionEvent e) {
				JOptionPane.showMessageDialog(JConsole.this, programAnnot.description());
				}
			});
		this.actionMap.put("launcher.jconsole.exit", new AbstractAction("Exit") {
			@Override
			public void actionPerformed(ActionEvent e) {
				doMenuQuit();
				}
			});
		this.actionMap.put("launcher.jconsole.run", new AbstractAction("Run") {
			@Override
			public void actionPerformed(ActionEvent e) {
				doMenuRun();
				}
			});
		this.actionMap.put("launcher.jconsole.stop", new AbstractAction("Stop") {
			@Override
			public void actionPerformed(ActionEvent e) {
				doMenuStop();
				}
			});

		this.actionMap.put("launcher.jconsole.insertpathr", new AbstractAction("Insert Path/R") {
			@Override
			public void actionPerformed(ActionEvent e) {
				doMenuInsertPath(true);
				}
			});
		this.actionMap.put("launcher.jconsole.insertpathw", new AbstractAction("Insert Path/W") {
			@Override
			public void actionPerformed(ActionEvent e) {
				doMenuInsertPath(false);
				}
			});
		
		JMenu jmenu=new JMenu("File");
		jmenu.add(this.actionMap.get("launcher.jconsole.about"));
		jmenu.add(this.actionMap.get("launcher.jconsole.exit"));
		menubar.add(jmenu);
		jmenu=new JMenu("Tool");
		jmenu.add(this.actionMap.get("launcher.jconsole.run"));
		jmenu.add(this.actionMap.get("launcher.jconsole.insertpathr"));
		jmenu.add(this.actionMap.get("launcher.jconsole.insertpathw"));
		menubar.add(jmenu);
		JPanel contentPane=new JPanel(new BorderLayout(5,5));
		final JToolBar toolbar=new JToolBar();
		contentPane.add(toolbar,BorderLayout.NORTH);
		toolbar.add(this.actionMap.get("launcher.jconsole.run"));
		toolbar.add(this.actionMap.get("launcher.jconsole.insertpathr"));
		toolbar.add(this.actionMap.get("launcher.jconsole.insertpathw"));
		
		this.textArea = new JTextArea(5, 40);
		
		pref= this.prefs.get(this.clazz.getName()+".text", null);
		if(pref!=null)
			{
			this.textArea.setText(pref);
			}
		
		contentPane.add(new JScrollPane(this.textArea),BorderLayout.CENTER);
		
		final JPanel bot=new JPanel(new FlowLayout(FlowLayout.TRAILING,5,5));
		contentPane.add(bot,BorderLayout.SOUTH);
		
		bot.add(new JButton(this.actionMap.get("launcher.jconsole.stop")));
		bot.add(new JButton(this.actionMap.get("launcher.jconsole.run")));
		
		setContentPane(contentPane);
		
		this.addWindowListener(new WindowAdapter() {
			@Override
			public void windowClosing(final WindowEvent e) {
				try
					{
					prefs.put(clazz.getName()+".text", textArea.getText());
					if(lastDir!=null) prefs.put(clazz.getName()+".lastDir",lastDir.getPath());
					prefs.flush();
					}
				catch(Exception err)
					{
					LOG.error(err);
					}
				}
			});
		}
	

	protected void doMenuInsertPath(boolean openDialog) {
		final JFileChooser fc=new JFileChooser(this.lastDir);
		if(openDialog)
			{
			if(fc.showOpenDialog(JConsole.this)!=JFileChooser.APPROVE_OPTION) return;
			}
		else
			{
			if(fc.showSaveDialog(JConsole.this)!=JFileChooser.APPROVE_OPTION) return;
			}
		final File f=fc.getSelectedFile();
		if(f==null) return;
		this.lastDir=f.getParentFile();
		this.textArea.insert(f.getPath(), this.textArea.getCaretPosition());
		}
	
	private synchronized void doMenuStop() {
		if(this.runner!=null)
			{
			try { this.runner.interrupt();}
			catch(Throwable err)  {}
			this.runner=null;
			}
		}

	protected synchronized void doMenuRun() 
		{
		if(this.runner!=null)
			{	
			JOptionPane.showMessageDialog(JConsole.this, "Program is already running....");
			return;
			}
		final List<String> args;
		try {
			args = this.getArguments();
			}
		catch(final IllegalArgumentException err)
			{
			LOG.error(err);
			JOptionPane.showMessageDialog(JConsole.this,
					String.valueOf(err.getMessage()),
					"Error",
					JOptionPane.ERROR_MESSAGE
					);
			return;
			}
		Runner run =new Runner();
		final Launcher instance;
		try {
			run.instance = this.clazz.newInstance();
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			JOptionPane.showMessageDialog(JConsole.this,
					String.valueOf(err.getMessage()),
					"Error creating a new instance of "+clazz.getName(),
					JOptionPane.ERROR_MESSAGE
					);
			return;
			}
		run.args = args.toArray(new String[args.size()]);
		this.runner=run;
		this.runner.start();
		
		}
	protected void doMenuQuit() 
		{
		setVisible(false);
		dispose();
		}
	protected List<String> getArguments() {
		final List<String> args = new ArrayList<>();
		final String s=this.textArea.getText();
		int i=0;
		while(i<s.length())
			{
			if(Character.isWhitespace(s.charAt(i))) {++i;continue;}
			if(s.charAt(i)=='\"' || s.charAt(i)=='\'')
				{
				char quote=s.charAt(i);
				i++;
				final StringBuilder sb=new StringBuilder();
				while(i< s.length())
					{
					char c = s.charAt(i);
					++i;
					if(c==quote) break;
					if(c=='\\')
						{
						if(i+1>=s.length())
							{	
							throw new IllegalArgumentException("Unclosed string after "+sb.toString());
							}
						c= s.charAt(i);
						switch(c)
							{
							case '\'': sb.append('\'');break;
							case '\"': sb.append('\"');break;
							case 'n': sb.append('\n');break;
							case 't': sb.append('\t');break;
							case '\\': sb.append('\\');break;
							default: throw new IllegalArgumentException("Unknown escape sequence after: "+sb.toString());
							}
						}
					else
						{
						sb.append(c);
						}
					}
				args.add(sb.toString());
				}
			else
				{
				final StringBuilder sb=new StringBuilder();
				while(i< s.length() && !Character.isWhitespace(s.charAt(i)))
					{
					sb.append(s.charAt(i));
					i++;
					}
				args.add(sb.toString());
				}
			}
		return args;
		}

	}



public Launcher()
	{
	final Class<?> clazz=Launcher.this.getClass();
	final Program programdesc=clazz.getAnnotation(Program.class);
	if(programdesc!=null)
		{
		this.programName = programdesc.name();
		}
	else
		{
		this.programName = getClass().getSimpleName();
		}	
	
	try {
		/** set locale http://seqanswers.com/forums/showthread.php?p=174020#post174020 */
		Locale.setDefault(Locale.US);
		}
	catch(final java.security.AccessControlException err) {
		System.err.println("Cannot set Locale to US for security reasons"+err.getMessage());
		}
	try {
		System.setProperty("file.encoding", "UTF-8");
		}
	catch(final java.security.AccessControlException err) {
		System.err.println("Cannot set file.encoding to UTF-8 for security reasons"+err.getMessage());
		}
	try {
	/* https://bugs.openjdk.java.net/browse/JDK-8028111 */
	System.setProperty("jdk.xml.entityExpansionLimit","0");
	}
	catch(final java.security.AccessControlException err) {
	}
	
	 final Map<Class, Class<? extends IStringConverter<?>>> MAP = new HashMap() {{
		    put(VcfWriterOnDemand.class, VcfWriterOnDemandConverter.class);
		    put(VariantContextWriter.class, VcfWriterOnDemandConverter.class);
		    put(Dimension.class,DimensionConverter.class);
		}};	
	this.jcommander.addConverterFactory(new IStringConverterFactory() {
			@Override
			public Class<? extends IStringConverter<?>> getConverter(Class forType) {		
				return MAP.get(forType);
				}
			});
	this.jcommander.addObject(this);	
	}

public String getProgramName()
	{
	return this.programName;
	}

public String getCompileDate()
	{
	if(this.compileDate==null)
		{
		this.compileDate="undefined";
		loadManifest();
		}
	return compileDate;
}

public String getGitHash()
{
if(this.gitHash==null)
	{
	this.gitHash="1.0";
	loadManifest();
	}
return this.gitHash;
}

public String getVersion()
{
return getGitHash();
}


private void loadManifest()
	{
	try
		{
		final Enumeration<URL> resources = getClass().getClassLoader()
				  .getResources("META-INF/MANIFEST.MF");//not '/META-INF'
		while (resources.hasMoreElements())
			{
			final URL url=resources.nextElement();
			InputStream in=url.openStream();
			if(in==null)
				{
				continue;
				}
			
			Manifest m=new Manifest(in);
			in.close();
			in=null;
			final java.util.jar.Attributes attrs=m.getMainAttributes();
			if(attrs==null)
				{
				continue;
				}
			String s =attrs.getValue("Git-Hash");
			if(s!=null && !s.isEmpty() && !s.contains("$")) //ant failed
				{
				this.gitHash=s;
				}
			s =attrs.getValue("Compile-Date");
			if(s!=null && !s.isEmpty()) //ant failed
				{
				this.compileDate=s;
				}
			}
		}	
	catch(Exception err)
		{
		
		}
	
}



protected JCommander getJCommander()
	{
	return this.jcommander;
	}

/** called AFTER argc/argv has been initialized */
protected int initialize() {
	return 0;
	}

/** called AFTER the work is done, before returning the status */
protected void cleanup() {
	
	}


public String getProgramCommandLine() {
	return String.join(" ",this.argcargv);
	}

protected Status parseArgs(final String args[])
	{
	this.argcargv = Arrays.asList(args);
	try
	  	{
		getJCommander().parse(args);
	  	}
	 catch(final com.beust.jcommander.ParameterException err) {
		stderr().println("There was an error in the input parameters.");
		stderr().println(err.getMessage());
		return Status.EXIT_FAILURE; 
	 	}
	 
	 if (this.print_help || this.print_markdown_help) return Status.PRINT_HELP;
	 if (this.print_version) return Status.PRINT_VERSION;
	 return Status.OK;
	}

protected VcfIterator openVcfIterator(final String inputNameOrNull) throws IOException {
	return VCFUtils.createVcfIterator(inputNameOrNull);
}

protected VariantContextWriter openVariantContextWriter(final File outorNull) throws IOException {
	return VCFUtils.createVariantContextWriter(outorNull);
}

protected InputStream openInputStream(final String inOrNull) throws IOException {
	return(inOrNull==null?
			stdin():
			IOUtils.openURIForReading(inOrNull)
			);
}

protected BufferedReader openBufferedReader(final String inOrNull) throws IOException {
	return(inOrNull==null?
			new BufferedReader(new InputStreamReader(stdin())):
			IOUtils.openURIForBufferedReading(inOrNull)
			);
	}


protected VCFHeader addMetaData(final VCFHeader header) 
	{
	return header;
	}

protected int doVcfToVcf(final String inputName,final VcfIterator iterin,final VariantContextWriter out){
	LOG.debug("using default doVcfToVcf ??");
	VCFUtils.copyHeaderAndVariantsTo(iterin, out);
	return 0;
	}
protected int doVcfToVcf(final String inputNameOrNull,final File outorNull){
	VcfIterator iterin=null;
	VariantContextWriter w=null;
	int ret=0;
	try {
		iterin = openVcfIterator(inputNameOrNull);
		w = openVariantContextWriter(outorNull);
		ret=doVcfToVcf(inputNameOrNull==null?"<STDIN>":inputNameOrNull,iterin,w);
		w.close();
		w=null;
		iterin.close();
		iterin=null;
		return ret;
		}
	catch(final Exception err)
		{
		LOG.error(err);
		return -1;
		}
	finally
		{
		CloserUtil.close(iterin);
		CloserUtil.close(w);
		}
	}

protected int doVcfToVcf(final List<String> inputs,final File outorNull) {
	return doVcfToVcf(oneFileOrNull(inputs),outorNull);
	}


protected String oneAndOnlyOneFile(final List<String> args) {
	switch(args.size())
		{
		case 1: return args.get(0);
		default: throw new JvarkitException.CommandLineError("Expected one and only one argument but got "+args.size());
		}
}


protected String oneFileOrNull(final List<String> args) {
	switch(args.size())
	{
	case 0: return null;
	case 1: return args.get(0);
	default: throw new JvarkitException.CommandLineError("Expected one or zero argument but got "+args.size());
	}
}

protected boolean evalJavaScriptBoolean(
		final javax.script.CompiledScript compiledScript,
		final javax.script.Bindings bindings) throws javax.script.ScriptException
		{
		Object result = compiledScript.eval(bindings);
		if(result==null) return false;
		if(result instanceof Boolean)
			{
			if(Boolean.FALSE.equals(result)) return false;
			}
		else if(result instanceof Number)
			{
			if(((Number)result).intValue()!=1) return false;
			}
		else
			{
			LOG.warn("Script returned something that is not a boolean or a number:"+result.getClass());
			 return false;
			}
		return true;
		}


/** reads a Bed file and convert it to a IntervalTreeMap<Boolean> */
protected IntervalTreeMap<Boolean> readBedFileAsBooleanIntervalTreeMap(final java.io.File file) throws java.io.IOException
	{
	java.io.BufferedReader r=null;
	try
		{
		final  IntervalTreeMap<Boolean> intervals = new IntervalTreeMap<Boolean>();
		r=com.github.lindenb.jvarkit.io.IOUtils.openFileForBufferedReading(file);
		final BedLineCodec bedCodec=new BedLineCodec();
		r.lines().
			filter(line->!(line.startsWith("#") ||  com.github.lindenb.jvarkit.util.bio.bed.BedLine.isBedHeader(line) ||  line.isEmpty())).
			map(line->bedCodec.decode(line)).
			filter(B->B!=null).
			map(B->B.toInterval()).
			filter(L->L.getStart()<L.getEnd()).
			forEach(L->{intervals.put(L,true); }
			);	
		return intervals;
		}
	finally
		{
		htsjdk.samtools.util.CloserUtil.close(r);
		}
	}


/** compile the javascript script. Can be either from JavascriptFile or JavascriptExpr */
protected javax.script.CompiledScript compileJavascript(
		final String jsExpression,
		final File jsFile
		) throws Exception
	{
	if( jsExpression!=null && jsFile!=null)
		{
		throw new RuntimeException("Both javascript expression and file defined.");
		}
	
	
	if(jsExpression==null && jsFile==null)
		{
		throw new RuntimeException("User error : Undefined script. Check your parameters.");
		}
		
	LOG.info("getting javascript manager");
	final javax.script.ScriptEngineManager manager = new javax.script.ScriptEngineManager();
	final javax.script.ScriptEngine engine = manager.getEngineByName("js");
	if(engine==null)
		{
		throw new RuntimeException("not available ScriptEngineManager: javascript. Use the SUN/Oracle JDK ?");
		}
	final javax.script.Compilable compilingEngine = (javax.script.Compilable)engine;
	if(jsFile!=null)
		{
		LOG.info("Compiling "+jsFile);
		java.io.FileReader r = null;
		try
			{
			r = new java.io.FileReader(jsFile);
			return compilingEngine.compile(r);
			}
		finally
			{
			htsjdk.samtools.util.CloserUtil.close(r);
			}
		}
	else if(jsExpression!=null)
		{
		LOG.info("Compiling "+jsExpression);
		return compilingEngine.compile(jsExpression);
		}
	else
		{
		throw new RuntimeException("illegal state");
		}
	}
/** END : JAVASCRIPT SECTION ************************************************/


public int doWork(final List<String> args)
	{
	return -1;
	}
public int instanceMain(final String args[]) {
	int ret=0;
	try 
		{
		final Status status = parseArgs(args);
		switch(status)
			{
			case EXIT_FAILURE: return -1;
			case EXIT_SUCCESS: return 0;
			case PRINT_HELP: getJCommander().usage(); return 0;
			case PRINT_VERSION: return 0;
			case OK:break;
			}
		
		try 
			{
			ret=initialize();
			if(ret!=0) return ret;
			}
		catch(final Throwable err)
			{
			LOG.severe(err.getMessage());
			return -1;
			}
		try 
			{
			ret=doWork(getFilenames());
			if(ret!=0) return ret;
			}
		catch(final Throwable err)
			{
			LOG.severe(err.getMessage());
			return -1;
			}
		}
	finally
		{
		cleanup();
		}
	return 0;
	}

public List<String> getFilenames() {
	return Collections.unmodifiableList(files);
	}
public List<File> getFiles() {
	return getFilenames().stream().
			map(S->new File(S)).
			collect(Collectors.toList());
	}

public PrintStream stdout() { return System.out;}
public PrintStream stderr() { return System.err;}
public InputStream stdin() { return System.in;}

/** open output (file or stdout) as PrintWriter */
protected java.io.PrintWriter openFileOrStdoutAsPrintWriter(File out) throws java.io.IOException
	{
	if(out!=null)
		{
		if(out.getName().endsWith(".gz"))
			{
			return new java.io.PrintWriter(this.openFileOrStdoutAsStream(out));
			}
		return new java.io.PrintWriter(out);
		}
	else
		{
		return new java.io.PrintWriter( stdout() );
		}
	}


/** open output (file or stdout) as PrintStream */
protected java.io.PrintStream openFileOrStdoutAsPrintStream(File out) throws java.io.IOException
	{
	if(out!=null)
		{
		if(out.getName().endsWith(".gz"))
			{
			final java.io.OutputStream os = this.openFileOrStdoutAsStream(out);
			if(os instanceof java.io.PrintStream) {
				return java.io.PrintStream.class.cast(os);
				}
			else
				{
				return new java.io.PrintStream(os);
				}
			}
		return new java.io.PrintStream(out);
		}
	else
		{
		return stdout();
		}
	}

/** open output (file or stdout) as OutputStream */
protected java.io.OutputStream openFileOrStdoutAsStream(final File out) throws java.io.IOException
	{
	if(out!=null)
		{
		return  IOUtils.openFileForWriting(out);
		}
	else
		{
		return stdout();
		}
	}

/** create a new SamReaderFactory */
protected htsjdk.samtools.SamReaderFactory createSamReaderFactory()
	{
	return  htsjdk.samtools.SamReaderFactory.makeDefault().validationStringency(htsjdk.samtools.ValidationStringency.LENIENT);
	}

/** open a new SAM reader; If inputName==null, it reads from stdin */
protected htsjdk.samtools.SamReader openSamReader(final String inputName)
	{
	final htsjdk.samtools.SamReaderFactory srf= this.createSamReaderFactory();
	if(inputName==null)
		{
		return srf.open(htsjdk.samtools.SamInputResource.of(stdin()));
		}
	else
		{
		return srf.open(htsjdk.samtools.SamInputResource.of(inputName));
		}
	}


protected static class StringToMd5 implements Function<String, String>
	{
	private final java.security.MessageDigest _md5;

	public StringToMd5() {
		try {
			_md5 = java.security.MessageDigest.getInstance("MD5");
		} catch (java.security.NoSuchAlgorithmException e) {
			throw new RuntimeException("MD5 algorithm not found", e);
		}}
		
	@Override
	public String apply(String in) {
		_md5.reset();
		_md5.update(in.getBytes());
		String s = new java.math.BigInteger(1, _md5.digest()).toString(16);
		if (s.length() != 32) {
			final String zeros = "00000000000000000000000000000000";
			s = zeros.substring(0, 32 - s.length()) + s;
		}
		return s;
		}
	}


/** just created to make a transition between XML and Jcommander. Remove in the future */
@Deprecated
protected int wrapException(final Object msg) 
	{
	LOG.error(msg);
	return -1;
	}

/** just created to make a transition with old programs. Remove in the future */
@Deprecated
protected String getMessageBundle(final String s){
	return String.valueOf(s);
	}

public void instanceMainWithExit( final String args[]) {
	System.exit( instanceMain(args) );
	}

}
