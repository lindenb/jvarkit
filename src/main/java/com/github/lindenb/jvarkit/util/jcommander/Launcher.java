package com.github.lindenb.jvarkit.util.jcommander;

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
import java.util.stream.Collectors;
import java.util.zip.Deflater;

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
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
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
		}
		super.usage(sb);
		
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
@Parameter(names = {"--version"}, help = true,description="print version and exits")
private boolean print_version = false;
@Parameter(description = "Files")
private List<String> files = new ArrayList<>();



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


public static class WritingBamArgs
	{
	@Parameter(names={"--xx"},description="Compression Level.",converter=CompressionConverter.class)
	public int compressionLevel=5;
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


public Launcher()
	{
	
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
	 final Map<Class, Class<? extends IStringConverter<?>>> MAP = new HashMap() {{
		    put(VcfWriterOnDemand.class, VcfWriterOnDemandConverter.class);
		    put(VariantContextWriter.class, VcfWriterOnDemandConverter.class);
		}};	
	this.jcommander.addConverterFactory(new IStringConverterFactory() {
			@Override
			public Class<? extends IStringConverter<?>> getConverter(Class forType) {		
				return MAP.get(forType);
				}
			});
	this.jcommander.addObject(this);	
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
	 
	 if (this.print_help) return Status.PRINT_HELP;
	 if (this.print_version) return Status.PRINT_VERSION;
	 return Status.OK;
	}

protected VcfIterator openVcfIterator(final String inputNameOrNull) throws IOException {
	return VCFUtils.createVcfIterator(inputNameOrNull);
}

protected VariantContextWriter openVariantContextWriter(final File outorNull) throws IOException {
	return VCFUtils.createVariantContextWriter(outorNull);
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

/** just created to make a transition between XML and Jcommander. Remove in the future */
@Deprecated
protected int wrapException(final Object msg) 
	{
	LOG.error(msg);
	return -1;
	}


public void instanceMainWithExit( final String args[]) {
	System.exit( instanceMain(args) );
	}

}
