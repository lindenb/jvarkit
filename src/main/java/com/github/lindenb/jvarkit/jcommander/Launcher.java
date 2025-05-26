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
package com.github.lindenb.jvarkit.jcommander;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.zip.Deflater;

import com.beust.jcommander.IStringConverter;
import com.beust.jcommander.IStringConverterFactory;
import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.bed.BedLineCodec;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.bio.samfilter.SamRecordFilterFactory;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.AbstractProgressLogger;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.ProgressLoggerInterface;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;



public class Launcher {
private static final Logger LOG=Logger.of( Launcher.class);
public static final String OPT_OUPUT_FILE_OR_STDOUT="Output file. Optional . Default: stdout";
public static final String INDEXED_FASTA_REFERENCE_DESCRIPTION="Indexed fasta Reference file. "+
		"This file must be indexed with samtools faidx and with picard/gatk CreateSequenceDictionary or samtools dict";
public static final String DICTIONARY_SOURCE="A SAM Sequence dictionary source: it can be a *.dict file, a fasta file indexed with 'picard CreateSequenceDictionary' or 'samtools dict', or any hts file containing a dictionary (VCF, BAM, CRAM, intervals...)";

public static final String CRAM_INDEXED_REFENCE="For reading/writing CRAM files. "+ INDEXED_FASTA_REFERENCE_DESCRIPTION;


/** description used when building custom URLs (e.g: VcfServer ) */
public static final String USER_CUSTOM_INTERVAL_URL_DESC="A custom URL for a web browser. The following words will be replaced by their values: ${CHROM}, ${START}, ${END}. "
		+ "For example for IGV that would be: 'http://localhost:60151/goto?locus=${CHROM}%3A${START}-${END}' (see http://software.broadinstitute.org/software/igv/book/export/html/189)";

protected static final int RETURN_OK=0;
public enum Status { OK, PRINT_HELP,PRINT_VERSION,EXIT_SUCCESS,EXIT_FAILURE};



@ParametersDelegate
private CmdUsageBuilder usageBuilder = null;

/** custom instance of jcommander, don't add same command twice. */
private class MyJCommander extends JCommander
	{
	/** when registering the option for jcommander, we take care of not visiting object twice */
	private Collection<Object> ojectsVisitedByJCommander = new ArrayList<>();

	
	
	@Override
	public void addCommand(String name, Object object, String... aliases) 
			{
			if(this.ojectsVisitedByJCommander.stream().anyMatch(O->O==object)) return;
			this.ojectsVisitedByJCommander.add(object);
			super.addCommand(name, object, aliases);
			}
	}
	





/** original arc/argv */
private List<String> argcargv=Collections.emptyList();
private final JCommander jcommander = new MyJCommander();

@Parameter(description = "Files")
private List<String> files = new ArrayList<>();

private String programName="";


public class WritingSortingCollection
	{
	@Parameter(names={"--maxRecordsInRam"},description="When writing  files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file  handles needed to sort a file, and increases the amount of RAM needed")
	public int maxRecordsInRam=50_000;
	
	@Parameter(names={"--tmpDir"},description= "tmp working directory. Default: java.io.tmpDir")
	private List<File> tmpDirs=new ArrayList<>();
	
	
	public WritingSortingCollection maxRecordsInRam(final int n)
		{
		this.maxRecordsInRam = n;
		return this;
		}
	public int getMaxRecordsInRam() { return this.maxRecordsInRam;}
	public List<File> getTmpDirectories() {
		final List<File> L= new ArrayList<>(this.tmpDirs);
		if(L.isEmpty() )
			{
			L.add(IOUtils.getDefaultTmpDir());
			}
		return L;
		}
	/** convert getTmpDirectories to Array of Path */
	public Path[] getTmpPaths() {
		return getTmpDirectories().stream().
				map(F->F.toPath()).
				toArray((i)->new Path[i]);
		}
	/** get One path */
	public Path getTmpPath() {
		return getTmpPaths()[0];
		}
	}	

public static enum WritingSamReaderType
	{
	BAM,SAM,CRAM
	}

public class WritingBamArgs
	{
	private Path referenceFile = null;
	private ProgressLoggerInterface progressLogger = null;
	
	@Parameter(names={"--bamcompression"},description="Compression Level. 0: no compression. 9: max compression;")
	public int compressionLevel=5;
	@Parameter(names={"--samoutputformat"},description="Sam output format.")
	public WritingSamReaderType samoutputformat = WritingSamReaderType.SAM;
	
	/** creates a SAMFileWriterFactory */
	public htsjdk.samtools.SAMFileWriterFactory createSAMFileWriterFactory() {
		final SAMFileWriterFactory sfw =  new SAMFileWriterFactory();
		int n= this.compressionLevel;
		if(n<0) n= Deflater.NO_COMPRESSION;
		if(n>9) n= Deflater.BEST_COMPRESSION;
		sfw.setCompressionLevel(n);
		return sfw;
		}
	
	public WritingBamArgs setProgressLogger(final Logger log) {
		if(log!=null) {
			return setProgressLogger(new AbstractProgressLogger(getProgramName(),"Writer",1_000_000)
				{
				@Override
				protected void log(String... message) {
					log.info(String.join(" ", message));
					}
					});
				}
		return this;
		}
	
	public WritingBamArgs setProgressLogger(final ProgressLoggerInterface progressLogger) {
		this.progressLogger = progressLogger;
		return this;
	}
	
	public WritingBamArgs setSamOutputFormat(WritingSamReaderType samoutputformat) {
		this.samoutputformat = samoutputformat;
		return this;
		}
	
	@Deprecated
	public WritingBamArgs setReferenceFile(final File referenceFile) {
		this.referenceFile = (referenceFile==null?null:referenceFile.toPath());
		return this;
		}
	
	public WritingBamArgs setReferencePath(final Path referenceFile) {
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
	/** return reference file as java.io.File. default implementation returns null */
	@Deprecated
	public File getReferenceFile() {
		final Path p=getReferencePath();
		return p==null?null:p.toFile();
		}
	/** return reference file as java.nio.Path. default implementation returns null */
	public Path getReferencePath() {
		return this.referenceFile;
		}
	
	/** get File extension for creating output file , including the dot */
	public String getFileExtension() {
			if(this.samoutputformat==null) return FileExtensions.SAM;
			switch(this.samoutputformat) {
				case BAM: return FileExtensions.BAM;// ".bam";
				case CRAM: return FileExtensions.CRAM;
				default: return FileExtensions.SAM;
			}
		}
	
	public SAMFileWriter openSamWriter(final Path outputFileOrNull,SAMFileHeader header,boolean presorted) {
		final htsjdk.samtools.SAMFileWriterFactory sfw= this.createSAMFileWriterFactory();
		
		final SAMFileWriter w;
		if(outputFileOrNull==null)
			{
			if( this.samoutputformat!=null &&
				this.samoutputformat.equals(WritingSamReaderType.BAM))
				{
				w = sfw.makeBAMWriter(header, presorted, stdout());
				}
			else if(this.samoutputformat!=null &&
					this.samoutputformat.equals(WritingSamReaderType.CRAM))
				{
				w = sfw.makeCRAMWriter(header,stdout(),getReferencePath());
				}
			else if(this.samoutputformat==null || this.samoutputformat.equals(WritingSamReaderType.SAM))
				{
				w = sfw.makeSAMWriter(header, presorted, stdout());
				}
			else
				{
				throw new IllegalStateException("Bad output format "+this.samoutputformat+" expected one of "+Arrays.toString(WritingSamReaderType.values()));
				}
			}
		else
			{
			w =  sfw.makeWriter(header, presorted, outputFileOrNull, getReferencePath());
			}
		if(this.progressLogger!=null) {
			w.setProgressLogger(this.progressLogger);
			}
		return w;
		}
	
	public SAMFileWriter openSAMFileWriter(File outputFileOrNull,SAMFileHeader header,boolean presorted) {
		return this.openSamWriter(outputFileOrNull==null?null:outputFileOrNull.toPath(),header,presorted);
		}
	}


@SuppressWarnings("unchecked")
public Launcher()
	{
	final Class<?> clazz=Launcher.this.getClass();
	this.usageBuilder = new CmdUsageBuilder(clazz);
	if(this.usageBuilder.hasProgram())
		{
		this.programName = this.usageBuilder.getProgram().name();
		}
	else
		{
		this.programName = getClass().getSimpleName();
		}	
	
	try {
		/** set locale http://seqanswers.com/forums/showthread.php?p=174020#post174020 */
		Locale.setDefault(Locale.US);
		}
	catch(final Throwable /*java.security.AccessControlException deprecated*/ err) {
		System.err.println("Cannot set Locale to US for security reasons"+err.getMessage());
		}
	try {
		System.setProperty("file.encoding", "UTF-8");
		}
	catch(final Throwable /*java.security.AccessControlException deprecated*/ err) {
		System.err.println("Cannot set file.encoding to UTF-8 for security reasons : "+err.getMessage());
		}
	try {
	/* https://bugs.openjdk.java.net/browse/JDK-8028111 */
	System.setProperty("jdk.xml.entityExpansionLimit","0");
	}
	catch(final Throwable /*java.security.AccessControlException deprecated*/ err) {
	}
	
	@SuppressWarnings({"rawtypes","serial"})
	final Map<Class, Class<? extends IStringConverter<?>>> MAP = new HashMap() {{
		    put(SamRecordFilter.class,SamRecordFilterFactory.class);
		}};	

	this.jcommander.addConverterFactory(new IStringConverterFactory() {
			@Override
			public Class<? extends IStringConverter<?>> getConverter(@SuppressWarnings("rawtypes") Class forType) {		
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
	return this.usageBuilder.getCompileDate();
	}

public String getGitHash()
	{
	return this.usageBuilder.getGitHash();
	}

public String getVersion()
	{
	return this.usageBuilder.getVersion();
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
	return String.join(" ",this.getRawArguments());
	}

/** get the original args on the command line, before jcommander processnig */
public List<String> getRawArguments() {
	return this.argcargv;
	}

protected Status parseArgs(final String args[])
	{
	this.argcargv = Arrays.asList(args);
	try
	  	{
		getJCommander().parse(args);
	  	}
	 catch(final com.beust.jcommander.ParameterException err) {
		stderr().print("There was an error in the input parameters: ");
		stderr().println(err.getMessage());
		return Status.EXIT_FAILURE;
	 	}
	
	 if (this.usageBuilder.shouldPrintUsage()) return Status.PRINT_HELP;
	 if (this.usageBuilder.print_version) return Status.PRINT_VERSION;
	 return Status.OK;
	}

protected VCFIterator openVCFIterator(final String inputNameOrNull) throws IOException {
	return VCFUtils.createVCFIterator(inputNameOrNull);
}

/**
 * dict will be using for sorting or indexing. May be null it no indexing.
*/
protected VariantContextWriter openVariantContextWriter(final File outorNull) throws IOException {
	if( outorNull == null)
		{
		return VCFUtils.createVariantContextWriterToOutputStream(stdout());
		}
	else
		{
		return VCFUtils.createVariantContextWriter(outorNull);
		}
	}


protected InputStream openInputStream(final String inOrNull) throws IOException {
	return(inOrNull==null || inOrNull.equals("-")?
			stdin():
			IOUtils.openURIForReading(inOrNull)
			);
}

protected BufferedReader openBufferedReader(final String inOrNull) throws IOException {
	return(inOrNull==null || inOrNull.equals("-")?
			new BufferedReader(new InputStreamReader(stdin())):
			IOUtils.openURIForBufferedReading(inOrNull)
			);
	}


protected VCFHeader addMetaData(final VCFHeader header) 
	{
	return header;
	}

protected int doVcfToVcf(final String inputName,final VCFIterator iterin,final VariantContextWriter out){
	LOG.debug("using default doVcfToVcf ??");
	VCFUtils.copyHeaderAndVariantsTo(iterin, out);
	return 0;
	}



protected int doVcfToVcf(final String inputNameOrNull,final File outorNull){
	try (VCFIterator iterin = openVCFIterator(inputNameOrNull)) {
		try(VariantContextWriter w = openVariantContextWriter(outorNull) ) {
			int ret=doVcfToVcf(inputNameOrNull==null?"<STDIN>":inputNameOrNull,iterin,w);
			return ret;
			}
		}
	catch(final Exception err)
		{
		LOG.error(err);
		return -1;
		}
	}


protected int doVcfToVcfPath(final String inputNameOrNull,final WritingVariantsDelegate delegate,final Path outorNull){
	try(VCFIterator iterin = openVCFIterator(inputNameOrNull)) {
		try(VariantContextWriter w = delegate.dictionary(iterin.getHeader()).open(outorNull)) {
			int ret=doVcfToVcf(inputNameOrNull==null?"<STDIN>":inputNameOrNull,iterin,w);
			return ret;
			}
		}
	catch(final Exception err)
		{
		LOG.error(err);
		return -1;
		}
	}


protected int doVcfToVcf(final List<String> inputs,final File outorNull) {
	return doVcfToVcf(oneFileOrNull(inputs),outorNull);
	}

protected int doVcfToVcfPath(final List<String> inputs,final WritingVariantsDelegate delegate, Path outorNull) {
	return doVcfToVcfPath(oneFileOrNull(inputs),delegate, outorNull);
	}


protected String oneAndOnlyOneFile(final List<String> args) {
	switch(args.size())
		{
		case 1: return args.get(0);
		default: throw new JvarkitException.CommandLineError("Expected one and only one argument but got "+args.size()+" : "+args.toString());
		}
}


protected String oneFileOrNull(final List<String> args) {
	switch(args.size())
	{
	case 0: return null;
	case 1: return args.get(0);
	default: throw new JvarkitException.CommandLineError("Expected one or zero argument but got "+args.size()+" : "+args.toString());
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
			filter(line->!(line.startsWith("#") ||  com.github.lindenb.jvarkit.bed.BedLine.isBedHeader(line) ||  line.isEmpty())).
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
		throw new JvarkitException.JavaScriptEngineNotFound();
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

/** instanceMain with Collection to make test-unit easier to write , just call `args.toArray` */
public final int instanceMain(final List<String> args) {
	return instanceMain(args.toArray(new String[args.size()]));
	}
	
public int instanceMain(final String args[]) {
	int ret=RETURN_OK;
	try 
		{
		final Status status = parseArgs(args);
		switch(status)
			{
			case EXIT_FAILURE: return -1;
			case EXIT_SUCCESS: return 0;
			case PRINT_HELP: this.usageBuilder.usage(getJCommander()); return 0;
			case PRINT_VERSION: System.out.println(getVersion());return 0;
			case OK:break;
			}
		
		try 
			{
			ret = initialize();
			if(ret!=RETURN_OK) return ret;
			}
		catch(final Throwable err)
			{
			LOG.severe(err.getMessage());
			return -1;
			}
		try 
			{
			ret=doWork(getFilenames());
			if(ret!=RETURN_OK) return ret;
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

private PrintStream _stdout = System.out;
public PrintStream stdout() { return _stdout;}
public PrintStream stdout(final PrintStream out) {PrintStream old=_stdout; this._stdout=out; return old;}
private PrintStream _stderr = System.err;

public PrintStream stderr() { return _stderr;}
public PrintStream stderr(final PrintStream out) {PrintStream old=_stderr; this._stderr=out; return old;}

private InputStream _stdin = System.in;
public InputStream stdin() { return _stdin;}
public InputStream stdin(final InputStream in) {InputStream old=_stdin; this._stdin=in; return old;}



/** open output (file or stdout) as PrintWriter */
protected java.io.PrintWriter openFileOrStdoutAsPrintWriter(final File out) throws java.io.IOException
	{
	return this.openPathOrStdoutAsPrintWriter(out==null?null:out.toPath());
	}

/** open output (path or stdout) as PrintWriter */
protected java.io.PrintWriter openPathOrStdoutAsPrintWriter(final Path out) throws java.io.IOException
	{
	if(out==null || out.toString().equals("-"))
		{
		return new java.io.PrintWriter( stdout() );
		}
	else
		{
		return IOUtils.openPathForPrintWriter(out);
		}
	}


/** open output (file or stdout if out is null ) as PrintStream */
protected java.io.PrintStream openFileOrStdoutAsPrintStream(final File out) throws java.io.IOException
	{
	return this.openPathOrStdoutAsPrintStream(out==null?null:out.toPath());
	}

/** open output (file or stdout if out is null ) as PrintStream */
protected java.io.PrintStream openPathOrStdoutAsPrintStream(final Path out) throws java.io.IOException
	{
	if(out==null || out.toString().equals("-"))
		{
		return stdout();
		}
	else
		{
		return new PrintStream(IOUtils.openPathForWriting(out));
		}
	}


/** open output (file or stdout) as OutputStream */
protected java.io.OutputStream openFileOrStdoutAsStream(final File out) throws java.io.IOException
	{
	return openPathOrStdoutAsStream(out==null?null:out.toPath());
	}

/** open output (file or stdout) as OutputStream */
protected java.io.OutputStream openPathOrStdoutAsStream(final Path out) throws java.io.IOException
	{
	if(out==null || out.toString().equals("-"))
		{
		return stdout();
		}
	else
		{
		return  IOUtils.openPathForWriting(out);
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



/** extract case controls in VCF header injected with VcfInjectPedigree */
protected java.util.Set<com.github.lindenb.jvarkit.util.Pedigree.Person> getCasesControlsInPedigree(final htsjdk.variant.vcf.VCFHeader header) {
	return new com.github.lindenb.jvarkit.util.Pedigree.CaseControlExtractor().extract(header);
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
	final int ret = instanceMain(args);

	if(ret!=RETURN_OK) {
		LOG.info(getProgramName()+" Exited with failure ("+ret+")");
		}
	System.exit( ret );
	}

/** check we're running under webstart http://stackoverflow.com/questions/216315*/
public static boolean isRunningJavaWebStart() {
    try {
      Class.forName("javax.jnlp.ServiceManager"); // or System.getProperty("javawebstart.version", null) != null;
      return true;
    	} 
    catch (ClassNotFoundException ex) {
       return false;
    	}
	}


/** build a cutsom URL from an interval. see {@link #USER_CUSTOM_INTERVAL_URL_DESC}. 
 * @return the new url of null if pattern is blank or interval is null */
public static String createUrlFromInterval(final String pattern,final Interval interval)
	{
	if(StringUtil.isBlank(pattern) || interval==null) return null;
	return pattern.
		replaceAll(Pattern.quote("${CHROM}"),StringUtils.escapeHttp(interval.getContig())).
		replaceAll(Pattern.quote("${START}"),String.valueOf(interval.getStart()).
		replaceAll(Pattern.quote("${END}"),String.valueOf(interval.getEnd())))
		;
	}

/** when using doVcfToVcfMultipleStream, this will be the name of the zip */
protected String getVcfMultipleStreamZipDirectory() {
	return getProgramName();
	}

/** can be called instead of doVcfVcf to hande a stream of multiple vcfs */
protected int doVcfToVcfMultipleStream(
		final String inputName,final File outputFile)
	{
	java.io.FileOutputStream fout = null;
	java.util.zip.ZipOutputStream zout = null;
	htsjdk.tribble.readers.LineIterator lineIter = null;
	java.io.PrintStream pw = null;
	htsjdk.variant.variantcontext.writer.VariantContextWriter vcw=null;
	try {
		if(inputName==null) {
			lineIter = com.github.lindenb.jvarkit.io.IOUtils.openStreamForLineIterator(stdin());
		} else {
			lineIter = com.github.lindenb.jvarkit.io.IOUtils.openURIForLineIterator(inputName);
		}
		if(!lineIter.hasNext()) {
			LOG.warn("No input found. Couldn't read any VCF header file");
		}
		
		if(outputFile!=null && outputFile.getName().endsWith(".zip")) {
			fout = new java.io.FileOutputStream(outputFile);
			zout = new java.util.zip.ZipOutputStream(fout);
		} else {
			pw = openFileOrStdoutAsPrintStream(outputFile);
		}
		
		int n_vcf=0;
		while(lineIter.hasNext()) {
		final String filename = String.format("%04d.vcf",(++n_vcf));
			
		/* create VCF iterator */
		final htsjdk.variant.vcf.VCFIterator in = com.github.lindenb.jvarkit.util.vcf.VCFUtils.createVCFIteratorFromLineIterator(lineIter, true);
		
			
		/* if zip: add new entry */
		if(zout!=null) {
			final java.util.zip.ZipEntry entry = new java.util.zip.ZipEntry(getVcfMultipleStreamZipDirectory()+"/"+filename);
			entry.setComment("Generated with "+getProgramName()+" v."+getVersion());
			zout.putNextEntry(entry);
			pw = new java.io.PrintStream(zout);
		}
		
		/* create VariantWriter */
		vcw = com.github.lindenb.jvarkit.util.vcf.VCFUtils.createVariantContextWriterToOutputStream(
		  com.github.lindenb.jvarkit.io.IOUtils.uncloseableOutputStream(pw));
		
		int errors = doVcfToVcf(filename, in, vcw);
		if(errors!=0) {
			return errors;
		}
		in.close();
		vcw.close();
		vcw=null;
		if(zout!=null)   {
			zout.closeEntry();
			pw=null;
		}
		}			
		
		if(zout!=null) {
			zout.finish();
			zout.flush(); zout=null;
			fout.close(); fout=null;
		}
		else {
			pw.flush();
			}
		return RETURN_OK;
	} catch(final Exception err) {
		LOG.error(err);
		return -1;
	} finally {
		htsjdk.samtools.util.CloserUtil.close(zout);
		htsjdk.samtools.util.CloserUtil.close(fout);
		htsjdk.samtools.util.CloserUtil.close(lineIter);
		htsjdk.samtools.util.CloserUtil.close(pw);
	}
	
}



}
