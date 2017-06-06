/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.util.jcommander;

import java.awt.Dimension;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.StringWriter;
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
import java.util.Random;
import java.util.function.Function;
import java.util.jar.Manifest;
import java.util.stream.Collectors;
import java.util.zip.Deflater;

import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.w3c.dom.Document;
import org.w3c.dom.Element;

import com.beust.jcommander.IStringConverter;
import com.beust.jcommander.IStringConverterFactory;
import com.beust.jcommander.IValueValidator;
import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParameterDescription;
import com.beust.jcommander.ParameterException;
import com.beust.jcommander.ParametersDelegate;
import com.beust.jcommander.converters.IntegerConverter;
import com.github.lindenb.jvarkit.annotproc.IncludeSourceInJar;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.bio.samfilter.SamFilterParser;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.semontology.Term;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;


@IncludeSourceInJar
public class Launcher {
private static final Logger LOG=Logger.build( Launcher.class).make();
public static final String OPT_OUPUT_FILE_OR_STDOUT="Output file. Optional . Default: stdout";
public static final String INDEXED_FASTA_REFERENCE_DESCRIPTION="Indexed fasta Reference file. "+
		"This file must be indexed with samtools faidx and with picard CreateSequenceDictionary";

protected static final int RETURN_OK=0;
public enum Status { OK, PRINT_HELP,PRINT_VERSION,EXIT_SUCCESS,EXIT_FAILURE};

/** need to decouple from Launcher for JXF applications that cannot extends 'Launcher' */
public static  class UsageBuider
	{
	/** git hash in the manifest */
	private String gitHash = null;
	/** compile date in the manifest */
	private String compileDate = null;
	/** main class */
	private Class<?> mainClass=Object.class;
	
	@Parameter(names = {"-h","--help"},description="print help and exit", help = true)
	public boolean print_help = false;
	@Parameter(names = {"--helpFormat"},description="What kind of help", help = true)
	public HelpFormat helpFormat = HelpFormat.usage;
	@Parameter(names = {"--version"}, help = true,description="print version and exit")
	public boolean print_version = false;

	private enum HelpFormat {usage,markdown,xml}
	
	public UsageBuider(Class<?> mainClass) {
		this.mainClass = mainClass;
	}
	
	public boolean shouldPrintUsage() {
		return this.print_help;
	}
	
	public void usage(final JCommander jc) {
		if(this.helpFormat.equals(HelpFormat.xml))
			{
			final Document dom;
			try {
				dom = this.xmlUsage(jc);
				}
			catch(final Exception err)
				{
				err.printStackTrace();
				System.err.println("An error occured. Cannot produce XML");
				return ;
				}

			try {
				final Transformer transformer = TransformerFactory.newInstance().newTransformer();

				final DOMSource source = new DOMSource(dom);
				final StreamResult result = new StreamResult(new StringWriter());

				if(!this.helpFormat.equals(HelpFormat.xml))
					{
					System.err.println("TODO");
					}
				else
					{
					transformer.setOutputProperty(OutputKeys.STANDALONE,"yes");
					transformer.transform(source, result);
					final String xmlString = result.getWriter().toString();
					System.out.println(xmlString);
					}
				}
			catch(final Exception err)
				{
				err.printStackTrace();
				System.err.println("An error occured. Cannot produce XML");
				}
			}
		else
			{
			System.out.println(getUsage(jc));
			}
		}

	
	public String getUsage(final JCommander jc) {
		final StringBuilder sb=new StringBuilder();
		this.usage(jc,sb);
		return sb.toString();
		}
	public String hyperlink(final String url)
		{
		return "["+url+"]("+url+")";
		}
	
	private void include(final StringBuilder sb,String className) {
	InputStream in=null;
	try {
		int dollar=className.indexOf('$');
		if(dollar!=-1) className=className.substring(0, dollar);
		className=className.replace('.', '/')+".java";
		in=getMainClass().getResourceAsStream("/"+className);
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
					if(!ok) LOG.warn("END_"+"DOC without BEGIN");
					ok=false;
					}
				else if(ok)
					{
					if(line.trim().startsWith("@@INCLUDE"))
						{
						int n=line.indexOf(" ");
						if(n==-1)  n=line.indexOf("\t");
						if(n==-1)  n=line.indexOf("=");
						if(n!=-1) {
							line=line.substring(n+1).trim();
							if(!line.isEmpty())
								{
								include(sb, line);
								}
							}
						}
					else
						{
						sb.append(line).append("\n");
						}
					}
				}
			r.close();
			if(ok) LOG.warn("BEGIN_"+"DOC without END");
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
	
	/** create usage as XML. No DTD. Might change in the future */
	private Document xmlUsage(final JCommander jc) throws Exception {
		final Class<?> clazz = getMainClass();			
		final Program programdesc = clazz.getAnnotation(Program.class);

		final Document dom = DocumentBuilderFactory.newInstance().
					newDocumentBuilder().
					newDocument();
		final Element root = dom.createElement("program");
		Element e1;
		dom.appendChild(root);
		if(programdesc!=null)
			{
			root.setAttribute("name", programdesc.name());
			e1 = dom.createElement("description");
			root.appendChild(e1);
			e1.appendChild(dom.createTextNode(programdesc.description()));
			for(final String kw : programdesc.keywords())
				{
				e1 = dom.createElement("keyword");
				root.appendChild(e1);
				e1.appendChild(dom.createTextNode(kw));
				}
			for(final Term t : programdesc.terms())
				{
				e1 = dom.createElement("term");
				root.appendChild(e1);
				e1.setAttribute("id", t.getAccession());
				e1.appendChild(dom.createTextNode(t.getLabel()));
				}
			for(final int bid : programdesc.biostars())
				{
				e1 = dom.createElement("biostar");
				root.appendChild(e1);
				e1.appendChild(dom.createTextNode("https://www.biostars.org/p/"+bid));
				}
			}
		e1 = dom.createElement("parameters");
		root.appendChild(e1);
		for(final ParameterDescription pd:jc.getParameters())
			{
			final Element e2= dom.createElement("parameter");
			e2.setAttribute("help",String.valueOf(pd.isHelp()));
			e1.appendChild(e2);
			for(final String name : pd.getParameter().names())
				{
				final Element e3 = dom.createElement("name");
				e2.appendChild(e3);
				e3.appendChild(dom.createTextNode(name));
				}
			final Element e3 = dom.createElement("description");
			e2.appendChild(e3);
			e3.appendChild(dom.createTextNode(pd.getDescription()));
			}
		
	return dom;
	}	
	
	
	public void usage(final JCommander jc,final StringBuilder sb) {

		final Class<?> clazz = getMainClass();
		
		final Program programdesc = clazz.getAnnotation(Program.class);
		
		if(programdesc!=null){
			jc.setProgramName(programdesc.name());
		} else
			{
			jc.setProgramName(clazz.getSimpleName());
			}
		if(this.helpFormat.equals(HelpFormat.markdown))  
			{
			sb.append("# "+clazz.getSimpleName()+"\n\n");
			
			if(programdesc!=null){
				
				sb.append(programdesc.description()).
					append("\n\n");
				
				if(!programdesc.deprecatedMsg().isEmpty())
					{
					sb.append("\n## DEPRECATED\n\n").
						append(programdesc.deprecatedMsg()).
						append("\n");
					}
				}
			}
		
		
		if(this.helpFormat.equals(HelpFormat.markdown)) sb.append("\n## Usage\n\n```\n");
		jc.usage(sb);
		if(this.helpFormat.equals(HelpFormat.markdown)) sb.append("\n```\n\n");

		if(programdesc!=null && this.helpFormat.equals(HelpFormat.markdown)){
			if(programdesc.keywords()!=null && programdesc.keywords().length>0) {
				sb.append("\n## Keywords\n\n");
				for(String sk:programdesc.keywords()) sb.append(" * "+sk+"\n");
				sb.append("\n\n");
			}
			if(programdesc.biostars()!=null && programdesc.biostars().length>0) {
				sb.append("\n## See also in Biostars\n\n");
				for(int postid:programdesc.biostars()) sb.append(" * "+hyperlink("https://www.biostars.org/p/"+postid)+"\n");
				sb.append("\n\n");
			}	
		}
		
		if(this.helpFormat.equals(HelpFormat.markdown))
			{
			final String progName=(programdesc==null?"software":programdesc.name());
			sb.append("## Compilation\n");
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
			sb.append(hyperlink("https://github.com/lindenb/jvarkit/tree/master/src/main/java/"+
				clazz.getName().replace('.','/')+".java\n"));
			sb.append("\n");
			
			sb.append("## Contribute\n");
			sb.append("\n");
			sb.append("- Issue Tracker: "+hyperlink("http://github.com/lindenb/jvarkit/issues")+"\n");
			sb.append("- Source Code: "+hyperlink("http://github.com/lindenb/jvarkit")+"\n");
			sb.append("\n");
			sb.append("## License\n");
			sb.append("\n");
			sb.append("The project is licensed under the MIT license.\n");
			sb.append("\n");
			sb.append("## Citing\n");
			sb.append("\n");
			sb.append("Should you cite **"+progName +"** ? "+hyperlink("https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md")+"\n");
			sb.append("\n");
			sb.append("The current reference is:\n");
			sb.append("\n");
			sb.append(hyperlink("http://dx.doi.org/10.6084/m9.figshare.1425030")+"\n");
			sb.append("\n");
			sb.append("> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.\n");
			sb.append("> "+hyperlink("http://dx.doi.org/10.6084/m9.figshare.1425030")+"\n");
			sb.append("\n");
			}
		if( this.helpFormat.equals(HelpFormat.markdown) ) {
			include(sb,clazz.getName());
			}
		}
	
	public Class<?> getMainClass()
		{
		return mainClass;
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
			final Enumeration<URL> resources = getMainClass().getClassLoader()
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

	}

@ParametersDelegate
private UsageBuider usageBuilder = null;

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
	@Parameter(names={"--xxx"},description="Compression Level.",converter=CompressionConverter.class)
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


public static class RandomConverter
implements IStringConverter<Random>
{
public static Random now() {
	return new Random(System.currentTimeMillis());
}

@Override
public Random convert(final String v) {	
	if(v==null || v.isEmpty() || v.equals("-1") || v.equals("now") || v.equals("timestamp"))
		{
		return now();
		}
	else
		{
		try {
			final long t = Long.parseLong(v);
			return new Random(t);
			}
		catch(NumberFormatException err)
			{
			throw new ParameterException(err);
			}
		}
	
	}
}

public static class IndexedFastaSequenceFileConverter
implements IStringConverter<IndexedFastaSequenceFile>
	{
	@Override
	public IndexedFastaSequenceFile convert(final String path) {
		if(path==null) return null;
		File faidx=new File(path);
		if(!faidx.exists()) throw new ParameterException("file doesn't exists: "+path);
		try {
			return new IndexedFastaSequenceFile(faidx);
		} catch (FileNotFoundException err) {
			throw new ParameterException(err);
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
	final Class<?> clazz=Launcher.this.getClass();
	this.usageBuilder = new UsageBuider(clazz);
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
	
	 @SuppressWarnings({"rawtypes","unchecked","serial"})
	final Map<Class, Class<? extends IStringConverter<?>>> MAP = new HashMap() {{
		    put(Dimension.class,DimensionConverter.class);
		    put(SamRecordFilter.class,SamFilterParser.StringConverter.class);
		    put(Random.class,RandomConverter.class);
		}};	
	this.jcommander.addConverterFactory(new IStringConverterFactory() {
			@SuppressWarnings("unchecked")
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
		stderr().println("There was an error in the input parameters.");
		stderr().println(err.getMessage());
		return Status.EXIT_FAILURE; 
	 	}
	 
	 if (this.usageBuilder.shouldPrintUsage()) return Status.PRINT_HELP;
	 if (this.usageBuilder.print_version) return Status.PRINT_VERSION;
	 return Status.OK;
	}

protected VcfIterator openVcfIterator(final String inputNameOrNull) throws IOException {
	return VCFUtils.createVcfIterator(inputNameOrNull);
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
	try {
		iterin = openVcfIterator(inputNameOrNull);
		w = openVariantContextWriter(outorNull);
		int ret=doVcfToVcf(inputNameOrNull==null?"<STDIN>":inputNameOrNull,iterin,w);
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
/** extract case controls in VCF header injected with VcfInjectPedigree */
protected java.util.Set<com.github.lindenb.jvarkit.util.Pedigree.Person> getCasesControlsInPedigree(final htsjdk.variant.vcf.VCFHeader header) {
	final com.github.lindenb.jvarkit.util.Pedigree pedigree = Pedigree.newParser().parse(header);
	if(pedigree.isEmpty())
		{
		throw new IllegalArgumentException("No pedigree found in header. use VcfInjectPedigree to add it");
		}
	if(!pedigree.verifyPersonsHaveUniqueNames()) {
		throw new IllegalArgumentException("I can't use this pedigree in VCF because two samples have the same ID.");
	}

	final java.util.Set<String> samplesNames= new java.util.HashSet<>(header.getSampleNamesInOrder());
	final java.util.Set<com.github.lindenb.jvarkit.util.Pedigree.Person> individuals = new java.util.HashSet<>(pedigree.getPersons());
	
	
	final java.util.Iterator<com.github.lindenb.jvarkit.util.Pedigree.Person> iter= individuals.iterator();
	while(iter.hasNext())
	{
		final com.github.lindenb.jvarkit.util.Pedigree.Person person = iter.next();
		if(!(samplesNames.contains(person.getId()) && (person.isAffected() || person.isUnaffected()))) {
			LOG.warn("Ignoring "+person+" because it is not present in VCF header or status is unknown");
			iter.remove();
		}
	}
	
	LOG.info("Individuals :"+individuals.size() +
		" affected :"+individuals.stream().filter(P->P.isAffected()).count() +
		" unaffected :"+individuals.stream().filter(P->P.isUnaffected()).count()
		);

	return java.util.Collections.unmodifiableSet( individuals );
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
		final com.github.lindenb.jvarkit.util.vcf.VcfIterator in = com.github.lindenb.jvarkit.util.vcf.VCFUtils.createVcfIteratorFromLineIterator(lineIter, true);
		
			
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
