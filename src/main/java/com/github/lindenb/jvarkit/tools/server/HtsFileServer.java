/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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
et
*/
package com.github.lindenb.jvarkit.tools.server;


import java.io.BufferedReader;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.nio.charset.Charset;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import org.eclipse.jetty.server.Server;
import org.eclipse.jetty.servlet.ServletContextHandler;
import org.eclipse.jetty.servlet.ServletHolder;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.IntervalParserFactory;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFCodec;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFLine;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.cram.ref.CRAMReferenceSource;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;
/**
BEGIN_DOC

## input

input is a set of indexed Vcf/Bam file or a file with the suffix `.list` containing the path to the files.
 
## Example

```
java -jar dist/htsfileserver.jar -R src/test/resources/rotavirus_rf.fa src/test/resources/S*.bam  src/test/resources/rotavirus_rf.*.vcf.gz
```

# Screenshot

![https://i.imgur.com/ObRsVxE.png](https://i.imgur.com/ObRsVxE.png)

END_DOC
 
 */

@Program(name="htsfileserver",
	description="Jetty Based http server serving Vcf,Bam,Tabix files.",
	creationDate="20200405",
	modificationDate="20200406",
	keywords={"vcf","bam","server","tabix"},
	biostars={430718},
	generate_doc=true
	)
public  class HtsFileServer extends Launcher {
	private static final String KEY_INTERVAL="interval";
	private static final String KEY_ACTION="action";
	private static final String KEY_FILEID="fid";
	private static final String VALUE_DUMP="dump";
	
	private static final Logger LOG = Logger.build(HtsFileServer.class).make();
	
	@Parameter(names="--port",description="server port.")
	private int serverPort = 8080;
	@Parameter(names= {"--gtf"},description="Optional GTF file. Will be used to retrieve an interval by gene name")
	private Path gtfFile = null;
	@Parameter(names= {"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path faidxRef = null;
	@Parameter(names= {"-G","--no-genotype"},description="remove genotypes from vcf")
	private boolean remove_genotype_vcf = false;

	
	/** base handler for any data that can be querid by interval */
	private abstract class AbstractInput {
		private final Path path;
		private final String md5;
		AbstractInput(final Path path) throws IOException{
			this.path = path;
			IOUtil.assertFileIsReadable(path);
			this.md5 = StringUtils.md5(path.getFileName().toString()).substring(0,10);
			}
		public String getMd5()
			{
			return md5;
			}
		public Path getPath()
			{
			return path;
			}
		abstract String getOutputName();
		abstract String getContentType();
		abstract void dump(final Locatable loc,final OutputStream os,final CRAMReferenceSource ref) throws IOException;
		}
	
	/** vcf implementation */
	private class VcfInput extends AbstractInput {
		private boolean has_genotypes;
		VcfInput(final Path path) throws IOException{
			super(path);
			// try open with index
			try(VCFReader ignore= VCFReaderFactory.makeDefault().open(path, true)) {
				has_genotypes = ignore.getHeader().hasGenotypingData();
				//nothing
				}
			}
		
		@Override
		String getOutputName()
			{
			String s= getPath().getFileName().toString();
			if(!s.endsWith(".gz")) s+=".gz";
			return s;
			}
		
		@Override
		String getContentType()
			{
			return "application/gzip";
			}
		
		@Override
		void dump(final Locatable loc, final OutputStream os,final CRAMReferenceSource ignore) throws IOException
			{
			if(loc==null && (!this.has_genotypes || !remove_genotype_vcf) && getPath().getFileName().toString().endsWith(FileExtensions.COMPRESSED_VCF)) {
				IOUtils.copyTo(getPath(), os);
				return;
			}
			final BlockCompressedOutputStream bos = new BlockCompressedOutputStream(os,(Path)null);
			final VariantContextWriterBuilder vcb = new VariantContextWriterBuilder();
			vcb.setOutputStream(bos);
			vcb.setCreateMD5(false);
			vcb.setReferenceDictionary(dictionary);
			vcb.clearOptions();
			try(VariantContextWriter w=vcb.build()) {
				try(VCFReader r = VCFReaderFactory.makeDefault().open(getPath(), true)) {
					final VCFHeader header = r.getHeader();
					if(remove_genotype_vcf && this.has_genotypes) {
						w.writeHeader(new VCFHeader(header.getMetaDataInInputOrder(),Collections.emptyList()));
						}
					else
						{
						w.writeHeader(header);
						}
					
					CloseableIterator<VariantContext> iter;
					if(loc==null) {
						iter = r.iterator();
						}
					else
						{
						iter = r.query(loc);
						}
					while(iter.hasNext()) {
						VariantContext ctx=iter.next();
						if(remove_genotype_vcf) {
							ctx =  new VariantContextBuilder(ctx).noGenotypes().make();
							}
						w.add(ctx);
						}
					iter.close();
					}
				}
			bos.flush();
			bos.close();
			}
		}
	/** implementation for Bam */
	private class BamInput extends AbstractInput {
		BamInput(final Path path) throws IOException{
			super(path);
			}
		@Override
		String getContentType()
			{
			return "application/bam";
			}
		@Override
		String getOutputName()
			{
			String s= getPath().getFileName().toString();
			if(!s.endsWith(FileExtensions.BAM)) s+=FileExtensions.BAM;
			return s;
			}
		@Override
		void dump(final Locatable loc, OutputStream os,final CRAMReferenceSource ref) throws IOException
			{
			if(loc==null && getPath().getFileName().toString().endsWith(FileExtensions.BAM)) {
				IOUtils.copyTo(getPath(), os);
				return;
				}
			
			final SAMFileWriterFactory vcb = new SAMFileWriterFactory()
					.setCreateIndex(false)
					.setCreateMd5File(false)
					;
			
			try(SamReader r = SamReaderFactory.makeDefault().
						validationStringency(ValidationStringency.SILENT).
						referenceSource(ref).
						open(getPath())) {
				final SAMFileHeader header = r.getFileHeader();
				try(SAMFileWriter w=vcb.makeBAMWriter(header, true, os)) {
					final CloseableIterator<SAMRecord> iter;
					if(loc==null) {
						iter = r.iterator();
						}
					else
						{
						iter = r.query(loc.getContig(), loc.getStart(), loc.getEnd(), false);
						}
					while(iter.hasNext()) {
						w.addAlignment(iter.next());
						}
					iter.close();
					}
				
				}
			}
		}
	
	/** implementation for tabix files */
	private class TabixInput extends AbstractInput {
		public TabixInput(final Path path) throws IOException{
			super(path);
			try(TabixReader tr=new TabixReader(path.toString())) {
				//do nothing
				}
			}
		@Override
		String getOutputName() {
			return getPath().getFileName().toString();
			}

		@Override
		String getContentType() {
			return "application/gzip";
		}

		@Override
		void dump(Locatable loc, OutputStream os, CRAMReferenceSource ref) throws IOException {
			if(loc==null) {
				IOUtils.copyTo(getPath(), os);
				return;
			}
			try(TabixReader tb = new TabixReader(getPath().toString())) {
				try(final BlockCompressedOutputStream bos = new BlockCompressedOutputStream(os,(Path)null)) {
					final PrintWriter pw = new PrintWriter(bos);
					final ContigNameConverter converter = ContigNameConverter.fromContigSet(tb.getChromosomes());
					final Optional<SimpleInterval> r= converter.convertToSimpleInterval(loc);
					if(r.isPresent()) {
						final Locatable loc2 = r.get();
						final TabixReader.Iterator iter = tb.query(loc2.getContig(),loc2.getStart(),loc2.getEnd());
						for(;;) {
							final String line = iter.next();
							if(line==null) break;
							pw.println(line);
							}
						}
					pw.flush();
					pw.close();
					}
				}
			}
		}
	
	/** implementation for tabix files */
	private class FastaInput extends AbstractInput {
		public FastaInput(final Path path) throws IOException{
			super(path);
			try(ReferenceSequenceFile ref= ReferenceSequenceFileFactory.getReferenceSequenceFile(path)) {
				SequenceDictionaryUtils.extractRequired(ref);
				}
			}
		@Override
		String getOutputName() {
			String s= getPath().getFileName().toString();
			if(!s.endsWith(".gz")) s+=".gz";
			return s;
			}

		@Override
		String getContentType() {
			return "application/gzip";
		}

		@Override
		void dump(Locatable loc, OutputStream os, CRAMReferenceSource ignore) throws IOException {
			if(loc==null) {
				if(getPath().getFileName().toString().endsWith(".gz")) {
					IOUtils.copyTo(getPath(), os);
					}
				else
					{
					try(final BlockCompressedOutputStream bos = new BlockCompressedOutputStream(os,(Path)null)) {
						IOUtils.copyTo(getPath(), bos);
						}
					}
				return;
				}
			try(ReferenceSequenceFile ref= ReferenceSequenceFileFactory.getReferenceSequenceFile(getPath())) {
				try(final BlockCompressedOutputStream bos = new BlockCompressedOutputStream(os,(Path)null)) {
					final PrintWriter pw = new PrintWriter(bos);
					final ContigNameConverter converter = ContigNameConverter.fromOneDictionary(ref.getSequenceDictionary());
					final Optional<SimpleInterval> r= converter.convertToSimpleInterval(loc);
					if(r.isPresent()) {
						final SimpleInterval loc2 = r.get();
						final ReferenceSequence refseq = ref.getSubsequenceAt(loc2.getContig(), loc2.getStart(), loc2.getEnd());
						pw.print(">"+loc2.toNiceString());
						for(int i=0;i< refseq.length();i++) {
							if(i%60==0) pw.println();
							pw.print((char)refseq.getBases()[i]);
							}
						pw.println();
						}
					pw.flush();
					pw.close();
					}
				}
			}
		}
	
	private final Map<String,AbstractInput> htsMap = new HashMap<>();
	private SAMSequenceDictionary dictionary = null;
	
	@SuppressWarnings("serial")
	private  class HtsFileServerServlet extends HttpServlet {
		@Override
		protected void doGet(final HttpServletRequest req,final  HttpServletResponse resp) throws ServletException, IOException {
			doPost(req, resp);
			}
		@Override
		protected void doPost(final HttpServletRequest request, final HttpServletResponse response) throws ServletException, IOException {
			final String action = request.getParameter(KEY_ACTION);
			if(VALUE_DUMP.equals(action))
				{
				dumpData(request,response);
				}
			else
				{
				printPage(request,response);
				}
			}
		}
		
	
	
	/** dumpData */
	private void dumpData(final HttpServletRequest request,final HttpServletResponse response)	throws IOException, ServletException
		{
		final String fileids[]  = request.getParameterValues(KEY_FILEID);
		final List<AbstractInput> selected = new ArrayList<>();
		if(fileids!=null) {
			for(final String fid: fileids) {
				final AbstractInput input = this.htsMap.get(fid);
				if(input!=null) selected.add(input);
			}
		}
		// default: add all
		if(selected.isEmpty()) {
			selected.addAll(this.htsMap.values());
			}
		Locatable loc = null;
		/* parse user interval */
		final String intervalstr = request.getParameter(KEY_INTERVAL);
		if(!StringUtils.isBlank(intervalstr)) {
			loc = IntervalParserFactory.newInstance().
					dictionary(this.dictionary).
					make().
					apply(intervalstr).
					orElse(null);
			if(loc==null && this.gtfFile!=null) {
				final ContigNameConverter cvt = ContigNameConverter.fromOneDictionary(this.dictionary);
				final String geneName=intervalstr.trim();
				final GTFCodec codec = new GTFCodec();
				try(BufferedReader br = IOUtils.openPathForBufferedReading(gtfFile)) {
					br.lines().map(line->{
						if(StringUtils.isBlank(line) ||  line.startsWith("#"))return null;
						final String tokens[]= CharSplitter.TAB.split(line);
						if(tokens.length<9 ) return null;
						if(!(tokens[2].equals("gene") || tokens[2].equals("transcript"))) return null;
						if(StringUtils.indexOfIgnoreCase(tokens[8],geneName)==-1) return null;
						final GTFLine gtfLine = codec.decode(line);
						if(gtfLine==null) return null;
						
						if(tokens[2].equals("gene") ) {
							if(!(geneName.equals(gtfLine.getAttribute("gene_id")) || geneName.equals(gtfLine.getAttribute("gene_name")))) return null;
							}
						else if(tokens[2].equals("transcript") ) {
							if(!(geneName.equals(gtfLine.getAttribute("transcript_id")))) return null;
							}
					
						final String ctg = cvt.apply(gtfLine.getContig());
						if(StringUtils.isBlank(ctg)) return null;
							return new SimpleInterval(ctg,gtfLine.getStart(),gtfLine.getEnd());
							}).
						filter(R->R!=null).
						findFirst().
						orElse(null);
					}
				}	
			if(loc==null) loc = new SimpleInterval("undefined_interval",1,1);
			}
		
		
		String prefix= StringUtils.now()+".";
		if(loc!=null) {
			prefix += loc.getContig()+"_"+loc.getStart()+"_"+loc.getEnd()+".";
			}
		
		final CRAMReferenceSource referenceSource = (
				selected.stream().anyMatch(I->I instanceof BamInput)?
				new ReferenceSource(faidxRef)
				: null
				);
		
		try(PrintStream out = new PrintStream(response.getOutputStream())) {
			final String charset = StringUtils.ifBlank(request.getCharacterEncoding(), "UTF-8");
			response.setCharacterEncoding(charset);

			if(selected.size()==1) {
				final AbstractInput first=selected.get(0);
				final String fname = prefix+first.getOutputName();
				response.addHeader("Content-Type",first.getContentType());
				response.addHeader("Content-Disposition","attachment; name=\""+ fname +"\"; filename=\""+ fname +"\"");
				response.setContentType("data/binary; charset="+charset.toLowerCase());

				first.dump(loc, out,referenceSource);
				}
			else
				{
				final String fname = prefix+HtsFileServer.class.getSimpleName()+".zip";
				response.addHeader("Content-Type","application/zip");
				response.addHeader("Content-Disposition","attachment; name=\""+ fname +"\"; filename=\""+ fname +"\"");

				response.setContentType("data/binary; charset="+charset.toLowerCase());
				final ZipOutputStream zout = new ZipOutputStream(out, Charset.forName(charset));
				zout.setLevel(0);
				for(AbstractInput input: selected) {
					final ZipEntry zipEntry= new ZipEntry(prefix+input.getOutputName());
					zout.putNextEntry(zipEntry);
					OutputStream uos = IOUtils.uncloseableOutputStream(zout);
					input.dump(loc, uos,referenceSource);
					zout.closeEntry();
					if(out.checkError()) break;
					}
				
				zout.finish();
				}
			out.flush();
			}
		
		}

	/** print HTML page */
	private void printPage(final HttpServletRequest request,final HttpServletResponse response)	throws IOException, ServletException
		{
		 
		 final String charset = StringUtils.ifBlank(request.getCharacterEncoding(), "UTF-8");
		 response.setContentType("text/html; charset="+charset.toLowerCase());
		 response.setCharacterEncoding(charset);
		 PrintWriter pw = response.getWriter();
		
		
		 
		 final String title = ("HtsFileServer "+SequenceDictionaryUtils.getBuildName(this.dictionary).orElse("")).trim();
		
		 
		 try {
			final XMLStreamWriter w=XMLOutputFactory.newFactory().createXMLStreamWriter(pw);
			w.writeStartElement("html");
			
			w.writeStartElement("head");
			
			w.writeEmptyElement("meta");
			w.writeAttribute("charset", charset);
			
			w.writeStartElement("title");
			w.writeCharacters(title);
			w.writeEndElement();//title
			
			w.writeStartElement("style");
			w.writeCharacters(
					"body {background-color:#f0f0f0;color:#070707;font-size:18px;}"
					+ "h1 {text-align:center;color:#070707;}"
					+ ".lbl {font-weight: bold;}"
					+ ".files {font-family:monospace;font-size:12px; font-style: normal;color:gray;}"
					+ ".headerform {background-color:lavender;text-align:center;font-size:14px;}"
					);
			w.writeEndElement();//title
			
			w.writeStartElement("script");
			w.writeCharacters("var htsFiles={");
			boolean first=true;
			for(final AbstractInput input: this.htsMap.values()) {
				if(!first) w.writeCharacters(",");
				w.writeCharacters("\""+input.getMd5()+"\":\"");
				w.writeCharacters(input.getPath().toString());
				w.writeCharacters("\"");
				first = false;
				}
			w.writeCharacters("};");
			w.writeCharacters("function selectBams() {"+
					"for(k in htsFiles) {var f=htsFiles[k]; if(f.endsWith(\""+ FileExtensions.BAM+"\") || f.endsWith(\""+ FileExtensions.CRAM+"\")) document.getElementById(k).checked=true;}"+
					"}\n");
			w.writeCharacters("function selectVcfs() {"+
					"for(k in htsFiles) {var f=htsFiles[k]; if(f.endsWith(\""+ FileExtensions.VCF+"\") || f.endsWith(\""+ FileExtensions.COMPRESSED_VCF+"\")) document.getElementById(k).checked=true;}"+
					"}\n");
			w.writeCharacters("function selectInvert() {"+
					"for(k in htsFiles) { var e=document.getElementById(k); e.checked=!e.checked;}"+
					"}\n");
			w.writeCharacters("function selectAll(b) {"+
				"for(k in htsFiles) { document.getElementById(k).checked=b;}"+
				"}\n");	
					
			w.writeEndElement();//script
			
			w.writeEndElement();//head
			
			w.writeStartElement("body");
			
			
			
			w.writeStartElement("h1");
			w.writeCharacters(title);
			w.writeEndElement();//h1
			
			w.writeComment("BEGIN FORM");
			w.writeStartElement("form");
			w.writeAttribute("method", "POST");
			w.writeAttribute("action", "/page");
			
			w.writeStartElement("div");
			w.writeAttribute("class", "headerform");

			
			w.writeEmptyElement("input");
			w.writeAttribute("name", KEY_ACTION);
			w.writeAttribute("type", "hidden");
			w.writeAttribute("value",VALUE_DUMP);

			w.writeStartElement("label");
			w.writeAttribute("class", "lbl");
			w.writeAttribute("for", KEY_INTERVAL);
			w.writeCharacters("Interval"+(this.gtfFile==null?"":" or Gene")+":");
			w.writeEndElement();//label
			
			w.writeEmptyElement("input");
			w.writeAttribute("name", KEY_INTERVAL);
			w.writeAttribute("id", KEY_INTERVAL);
			w.writeAttribute("type", "text");
			w.writeAttribute("value","");
			w.writeAttribute("placeholder","chrom:start-end");
			
			w.writeEmptyElement("input");
			w.writeAttribute("type", "Submit");
			w.writeAttribute("class", "btn");
			w.writeAttribute("value", "Go");

				
			w.writeEmptyElement("br");

			
			
			w.writeStartElement("button");
			w.writeAttribute("class", "btn");
			w.writeAttribute("onclick", "selectBams(); return false;");
			w.writeCharacters("Select Bams");
			w.writeEndElement();//button
			w.writeCharacters(" ");
			w.writeStartElement("button");
			w.writeAttribute("class", "btn");
			w.writeAttribute("onclick", "selectVcfs();return false;");
			w.writeCharacters("Select Vcfs");
			w.writeEndElement();//button
			w.writeCharacters(" ");
			w.writeStartElement("button");
			w.writeAttribute("class", "btn");
			w.writeAttribute("onclick", "selectInvert();return false;");
			w.writeCharacters("Invert Selection");
			w.writeEndElement();//button
			w.writeCharacters(" ");
			w.writeStartElement("button");
			w.writeAttribute("class", "btn");
			w.writeAttribute("onclick", "selectAll(true);return false;");
			w.writeCharacters("Select All");
			w.writeEndElement();//button
			w.writeCharacters(" ");
			w.writeStartElement("button");
			w.writeAttribute("class", "btn");
			w.writeAttribute("onclick", "selectAll(false);return false;");
			w.writeCharacters("Select None");
			w.writeEndElement();//button

			w.writeEndElement();//div
			
			
			w.writeStartElement("div");
			w.writeAttribute("class", "files");
			w.writeStartElement("ul");
			
			for( AbstractInput htsFile: this.htsMap.values()) {
				w.writeStartElement("li");
			
				w.writeEmptyElement("input");
				w.writeAttribute("type", "checkbox");
				w.writeAttribute("name",KEY_FILEID);
				w.writeAttribute("id",htsFile.getMd5());
				w.writeAttribute("value",htsFile.getMd5());
				
				w.writeStartElement("label");
				w.writeAttribute("for",htsFile.getMd5());
				w.writeCharacters(htsFile.getPath().toString());
				w.writeEndElement();//label
				
				w.writeEndElement(); //lli
			}
			
			w.writeEndElement();//ul
			w.writeEndElement();//div
			
			w.writeEndElement();//form
			
			w.writeComment("END FORM");
			
			w.writeEmptyElement("hr");
			w.writeStartElement("div");
			w.writeCharacters("Author: Pierre Lindenbaum. ");
			w.writeCharacters(JVarkitVersion.getInstance().getLabel());
			w.writeEndElement();
			
			
			w.writeEndElement();//body
			w.writeEndElement();//html
			w.flush();
			w.close();
		 	}
		 catch(XMLStreamException err) { LOG.warn(err);throw new IOException(err);}
		 finally { pw.close(); }
		}
		
	
	
	
	@Override
	public int doWork(final List<String> args) {
	
	
		try {
			this.dictionary  = SequenceDictionaryUtils.extractRequired(this.faidxRef);
			
			final Set<String> filenames = new HashSet<>();
			for(final Path path: IOUtils.unrollPaths(args)) {
				final String fn = path.getFileName().toString();
				if(!filenames.add(fn)) {
					LOG.error("duplicate filename "+fn);
					return -1;
					}
				final AbstractInput input;
				final SAMSequenceDictionary dict;
				if(fn.endsWith(FileExtensions.BAM) || fn.endsWith(FileExtensions.CRAM)) {
					dict = SequenceDictionaryUtils.extractRequired(path);
					input = new BamInput(path);
					}
				else if(fn.endsWith(FileExtensions.VCF) || fn.endsWith(FileExtensions.COMPRESSED_VCF)) {
					dict = SequenceDictionaryUtils.extractRequired(path);
					input = new VcfInput(path);
					}
				else if(FileExtensions.FASTA.stream().anyMatch(X->fn.endsWith(X))) {
					dict = SequenceDictionaryUtils.extractRequired(path);
					input = new FastaInput(path);
					}
				else if(fn.endsWith(".gz")) {
					input = new TabixInput(path);
					dict = null;
					}
				else
					{
					LOG.error("unsupported format "+path);
					return -1;
					}
				
				if(dict!=null) SequenceUtil.assertSequenceDictionariesEqual(this.dictionary, dict);
				
				if(this.htsMap.containsKey(input.getMd5())) {
					LOG.error("duplicate key "+input.getMd5()+" "+input.getPath());
					return -1;
					}
				this.htsMap.put(input.getMd5(), input);
				}
			if(this.htsMap.isEmpty()) {
				LOG.error("No VCF/BAM/Tabix defined.");
				return -1;
				}
			
			
			final Server server = new Server(this.serverPort);
			
			final ServletContextHandler context = new ServletContextHandler();
	        context.addServlet(new ServletHolder(new HtsFileServerServlet()),"/*");
	        context.setContextPath("/");
	        context.setResourceBase(".");
	        server.setHandler(context);
			
		    
		    
		    LOG.info("Starting server "+getProgramName()+" on port "+this.serverPort);
		    server.start();
		    LOG.info("Server started. Press Ctrl-C to stop. Check your proxy settings ."
		    		+ " Open a web browser at http://localhost:"+this.serverPort+"/htsserver .");
		    server.join();
		    return 0;
			}
		catch (final Throwable err) {
			LOG.error(err);
			return -1;
			}
		
		}	


public static void main(final String[] args) throws Exception{
    new HtsFileServer().instanceMainWithExit(args);
	}

}
