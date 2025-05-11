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
package com.github.lindenb.jvarkit.tools.labguru;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.URI;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import org.apache.commons.compress.archivers.ArchiveEntry;
import org.apache.commons.compress.archivers.tar.TarArchiveInputStream;
import org.apache.http.HttpResponse;
import org.apache.http.client.HttpClient;
import org.apache.http.client.config.CookieSpecs;
import org.apache.http.client.config.RequestConfig;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.client.utils.URIBuilder;
import org.apache.http.impl.client.HttpClients;
import org.apache.http.util.EntityUtils;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.gatk.GATKConstants;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.variant.vcf.VcfHeaderExtractor;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

/**
BEGIN_DOC

Motivation: scan the files stored in labguru


usage:

```
jvarkit -jar jvarkit.jar scanlabguru -A api.txt <projectid1> <projectid2> <projectid3>
```



END_DOC
 */
@Program(
		name="scanlabguru",
		description="scan the files stored in labguru",
		modificationDate="20240325",
		creationDate="20240325",
		keywords={"labguru","lims","vcf","sam"}
		)
public class ScanLabGuru extends Launcher {
private static final Logger LOG = Logger.of(ScanLabGuru.class);


@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
private Path ouput = null;

@Parameter(names={"--base"},description="Labguru base URI.")
private String baseUri = "https://cle.inserm.fr";

@Parameter(names={"--api-key","-A"},description="API key (can be a string or an existing file",required = true)
private String apiKeySource="";


abstract class JsonWrapper {
	protected final JsonObject root;
	
	protected JsonWrapper(final JsonElement root) {
		this.root = root.getAsJsonObject();
		}

	public String getId() {
		return root.get("id").getAsString();
		}
	public String getTitle() {
		return root.get("title").getAsString();
		}
	
	public abstract String getTypeName();
	JsonObject fetch() {
		return  parseJson("/"+getTypeName()+"/"+getId(), null);
		}
	
	protected <T extends JsonWrapper> List<T> getItems(final String itemNames,final Function<JsonElement, T> mapper) {
		final JsonArray L=fetch().get(itemNames).getAsJsonArray();
		final List<T> result = new ArrayList<>(L.size());
		for(int i=0;i< L.size();i++) {
			result.add(mapper.apply(L.get(i)));
			}
		return result;
		}
	
	}
class Project extends JsonWrapper {
	Project(JsonElement root) {
		super(root);
		}
	
	Project(String id) {
		super(new JsonObject());
		super.root.addProperty("id", id);
		super.root.addProperty("title", id);
		}
	
	public String getTypeName() {
		return "projects";
		}
	public List<Milestone> getMilestones() {
		return getItems("milestones",J->new Milestone(J));
		}

	}
class Milestone extends JsonWrapper {
	Milestone(JsonElement root) {
		super(root);
		}
	public String getTypeName() {
		return "milestones";
		}
	public List<Experiment> getExperiments() {
		return getItems("experiments",J->new Experiment(J));
		}
	
	}
class Experiment extends JsonWrapper {
	Experiment(JsonElement root) {
		super(root);
		}
	public String getTypeName() {
		return "experiments";
		}
	public List<Attachment> getAttachments() {
		return getItems("attachments",J->new Attachment(J));
		}
	}


class Attachment extends JsonWrapper {
	Attachment(JsonElement root) {
		super(root);
		}
	public String getTypeName() {
		return "attachments";
		}
	public String getFilename() {
		return root.get("filename").getAsString();
		}
	
	private byte[] readBuffer(InputStream in) throws IOException {
		final int BUFFER_SIZE = 2_000_000;
		byte[] buffer = new byte[BUFFER_SIZE];
		int nRead=0;
		while(nRead<buffer.length) {
			int n= in.read(buffer,nRead,buffer.length-nRead);
			if(n<0) break;
			nRead+=n;
			}
		return buffer;
	}

	private void  vcf(XMLStreamWriter w,InputStream in,String filename) throws IOException,XMLStreamException {
		boolean debug=false;
		try(ByteArrayInputStream bis = new ByteArrayInputStream(readBuffer(in))) {
			final VCFHeader header=VcfHeaderExtractor.decode(bis);
			if(header!=null) {
				w.writeStartElement("vcf");
				w.writeAttribute("path", filename);
				final Set<String> flags = new HashSet<>();
				for(VCFInfoHeaderLine h:header.getInfoHeaderLines()) {
					if(h.getID().equals(VCFConstants.SVTYPE)) flags.add("SVTYPE");
					else if(h.getID().equals("SVLEN")) flags.add("SVLEN");
					}
				if(header.getInfoHeaderLine(GATKConstants.MQRankSum_KEY)!=null && 
					header.getOtherHeaderLines().stream().anyMatch(H->H.getValue().contains("GATKCommandLine"))) {
					flags.add("GATK");
					}
				if(header.getOtherHeaderLines().stream().anyMatch(V->V.getValue().contains("configManta.py"))) {
					flags.add("MANTA");
					}
				
				if(header.getFilterHeaderLine("FailDellyFilter")!=null || header.getFormatHeaderLine("RDCN")!=null) {
						flags.add("DELLY");
						}
				
				for(final String flag: flags) {
					w.writeStartElement("attribute");
					w.writeCharacters(flag);
					w.writeEndElement();
					}
				
				if(header.hasGenotypingData()) {
					w.writeStartElement("samples");
					for(final String sn : header.getGenotypeSamples()) {
						w.writeStartElement("sample");
						w.writeCharacters(sn);
						w.writeEndElement();
						}
					w.writeEndElement();
					}
				w.writeEndElement();
				
				if(debug) {
					VariantContextWriter wc=VCFUtils.createVariantContextWriterToOutputStream(System.err);
					wc.writeHeader(header);
					}
				}
			}
		catch(Throwable err3) {
			LOG.error("Cannot decode "+filename);
			err3.printStackTrace();
			}
		}

	
	private void sam(XMLStreamWriter w,InputStream in,String filename) throws IOException,XMLStreamException {
		try(ByteArrayInputStream bis = new ByteArrayInputStream(readBuffer(in))) {
			try {
				final SamReaderFactory srf = SamReaderFactory.makeDefault();
				srf.validationStringency(ValidationStringency.SILENT);
				try(SamReader sr=srf.open(SamInputResource.of(bis))) {
					SAMFileHeader header= sr.getFileHeader();
					w.writeStartElement("bam");
					w.writeAttribute("path", filename);
					w.writeStartElement("samples");
					for(final String sn : header.getReadGroups().stream().map(RG->RG.getSample()).filter(S->!StringUtils.isBlank(S)).collect(Collectors.toSet())) {
						w.writeStartElement("sample");
						w.writeCharacters(sn);
						w.writeEndElement();
						}
					w.writeEndElement();
						
					w.writeEndElement();
					}
				}
			catch(Throwable err3) {
				LOG.error("Cannot decode "+filename);
				err3.printStackTrace();
				}
			}
		}
	
	private void  byFileNameName(XMLStreamWriter w, final String filename,InputStream in) throws IOException,XMLStreamException {
		if(FileExtensions.VCF_LIST.stream().anyMatch(EX->filename.endsWith(EX))) {
			vcf(w,in,filename);
			}
		else if(filename.endsWith(FileExtensions.BAM) || filename.endsWith(FileExtensions.SAM)) {
			sam(w,in,filename);
			}
		}

	
	
	public void  scan(XMLStreamWriter w) throws XMLStreamException {
		try {
			final String filename = getFilename();
			boolean is_vcf = FileExtensions.VCF_LIST.stream().anyMatch(EX->filename.endsWith(EX));
			boolean is_zip = filename.endsWith(".zip");
			boolean is_tar = filename.endsWith(".tar") || filename.endsWith(".tar.gz");;
			boolean parseable= is_vcf || is_zip || is_tar;
			
			
			
			w.writeStartElement("attachment");
			w.writeAttribute("id", this.getId());
			w.writeAttribute("filename",filename);
			

			if(parseable) {
				final URIBuilder builder = new URIBuilder(getBaseApi()+"/"+getTypeName()+"/download");
				builder.setParameter("annotated","false");
				builder.setParameter("id",getId());
				builder.setParameter("token",api_key);
				final URI uri= builder.build();
				LOG.info(uri);
				final HttpGet httpGet = new HttpGet(uri);
				final HttpResponse  response = httpClient.execute(httpGet);
				try(InputStream in = response.getEntity().getContent()) {
					if(filename.endsWith(".tar") || filename.endsWith(".tar.gz")) {
						w.writeStartElement("tar");
						final TarArchiveInputStream tarin = new TarArchiveInputStream(IOUtils.mayBeGzippedInputStream(in));
						ArchiveEntry entry = null;
			            while ( (entry = tarin.getNextEntry()) != null ) 
			            	{
			            	final String entryName = entry.getName();
			            	w.writeStartElement("entry");
			            	w.writeAttribute("filename",entryName);
			            	byFileNameName(w,entryName,tarin);
			            	w.writeEndElement();
			            	}
						w.writeEndElement();
						}
					else if(filename.endsWith(".zip")) {
						w.writeStartElement("zip");
						try(ZipInputStream zis = new ZipInputStream(in))  {
							ZipEntry entry = null;
				            while ( (entry = zis.getNextEntry()) != null ) 
				            	{
				            	final String entryName = entry.getName();
				            	w.writeStartElement("entry");
				            	w.writeAttribute("filename",entryName);
				            	byFileNameName(w,entryName,zis);
				            	zis.closeEntry();
				            	w.writeEndElement();
				            	}
							}
						w.writeEndElement();
						}
					else
						{
						byFileNameName(w, filename, in);
						}
					}
				catch(Throwable err2) {
					LOG.error(err2);
					}
				}
			else
				{
				w.writeComment("not parseable");
				}
			w.writeEndElement();
			LOG.info("done "+getFilename());
			}
		catch(Throwable err) {
			LOG.error(err);
			}
		finally 
			{
			
			}
		}
	}


private String api_key = "";
private HttpClient httpClient;
private String getBaseApi() {
	return baseUri+"/api/v1";
	}

private JsonObject parseJson(String query,Map<String, String> params) {
	try {
		final JsonParser parser = new JsonParser();
		final URIBuilder builder = new URIBuilder(getBaseApi()+query);
		builder.setParameter("token",api_key);
		if(params!=null) {
			for(String k:params.keySet()) {
				builder.setParameter(k,String.valueOf(params.get(k)));
				}
			}
		final URI uri= builder.build();
		LOG.info(uri);
		final HttpGet httpGet = new HttpGet(uri);
		final HttpResponse  response = httpClient.execute(httpGet);
		return parser.parse(EntityUtils.toString(response.getEntity())).getAsJsonObject();
    	}
	catch (final Throwable err) {
        throw new RuntimeIOException(err);
    	}
	}
@Override
public int doWork(final List<String> args) {
	try {
		try {
			final Path p=Paths.get(this.apiKeySource);
			if(Files.exists(p) && Files.isReadable(p) && Files.isRegularFile(p)) {
				this.api_key = IOUtils.slurpPath(p).trim();
				}
			}
		catch(Throwable err) {
			}
		if(StringUtils.isBlank(this.api_key)) {
			this.api_key = this.apiKeySource;
			}
		if(StringUtils.isBlank(this.api_key)) {
			LOG.error("API_KEY is blank");
			return -1;
			}
		
		this.httpClient =  HttpClients.custom()
			    .setDefaultRequestConfig(RequestConfig.custom()
			            .setCookieSpec(CookieSpecs.STANDARD)
			            .build())
			        .build();

		
		try(OutputStream os=super.openPathOrStdoutAsStream(this.ouput)) {
			final XMLOutputFactory xof = XMLOutputFactory.newFactory();
			final XMLStreamWriter w=xof.createXMLStreamWriter(os);
			w.writeStartDocument("UTF-8", "1.0");
			w.writeStartElement("labguru");
			for(String projectId : args.stream().
						flatMap(S->Arrays.stream(S.split("[,; ]"))).
						filter(S->!StringUtils.isBlank(S)).
						collect(Collectors.toSet())) {
				w.writeStartElement("project");
				w.writeAttribute("id", projectId);
				Project p = new Project(projectId);
				for(Milestone m :p.getMilestones() ) {
					w.writeStartElement("milestone");
					w.writeAttribute("id", m.getId());
					w.writeAttribute("title", m.getTitle());

					for(Experiment ex :m.getExperiments() ) {
						w.writeStartElement("experiment");
						w.writeAttribute("id", ex.getId());
						w.writeAttribute("title", ex.getTitle());

						
						for(Attachment attach :ex.getAttachments() ) {
							attach.scan(w);
							}
						w.writeEndElement();
						}
					w.writeEndElement();
					}
				w.writeEndElement();
				}
			w.writeEndElement();
			os.flush();
			}
		return 0;
		}
	catch(Throwable err) {
		LOG.error(err);
		return -1;
		}
	finally 
		{
		}
	}

public static void main(String[] args) {
	new ScanLabGuru().instanceMainWithExit(args);
}

}
