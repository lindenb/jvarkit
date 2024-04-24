/*
The MIT License (MIT)

Copyright (c) 2024 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.ngsfiles;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.illumina.FastQName;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.variant.vcf.VcfHeaderExtractor;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.vcf.VCFHeader;

/**
BEGIN_DOC

## Example

```bash
find /projects/align01/ -type f |\
  java -jar dist/ngsfilessummary.jar 

SAMPLE1	BAM	/projects/align01/Samples/SAMPLE1/BAM/SAMPLE1_final.bam	321262321	Wed Jun 26 10:30:07 CEST 2013
SAMPLE1	FASTQ	/project/align01/fastq/SAMPLE1/SAMPLE1_CGATGT_L008_R1_002.fastq.gz	35828879	Fri Oct 18 16:15:58 CEST 2013
SAMPLE1	VCF	/projects/align01/Samples/SAMPLE1/VCF/SAMPLE1_variations.freebayes.vcf.gz	184191	Mon Jun 17 14:47:22 CEST 2013
SAMPLE1	VCF	/projects/align01/Samples/SAMPLE1/VCF/SAMPLE1_variations.gatk.vcf.gz	113341	Mon Jun 17 11:57:19 CEST 2013
SAMPLE1	VCF	/projects/align01/Samples/SAMPLE1/VCF/SAMPLE1_variations.samtools.vcf.gz	57518	Mon Jun 17 11:58:49 CEST 2013
SAMPLE2	BAM	/projects/align01/Samples/SAMPLE2/BAM/SAMPLE2_final.bam	286100773	Wed Jun 26 10:47:09 CEST 2013
SAMPLE2	FASTQ	/project/align01/fastq/SAMPLE2/SAMPLE2_CGATGT_L008_R1_002.fastq.gz	356828879	Fri Oct 18 16:15:58 CEST 2013
SAMPLE2	VCF	/projects/align01/Samples/SAMPLE2/VCF/SAMPLE2_variations.freebayes.vcf.gz	172970	Mon Jun 17 14:45:51 CEST 2013
SAMPLE2	VCF	/projects/align01/Samples/SAMPLE2/VCF/SAMPLE2_variations.gatk.vcf.gz	106390	Mon Jun 17 11:57:19 CEST 2013
SAMPLE2	VCF	/projects/align01/Samples/SAMPLE2/VCF/SAMPLE2_variations.samtools.vcf.gz	52709	Mon Jun 17 11:58:04 CEST 2013
```
END_DOC

 */
@Program(name="ngsfilessummary",
	description="Scan folders and generate a summary of the files (SAMPLE/BAM SAMPLE/VCF etc..). Useful to get a summary of your samples.",
	keywords= {"sam","bam","vcf","util"},
	modificationDate="20240324",
	creationDate = "20140430",
	jvarkit_amalgamion = true
	)
public class NgsFilesSummary extends Launcher
	{
	private enum Format {tsv,xml};
	private static final Logger LOG = Logger.build(NgsFilesSummary.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-header","--header"},description="[20180725]print header")
	private boolean show_header= false;
	@Parameter(names={"-R","--reference"},description="[20190905]restrict to thoses reference(s). Also is used to read CRAM files. A file with the '.list' suffix is interpreted as a list of paths to fasta REF.")
	private List<String> faidxPathList = new ArrayList<>();
	@Parameter(names={"-i","--indexed"},description="[20190905]VCF or BAM must be indexed")
	private boolean must_be_indexed=false;
	@Parameter(names={"-p","--partition"},description="For BAM files: "+SAMRecordPartition.OPT_DESC)
	private SAMRecordPartition partition = SAMRecordPartition.sample;
	@Parameter(names={"--format"},description="output format")
	private Format outputFormat = Format.tsv;

	
	private final Map<Path,SAMSequenceDictionary> ref2dict = new HashMap<>();
	
	private abstract class Reporter implements AutoCloseable {
		abstract void print(final String sample,final String type,final Path f,final String index, final Optional<Path> reference);
		@Override 
		public abstract void close();
		}
	
	private class TsvReporter extends Reporter {
		private final PrintWriter printWriter;
		
		TsvReporter(final PrintWriter pw) {
			this.printWriter = pw;
			if(show_header)
				{
				this.printWriter.print("#"+partition.name());
				this.printWriter.print('\t');
				this.printWriter.print("TYPE");
				this.printWriter.print('\t');
				this.printWriter.print("FILE");
				this.printWriter.print('\t');
				this.printWriter.print("INDEXED");
				this.printWriter.print('\t');
				this.printWriter.print("FILE_SIZE");
				this.printWriter.print('\t');
				this.printWriter.print("DATE");
				this.printWriter.print('\t');
				this.printWriter.print("REF");
				this.printWriter.println();
				}
		}
		
		@Override 
		void print(final String sample,final String type,final Path f,final String index, final Optional<Path> reference) {
			this.printWriter.print(StringUtils.isBlank(sample)?".":sample);
			this.printWriter.print('\t');
			this.printWriter.print(type);
			this.printWriter.print('\t');
			this.printWriter.print(f.toAbsolutePath());
			this.printWriter.print('\t');
			this.printWriter.print(index);

			if(Files.isRegularFile(f))
				{
				long size;
				Date modif;
				try {
					size=Files.size(f);
					}
				catch(IOException err) {
					size=-1L;
					}
				try {
					modif=new Date(Files.getLastModifiedTime(f).toMillis());
					}
				catch(final IOException err) {
					modif = null;
					}
				
				this.printWriter.print('\t');
				this.printWriter.print(size);
				this.printWriter.print('\t');
				this.printWriter.print(modif);
				}
			this.printWriter.print('\t');
			this.printWriter.print(reference.isPresent()?reference.get().toString():".");

			this.printWriter.println();			
			}
		@Override 
		public void close() {
			this.printWriter.flush();
			this.printWriter.close();
			}
		
		}

	private class XmlReporter extends Reporter {
		final OutputStream os;
		final XMLStreamWriter w;
		Path prev=null;
		XmlReporter(OutputStream os) {
			this.os = os;
			XMLOutputFactory xof= XMLOutputFactory.newInstance();
			try {
				this.w=xof.createXMLStreamWriter(os);
				this.w.writeStartDocument("UTF-8", "1.0");
				this.w.writeStartElement("summary");
			} catch (XMLStreamException e) {
				throw new RuntimeIOException(e);
				}
			}
		private void write(String tag,Object o) throws XMLStreamException {
			if(o==null) return;
			this.w.writeStartElement(tag);
			this.w.writeCharacters(String.valueOf(o));
			this.w.writeEndElement();
			}
		@Override
		void print(String sample, String type, Path path, String index, Optional<Path> reference) {
			try {
				if(prev==null || !this.prev.equals(path) ) {
					if(prev!=null) this.w.writeEndElement();
					this.w.writeStartElement(type);
					this.write("path",path.toString());
					if(reference.isPresent()) {
						this.write("reference",reference.get().toString());
						}
					if(Files.isRegularFile(path))
						{
						try {
							long size=Files.size(path);
							this.write("size",size);
							}
						catch(IOException err) {
							}
						try {
							Date modif=new Date(Files.getLastModifiedTime(path).toMillis());
							this.write("modified",modif);
							}
						catch(final IOException err) {
							}
						}
					}
				this.write("sample",sample);
				this.prev=path;
				}  
			catch (final XMLStreamException e) {
				throw new RuntimeIOException(e);
				}
			}
		@Override
		public void close() {
			try { 
				if(prev!=null) this.w.writeEndElement();
				this.w.writeEndElement();
				this.w.writeEndDocument();
				this.w.close(); } catch (XMLStreamException e) {
				throw new RuntimeIOException(e);
				}
			try { this.os.close(); } catch (IOException e) {
				throw new RuntimeIOException(e);
				}
			}
		
	}
	
    public NgsFilesSummary()
    	{
    	}		
    
  
    private Optional<Path> getKnownDict(final SAMSequenceDictionary dict) {
    	if(this.ref2dict.isEmpty()) return Optional.empty();
    	if(dict==null) return Optional.empty();
    	return this.ref2dict.entrySet().stream().
    			filter(KV->SequenceUtil.areSequenceDictionariesEqual(KV.getValue(),dict)).
    			map(KV->KV.getKey()).
    			findFirst();
    }
  
    	
    
    private void readBam(final Reporter w,final Path f)
    	{
		final SamReaderFactory srf = super.createSamReaderFactory();
		srf.validationStringency(ValidationStringency.SILENT);
    	try(SamReader r= srf.open(f)) {
			boolean indexed = r.hasIndex();
    		if(this.must_be_indexed && !indexed) return;

			
			final SAMFileHeader h=r.getFileHeader();
			final Optional<Path> faidx = getKnownDict(h.getSequenceDictionary());
			
    		final Set<String> names = this.partition.getPartitions(h);
    		
			if(!names.isEmpty())
				{
				for(final String sample:names)
					{
					w.print(sample,"BAM", f,String.valueOf(r.hasIndex()), faidx);
					}
				}
			else
				{
				w.print(null,"BAM", f,String.valueOf(r.hasIndex()),faidx);
				}
			} 
    	catch (final Exception e)
    		{
    		LOG.warning(e);
			}
    	}
   
    private void readFastq(final Reporter w,final Path f)
		{
    	//File parent=f.getParentFile();
    	//if(parent==null || super.VERBOSITY==Log.LogLevel.) return;
    	final FastQName fq=FastQName.parse(f.toFile());
    	
		
		if(!fq.isValid())
			{
			//bad name
			return;
			}
		
    	w.print(fq.getSample(),"FASTQ", f,".",Optional.empty());
		}
    
    private void readVCF(final Reporter w,final Path f)
		{
		final boolean indexed= VCFUtils.isTabixVcfPath(f) ||
				VCFUtils.isTribbleVcfPath(f) ||
				VCFUtils.isBcfIndexedPath(f)
				;
		if(this.must_be_indexed && !indexed) return;
		final VCFHeader header;
		try {
			header = VcfHeaderExtractor.decode(f);
			}
		catch(Throwable err) {
			LOG.error(err);
			return;
			}
		final Optional<Path> faidx = getKnownDict(header.getSequenceDictionary());

		final List<String> sns = header.getSampleNamesInOrder();
    	
		if(sns.isEmpty())
			{
			w.print(null,"VCF", f,String.valueOf(indexed),faidx);
			}
		else
    		{
        	for(final String sample:sns)
	        	{
	        	w.print(sample,"VCF", f,String.valueOf(indexed),faidx);
	    		}
    		}
		}

    private void scan(final Reporter reporter, final BufferedReader in) throws IOException
    	{
    	in.lines().
    		filter(L->!(StringUtil.isBlank(L) || L.startsWith("#"))).
    		map(L->Paths.get(L)).
    		filter(L->Files.isRegularFile(L)).
    		filter(L->Files.isReadable(L)).
    		forEach(L->scan(reporter,L));
    	}
    
    private void scan(final Reporter w,final Path f)
		{
		if(f==null) return;
		if(!Files.isRegularFile(f)) return;		
		if(!Files.isReadable(f)) return;		
		
		final String name=f.getFileName().toString();
		if(name.endsWith(FileExtensions.SAM)  && !this.must_be_indexed)
			{
			readBam(w,f);
			}
		else if(name.endsWith(FileExtensions.BAM) || name.endsWith(FileExtensions.CRAM))
			{
			readBam(w,f);
			}
		else if(FileExtensions.VCF_LIST.stream().anyMatch(E->name.endsWith(E)))
			{
			readVCF(w,f);
			}
		else if( (name.endsWith(".fastq") || name.endsWith(".fastq.gz") ||
				name.endsWith(".fq") || name.endsWith(".fq.gz")))
			{
			readFastq(w,f);
			}
		}

    @Override
	public int doWork(final List<String> args) {
    	Reporter  w = null;
    	try
			{
			for(Path faix:IOUtils.unrollPaths(this.faidxPathList)) {
				this.ref2dict.put(faix,SequenceDictionaryUtils.extractRequired(faix));
				}
			
			w =  this.outputFormat.equals(Format.tsv)?
				new TsvReporter(super.openPathOrStdoutAsPrintWriter(this.outputFile)):
				new XmlReporter(super.openPathOrStdoutAsStream(this.outputFile))
				;
			
			
			if(args.isEmpty())
				{
				try(BufferedReader r = new BufferedReader(new InputStreamReader(stdin()))) {
					scan(w,r);
					}
				}
			else
				{
				for(final String filename:args)
					{
					try(final BufferedReader r=IOUtils.openURIForBufferedReading(filename)) {
						scan(w,r);
						}
					}
				}
			w.close();
			w= null;
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
    	finally {
    		if(w!=null) try {w.close();} catch(Throwable err) {}
    		}
		}
	
	public static void main(final String[] args) {
		new NgsFilesSummary().instanceMainWithExit(args);
		}

}
