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
package com.github.lindenb.jvarkit.tools.htsvelocity;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringReader;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;

import org.apache.velocity.VelocityContext;
import org.apache.velocity.app.VelocityEngine;
import org.apache.velocity.context.Context;
import org.apache.velocity.runtime.RuntimeConstants;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.json.FromJson;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.samtools.util.IntervalParser;
import com.google.gson.stream.JsonReader;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;

/**
BEGIN_DOC

todo

END_DOC
*/

@Program(
	name="htsfreemarker",
	description="Apply apache velocity to VCF/BAM/JSON files.",
	keywords={"template","velocuty","vcf","bam"},
	creationDate = "20230616",
	modificationDate = "20230616",
	jvarkit_amalgamion = true
	)
public class HtsVelocity extends Launcher {
	private static final Logger LOG = Logger.of(HtsVelocity.class);
	@Parameter(names={"--output","-o"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"--reference","-R"},description="For reading CRAM." + INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private Path reference = null;
	@Parameter(names={"--json"},description="<name> <json-file> . Add this pair to the freemarker context.",arity = 2)
	private List<String> input_json_file = new ArrayList<>();
	@Parameter(names={"--json-string"},description="<name> <json-string> . Add this pair to the freemarker context.",arity = 2)
	private List<String> input_json_string = new ArrayList<>();
	@Parameter(names={"--string"},description="<name> <string>. Add this pair to the freemarker context.",arity = 2)
	private List<String> input_string = new ArrayList<>();
	@Parameter(names={"--vcf"},description="<name> <vcf-file>. Add this VCF to the freemarker context. Object created is [header:object,variants:list]",arity = 2)
	private List<String> input_vcf_files = new ArrayList<>();
	@Parameter(names={"--bam","--sam"},description="<name> <bam-file>. Add this Bam to the freemarker context. Object created is [header:object,reads:list]",arity = 2)
	private List<String> input_bam_files = new ArrayList<>();
	@Parameter(names={"--interval"},description="Restrict BAM/VCF/... file(s) to that interval. Files must be indexed.")
	private String interval_str = null;

	
	private Object remapBam(final String filename) {
		final Map<String,Object> hash = new LinkedHashMap<>();
		final SamReaderFactory srf = super.createSamReaderFactory();
		
		if(reference!=null) srf.referenceSequence(reference);
		
		try(SamReader reader = srf.open(Paths.get(filename))) {
			final SAMFileHeader header=reader.getFileHeader();
			hash.put("header", header);
			CloseableIterator<SAMRecord> iter;
			if(StringUtils.isBlank(this.interval_str))
					{
					iter = reader.iterator();
					}
			else
					{
					final Locatable loc = new IntervalParser(header.getSequenceDictionary()).apply(this.interval_str).get();
					iter = reader.query(loc.getContig(), loc.getStart(), loc.getEnd(), false);
					}
			hash.put("reads", iter.stream().collect(Collectors.toList()));
			iter.close();
			}
		catch(final IOException err) {
			throw new RuntimeIOException(err);
			}
		return hash;
		}

	
	private Object remapVcf(final String filename) {
		final Map<String,Object> hash = new LinkedHashMap<>();
		try(VCFReader reader = new VCFFileReader(Paths.get(filename),!StringUtils.isBlank(this.interval_str))) {
			final VCFHeader header = reader.getHeader();
			hash.put("header", header);
			try(CloseableIterator<VariantContext> iter=StringUtils.isBlank(this.interval_str)?
					reader.iterator():
					reader.query(new IntervalParser(header.getSequenceDictionary()).apply(this.interval_str).orElseThrow(IntervalParser.exception("Cannot open VCF for this interval "+ this.interval_str)))) {
				hash.put("variants", iter.stream().collect(Collectors.toList()));
				}
			}
		catch(final IOException err) {
			throw new RuntimeIOException(err);
			}
		return hash;
		}
	
	private Object remapJson(JsonReader reader) throws IOException {
        return new FromJson().parse(reader);
		}
	
	private void insertPairs(final List<String> pairs,final Map<String, Object> dest, final Function<String,Object> to_object) {
		for(int i=0;i< pairs.size();i+=2) {
			final String key = pairs.get(i);
			if(dest.containsKey(key)) {
				LOG.warn("Freemarker context for key \""+key+"\" was defined twice.");
				}
			final Object v = to_object.apply(pairs.get(i+1));
			dest.put(key, v);
			}
		}
	
	@Override
	public int doWork(List<String> args) {
		try {
			
			final String templateName = oneAndOnlyOneFile(args);
			
			final Context context = new VelocityContext();
			final File file=new File(templateName);
            IOUtil.assertFileIsReadable(file);
            
            final VelocityEngine ve = new VelocityEngine();
            ve.setProperty(RuntimeConstants.RESOURCE_LOADER, "file");
            ve.setProperty("file.resource.loader.class","org.apache.velocity.runtime.resource.loader.FileResourceLoader");
            //ve.setProperty("file.resource.loader.path",Collections.singletonList(file.getParent()==null?new File("."):file.getParent()));
            ve.setProperty("file.resource.loader.path",file.getParent()==null?".":file.getParent().toString());
           
			
			
	        final Map<String, Object> hashParams = new HashMap<String, Object>();
			insertPairs(input_json_file,hashParams, S->{
				try(JsonReader jsr = new JsonReader(new FileReader(S))) {
					return remapJson(jsr);
					}
				catch(IOException err) {
					throw new RuntimeIOException(err);
					}
				});
			insertPairs(input_json_string,hashParams, S->{
				try(JsonReader jsr = new JsonReader(new StringReader(S))) {
					return remapJson(jsr);
					}
				catch(IOException err) {
					throw new RuntimeIOException(err);
					}
				});
			insertPairs(input_string,hashParams, S->S);
			insertPairs(this.input_vcf_files,hashParams, S->remapVcf(S));
			insertPairs(this.input_bam_files,hashParams, S->remapBam(S));

		    
		    try(PrintWriter out = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
		    	 final org.apache.velocity.Template template = ve.getTemplate(file.getName());
                 template.merge( context, out);
			    out.flush();
			    }
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
	    
		return 0;
		}
	
	public static void main(String[] args) {
		new HtsVelocity().instanceMainWithExit(args);
	}

}
