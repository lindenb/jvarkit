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
package com.github.lindenb.jvarkit.tools.dict2xml;

import java.io.BufferedReader;
import java.io.OutputStream;
import java.nio.file.Path;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.dict.DictionaryXmlSerializer;
import com.github.lindenb.jvarkit.dict.SequenceDictionaryExtractor;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.SequenceUtil;

/**
BEGIN_DOC

## Motivation

extract SAM Sequence dictionaries from SAM/BAM/FASTA/VCF files and convert them to XML
Then we can use XSLT to generate code...

## Example

```
$ java -jar dist/jvarkit.jar dict2xml ~/src/jvarkit-git/src/test/resources/*.bam | xmllint --format -
chrom  start  end        path                                                                         tid  buildName  AS  M5  SP  UR
chr1   0      248956422  /home/lindenb/src/jvarkit-git/src/test/resources/ENCFF331CGL.rnaseq.b38.bam  0    GRCh38     .   .   .   .
chr2   0      242193529  /home/lindenb/src/jvarkit-git/src/test/resources/ENCFF331CGL.rnaseq.b38.bam  1    GRCh38     .   .   .   .
chr3   0      198295559  /home/lindenb/src/jvarkit-git/src/test/resources/ENCFF331CGL.rnaseq.b38.bam  2    GRCh38     .   .   .   .
chr4   0      190214555  /home/lindenb/src/jvarkit-git/src/test/resources/ENCFF331CGL.rnaseq.b38.bam  3    GRCh38     .   .   .   .
chr5   0      181538259  /home/lindenb/src/jvarkit-git/src/test/resources/ENCFF331CGL.rnaseq.b38.bam  4    GRCh38     .   .   .   .
chr6   0      170805979  /home/lindenb/src/jvarkit-git/src/test/resources/ENCFF331CGL.rnaseq.b38.bam  5    GRCh38     .   .   .   .
chr7   0      159345973  /home/lindenb/src/jvarkit-git/src/test/resources/ENCFF331CGL.rnaseq.b38.bam  6    GRCh38     .   .   .   .
chr8   0      145138636  /home/lindenb/src/jvarkit-git/src/test/resources/ENCFF331CGL.rnaseq.b38.bam  7    GRCh38     .   .   .   .
chr9   0      138394717  /home/lindenb/src/jvarkit-git/src/test/resources/ENCFF331CGL.rnaseq.b38.bam  8    GRCh38     .   .   .   .
```

```
$ find src/test/resources/ -type f -name "*.vcf.gz" | java -jar dist/jvarkit.jar dict2xml | xmllint --format -

```

END_DOC
 */
@Program(name="dict2xml",
	description="convert a SAM dictionary from vcf,sam,bam,dict, etc.. to XML.",
	keywords={"dict","xml","sam","bam","vcf"},
	creationDate="20240824",
	modificationDate="20240824",
	jvarkit_amalgamion =  true
	)
public class DictToXml extends Launcher {
	private static Logger LOG=Logger.build(DictToXml.class).make();

	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile=null;
	@Parameter(names={"--ignore-errors"},description="ignore errors, skip files that don't have a dictionary")
	private boolean skip_error = false;
	@Parameter(names={"--duplicate"},description="keep duplicates (default behavior is to keep one dictionary")
	private boolean duplicate = false;

	
	private SAMSequenceDictionary extract(SequenceDictionaryExtractor extractor,final String pathOrUrl) {
		if(!skip_error) {
			return extractor.extractRequiredDictionary(pathOrUrl);
			}
		else
			{
			return extractor.extractDictionary(pathOrUrl).orElse(null);
			}
		}
	
	@Override
	public int doWork(final List<String> args) {
		try {
			final SequenceDictionaryExtractor extractor=new SequenceDictionaryExtractor();
			final List<String> sources;
			
			if(args.isEmpty()) {
				try(BufferedReader br = IOUtils.openStdinForBufferedReader()) {
					sources = br.lines().collect(Collectors.toList());
					}
				}
			else
				{
				sources = IOUtils.unrollStrings(args);
				}
			final List<Map.Entry<SAMSequenceDictionary,String>> dictionaries=new ArrayList<>(sources.size());
			for(final String source: sources) {
				final SAMSequenceDictionary dict=  extract(extractor, source);
				if(dict==null) continue;
				if(!this.duplicate && dictionaries.stream().anyMatch(D->SequenceUtil.areSequenceDictionariesEqual(D.getKey(), dict))) {
					continue;
					}
				dictionaries.add(new AbstractMap.SimpleEntry<>(dict,source));
				}
			
			final XMLOutputFactory xof=XMLOutputFactory.newInstance();
			try (OutputStream os = super.openPathOrStdoutAsStream(outputFile)) {
				final XMLStreamWriter w = xof.createXMLStreamWriter(os);
				w.writeStartDocument("UTF-8", "1.0");
				w.writeStartElement("dictionaries");
				
				for(final Map.Entry<SAMSequenceDictionary,String> kv: dictionaries) {
					final DictionaryXmlSerializer serializer=new  DictionaryXmlSerializer() {
						@Override
						protected void writeOtherAttributes(XMLStreamWriter w, SAMSequenceDictionary dict)
								throws XMLStreamException {
							w.writeAttribute("source", kv.getValue());
							}
						};
					serializer.writeDictionary(w, kv.getKey());
					}
				
				w.writeEndDocument();		
				w.flush();
				w.close();
				os.flush();
				}
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	public static void main(final String[] args) {
		new DictToXml().instanceMainWithExit(args);
	}

}
