/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.log.Logger;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.SequenceUtil;

/**
BEGIN_DOC

## Motivation

extract SAM Sequence dictionaries from SAM/BAM/FASTA/VCF files and convert them to XML
Then we can use XSLT to generate code...

## Example

```
$ java -jar dist/jvarkit.jar dict2xml "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz" | xmllint --format - | head
<?xml version="1.0" encoding="UTF-8"?>
<dictionaries>
  <dictionary md5="c4e11bf85a6ab2f944340f409c751f32" length="3137454505" count="86" source="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz" build="GRCh37">
    <sequence name="1" length="249250621" index="0" offset="0"/>
    <sequence name="2" length="243199373" index="1" offset="249250621"/>
    <sequence name="3" length="198022430" index="2" offset="492449994"/>
    <sequence name="4" length="191154276" index="3" offset="690472424"/>
    <sequence name="5" length="180915260" index="4" offset="881626700"/>
    <sequence name="6" length="171115067" index="5" offset="1062541960"/>
    <sequence name="7" length="159138663" index="6" offset="1233657027"/>
(...)
```

```
$ find src/test/resources/ -type f -name "*.vcf.gz" | java -jar dist/jvarkit.jar dict2xml | xmllint --format - 
<?xml version="1.0" encoding="UTF-8"?>
<dictionaries>
  <dictionary md5="4677ece43eea2b029d0d33fe130ea6c7" length="3137454505" count="86" source="src/test/resources/roxan
.hs37d5.csq.vcf.gz" build="GRCh37">
    <sequence name="chr1" length="249250621" index="0" offset="0"/>
(...)
    <sequence name="hs37d5" length="35477943" index="85" offset="3101976562"/>
  </dictionary>
  <dictionary md5="bd7e0928fc3c810e48fafc53a4222ed5" length="18490" count="11" source="src/test/resources/S4.vcf.gz"
>
    <sequence name="RF01" length="3302" index="0" offset="0"/>
(...)
    <sequence name="RF10" length="751" index="9" offset="17073"/>
    <sequence name="RF11" length="666" index="10" offset="17824"/>
  </dictionary>
  <dictionary md5="9a5c58c2c91e731135b27ed14974523a" length="3101976562" count="85" source="src/test/resources/gnoma
d.exomes.r2.0.1.sites.vcf.gz" build="GRCh37">
    <sequence name="1" length="249250621" index="0" offset="0"/>
(...)
    <sequence name="NC_007605" length="171823" index="84" offset="3101804739"/>
  </dictionary>
  <dictionary md5="635de5cb51973d45844fa713ac0b7719" length="85" count="2" source="src/test/resources/toy.vcf.gz">
    <sequence name="ref" length="45" index="0" offset="0"/>
    <sequence name="ref2" length="40" index="1" offset="45"/>
  </dictionary>
  <dictionary md5="f8d942cb3fc6ebef618a0b0ba3f4ef99" length="3095677412" count="24" source="src/test/resources/gnoma
d_v2_sv.sites.vcf.gz" build="GRCh37">
    <sequence name="1" length="249250621" index="0" offset="0"/>
    <sequence name="2" length="243199373" index="1" offset="249250621"/>
    <sequence name="3" length="198022430" index="2" offset="492449994"/>
    <sequence name="4" length="191154276" index="3" offset="690472424"/>
(...)
    <sequence name="Y" length="59373566" index="23" offset="3036303846"/>
  </dictionary>
</dictionaries>

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
	private static Logger LOG=Logger.of(DictToXml.class);

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
