package com.github.lindenb.jvarkit.tools.bgen.summary;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bgen.BGenHeader;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.BinaryCodec;

/**
BEGIN_DOC

```
$ java -jar dist/jvarkit.jar bgensummary ~/jeter.bgen | xmllint --format -
<?xml version="1.0" encoding="UTF-8"?>
<bgen-summary>
  <bgen-file filename="/home/lindenb/jeter.bgen" compression="e_ZstdCompression" n-variants="109745123" n-samples="490541" layout="e_Layout2" snps-offset="8339225" anonymous="true"/>
</bgen-summary>
```

END_DOC
 */
@Program(name="bgensummary",
description="bgen file summary",
keywords={"bgen"},
creationDate="20250507",
modificationDate="20250507",
jvarkit_amalgamion =  true
)
public class BGenSummary extends Launcher {
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile=null;

	private static final Logger LOG = Logger.build(BGenSummary.class).make();
	
	
	private void summary(XMLStreamWriter w,InputStream in,String path) throws IOException,XMLStreamException {
		final BinaryCodec binaryCodec = new BinaryCodec(in);
		final long snps_offset =  binaryCodec.readUInt();
		final BGenHeader header= BGenHeader.of(binaryCodec);
		
		w.writeStartElement("bgen-file");
		w.writeAttribute("filename", path);
		w.writeAttribute("compression", header.getCompression().name());
		w.writeAttribute("n-variants", String.valueOf(header.getNVariants()));
		w.writeAttribute("n-samples", String.valueOf(header.getNSamples()));
		w.writeAttribute("layout",header.getLayout().name());
		w.writeAttribute("snps-offset",String.valueOf(snps_offset));
		w.writeAttribute("anonymous",String.valueOf(header.hasAnonymousSamples()));
		if(!header.hasAnonymousSamples()) {
			w.writeStartElement("samples");
			w.writeAttribute("n-samples", String.valueOf(header.getNSamples()));
			for(int i=0;i< header.getNSamples();++i) {
				w.writeStartElement("samples");
				w.writeAttribute("idx", String.valueOf(i));
				w.writeCharacters(header.getSamples().get(i));
				w.writeEndElement();
				}
			w.writeEndElement();
			}
		
		w.writeEndElement();
		}
	@Override
	public int doWork(List<String> args) {
		try {
			XMLOutputFactory xof=XMLOutputFactory.newInstance();
			try(OutputStream os=super.openPathOrStdoutAsStream(outputFile)) {
				final XMLStreamWriter w= xof.createXMLStreamWriter(os,"UTF-8");
				w.writeStartDocument("UTF-8", "1.0");
				w.writeStartElement("bgen-summary");
				List<Path> inputs = IOUtils.unrollPaths(args);
				if(inputs.isEmpty()) {
					summary(w,stdin(),"<stdin>");
					}
				else
					{
					for(Path path : inputs) {
						try(InputStream in = Files.newInputStream(path) ) {
							summary(w,in,path.toString());
							}
						}
					}
				w.writeEndElement();
				w.writeEndDocument();
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
		new BGenSummary().instanceMainWithExit(args);
	}

}
