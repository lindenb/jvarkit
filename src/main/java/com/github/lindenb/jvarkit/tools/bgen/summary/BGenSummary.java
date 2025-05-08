package com.github.lindenb.jvarkit.tools.bgen.summary;

import java.io.IOException;
import java.io.OutputStream;
import java.nio.file.Path;
import java.util.List;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bgen.BGenHeader;
import com.github.lindenb.jvarkit.bgen.BGenReader;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;


/**
BEGIN_DOC

## Example

```
$ java -jar dist/jvarkit.jar bgensummary *.bgen | xmllint --format -
$ java -jar dist/jvarkit.jar bgensummary ~/src/jbgen/src/test/resources/htsjdk/tribble/bgen/*bgen | xmllint --format -
<?xml version="1.0" encoding="UTF-8"?>
<bgen-summary>
  <bgen-file filename=complex.bgen" compression="e_ZlibCompression" n-variants="10" n-samples="4" layout="e_Layout2" snps-offset="68" anonymous="false">
    <samples n-samples="4">
      <samples idx="0">sample_0</samples>
      <samples idx="1">sample_1</samples>
      <samples idx="2">sample_2</samples>
      <samples idx="3">sample_3</samples>
    </samples>
  </bgen-file>
  <bgen-file filename="example.v11.bgen" compression="e_ZlibCompression" n-variants="199" n-samples="500" layout="e_Layout1" snps-offset="5920" anonymous="true"/>
  <bgen-file filename="haplotypes.bgen" compression="e_ZlibCompression" n-variants="4" n-samples="4" layout="e_Layout2" snps-offset="68" anonymous="false">
    <samples n-samples="4">
      <samples idx="0">sample_0</samples>
      <samples idx="1">sample_1</samples>
      <samples idx="2">sample_2</samples>
      <samples idx="3">sample_3</samples>
    </samples>
  </bgen-file>
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
	@Parameter(names={"-n"},description="print the summary for 'n' variants")
	private long n_variants=10;
	@Parameter(names={"-g"},description="print the genotypes for 'g' genotypes")
	private long n_genotypes=3;

	private static final Logger LOG = Logger.build(BGenSummary.class).make();
	
	
	private void summary(XMLStreamWriter w,BGenReader reader,String path) throws IOException,XMLStreamException {
		final BGenHeader header= reader.getHeader();
		
		w.writeStartElement("bgen-file");
		w.writeAttribute("filename", path);
		w.writeAttribute("compression", header.getCompression().name());
		w.writeAttribute("n-variants", String.valueOf(header.getNVariants()));
		w.writeAttribute("n-samples", String.valueOf(header.getNSamples()));
		w.writeAttribute("layout",header.getLayout().name());
		w.writeAttribute("snps-offset",String.valueOf(reader.getSnpsOffset()));
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
		if(this.n_variants>0) {
			w.writeStartElement("variants");
			for(long i=0;i< this.n_variants;i++) {
				BGenReader.Variant ctx = reader.readVariant();
				if(ctx==null) break;
				if(n_genotypes==0) {
					reader.skipGenotypes();
					}
				else
					{
					ctx = reader.readGenotypes();
					}
				w.writeStartElement("variant");
				w.writeAttribute("idx", String.valueOf(i));
				w.writeAttribute("chrom", ctx.getContig());
				w.writeAttribute("pos", String.valueOf(ctx.getPosition()));
				if(!StringUtils.isBlank(ctx.getId())) {
					w.writeAttribute("id", String.valueOf(ctx.getId()));
					}
				if(!StringUtils.isBlank(ctx.getRsId())) {
					w.writeAttribute("rsid", String.valueOf(ctx.getRsId()));
					}
				w.writeStartElement("alleles");
				w.writeAttribute("count", String.valueOf(ctx.getNAlleles()));
				for(int k=0;k< ctx.getNAlleles();++k) {
					w.writeStartElement("allele");
					w.writeAttribute("idx", String.valueOf(k));
					w.writeCharacters(ctx.getAllele(k));
					w.writeEndElement();
					}
				w.writeEndElement();
				
				if(n_genotypes>0) {
					w.writeStartElement("genotypes");
					for(int k=0;k< ctx.getGenotypes().size() && k< n_genotypes;++k) {
						final BGenReader.Genotype gt= ctx.getGenotype(k);
						w.writeStartElement("genotype");
						w.writeAttribute("idx", String.valueOf(k));
						w.writeAttribute("ploidy", String.valueOf(gt.getPloidy()));
						w.writeAttribute("phased", String.valueOf(gt.isPhased()));
						w.writeAttribute("missing", String.valueOf(gt.isMissing()));
						if(!gt.isMissing()) {
							for(double v: gt.getProbs()) {
								w.writeStartElement("value");
								w.writeCharacters(String.valueOf(v));
								w.writeEndElement();
								}
							}
						w.writeEndElement();
						}
					w.writeEndElement();
					}
				
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
				final List<String> inputs = IOUtils.unrollStrings(args);
				if(inputs.isEmpty()) {
					try(BGenReader r=new BGenReader(System.in)) {
						summary(w,r,"<stdin>");
						}
					}
				else
					{
					for(String path : inputs) {
						try(BGenReader r=new BGenReader(path)) {
							summary(w,r,path);
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
