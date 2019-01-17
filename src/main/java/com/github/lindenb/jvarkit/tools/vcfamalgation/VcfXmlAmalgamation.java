/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.vcfamalgation;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Unmarshaller;
import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.tools.burden.VcfBurdenFilterExac;
import com.github.lindenb.jvarkit.tools.burden.VcfBurdenFisherH;
import com.github.lindenb.jvarkit.tools.burden.VcfInjectPedigree;
import com.github.lindenb.jvarkit.tools.burden.VcfMoveFiltersToInfo;
import com.github.lindenb.jvarkit.tools.misc.VcfHead;
import com.github.lindenb.jvarkit.tools.misc.VcfSetSequenceDictionary;
import com.github.lindenb.jvarkit.tools.misc.VcfTail;
import com.github.lindenb.jvarkit.tools.vcfbigwig.VCFBigWig;
import com.github.lindenb.jvarkit.tools.vcffilterjs.VcfFilterJdk;
import com.github.lindenb.jvarkit.tools.vcffilterso.VcfFilterSequenceOntology;
import com.github.lindenb.jvarkit.tools.vcftrios.VCFTrios;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.DelegateVariantContextWriter;
import com.github.lindenb.jvarkit.util.vcf.VariantContextWriterFactory;
import htsjdk.variant.vcf.VCFIterator;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
/**
BEGIN_DOC

## Example


END_DOC
**/
@Program(
		name="vcfamalgamation",
		description="Builds a complex VCF filtering engine using the Java Architecture for XML Binding  (JAXB) and the filters already defined in jvarkit.",
		keywords={"vcf","jaxb","xml","filter"},
		generate_doc=false
		)
public class VcfXmlAmalgamation extends Launcher {
	private static final Logger LOG = Logger.build(VcfXmlAmalgamation.class).make();

	private static final  Class<?> SUPPORTED_CLASSES[]=new Class[]{
			VCFTrios.CtxWriterFactory.class,
			VcfSetSequenceDictionary.CtxWriterFactory.class,
			VcfHead.CtxWriterFactory.class,
			VcfTail.CtxWriterFactory.class,
			VcfFilterSequenceOntology.CtxWriterFactory.class,
			VcfInjectPedigree.CtxWriterFactory.class,
			VcfMoveFiltersToInfo.CtxWriterFactory.class,
			VCFBigWig.CtxWriterFactory.class,
			VcfBurdenFilterExac.CtxWriterFactory.class,
			VcfBurdenFisherH.CtxWriterFactory.class,
			VcfFilterJdk.CtxWriterFactory.class
		};
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-list","--list"},description="list supported jvarkit classes and exit",help=true)
	private boolean listSupportedClasses = false;

	@ParametersDelegate
	private CtxWriterFactory component = new CtxWriterFactory();

	
	private interface PipelineInput
		{
		}
	private interface PipelineOutput
		{
		}
	

	
	
	public static class VcfReformater implements PipelineOutput
		{
		protected PrintWriter out;
		protected void println() { this.out.println();}
		protected void print(final Object o) { this.out.print(o);}
		protected void println(final Object o) { this.out.println(o);}
		protected void add(final VariantContext ctx) {}
		}
	
	//TODO
	@XmlRootElement(name="vcf-tee")
	public static class VcfTeeWriterFactory implements VariantContextWriterFactory {
		private class VcfTeeWriter extends DelegateVariantContextWriter
			{
			private final VariantContextWriterFactory teeFactory;
			VcfTeeWriter(final VariantContextWriter delegate,final VariantContextWriterFactory teeFactory) {
				super(delegate);
				this.teeFactory =  teeFactory;
				}
			
			
			}
		@Override
		public VariantContextWriter open(final VariantContextWriter delegate) {
			return new VcfTeeWriter(delegate,null);//TODO
			}
		}
	
	@XmlRootElement(name="vcf-xml-amalgation")
	@XmlAccessorType(XmlAccessType.FIELD)
	public static class CtxWriterFactory implements VariantContextWriterFactory {
		
		@Parameter(names={"-c","--xml","--config"},description=
				"XML config file. XMl root name will be ignored. Children are xml nodes",
				required=true
				)
		@XmlElement(name="xml")
		private File xmlConfigFile = null;
		
		public void setXmlConfigFile(final File xmlConfigFile) {
			this.xmlConfigFile = xmlConfigFile;
			}
		public File getXmlConfigFile() {
			return xmlConfigFile;
		}
		
		@XmlTransient
		private final List<VariantContextWriterFactory> variantContextWriterFactories = new ArrayList<>();
		
		@Override
		public int initialize() {
			try {
				IOUtil.assertFileIsReadable(this.xmlConfigFile);
				final DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
				final DocumentBuilder db = dbf.newDocumentBuilder();
				final Document dom = db.parse(this.xmlConfigFile);
					
				final Element root=  dom.getDocumentElement();
				if(root==null) {
					LOG.error("not root in xml document.");
					return -1;
					}
				
				
				final JAXBContext jxbctx = JAXBContext.newInstance(
						SUPPORTED_CLASSES
						);
				final Unmarshaller unmarshaller = jxbctx.createUnmarshaller();
				for(Node c = root.getFirstChild();
						c!=null;
						c = c.getNextSibling()
						)
					{
					if(c.getNodeType()!=Node.ELEMENT_NODE) continue;
					final Object unmarshalled;
					try {
						unmarshalled = unmarshaller.unmarshal(c);
						}
					catch(final JAXBException err2) {
						LOG.error("Cannot convert node <"+c.getNodeName()+">. supported Classes are :\n"+
								Arrays.stream(SUPPORTED_CLASSES).map(C->C.getName()).collect(Collectors.joining(";\n" ))
								,err2);
						throw err2;
						}
					if( unmarshalled == null)
						{
						LOG.error("cannot  unmarshall <"+c.getNodeName()+">");
						return -1;
						}
					if(!(unmarshalled instanceof VariantContextWriterFactory))
						{
						LOG.error("object  <"+c.getNodeName()+"> is not an instance of "+ 
								VariantContextWriterFactory.class.getName()+" but "+ unmarshalled.getClass());
						return -1;
						}
					final VariantContextWriterFactory vcwf = VariantContextWriterFactory.class.cast(unmarshalled);
					this.variantContextWriterFactories.add(vcwf);
					}
				for(final VariantContextWriterFactory vcwf: this.variantContextWriterFactories)
					{
					if(vcwf.initialize()!=0) {
						LOG.info("initialize failed");
						}
					}
				return 0;
				}
			catch(final Exception err) {
				LOG.error(err);
				return -1;
				}
			}
		
		@Override
		public VariantContextWriter open(final VariantContextWriter delegate) {
			VariantContextWriter base = delegate;
			try {
				for(int i=this.variantContextWriterFactories.size()-1;i>=0;i--) {
					base = this.variantContextWriterFactories.get(i).open(base);
					}
				}
			catch(final Exception err)
				{
				LOG.error(err);
				}
			return base;
			}
		@Override
		public void close() throws IOException {
			for(final VariantContextWriterFactory vcwf: this.variantContextWriterFactories)
				{
				CloserUtil.close(vcwf);
				}
			this.variantContextWriterFactories.clear();
			}
		}
	
	@Override
	public int doVcfToVcf(final String inputName, final VCFIterator r, final VariantContextWriter delegate)
		{
		final VariantContextWriter out = this.component.open(delegate);
		final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(r.getHeader()).logger(LOG);
		out.writeHeader(r.getHeader());
		while(r.hasNext())
			{
			out.add(progress.watch(r.next()));
			}
		out.close();
		progress.finish();
		return 0;
		}

	
	@Override
	public int doWork(final List<String> args) {
		try {
			if(this.listSupportedClasses)
				{
				final PrintWriter pw = super.openFileOrStdoutAsPrintWriter(outputFile);
				for(final Class<?> C: SUPPORTED_CLASSES)
					{
					pw.println(C.getName());
					}
				pw.flush();
				pw.close();
				return 0;
				}
			
			if(this.component.initialize()!=0) return -1;
			return doVcfToVcf(args, this.outputFile);
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.component);
			}
		}

	
	public static void main(final String[] args) {
		new VcfXmlAmalgamation().instanceMainWithExit(args);
	}
	
}
