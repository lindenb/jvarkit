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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.util.ArrayList;
import java.util.Base64;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import javax.xml.namespace.NamespaceContext;
import javax.xml.namespace.QName;
import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathExpression;
import javax.xml.xpath.XPathFactory;
import javax.xml.xpath.XPathVariableResolver;

import org.xml.sax.InputSource;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import htsjdk.variant.vcf.VCFIterator;

@Program(name="vcffilterxpath",
	description="Filter a VCF with a XPATH expression on a INFO tag containing a base64 encodede xml document",
	keywords={"vcf","xml","xpath"}
	)
public class VcfFilterXPath
	extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfFilterXPath.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	/** the INFO tag to use in the VCF input */
	@Parameter(names="-T",description=" (info tag) INFO tag containing a base64-encoded XML document")
	private String infoTag=null;
	/** user xpath expression */
	@Parameter(names="-x",description="(xpath) XPath expression")

	private String xpathExpression=null;
	/** compiled XPath expression */
	private XPathExpression xpathExpr=null;
	/** namespace mapping for xpath object */
	private final Map<String, String> prefix2uri = new HashMap<String, String>();
	@Parameter(names="-n",description="prefix=uri (add this namespace mapping to xpath context)")
	private List<String> __prefix2uri = new ArrayList<String>();

	
	/** variable mapping for xpath object */
	private final Map<QName,Object> xpathVariableMap = new HashMap<QName, Object>();
	
	public VcfFilterXPath()
		{
		}
	
	public void setInfoTag(String infoTag) {
		this.infoTag = infoTag;
		}
	public void setXpathExpression(String xpathExpression) {
		this.xpathExpression = xpathExpression;
		}
	
	@Override
	protected int doVcfToVcf(String inputName, VCFIterator in, VariantContextWriter out) {
		try {
			//TODO in jdk8 replace with http://docs.oracle.com/javase/8/docs/api/java/util/Base64.html
			VCFHeader header=in.getHeader();
			VCFInfoHeaderLine infoHeader = header.getInfoHeaderLine(this.infoTag);
			if(infoHeader==null)
				{
				LOG.warning("No INFO header line for "+this.infoTag+" in "+inputName);
				}
			else if(!(infoHeader.getCountType()==VCFHeaderLineCount.INTEGER &&
					 infoHeader.getCount()==1 &&
					 infoHeader.getType()==VCFHeaderLineType.String))
				{
				LOG.warning("Bad definition of INFO header line for "+this.infoTag+" in "+inputName+" expected one 'string' got "+infoHeader);
				infoHeader=null;
				}
			
			final VCFHeader h2=new VCFHeader(header);
			h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
			h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
			JVarkitVersion.getInstance().addMetaData(getClass().getSimpleName(), h2);			
			
			final SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
			
			out.writeHeader(h2);
			while(in.hasNext() )
				{	
				VariantContext ctx = progess.watch(in.next());
				if(infoHeader==null)//no tag in header
					{
					out.add(ctx);
					continue;
					}
				
				Object o=ctx.getAttribute(this.infoTag);
				if(o==null)
					{
					out.add(ctx);
					continue;
					}
				StringBuilder base64=new StringBuilder(o.toString());
				while(base64.length()%4!=0) base64.append('=');
				ByteArrayInputStream xmlBytes =  new ByteArrayInputStream(Base64.getDecoder().decode(base64.toString()));
				InputSource inputSource=new InputSource(xmlBytes);
				xpathVariableMap.put(new QName("chrom"), ctx.getContig());
				xpathVariableMap.put(new QName("start"), ctx.getStart());
				xpathVariableMap.put(new QName("end"), ctx.getEnd());
				xpathVariableMap.put(new QName("id"), ctx.hasID()?ctx.getID():".");
				boolean accept=(Boolean)xpathExpr.evaluate(inputSource, XPathConstants.BOOLEAN);
				if(accept)
					{
					out.add(ctx);
					}
				if(out.checkError()) break;
				}
			
			progess.finish();
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			xpathVariableMap.clear();
			}
		}
	
	@Override
	public int doWork(List<String> args) {
		if(this.infoTag==null || this.infoTag.isEmpty())
			{
			LOG.error("Info Tag is undefined");
			return -1;
			}
		if(this.xpathExpression==null || this.xpathExpression.isEmpty())
			{
			LOG.error("XPath Expression is undefined");
			return -1;
			}
		
		for(final String s:__prefix2uri) {
			int eq=s.indexOf('=');
			if(eq<=0)
				{
				LOG.error("'=' missing in "+s);
				return -1;
				}
			
			this.prefix2uri.put(s.substring(0,eq),s.substring(eq+1));
			}

	try {
		XPathFactory xpf = XPathFactory.newInstance();
		XPath xpath= xpf.newXPath();
		xpath.setXPathVariableResolver(new XPathVariableResolver()
			{
			@Override
			public Object resolveVariable(QName qname)
				{
				return xpathVariableMap.get(qname);
				}
			});
		xpath.setNamespaceContext(new NamespaceContext()
			{
			@Override
			public Iterator<String> getPrefixes(String namespaceURI)
				{
				final List<String> L=new ArrayList<>();
				for(final String pfx:prefix2uri.keySet())
					{
					if(prefix2uri.get(pfx).equals(namespaceURI))
						{
						L.add(pfx);
						}
					}

				return L.iterator();
				}
			
			@Override
			public String getPrefix(final String namespaceURI)
				{
				for(String pfx:prefix2uri.keySet())
					{
					if(prefix2uri.get(pfx).equals(namespaceURI))
						{
						return pfx;
						}
					}
				return null;
				}
			
			@Override
			public String getNamespaceURI(String prefix)
				{
				return prefix2uri.get(prefix);
				}
			});
		this.xpathExpr=xpath.compile(this.xpathExpression);
		return doVcfToVcf(args, outputFile);
		}
	catch(Exception er) {
		LOG.error(er);
		return -1;
		}
	
	}
		
		
		
	public static void main(String[] args)
		{
		new VcfFilterXPath().instanceMainWithExit(args);
		}
	}
