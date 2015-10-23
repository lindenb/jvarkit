/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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


History:
* 2015 : creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import javax.xml.namespace.NamespaceContext;
import javax.xml.namespace.QName;
import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathException;
import javax.xml.xpath.XPathExpression;
import javax.xml.xpath.XPathFactory;
import javax.xml.xpath.XPathVariableResolver;

import org.xml.sax.InputSource;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;


public class VcfFilterXPath
	extends AbstractVcfFilterXPath
	{
	
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(VcfFilterXPath.class);

	@Override
	public  Command createCommand() {
			return new MyCommand();
		}
		 
	public  static class MyCommand extends AbstractVcfFilterXPath.AbstractVcfFilterXPathCommand
		{		

	private XPathExpression xpathExpr = null;
	/** namespace mapping for xpath object */
	private final Map<String, String> prefix2uri = new HashMap<String, String>();
	/** variable mapping for xpath object */
	private final Map<QName,Object> xpathVariableMap = new HashMap<QName, Object>();

	@Override
		protected Collection<Throwable> doVcfToVcf(String inpuSource,
				VcfIterator in, VariantContextWriter out) throws IOException {
		try {
			final sun.misc.BASE64Decoder base64Decoder = new sun.misc.BASE64Decoder();
			VCFHeader header=in.getHeader();
			VCFInfoHeaderLine infoHeader = header.getInfoHeaderLine(this.infoTag);
			if(infoHeader==null)
				{
				LOG.warn("No INFO header line for "+this.infoTag+" in "+inpuSource);
				}
			else if(!(infoHeader.getCountType()==VCFHeaderLineCount.INTEGER &&
					 infoHeader.getCount()==1 &&
					 infoHeader.getType()==VCFHeaderLineType.String))
				{
				LOG.warn("Bad definition of INFO header line for "+this.infoTag+" in "+inpuSource+" expected one 'string' got "+infoHeader);
				infoHeader=null;
				}
			
			VCFHeader h2=new VCFHeader(header);
			addMetaData(header);
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
				ByteArrayInputStream xmlBytes =  new ByteArrayInputStream(base64Decoder.decodeBuffer(base64.toString()));
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
			return RETURN_OK;
			}
		catch(XPathException err)
			{
			throw new IOException(err);
			}
		finally
			{
			xpathVariableMap.clear();
			}
		}
		
	@Override
	public Collection<Throwable> initializeKnime() {
		if(this.infoTag==null || this.infoTag.isEmpty())
			{
			return wrapException("Info Tag is undefined");
			}
		if(this.xpathExpression==null || this.xpathExpression.isEmpty())
			{
			return wrapException("XPath Expression is undefined");
			}
		try {
			XPathFactory xpf = XPathFactory.newInstance();
			XPath xpath= xpf.newXPath();
			xpath.setXPathVariableResolver(new XPathVariableResolver() {
				
				@Override
				public Object resolveVariable(QName variableName) {
					return xpathVariableMap.get(variableName);
					}
				});
			xpath.setNamespaceContext(new NamespaceContext()
				{
				@Override
				public Iterator<? extends Object> getPrefixes(String namespaceURI)
					{
					List<String> L=new ArrayList<>();
					for(String pfx:prefix2uri.keySet())
						{
						if(prefix2uri.get(pfx).equals(namespaceURI))
							{
							L.add(pfx);
							}
						}

					return L.iterator();
					}
				
				@Override
				public String getPrefix(String namespaceURI)
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
		} catch (Exception e) {
			return wrapException(e);
			}
		
		for(final String s: super.namespaceStrList)
			{
			int eq=s.indexOf('=');
			if(eq<=0)
				{
				return wrapException("'=' missing in "+s);
				}
			this.prefix2uri.put(s.substring(0,eq),s.substring(eq+1));
			}
		
		return super.initializeKnime();
		}
	
	@Override
		protected Collection<Throwable> call(String inputName) throws Exception {
			return doVcfToVcf(inputName);
			}

	}
		
	public static void main(String[] args)
		{
		new VcfFilterXPath().instanceMainWithExit(args);
		}
	}
