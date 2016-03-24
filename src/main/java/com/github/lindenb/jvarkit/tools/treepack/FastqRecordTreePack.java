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
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.treepack;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.Collection;
import java.util.List;

import javax.script.CompiledScript;
import javax.script.ScriptException;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.w3c.dom.Attr;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.FastqReader;
import com.github.lindenb.jvarkit.util.picard.FourLinesFastqReader;

import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.CloserUtil;

public class FastqRecordTreePack extends AbstractFastqRecordTreePack
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(FastqRecordTreePack.class);

	public FastqRecordTreePack()
		{
		
		}
	
	
	 private void scan(FastqReader iter)
		 {
		 while(iter.hasNext())
			 {
			 final FastqRecord rec=iter.next();
			 super.bindings.put("record", rec);
			 super.nodeFactoryChain.watch(rootNode,rec);
			 }
		 }

	private void parseConfigFile() throws IOException{
		if(super.configFile==null || !super.configFile.exists()) {
			throw new IOException("Undefined config file option -"+OPTION_CONFIGFILE);
		}
		try {
			LOG.info("getting javascript manager");
			final javax.script.ScriptEngineManager manager = new javax.script.ScriptEngineManager();
			final javax.script.ScriptEngine engine = manager.getEngineByName("js");
			if(engine==null)
				{
				throw new RuntimeException("not available ScriptEngineManager: javascript. Use the SUN/Oracle JDK ?");
				}
			final javax.script.Compilable compilingEngine = (javax.script.Compilable)engine;
			
			final DocumentBuilderFactory dbf=DocumentBuilderFactory.newInstance();
			final DocumentBuilder db=dbf.newDocumentBuilder();
			LOG.info("parsing "+configFile);
			final Document dom = db.parse(super.configFile);
			final Element root=dom.getDocumentElement();
			if(root==null) throw new RuntimeException("not root node in "+this.configFile);
			if(!root.getTagName().equals("treepack"))
				{
				throw new IOException("Bad root in "+configFile+" expected root node <treepack> but got "+root.getTagName());
				}
			Attr att=root.getAttributeNode("width");
			if(att!=null) {
				super.viewRect.width=Math.max(100, Integer.parseInt(att.getValue()));
			}
			att=root.getAttributeNode("height");
			if(att!=null) {
					super.viewRect.height=Math.max(100, Integer.parseInt(att.getValue()));
				}
			for(org.w3c.dom.Node c=root.getFirstChild();c!=null;c=c.getNextSibling()) {
				if(c.getNodeType()!=Node.ELEMENT_NODE) continue;
				Element e1=Element.class.cast(c);
				if(e1.getTagName().equals("node")) {
				att= e1.getAttributeNode("name");
				if(att==null) throw new IOException("missing attribute 'name' in element " +e1.getTagName()+" in "+super.configFile);
				final String name=att.getValue().trim();
				if(name.isEmpty())  throw new IOException("empty attribute 'name' in element " +e1.getTagName()+" in "+super.configFile);
				
				final String content=e1.getTextContent();
				if(content==null || content.trim().isEmpty())  throw new IOException("empty text content under element " +e1.getTagName()+" in "+super.configFile);
				CompiledScript compiled =null;
				try {
					compiled=compilingEngine.compile(content);
					}
				catch(ScriptException err) {
					throw new IOException("Cannot compile node "+e1.getTagName(), err);
				}

				final JsNodeFactory js=new JsNodeFactory(name,compiled); 
				super.nodeFactoryChain.append(js);
				}
				else {
					throw new IOException("Unknown node under root "+e1);
				}
			}
			
			} 
		catch(IOException err) {
			throw err;
		}
		catch(Exception err) {
			throw new IOException(err);
		}
		finally {
			}
		}

	
	@Override
	public Collection<Throwable> call() throws Exception {
		
		setDimension(super.dimensionStr);

		
		FastqReader fqr=null;
		final List<String> args = getInputFiles();
		try
			{
			parseConfigFile();
			if(super.nodeFactoryChain.next==null) {
			return wrapException("no path defined");
			}
			
			if(args.isEmpty())
				{
				LOG.info("Reading stdin");
				fqr=new FourLinesFastqReader(stdin());
				scan(fqr);
				CloserUtil.close(fqr);
				LOG.info("Done stdin");
				}
			else
				{
				for(final String filename:args)
					{
					InputStream in=null;
					LOG.info("Reading "+filename);
					if(IOUtils.isRemoteURI(filename))
						{
						in=IOUtils.openURIForReading(filename);
						fqr=new FourLinesFastqReader(in);
						}
					else
						{
						fqr=new FourLinesFastqReader(new File(filename));
						}
					scan(fqr);
					LOG.info("Done "+filename);
					CloserUtil.close(fqr);
					CloserUtil.close(in);
					}
				}
			
			this.layout();
			this.svg();
			return RETURN_OK;
			}
		catch (Exception e)
			{
			return wrapException(e);
			}
		finally
			{
			CloserUtil.close(fqr);
			}
		}
	
	public static void main(String[] args)
		{
		new FastqRecordTreePack().instanceMainWithExit(args);
		}

	}
