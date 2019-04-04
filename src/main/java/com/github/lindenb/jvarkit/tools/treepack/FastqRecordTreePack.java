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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.treepack;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
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

/**

BEGIN_DOC




### Synopsis




```
$ java -jar dist/fastqrecordtreepack.jar -c config.xml (stdin|fq1.gz fq2.gz ...) gt; out.svg
```





### XML config


XML root is <treepack>. children is '<node>' .
A '<node>' has an attribute 'name'. The text content of the <node> will be evaluated as a javascript expression with the embedded javascript engine.
The javascript engine injects record a https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/fastq/FastqRecord.html and
header a https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMFileHeader.html.



### Example



```
$ cat config.xml

<?xml version="1.0"?>
<treepack>
	<node name="length">record.length()</node>
	<node name="firstBase">(record.length()&gt;0?record.getReadString().charAt(0):null)</node>
</treepack>

```

s


```
$ curl -s "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data/NA21144/sequence_read/ERR047877.filt.fastq.gz" |\
   gunzip -c | java -jar dist/fastqrecordtreepack.jar -c) config.xml  > out.svg

```



![img](https://pbs.twimg.com/media/Bem-_tVCEAA9uT1.jpg:large)



### See also



 *  VcfTreePack
 *  http://www.cs.umd.edu/hcil/treemap-history/






END_DOC
*/
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
/**

BEGIN_DOC

### XML config

XML root is <treepack>. children is '<node>' .
A '<node>' has an attribute 'name'. The text content of the <node> will be evaluated as a javascript expression with the embedded javascript engine.
The javascript engine injects **record** a [https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/fastq/FastqRecord.html](https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/fastq/FastqRecord.html) and
**header** a [https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMFileHeader.html](https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMFileHeader.html).


### Example

```
$ cat config.xml

<?xml version="1.0"?>
<treepack>
	<node name="length">record.length()</node>
	<node name="firstBase">(record.length()&gt;0?record.getReadString().charAt(0):null)</node>
</treepack>

```
s

```
$ curl -s "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data/NA21144/sequence_read/ERR047877.filt.fastq.gz" |\
   gunzip -c | java -jar dist/fastqrecordtreepack.jar -c) config.xml  > out.svg

```


![https://pbs.twimg.com/media/Bem-_tVCEAA9uT1.jpg:large](https://pbs.twimg.com/media/Bem-_tVCEAA9uT1.jpg:large)



END_DOC

 */
@Program(name="fastqrecordtreepack",description="Create a TreeMap from one or more Fastq file. Ouput is a SVG file")
public class FastqRecordTreePack extends AbstractTreePackCommandLine
	{
	private static final Logger LOG = Logger.build(FastqRecordTreePack.class).make();


	@Parameter(names={"-c","--config"},description="XML config file")
	private File configFile = null;

	@Parameter(names={"-x","--dimension"},description="dimension of the output rectangle")
	private String dimensionStr = "1000x1000";


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
		if(this.configFile==null || !this.configFile.exists()) {
			throw new IOException("Undefined config file option");
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
			final Document dom = db.parse(this.configFile);
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
				if(att==null) throw new IOException("missing attribute 'name' in element " +e1.getTagName()+" in "+this.configFile);
				final String name=att.getValue().trim();
				if(name.isEmpty())  throw new IOException("empty attribute 'name' in element " +e1.getTagName()+" in "+this.configFile);
				
				final String content=e1.getTextContent();
				if(content==null || content.trim().isEmpty())  throw new IOException("empty text content under element " +e1.getTagName()+" in "+this.configFile);
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
	public int doWork(List<String> args) {
		setDimension(this.dimensionStr);

		
		FastqReader fqr=null;
		try
			{
			parseConfigFile();
			if(super.nodeFactoryChain.next==null) {
				LOG.error("no path defined");
				return -1;
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
			this.svg(this.outputFile);
			return RETURN_OK;
			}
		catch (Exception e)
			{
			LOG.error(e);
			return -1;
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
