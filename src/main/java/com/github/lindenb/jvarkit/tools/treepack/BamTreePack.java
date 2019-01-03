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
package com.github.lindenb.jvarkit.tools.treepack;

import java.io.File;
import java.io.IOException;
import java.util.List;

import javax.script.CompiledScript;
import javax.script.ScriptException;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.w3c.dom.Attr;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;

import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.util.CloserUtil;

/**

BEGIN_DOC




### Synopsis




```
$ $ java -jar dist/bamtreemap.jar -c config.xml (stdin|bam1 bam2 ....) > out.svg
```





### XML config


XML root is <treepack>. children is '<node>' .
A '<node>' has an attribute 'name'. The text content of the <node> will be evaluated as a javascript expression with the embedded javascript engine.
The javascript engine injects record a https://github.com/samtools/htsjdk/blob/master/src/java/htsjdk/samtools/SAMRecord.java and
header a https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMFileHeader.html.





### Example




```
$ cat config.xml

<?xml version="1.0"?>
<treepack>
  <node name="chr">(record.getReadUnmappedFlag()?"UNMAPPED":record.getContig())</node>
  <node name="mapq">(record.getReadUnmappedFlag()?null:record.getMappingQuality())</node>
</treepack>


```





```
$ java  -jar dist/bamtreepack.jar -c config.xml  \
  "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/NA12340/alignment/NA12340.mapped.ILLUMINA.bwa.CEU.low_coverage.20101123.bam" \
  "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/NA12273/alignment/NA12273.mapped.ILLUMINA.bwa.CEU.low_coverage.20101123.bam" > out.svg

```



![img](https://pbs.twimg.com/media/BeaCYQgCEAAZxio.png:large)



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
The javascript engine injects **record** a [https://github.com/samtools/htsjdk/blob/master/src/java/htsjdk/samtools/SAMRecord.java](https://github.com/samtools/htsjdk/blob/master/src/java/htsjdk/samtools/SAMRecord.java) and
**header** a [https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMFileHeader.html](https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMFileHeader.html).




### Example


```
$ cat config.xml

<?xml version="1.0"?>
<treepack>
  <node name="chr">(record.getReadUnmappedFlag()?"UNMAPPED":record.getContig())</node>
  <node name="mapq">(record.getReadUnmappedFlag()?null:record.getMappingQuality())</node>
</treepack>


```



```
$ java  -jar dist/bamtreepack.jar -c config.xml  \
  "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/NA12340/alignment/NA12340.mapped.ILLUMINA.bwa.CEU.low_coverage.20101123.bam" \
  "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/NA12273/alignment/NA12273.mapped.ILLUMINA.bwa.CEU.low_coverage.20101123.bam" > out.svg

```


![https://pbs.twimg.com/media/BeaCYQgCEAAZxio.png:large](https://pbs.twimg.com/media/BeaCYQgCEAAZxio.png:large)


END_DOC

*/
@Program(name="bamtreepack",description="Create a TreeMap from one or more SAM/BAM file. Ouput is a SVG file.",
keywords={"bam","treepack"})
public class BamTreePack extends AbstractTreePackCommandLine
	{
	private static final Logger LOG = Logger.build(BamTreePack.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;


	@Parameter(names={"-c","--config"},description="XML config file")
	private File configFile = null;

	@Parameter(names={"-x","--dimension"},description="dimension of the output rectangle")
	private String dimensionStr = "1000x1000";

	
	public BamTreePack()
		{
		
		}
	
	
	 private void scan(final SamReader sfr)
		 {
		 super.bindings.put("header", sfr.getFileHeader());
		 final SAMRecordIterator iter=sfr.iterator();
		 final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(sfr.getFileHeader());
		 while(iter.hasNext())
			 {
			 final SAMRecord rec= progress.watch(iter.next());
			 super.bindings.put("record", rec);
			 super.nodeFactoryChain.watch(rootNode, rec);
			 }
		 progress.finish();
		 }
	  
	
		private void parseConfigFile() throws IOException{
			if(this.configFile==null || !this.configFile.exists()) {
				throw new IOException("Undefined config file ");
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
		
		
		SamReader in=null;
		try
			{
			parseConfigFile();
			if(super.nodeFactoryChain.next==null)
				{
				LOG.error("no path defined");
				return -1;
				}

			
			final SamReaderFactory srf = super.createSamReaderFactory();
			if(args.isEmpty())
				{
				LOG.info("Reading stdin");
				in = srf.open(SamInputResource.of(stdin()));
				scan(in);
				CloserUtil.close(in);in=null;
				LOG.info("Done stdin");
				}
			else
				{
				for(final String filename:args)
					{
					LOG.info("Reading "+filename);
					in= srf.open(SamInputResource.of(filename));
					scan(in);
					LOG.info("Done "+filename);
					CloserUtil.close(in);
					}
				}
			this.layout();
			this.svg(this.outputFile);
			return RETURN_OK;
			}
		catch (final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(in);
			}
		}

	
	public static void main(String[] args)
		{
		new BamTreePack().instanceMainWithExit(args);

		}

	}
