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
import java.util.List;

import javax.script.CompiledScript;
import javax.script.ScriptException;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.w3c.dom.Attr;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import htsjdk.variant.vcf.VCFIterator;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;


/**

BEGIN_DOC




### Synopsis




```
$ java -jar dist/vcftreepack.jar -c config.xml (stdin|vcf1 vcf2 vcf3) > out.svg
```





### XML config


XML root is <treepack>. children are '<node>' or '<sample>'.
A '<node>' has an attribute 'name'. The text content of the <node> will be evaluated as a javascript expression with the embedded javascript engine.
The javascript engine injects :

 *  variant a https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/VariantContext.html
 *  header a https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/vcf/VCFHeader.html.
 *  style a java.util.concurrent.atomic.AtomicReference<String> that, if it is not empty, will be used as the SVG/CSS style of the current node.


A special node '<sample>' can be called once. <sample> will generate a new node splitting by sample name and injects a new variable 'genotype' a https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/Genotype.html .
A special node '<alt>' can be called once. <sample> will generate a new node splitting by ALT allele name and injects a new variable 'alt' a https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/Allele.html .



### Example




```
$ cat config.xml

<?xml version="1.0"?>
<treepack>
  <samples/>
  <node name="ref">variant.getReference()</node>
  <node name="alt">variant.getAlternateAlleles()</node>
</treepack>

```






```

<?xml version="1.0"?>
<treepack>
  <samples/>
  <node name="indels">
  function fun() {
	  if(variant.isSNP()) return null;
	  if(!genotype.isCalled() || genotype.isHomRef()) return null;
	  return "INDEL";
	  }
fun(); 
  </node>
</treepack>

$ echo 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT |\
tr " " "\n" |\
awk '{printf("ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/ByPopulation/YRI-1412-%s.vcf.gz\n",$1);}' |\
xargs java -jar dist/vcftreepack.jar -c config.xml -x 5000x5000 > output.svg

```



![img](https://pbs.twimg.com/media/BeR15u-CAAAwE9Z.png:large)



### See also



 *  BamTreePack
 *  FastqRecordTreePack
 *  http://www.cs.umd.edu/hcil/treemap-history/






END_DOC
*/


@Program(name="vcftreepack",description="Create a TreeMap from one or more VCF. Ouput is a SVG file.")
public class VcfTreePack extends  AbstractTreePackCommandLine
	{
	
	private static final Logger LOG = Logger.build(VcfTreePack.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;


	@Parameter(names={"-c","--config"},description="XML config file")
	private File configFile = null;

	@Parameter(names={"-x","--dimension"},description="dimension of the output rectangle")
	private String dimensionStr = "1000x1000";

	
	
	private class SampleNodeFactory extends NodeFactory {
		@Override
		public String getName() {
			return "sample";
			}
		
		@Override
		public void watch(TreeNode parentNode, Object o) {
			if(parentNode.factory!=this.prev) throw new IllegalStateException();

			if(!(o instanceof VariantContext)) throw new IllegalStateException();
			final VariantContext ctx=VariantContext.class.cast(o);
			for(final String key:ctx.getSampleNames()) {
				TreeNode child = increment(parentNode, key);
				final Genotype g= ctx.getGenotype(key);
				bindings.put("genotype", g);
				if(next!=null) next.watch(child, ctx);
			}
			bindings.remove("genotype");
		}
	}

	private class AltNodeFactory extends NodeFactory {
		@Override
		public String getName() {
			return "alt";
			}
		
		@Override
		public void watch(TreeNode parentNode, Object o) {
			if(parentNode.factory!=this.prev) throw new IllegalStateException();

			if(!(o instanceof VariantContext)) throw new IllegalStateException();
			final VariantContext ctx=VariantContext.class.cast(o);
			for(final Allele alt:ctx.getAlleles()) {
				final String key=alt.getDisplayString();
				TreeNode child = increment(parentNode, key);
				bindings.put("alt", alt);
				if(next!=null) next.watch(child,ctx);
			}
			bindings.remove("alt");
		}
	}


	void scan(final VCFIterator iter)
		 {
		 final VCFHeader header=iter.getHeader();
		 super.bindings.put("header", header);
		 final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
		 while(iter.hasNext())
			 {
			 final VariantContext ctx=progress.watch(iter.next());
			 bindings.put("variant",ctx);
			 super.nodeFactoryChain.watch(super.rootNode,ctx);
			 }
		 progress.finish();
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
			boolean found_sample=false;
			boolean found_alt=false;
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
				else if(!found_sample && (e1.getTagName().equals("samples") || e1.equals("genotypes"))) {
					found_sample=true;
					super.nodeFactoryChain.append(new SampleNodeFactory());
				} else if(!found_alt && (e1.getTagName().equals("alts") || e1.equals("alternate"))) {
					found_alt=true;
					super.nodeFactoryChain.append(new AltNodeFactory());
				} else
				{
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
		
		
		
		VCFIterator in=null;
 		try
			{
			parseConfigFile();
			
			if(super.nodeFactoryChain.next==null)
				{
				LOG.error("no path defined");
				return -1;
				}
			
			if(args.isEmpty())
				{
				LOG.info("Reading stdin");
				in=VCFUtils.createVCFIteratorFromStream(stdin());
				scan(in);
				CloserUtil.close(in);
				}
			else
				{
				for(final String f: args)
					{
					in=VCFUtils.createVCFIterator(f);
					scan(in);
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

	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new VcfTreePack().instanceMainWithExit(args);
		}

	}
