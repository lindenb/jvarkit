package com.github.lindenb.jvarkit.util.jcommander;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.Arrays;
import java.util.stream.Collectors;

import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathFactory;

import org.w3c.dom.Document;
import org.w3c.dom.DocumentFragment;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParameterDescription;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloserUtil;



/** need to decouple from Launcher for JXF applications that cannot extends 'Launcher' */
public  class CmdUsageBuilder
	{
	static final Logger LOG=Logger.build( CmdUsageBuilder.class).make();

	/** main class */
	private final Class<?> _mainClass;
	
	@Parameter(names = {"-h","--help"},description="print help and exit", help = true)
	public boolean print_help = false;
	@Parameter(names = {"--helpFormat"},description="What kind of help. One of [usage,markdown,xml].", help = true)
	public String helpFormatStr = "usage";
	@Parameter(names = {"--version"}, help = true,description="print version and exit")
	public boolean print_version = false;

	public CmdUsageBuilder(final Class<?> mainClass)
		{
		this._mainClass = mainClass;
		}
	
	/** return Program annotation associated to main class. May be null*/
	public Program getProgram() {
		return getMainClass().getAnnotation(Program.class);
		}
	
	/** return true if main class has Program annotation */
	public boolean hasProgram() {
		return getProgram()!=null;
		}
	
	public boolean shouldPrintUsage() {
		return this.print_help;
		}
	
	
	public void usage(final JCommander jc) {
			/* usage is XML ******************** */
			if(this.helpFormatStr.equalsIgnoreCase("xml"))
				{
				final Document dom;
				try {
					dom = this.xmlUsage(jc);
					}
				catch(final Exception err)
					{
					err.printStackTrace();
					LOG.warn("An error occured. Cannot produce usage as XML.");
					return ;
					}
	
				try {
					final Transformer transformer = TransformerFactory.newInstance().newTransformer();
					final DOMSource source = new DOMSource(dom);
					final StreamResult result = new StreamResult(new StringWriter());
					transformer.setOutputProperty(OutputKeys.STANDALONE,"yes");
					transformer.transform(source, result);
					final String xmlString = result.getWriter().toString();
					System.out.println(xmlString);
					}
				catch(final Exception err)
					{
					err.printStackTrace();
					LOG.warn("An error occured. Cannot produce usage as XML/XSLT.");
					}
				}
			else if(this.helpFormatStr.equalsIgnoreCase("markdown"))
				{
				final StringBuilder sb=new StringBuilder();
				this.markdown(jc,sb);
				System.out.println( sb );
				}
			else if(this.helpFormatStr.equals("make-doc")) {
				this.updateDoc(jc);
				}
			else
				{
				final StringBuilder sb=new StringBuilder();
				this.standard(jc,sb);
				System.out.println( sb );
				}
			}

		
		public String getUsage(final JCommander jc) {
			final StringBuilder sb=new StringBuilder();
			this.standard(jc,sb);
			return sb.toString();
			}
		public String hyperlink(final String url)
			{
			return "["+url+"]("+url+")";
			}
		
		private void include(final StringBuilder sb,String className) {
		InputStream in=null;
		try {
			int dollar=className.indexOf('$');
			if(dollar!=-1) className=className.substring(0, dollar);
			className=className.replace('.', '/')+".java";
			in=getMainClass().getResourceAsStream("/"+className);
			
			
			
			if(in!=null){
				BufferedReader r=new BufferedReader(new InputStreamReader(in));
				String line;
				boolean ok=false;
				while((line=r.readLine())!=null)
					{
					if(line.contains("BEGIN"+"_DOC"))
						{
						ok=true;
						}
					else if(line.contains("END"+"_DOC") && !line.contains(".END_DOCUMENT"))
						{
						if(!ok) LOG.warn("END_"+"DOC without BEGIN"+"_DOC "+line.trim());
						ok=false;
						}
					else if(ok)
						{
						if(line.trim().startsWith("@@INCLUDE"))
							{
							int n=line.indexOf(" ");
							if(n==-1)  n=line.indexOf("\t");
							if(n==-1)  n=line.indexOf("=");
							if(n!=-1) {
								line=line.substring(n+1).trim();
								if(!line.isEmpty())
									{
									include(sb, line);
									}
								}
							}
						else
							{
							sb.append(line).append("\n");
							}
						}
					}
				r.close();
				if(ok) LOG.warn("BEGIN_"+"DOC without END");
				}
			else
				{
				LOG.debug("cannot find java code for "+className);
				}
			}
		catch(final Exception err) {
			
			}
		finally
			{
			CloserUtil.close(in);
			}
		}
		
		/** create usage as XML. No DTD. Might change in the future */
		private Document xmlUsage(final JCommander jc) throws Exception {

			final Document dom = DocumentBuilderFactory.newInstance().
						newDocumentBuilder().
						newDocument();
			final Element root = dom.createElement("program");
			Element e1;
			dom.appendChild(root);
			if(hasProgram())
				{
				root.setAttribute("name", getProgram().name());
				e1 = dom.createElement("description");
				root.appendChild(e1);
				e1.appendChild(dom.createTextNode(getProgram().description()));
				for(final String kw : getProgram().keywords())
					{
					e1 = dom.createElement("keyword");
					root.appendChild(e1);
					e1.appendChild(dom.createTextNode(kw));
					}
				
				for(final int bid : getProgram().biostars())
					{
					e1 = dom.createElement("biostar");
					root.appendChild(e1);
					e1.appendChild(dom.createTextNode("https://www.biostars.org/p/"+bid));
					}
				}
			e1 = dom.createElement("parameters");
			root.appendChild(e1);
			for(final ParameterDescription pd:jc.getParameters())
				{
				final Element e2= dom.createElement("parameter");
				e2.setAttribute("help",String.valueOf(pd.isHelp()));
				e1.appendChild(e2);
				for(final String name : pd.getParameter().names())
					{
					final Element e3 = dom.createElement("name");
					e2.appendChild(e3);
					e3.appendChild(dom.createTextNode(name));
					}
				final Element e3 = dom.createElement("description");
				e2.appendChild(e3);
				e3.appendChild(dom.createTextNode(pd.getDescription()));
				}
			
		return dom;
		}	
		
		
		
		
		private void markdown(final JCommander jc,final StringBuilder sb) {
			if(hasProgram()){
				jc.setProgramName(getProgram().name());
				} 
			else
				{
				jc.setProgramName(getMainClass().getSimpleName());
				}

			
			sb.append("# "+getMainClass().getSimpleName()+"\n\n");
			
			sb.append("![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)").append("\n\n");
			
			if(hasProgram())	{
				sb.append(getProgram().description()).
					append("\n\n");
				
				if(!StringUtils.isBlank(getProgram().deprecatedMsg()))
					{
					sb.append("\n## DEPRECATED\n\n").
						append(getProgram().deprecatedMsg()).
						append("\n");
					}
				}
				
			
			
			sb.append("\n## Usage\n\n```\n");
			jc.usage(sb);
			sb.append("\n```\n\n");

			if(hasProgram()){
				if(getProgram().keywords()!=null && getProgram().keywords().length>0) {
					sb.append("\n## Keywords\n\n");
					for(final String sk:getProgram().keywords()) sb.append(" * "+sk+"\n");
					sb.append("\n\n");
				}
				if(getProgram().biostars()!=null && getProgram().biostars().length>0) {
					sb.append("\n## See also in Biostars\n\n");
					for(final int postid:getProgram().biostars()) sb.append(" * "+hyperlink("https://www.biostars.org/p/"+postid)+"\n");
					sb.append("\n\n");
				}	
			}
			
			
				{
				final String progName=(!hasProgram()?"software":getProgram().name());
				sb.append("## Compilation\n");
				sb.append("\n");
				sb.append("### Requirements / Dependencies\n");
				sb.append("\n");
				sb.append("* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )\n");
				sb.append("\n");
				sb.append("\n");
				sb.append("### Download and Compile\n");
				sb.append("\n");
				sb.append("```bash\n");
				sb.append("$ git clone \"https://github.com/lindenb/jvarkit.git\"\n");
				sb.append("$ cd jvarkit\n");
				sb.append("$ ./gradlew "+progName+"\n");
				sb.append("```\n");
				sb.append("\n");
				sb.append("The java jar file will be installed in the `dist` directory.\n");
				sb.append("\n");
				
				if(hasProgram() && !StringUtils.isBlank(getProgram().creationDate()))
					{
					sb.append("\n## Creation Date\n\n").
						append(getProgram().creationDate()).
						append("\n\n");
					}


				sb.append("## Source code \n\n");
				
				sb.append(hyperlink("https://github.com/lindenb/jvarkit/tree/master/src/main/java/"+
					getMainClass().getName().replace('.','/')+".java")+"\n");
				sb.append("\n");
				
				final File unitTestFile = new File("src/test/java/"+getMainClass().getName().replace('.','/')+"Test.java");
				if(unitTestFile.exists()) {
					sb.append("### Unit Tests\n\n");
					sb.append(hyperlink("https://github.com/lindenb/jvarkit/tree/master/src/test/java/"+
							getMainClass().getName().replace('.','/')+"Test.java")+"\n");
					sb.append("\n");
				}
				sb.append("\n");
				sb.append("## Contribute\n");
				sb.append("\n");
				sb.append("- Issue Tracker: "+hyperlink("http://github.com/lindenb/jvarkit/issues")+"\n");
				sb.append("- Source Code: "+hyperlink("http://github.com/lindenb/jvarkit")+"\n");
				sb.append("\n");
				sb.append("## License\n");
				sb.append("\n");
				sb.append("The project is licensed under the MIT license.\n");
				sb.append("\n");
				sb.append("## Citing\n");
				sb.append("\n");
				sb.append("Should you cite **"+progName +"** ? "+hyperlink("https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md")+"\n");
				sb.append("\n");
				sb.append("The current reference is:\n");
				sb.append("\n");
				
				
				final String references[]= (!hasProgram()?new String[0]:getProgram().references());
				
				
				if(references==null || references.length==0) {
					sb.append(hyperlink("http://dx.doi.org/10.6084/m9.figshare.1425030")+"\n");
					sb.append("\n");
					sb.append("> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.\n");
					sb.append("> "+hyperlink("http://dx.doi.org/10.6084/m9.figshare.1425030")+"\n");
					sb.append("\n");
					}
				else
					{
					for(final String ref:references) {
						sb.append(" * ");
						sb.append(ref);
						sb.append("\n");
						}
					sb.append("\n");
					}
				}
			
			include(sb,getMainClass().getName());
			}
		
		private void standard(final JCommander jc,final StringBuilder sb) {

			if(hasProgram()){
				jc.setProgramName(getProgram().name());
			} else
				{
				jc.setProgramName(getMainClass().getSimpleName());
				}
			jc.usage(sb);
			}
		
		public Class<?> getMainClass()
			{
			return this._mainClass;
			}
		
		public String getCompileDate()
			{
			return JVarkitVersion.getInstance().getCompilationDate();
			}

		public String getGitHash()
			{
			return JVarkitVersion.getInstance().getGitHash();
			}
		
		public String getVersion()
			{
			return getGitHash();
			}
		
		private void updateDoc(final JCommander jc) {
			if(!hasProgram()) return;
			if(!getProgram().generate_doc()) return;
			final String doc_dir = System.getProperty("jvarkit.doc.dir", "");
			if(StringUtils.isBlank(doc_dir)) return ;
			final File thisDir = new File(doc_dir);
			if(!thisDir.exists()) return;
			if(!thisDir.isDirectory()) return;
			final File index_html =  new File(thisDir,"index.html");
			if(index_html.exists()) {
				try 
					{
					final Document dom = DocumentBuilderFactory.
							newInstance().
							newDocumentBuilder().
							parse(index_html);
					Element tr= dom.createElement("tr");
					tr.setAttribute("id", getProgram().name());
					// name
					Element td = dom.createElement("th");
					tr.appendChild(td);
					Element a= dom.createElement("a");
					a.setAttribute("href",getMainClass().getSimpleName()+ ".html");
					a.setAttribute("title",getMainClass().getSimpleName().toString());
					td.appendChild(a);
					a.appendChild(dom.createTextNode( getProgram().name()));
					//desc
					td = dom.createElement("td");
					tr.appendChild(td);
					td.appendChild(dom.createTextNode(getProgram().description()));
					//keywords
					td = dom.createElement("td");
					tr.appendChild(td);
					td.appendChild(dom.createTextNode(Arrays.asList(getProgram().keywords()).stream().collect(Collectors.joining(" "))));
					//terms (deprecated)
					td = dom.createElement("td");
					tr.appendChild(td);
					td.appendChild(dom.createTextNode(""));
					//misc
					td = dom.createElement("td");
					tr.appendChild(td);
					final DocumentFragment misc=dom.createDocumentFragment();
					if(!StringUtils.isBlank(getProgram().deprecatedMsg())) {
						Element span = dom.createElement("span");
						span.setAttribute("style", "color:orange;");
						span.setAttribute("title", "deprecated");
						span.appendChild(dom.createTextNode(getProgram().deprecatedMsg()));
						misc.appendChild(span);
						misc.appendChild(dom.createTextNode(". "));
						}
					td.appendChild(misc);
					
					
					final XPath xpath = XPathFactory.newInstance().newXPath();
					final NodeList nodeList=(NodeList)xpath.evaluate("//table[1]/tbody/tr[@id='"+getProgram().name()+"']", dom, XPathConstants.NODESET);
					for(int i=0;i< nodeList.getLength();++i) {
						Node oldTr= nodeList.item(i);
						if(tr!=null) {
							oldTr.getParentNode().replaceChild(tr, oldTr);
							tr=null;
							}
						else
							{
							oldTr.getParentNode().removeChild(oldTr);
							}
						}
					if(tr!=null) {
						Node tbody = (Node)xpath.evaluate("//table[1]/tbody",dom,XPathConstants.NODE);
						if( tbody != null ) {
							tbody.appendChild(tr);
							tbody.appendChild(dom.createTextNode("\n"));
							tr=null;
							}
						}
					if(tr!=null)
						{
						LOG.warn("Cannot insert new doc in "+index_html);
						}
					else
						{
						final Transformer transformer = TransformerFactory.newInstance().newTransformer();
						transformer.setOutputProperty(OutputKeys.INDENT, "no");
						transformer.setOutputProperty(OutputKeys.METHOD, "xml");
						transformer.setOutputProperty(OutputKeys.OMIT_XML_DECLARATION, "yes");
						transformer.transform(new DOMSource(dom), new StreamResult(index_html));
						}
					}
				catch(final Exception err)
					{
					LOG.warn(err);
					}
				
				final File mdFile = new File(thisDir,getMainClass().getSimpleName()+".md");
				final StringBuilder mdBuilder=new StringBuilder();
				this.markdown(jc, mdBuilder);
				try (final PrintWriter pw=new PrintWriter(mdFile)) {
					pw.print(mdBuilder);
					pw.flush();
					}
				catch(final IOException err) {
					LOG.warning(err);
					}
				} // END UPDATE index
			else
				{
				LOG.warn("Cannot get "+index_html);
				}
			
			}
		}
		

