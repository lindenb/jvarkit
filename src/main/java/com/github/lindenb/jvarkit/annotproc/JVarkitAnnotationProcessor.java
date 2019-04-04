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
package com.github.lindenb.jvarkit.annotproc;

import java.io.File;
import java.io.FileReader;
import java.io.Reader;
import java.io.Writer;
import java.util.Arrays;
import java.util.Set;
import java.util.stream.Collectors;

import javax.annotation.processing.AbstractProcessor;
import javax.annotation.processing.Filer;
import javax.annotation.processing.ProcessingEnvironment;
import javax.annotation.processing.RoundEnvironment;
import javax.annotation.processing.SupportedSourceVersion;
import javax.lang.model.SourceVersion;
import javax.lang.model.element.ElementKind;
import javax.lang.model.element.TypeElement;
import javax.lang.model.type.DeclaredType;
import javax.lang.model.type.TypeKind;
import javax.lang.model.type.TypeMirror;
import javax.tools.FileObject;
import javax.tools.StandardLocation;
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
import org.w3c.dom.NodeList;
import org.w3c.dom.Node;

import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

@SupportedSourceVersion(SourceVersion.RELEASE_11)
@Deprecated
public class JVarkitAnnotationProcessor extends AbstractProcessor{
	private static final Logger LOG = Logger.build(JVarkitAnnotationProcessor.class).make();
	
	@Override
	public synchronized void init(final ProcessingEnvironment processingEnv) {
	   super.init(processingEnv);
	}
	
	
	
	@Override
	public Set<String> getSupportedAnnotationTypes() {
		return Arrays.asList(
				IncludeSourceInJar.class,
				Program.class
				).stream().
			map(C->C.getName()).
			collect(Collectors.toSet())
			;
		}
	
	private void copySource(final javax.lang.model.element.Element element)
		{
		
		if(element==null || element.getKind()!=ElementKind.CLASS) return;
		final String thisDir = System.getProperty("jvarkit.this.dir");	
		if(thisDir==null || thisDir.isEmpty()) {
			System.err.println("jvarkit.this.dir is not defined");
			return ;
			}
		Writer writer = null;
		Reader reader = null;
		try 
			{
			do
				{
				String className= element.toString();
				if(className==null ||  className.isEmpty() || className.equals("java.lang.Object") ) return;
				int dollar  = className.indexOf('$');
				if(dollar!=-1) className = className.substring(0,dollar);
				File javaFile = new File(thisDir+"/src/main/java/"+className.replace('.','/') +".java");
				if(!javaFile.exists()) break;
				final Filer filer = super.processingEnv.getFiler();					
				final String packageName;
				final String fileName;
				int dot= className.lastIndexOf('.');
				if(dot==-1)
					{
					packageName = "";
					fileName = className;
					}
				else
					{
					packageName = className.substring(0,dot);
					fileName = className.substring(dot+1);
					}
				final FileObject fo=filer.createResource(StandardLocation.CLASS_OUTPUT,
						packageName,
						fileName+".java");
				
				if( new File(fo.getName()).exists()) break;
				
				System.err.println("## Copying "+ javaFile+ " -> "+ fo.getName());
				reader = new FileReader(javaFile);
				writer=fo.openWriter();
				final char  buffer[]=new char[1024];
				int nRead;
				while((nRead=reader.read(buffer))!=-1)
					{
					writer.write(buffer,0,nRead);
					}
				writer.close();writer=null;
				reader.close();reader=null;
				} while(false);
			}
		catch(final Exception err)
			{
			System.err.println(err.getMessage());
			}
		finally
			{
			if(writer!=null) try {writer.close();}catch(Exception err) {}
			if(reader!=null) try {reader.close();}catch(Exception err) {}
			}
		final TypeMirror superclass = ((TypeElement) element).getSuperclass();
		if(superclass==null || superclass.getKind() == TypeKind.NONE) {
			return ;
			}
		final DeclaredType kind = (DeclaredType) superclass;
		copySource(kind.asElement());
		}
	
	@Override
	public boolean process(final Set<? extends TypeElement> annotations,
		final RoundEnvironment roundEnv
		) {
		final String mainClass = System.getProperty("jvarkit.main.class");
		final String thisDir = System.getProperty("jvarkit.this.dir");		
				
		
		roundEnv.getElementsAnnotatedWith(IncludeSourceInJar.class).stream().
			filter(E->E.getKind()==ElementKind.CLASS).
			filter(E-> E.getAnnotation(IncludeSourceInJar.class) !=null).
			forEach(E->{					
				if( thisDir==null || thisDir.isEmpty()) return;
				copySource(E);
				});
		
		/* find if roundEnv contains main class annotated with 'Program' annotation
		 * if true: generate a file that will tell Make to compile the markdown
		 * documentation
		 */
		roundEnv.getElementsAnnotatedWith(Program.class).stream().
			filter(E->E.getKind()==ElementKind.CLASS).
			filter(E->{final Program prog=E.getAnnotation(Program.class); return prog!=null && prog.generate_doc();}).
			forEach(E->{
				copySource(E);
				final String className=E.toString();
				if(mainClass==null ) return;
				if(!mainClass.equals(className)) return;
				final Program prog=E.getAnnotation(Program.class);
				if(prog==null || !prog.generate_doc()) return;

				try 
					{
					final Filer filer = super.processingEnv.getFiler();
					FileObject fo=filer.createResource(StandardLocation.CLASS_OUTPUT, "", "markdown.flag" );
					fo.openWriter().append(String.valueOf(mainClass)).close();
					}
				catch(final Exception err)
					{
					LOG.warn(err);
					}
				if(thisDir!=null) {
					final File index_html =  new File(thisDir,"docs/index.html");
					if(index_html.exists()) {
						try 
							{
							final Document dom = DocumentBuilderFactory.
									newInstance().
									newDocumentBuilder().
									parse(index_html);
							Element tr= dom.createElement("tr");
							tr.setAttribute("id", prog.name());
							// name
							Element td = dom.createElement("th");
							tr.appendChild(td);
							Element a= dom.createElement("a");
							a.setAttribute("href",E.getSimpleName()+ ".html");
							a.setAttribute("title",E.getSimpleName().toString());
							td.appendChild(a);
							a.appendChild(dom.createTextNode(prog.name()));
							//desc
							td = dom.createElement("td");
							tr.appendChild(td);
							td.appendChild(dom.createTextNode(prog.description()));
							//keywords
							td = dom.createElement("td");
							tr.appendChild(td);
							td.appendChild(dom.createTextNode(Arrays.asList(prog.keywords()).stream().collect(Collectors.joining(" "))));
							//terms (deprecated)
							td = dom.createElement("td");
							tr.appendChild(td);
							td.appendChild(dom.createTextNode(""));
							//misc
							td = dom.createElement("td");
							tr.appendChild(td);
							final DocumentFragment misc=dom.createDocumentFragment();
							if( prog.deprecatedMsg()!=null && !prog.deprecatedMsg().isEmpty()) {
								Element span = dom.createElement("span");
								span.setAttribute("style", "color:orange;");
								span.setAttribute("title", "deprecated");
								span.appendChild(dom.createTextNode( prog.deprecatedMsg()));
								misc.appendChild(span);
								misc.appendChild(dom.createTextNode(". "));
								}
							td.appendChild(misc);
							
							
							final XPath xpath = XPathFactory.newInstance().newXPath();
							final NodeList nodeList=(NodeList)xpath.evaluate("//table[1]/tbody/tr[@id='"+prog.name()+"']", dom, XPathConstants.NODESET);
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
						}
					else
						{
						LOG.warn("Cannot get "+index_html);
						}
					}
				});
		return true;
	}
	
	

}
