/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.nextflow.nf2xml;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.StringReader;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.testng.annotations.Parameters;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.nextflow.parser.NextflowParser;
import com.github.lindenb.jvarkit.nextflow.parser.ParseException;


public class NextflowToXml extends Launcher {
	@Parameter(names= {"-o","--out"},description = OPT_OUPUT_FILE_OR_STDOUT )
	private File xmlFileOut=null;
	
	public NextflowToXml() {
	}
	
	private void parseNextflowScript(Document dom,Element root,Path p) throws IOException,ParseException {
		final Document dom2= NextflowParser.parse(p);
		final Element root2 = dom2.getDocumentElement();
		for(Node c=root2.getFirstChild();c!=null;c=c.getNextSibling()) {
			if(c.getNodeType()!=Node.ELEMENT_NODE) continue;
			root.appendChild(dom.importNode(c, true));
			}
		}
	
	@Override
	public int doWork(final List<String> args) {
		try {
			final String input = oneFileOrNull(args);
			final DocumentBuilder db = DocumentBuilderFactory.newInstance().newDocumentBuilder();
			final Document dom = db.newDocument();
			final Element process=dom.createElement("process");
			dom.appendChild(process);
			final Element description = dom.createElement("description");
			String about="";
			process.appendChild(description);
			
			final Element inputE = dom.createElement("input");
			process.appendChild(inputE);

			final Element outputE = dom.createElement("output");
			process.appendChild(outputE);

			final Element tools = dom.createElement("tools");
			process.appendChild(tools);

			final Element creation = dom.createElement("creation");
			process.appendChild(creation);
			final Element modification = dom.createElement("modification");
			process.appendChild(modification);

			
			try(BufferedReader br=StringUtils.isBlank(input)? IOUtils.openStreamForBufferedReader(stdin()):IOUtils.openURIForBufferedReading(input)) {
				List<String> lines =  br.lines().
						map(L->L.replaceAll("[ \t]+"," ").trim()).
						filter(L->!StringUtils.isBlank(L)).
						collect(Collectors.toList());
				String procline=lines.stream().filter(L->L.startsWith("process ")).findFirst().orElseThrow();
				procline=procline.substring(8);
				about= "description missing for "+procline;
				process.setAttribute("name", procline);
				
				for(int i=0;i< lines.size();i++) {
					if(lines.get(i).equals("input:")) {
						for(int j=i+1;j< lines.size();j++) {
							if(lines.get(j).matches("[a-z]+:")) {
								break;	
								}
							try(StringReader r=new StringReader(lines.get(j))) {
								final Element item=new NextflowParser(r,dom).input_block_item();
								inputE.appendChild(item);
								}
							}
						break;
						}
					}
				
				
				for(int i=0;i< lines.size();i++) {
					if(lines.get(i).equals("output:")) {
						for(int j=i+1;j< lines.size();j++) {
							if(lines.get(j).matches("[a-z]+:")) {
								break;	
								}
							try(StringReader r=new StringReader(lines.get(j))) {
								final Element item=new NextflowParser(r,dom).output_block_item();
								
								inputE.appendChild(item);
								}
							}
						break;
						}
					}
				
				}
			
			
			description.appendChild(dom.createTextNode(about));
			
			final StreamResult out;
			if(this.xmlFileOut==null) {
				out = new StreamResult(stdout());
				}
			else
				{
				out = new StreamResult(this.xmlFileOut);
				}
			TransformerFactory.newInstance().newTransformer().transform(new DOMSource(dom),out);
			return 0;
			}
		catch(Throwable err) {
			err.printStackTrace();
			return -1;
			}
		}
	
	public static void main(final String[] args) {
		new NextflowToXml().instanceMainWithExit(args);
	}
}
