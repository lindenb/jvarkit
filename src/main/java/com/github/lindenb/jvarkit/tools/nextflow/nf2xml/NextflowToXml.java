package com.github.lindenb.jvarkit.tools.nextflow.nf2xml;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.nextflow.parser.NextflowParser;
import com.github.lindenb.jvarkit.nextflow.parser.ParseException;


public class NextflowToXml extends Launcher {

	public NextflowToXml() {
	}
	
	private void parseNextflowScript(Document dom,Element root,Path p) throws IOException,ParseException {
		Document dom2= NextflowParser.parse(p);
		Element root2=dom2.getDocumentElement();
		for(Node c=root2.getFirstChild();c!=null;c=c.getNextSibling()) {
			if(c.getNodeType()!=Node.ELEMENT_NODE) continue;
			root.appendChild(dom.importNode(c, true));
			}
		}
	
	@Override
	public int doWork(List<String> args) {
		try {
			final DocumentBuilder db = DocumentBuilderFactory.newInstance().newDocumentBuilder();
			final Document dom = db.newDocument();
			final Element root=dom.createElement("nextflow");
			dom.appendChild(root);
			if(args.isEmpty() || (args.size()==1 && args.get(0).equals("-"))) {
				try(BufferedReader br=IOUtils.openStreamForBufferedReader(stdin())) {
					for(;;) {
						String line = br.readLine();
						if(line==null) break;
						if(StringUtils.isBlank(line) || line.startsWith("#")) continue;
						Path script= Paths.get(line);
						parseNextflowScript(dom,root,script);
						}
					}
				}
			else
				{
				for(String line: args) {
					Path script= Paths.get(line);
					parseNextflowScript(dom,root,script);
					}
				}
			TransformerFactory.newInstance().newTransformer().transform(new DOMSource(dom),new StreamResult(stdout()));
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
