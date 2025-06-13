package com.github.lindenb.jvarkit.tools.nfcore;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.Document;
import org.w3c.dom.Element;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;

import htsjdk.samtools.util.IOUtil;

public class NFProcToMeta extends Launcher {
	private static final Logger LOG = Logger.of(NFProcToMeta.class);
	private void parse(String inputName,BufferedReader br) throws IOException, ParserConfigurationException{
		DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
		dbf.setNamespaceAware(false);
		dbf.setValidating(false);
		Document dom = dbf.newDocumentBuilder().newDocument();
		Element root= dom.createElement("nextflow");
		dom.appendChild(root);
		
		for(;;) {
			String line = br.readLine();
			if(line==null) break;
			Element procE = NextflowProcessParser.processDeclaration(dom, line);
			if(procE!=null) {
				dom.appendChild(procE);

				
				Element inputE= dom.createElement("inputs");
				procE.appendChild(inputE);
				
				Element outputE= dom.createElement("outputs");
				procE.appendChild(outputE);
				
				for(;;) {
					line = br.readLine();
					if(line==null) break;
					if(line.trim().equals("inputs:")) {
						
						}
					else if(line.trim().equals("outputs:")) {
						
						}
					}
				
				}
		}
		
		
	
		//return dom;
		}
	
	@Override
	public int doWork(List<String> args) {
		try {
			final String input = oneFileOrNull(args);
			if(!StringUtils.isBlank(input)) {
				final Path path  = Paths.get(input);
				IOUtil.assertFileIsReadable(path);
				if(!path.getFileName().toString().endsWith(".nf")) {
					LOG.error(input+": should end with '.nf'");
					return -1;
					}
				if(!path.getFileName().toString().equals("main.nf")) {
					LOG.warn("path should be main.nf not "+path.getFileName());
					}
				try(BufferedReader br =  IOUtils.openPathForBufferedReading(path)) {
					parse(input,br);

					}
				}
			else
				{
				try(BufferedReader br =  IOUtils.openStreamForBufferedReader(stdin())) {
					parse("<stdin>",br);
					}
				}
			return -1;
			}
		catch(Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	public static void main(String[] args) {
		new NFProcToMeta().instanceMain(args); 
	}

}
