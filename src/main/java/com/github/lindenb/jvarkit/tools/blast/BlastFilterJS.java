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
package com.github.lindenb.jvarkit.tools.blast;

import gov.nih.nlm.ncbi.blast.Hit;
import gov.nih.nlm.ncbi.blast.ObjectFactory;
import java.io.File;
import java.io.PrintWriter;
import java.util.List;

import javax.script.CompiledScript;
import javax.script.SimpleBindings;
import javax.xml.stream.XMLEventFactory;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLEventWriter;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;
import javax.xml.transform.stream.StreamSource;
import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBElement;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

/**

BEGIN_DOC




### Examples

Filter Hit having <Hit_len> <1500


```

$ java -jar dist/blastfilterjs.jar blastn.xml  -e 'parseInt(hit.getHitLen())<1500' 2> /dev/null |\
  xmllint --format - | grep "Hit_len"
  
   <Hit_len>1492</Hit_len>
   <Hit_len>1488</Hit_len>
   <Hit_len>1477</Hit_len>
   <Hit_len>1452</Hit_len>
   <Hit_len>1430</Hit_len>
   <Hit_len>1064</Hit_len>
   <Hit_len>1283</Hit_len>
   <Hit_len>1052</Hit_len>
   <Hit_len>1272</Hit_len>
   <Hit_len>693</Hit_len>
     

```


keep hsp having 100 > Hsp_align-len <= 200 


```

$ cat filter.js

// keep hsp having 100>= Hsp_align-len <= 200 
function rmhsps()
	{
	var hsps = hit.getHitHsps().getHsp();
	var i=0;
	while(i< hsps.size())
		{
		var hsp = hsps.get(i);
		var hsplen = parseInt(hsp.getHspAlignLen());
		
		if( hsplen < 100 || hsplen > 300 )
			{
			hsps.remove(i);
			}
		else
			{
			i++;
			}
		}
	return true;
	}
rmhsps();

```






```

$ java -jar dist/blastfilterjs.jar -f filter.js blastn.xml 2> /dev/null |\
	xmllint --format - | grep -F 'Hsp_align-len'

	 <Hsp_align-len>289</Hsp_align-len>
	 <Hsp_align-len>291</Hsp_align-len>
	 <Hsp_align-len>197</Hsp_align-len>
	 <Hsp_align-len>227</Hsp_align-len>

```





### See also


 *  






END_DOC
*/


/**
 * BlastFilterJS
 *
 */
@Program(name="blastfilterjs",
	description="Filters a BlastOutput with a javascript expression. The script injects each <Hit> as the variable 'blasthit'. The user script should return 'true' to keep the hit.",
	keywords={"blast","js","javascript","filter"}
	)
public class BlastFilterJS
	extends Launcher
	{
	private static final Logger LOG = Logger.build(BlastFilterJS.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-e","--expression"},description=" (js expression). Optional.")
	private String scriptExpr=null;
	@Parameter(names={"-f","--script"},description=" (js file). Optional.")
	private File scriptFile=null;


	@SuppressWarnings("unused")
	private static ObjectFactory _fool_javac=null;
	
	public BlastFilterJS()
		{
		
		}
	@Override
	public int doWork(List<String> args) {
	final CompiledScript compiledScript;
		Unmarshaller unmarshaller;
		Marshaller marshaller;
		try
			{
			compiledScript = super.compileJavascript(scriptExpr,scriptFile);
			
			
			JAXBContext jc = JAXBContext.newInstance("gov.nih.nlm.ncbi.blast");
			
			unmarshaller =jc.createUnmarshaller();
			marshaller =jc.createMarshaller();

			
			XMLInputFactory xmlInputFactory=XMLInputFactory.newFactory();
			xmlInputFactory.setProperty(XMLInputFactory.IS_NAMESPACE_AWARE, Boolean.TRUE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_COALESCING, Boolean.TRUE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_REPLACING_ENTITY_REFERENCES, Boolean.TRUE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_SUPPORTING_EXTERNAL_ENTITIES, Boolean.FALSE);
			marshaller.setProperty(Marshaller.JAXB_FRAGMENT, true);
			marshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, true);
			
			PrintWriter pw=openFileOrStdoutAsPrintWriter(outputFile);
			XMLOutputFactory xof=XMLOutputFactory.newFactory();
			xof.setProperty(XMLOutputFactory.IS_REPAIRING_NAMESPACES, Boolean.FALSE);
			XMLEventWriter w=xof.createXMLEventWriter(pw);
			
			
			final String inputName=oneFileOrNull(args);
			StreamSource src=null;
			if(inputName==null)
				{
				LOG.info("Reading stdin");
				src=new StreamSource(stdin());
				}
			else
				{
				LOG.info("Reading file "+inputName);
				src=new StreamSource(new File(inputName));
				}
			
			XMLEventReader r=xmlInputFactory.createXMLEventReader(src);
			
			XMLEventFactory eventFactory=XMLEventFactory.newFactory();
			
			SimpleBindings bindings=new SimpleBindings();
			while(r.hasNext())
				{
				XMLEvent evt=r.peek();
				switch(evt.getEventType())
					{
					case XMLEvent.START_ELEMENT:
						{
						StartElement sE= evt.asStartElement();
						Hit hit=null;
						JAXBElement<Hit> jaxbElement=null;
						if(sE.getName().getLocalPart().equals("Hit") )
							{
							jaxbElement= unmarshaller.unmarshal(r,Hit.class);
							hit= jaxbElement.getValue();
							}
						else
							{
							w.add(r.nextEvent());
							break;
							}
						
						if(hit!=null)
							{
							
							bindings.put("hit", hit);
							
							
							boolean accept = super.evalJavaScriptBoolean(compiledScript, bindings);
							
							if(accept)
								{
								marshaller.marshal(
										jaxbElement,
										w);
								w.add(eventFactory.createCharacters("\n"));
								}
							}
						
						break;
						}
					case XMLEvent.SPACE:break;
					default:
						{
						w.add(r.nextEvent());
						break;
						}
					}
				r.close();
				}
			w.flush();
			w.close();
			pw.flush();
			pw.close();
			return RETURN_OK;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			
			}
		}
	
	public static void main(String[] args) {
		new BlastFilterJS().instanceMainWithExit(args);
	}
}
