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
* 2015 creation copied from http://plindenbaum.blogspot.fr/2010/11/blastxmlannotations.html

*/
package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.io.FileOutputStream;
import java.util.ArrayList;
import java.util.List;

import javax.xml.stream.XMLEventFactory;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLEventWriter;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.events.XMLEvent;
import javax.xml.transform.stream.StreamSource;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloserUtil;

/**

BEGIN_DOC

## Example

```bash
$  java -jar dist-1.139/biostar165777.jar -o out__SPLIT__.xml -T Hit -N 5 ~/blastn.xml

$ ls -la ~/blastn.xml out*.xml
-rw-rw-r-- 1 lindenb lindenb 422606 nov.  14 12:47 /home/lindenb/blastn.xml
-rw-rw-r-- 1 lindenb lindenb  86319 nov.  14 16:17 out001.xml
-rw-rw-r-- 1 lindenb lindenb  83570 nov.  14 16:17 out002.xml
-rw-rw-r-- 1 lindenb lindenb  85096 nov.  14 16:17 out003.xml
-rw-rw-r-- 1 lindenb lindenb  88297 nov.  14 16:17 out004.xml
-rw-rw-r-- 1 lindenb lindenb  87123 nov.  14 16:17 out005.xml


$ grep -cF "<Hit>" ~/blastn.xml out*.xml
/home/lindenb/blastn.xml:100
out001.xml:20
out002.xml:20
out003.xml:20
out004.xml:20
out005.xml:20

```
END_DOC
 
 */

@Program(name="biostar165777",
	description="Split a XML file",
	keywords= {"xml"}
	)
public class Biostar165777 extends Launcher
	{

	private static final Logger LOG = Logger.build(Biostar165777.class).make();
	private static final String SPLIT_TOKEN="__SPLIT__";


	@Parameter(names={"-o","--output"},description="Output file. Must contains "+SPLIT_TOKEN ,required=true)
	private File outputFile = null;


	@Parameter(names={"-N","--count"},description="Number of files to be created")
	private int count = 100 ;

	@Parameter(names={"-T","--tag"},description="XML tag to be split.e.g 'Hit' in blast")
	private String tagName = null;

	
	public Biostar165777() {
		}
	
	
	@Override
	public int doWork(final List<String> args) {
		if(this.outputFile==null)
			{
			LOG.error("output file must be defined.");
			return -1;
			}
		if(!this.outputFile.getName().contains(SPLIT_TOKEN))
			{
			LOG.error("output file "+outputFile+" should contains the word "+SPLIT_TOKEN);
			return -1;
			}
		if(this.tagName==null)
			{
			LOG.error("Tag name was not defined");
			return -1;
			}
		if(this.count<1)
			{	
			LOG.error("bad count "+this.count);
			return -1;
			}
		XMLEventReader r = null;
		try {
			final List<FileOutputStream> fwriters = new ArrayList<>();
			final List<XMLEventWriter> writers = new ArrayList<>();
			final XMLOutputFactory xof=XMLOutputFactory.newFactory();
			for(int i=0;i< this.count;++i)
				{
				final File xmlout =new File(this.outputFile.getParentFile(),this.outputFile.getName().
						replaceAll(SPLIT_TOKEN, String.format("%03d", (i+1))));
				final FileOutputStream fos = new FileOutputStream(xmlout);
				fwriters.add(fos);
				writers.add(xof.createXMLEventWriter(fos,"UTF-8"));
				}
			
			//READ Whole XML file
			
			final XMLInputFactory xif = XMLInputFactory.newFactory();
			final String inputFile = oneFileOrNull(args);
			if(inputFile==null)
				{
				r = xif.createXMLEventReader(new StreamSource(stdin())); 
				}
			else
				{
				r = xif.createXMLEventReader(new StreamSource(new File(inputFile))); 
				}
			int numberOfItems = 0;
			final XMLEventFactory xef = XMLEventFactory.newFactory();
			while(r.hasNext())
				{
				
				XMLEvent evt = r.nextEvent();
				if(evt.isStartElement() && evt.asStartElement().getName().getLocalPart().equals(this.tagName))
					{
					++numberOfItems;
					XMLEventWriter w = writers.get(numberOfItems%this.count);
					int depth=1;
					w.add(evt);
					while(r.hasNext())
						{
						evt=r.nextEvent();
						if(evt.isStartElement())
							{
							++depth;
							}
						else if(evt.isEndElement())
							{
							--depth;
							}
						w.add(evt);
						if(depth==0) break;
						}
					while(r.hasNext())
						{	
						evt=r.peek();
						if(!evt.isCharacters()) break;
						if(!evt.asCharacters().getData().trim().isEmpty())
							{
							break;
							}		

						r.nextEvent();
						}
					continue;
					}
				if(evt.isStartDocument())
					{
					evt = xef.createStartDocument("UTF-8", "1.0");
					}
				for(final XMLEventWriter w:writers) w.add(evt);
				}
			r.close();
			for(final XMLEventWriter w: writers)
				{
				w.flush();
				CloserUtil.close(w);
				}
			for(final FileOutputStream w:fwriters)
				{	
				CloserUtil.close(w);
				}
			return RETURN_OK;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		}
	
	public static void main(final String[] args) {
		new Biostar165777().instanceMainWithExit(args);
		}
	}
