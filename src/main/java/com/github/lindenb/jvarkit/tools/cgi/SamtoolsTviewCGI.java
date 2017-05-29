package com.github.lindenb.jvarkit.tools.cgi;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloserUtil;

/**
##Motivation

CGI/Web based version of samtools tview


##Screenshot
![webtview.jpg](https://github.com/lindenb/jvarkit/blob/master/doc/webtview.jpg?raw=true)


##Compilation/Installation

### Create a Preference file.
create a java preferences file with the following key/example:

```xml
<?xml version="1.0"?>
<!DOCTYPE preferences SYSTEM "http://java.sun.com/dtd/preferences.dtd">
<preferences EXTERNAL_XML_VERSION="1.0">
  <root type="user">
    <map>
      <entry key="samtools.path" value="/commun/packages/samtools-0.1.19/samtools"/>
      <entry key="default.reference.path" value="/commun/reference/human_g1k_v37.fasta"/>
    </map>
  </root>
</preferences>
```

On my server, I wrote that file in the CGI-BIN folder of apache:

```
/var/www/cgi-bin/prefs.xml
```

### Compile the program

See also [[Compilation]].

Important : the JAR/library files of picard should be visible from the CGI-bin. The path defined in the build.properties should be absolute.

```bash
$ make tview.cgi
```
It creates a jar file (tviewweb.jar) and an executable shell script (tviewweb.cgi) in the 'dist' folder. Both should be placed in the cgi-bin script of your server:
```bash
$ sudo mv dist/tviewweb* /var/www/cgi-bin/
```

If needed, edit the script tviewweb.cgi and change the JVM property `Dprefs.file.xml=` to the correct place of your xml preference file.
```
java (...) -Dprefs.file.xml=/var/www/cgi-bin/prefs.xml (...)
````

 */

@Program(name="",description="CGI/Web based version of samtools tview")
public class SamtoolsTviewCGI extends AbstractCGICallApp
	{
	private static final Logger LOG=Logger.build(SamtoolsTviewCGI.class).make();

	private SamtoolsTviewCGI()
		{
		}

	/*
	private File convertFile(File f) throws IOException
		{
		if(f==null) return f;
		if(f.isAbsolute()) return f;
		String base=getPreferences().get("base.bam.directory",null);
		if(base==null) return f;
		return new File(base,f.getPath());
		}*/

	private void goButton(
			XMLStreamWriter w,
			String bamStr,
			String label,
			int shift,
			String chrom,
			int position1
			) throws XMLStreamException,IOException
		{
		if(position1+shift<1)
			{
			position1=1;
			shift=0;
			}
		if(chrom==null || bamStr==null || position1+shift<1) return;
		
		w.writeStartElement("form");
		w.writeAttribute("method", "POST");
		w.writeAttribute("class", "formbutton");

		w.writeEmptyElement("input");
		w.writeAttribute("name", "bam");
		w.writeAttribute("type", "hidden");
		w.writeAttribute("value", bamStr);
		
		w.writeEmptyElement("input");
		w.writeAttribute("name", "pos");
		w.writeAttribute("type", "hidden");
		w.writeAttribute("value",(chrom+":"+(position1+shift)));
		
		
		
		
		w.writeStartElement("button");
		w.writeAttribute("type", "submit");
		w.writeAttribute("title", chrom+":"+(position1+shift));
		w.writeEntityRef(label);
		w.writeEndElement();//button
		
		
		w.writeEndElement();//form
		}
		
	
	private void printForm(
			XMLStreamWriter w,
			String bamStr,
			String posStr
			) throws XMLStreamException
		{
		w.writeEmptyElement("hr");
		w.writeStartElement("form");
		w.writeAttribute("method", "POST");
		
		w.writeStartElement("label");
		w.writeAttribute("for", "bam");
		w.writeCharacters("Bam File(s)");
		w.writeEndElement();
		w.writeEmptyElement("br");
		
		w.writeStartElement("textarea");
		w.writeAttribute("width", "100%");		
		w.writeAttribute("name", "bam");		
		w.writeAttribute("id", "bam");		
		
		w.writeAttribute("placeholder", "/path/to/file.bam on server");
		w.writeAttribute("required", "true");
		w.writeCharacters(bamStr==null?"":bamStr.trim());		

		w.writeEndElement();
		w.writeEmptyElement("br");
		
		
		w.writeStartElement("label");
		w.writeAttribute("for", "pos");
		w.writeCharacters("Position");
		w.writeEndElement();
		
		w.writeEmptyElement("input");
		w.writeAttribute("type", "text");
		w.writeAttribute("name", "pos");		
		w.writeAttribute("id", "pos");		
		w.writeAttribute("value", posStr==null?"":posStr.trim());		
		w.writeAttribute("placeholder", "chrom:pos");
		
		w.writeEmptyElement("br");

		w.writeEmptyElement("input");
		w.writeAttribute("type", "submit");
		w.writeEmptyElement("br");
		
		
		w.writeEndElement();//form
		
		}
	
	@Override
	protected void doCGI()
		{
		setMimeHeaderPrinted(true);
		System.out.print("Content-type: text/html;charset=utf-8\n");
		System.out.println();
		System.out.flush();
		
		
		
		XMLStreamWriter w=null;
		try
			{
			/* parse CHROM and position */
			String chrom=null;
			int position1=-1;
			String posStr=this.getString("pos");
			if(posStr!=null)
				{
				posStr=posStr.trim();
				int colon=posStr.indexOf(':');
				if(colon>0 && colon+1 < posStr.length())
					{
					chrom=posStr.substring(0,colon);
					try {
						position1=Integer.parseInt(posStr.substring(colon+1));
					} catch (Exception e)
						{
						chrom=null;
						position1=-1;
						}
					}
				
				}

			
			XMLOutputFactory xof=XMLOutputFactory.newFactory();
			w=xof.createXMLStreamWriter(System.out,"UTF-8");
			w.writeStartElement("html");
			w.writeStartElement("head");
			
			w.writeStartElement("style");
			w.writeAttribute("type", "text/css");
			w.writeCharacters(
			"label { text-align:right; margin:5px;}\n"+
			"button {font-size:200%;min-width:100px;border: 1px solid; background-image:-moz-linear-gradient( top, gray, lightgray );margin:5px;}\n"+
			"button:hover {background-image:-moz-linear-gradient( top, lightgray, gray );}\n"+
			".code {border:1px solid black;font-family:monospace;font-size:14pt;color:white;background-color:black;max-width:100%;max-height:400px;overflow:auto;padding:20px;}\n"+
			".bigtitle {text-align:center;padding:10px;text-shadow: 3px 3px 4px gray; font-size:200%;}\n"+
			"textarea{box-sizing: border-box;width: 100%;}\n"+
			".formbutton{display:inline;}\n"
			);			
			w.writeEndElement();//style
			
			
			
			if(posStr!=null)
				{
				w.writeStartElement("title");
				w.writeCharacters(posStr);
				w.writeEndElement();//title
				}

			w.writeEndElement();//head
			w.writeStartElement("body");
			w.flush();
	
			if(posStr!=null)
				{
				w.writeStartElement("h1");
				w.writeAttribute("class","bigtitle");
				w.writeCharacters(posStr);
				w.writeEndElement();//title
				}

			
			
			if(posStr!=null && posStr.length()>0 && (chrom==null || position1<=0))
				{
				w.writeCharacters("BAD POS: "+posStr);
				w.writeEmptyElement("br");
				}

				
			
			try
				{
				getPreferences();
				}
			catch(Exception err)
				{
				writeHTMLException(w, err);
				w.flush();
				}
			
			
		
			/* get samtools loc */
			File samtools=null;
			try
				{
				samtools=new File(getPreferences().get("samtools.path", "samtools"));
				if(!samtools.exists())
					{
					w.writeCharacters("Cannot find samtools: "+samtools);
					w.writeEmptyElement("br");
					}
				}
			catch(Exception err)
				{
				writeHTMLException(w, err);
				}

		
		

			
			
			String bamString=this.getString("bam");
			String bamList[]=new String[0];

			
			if(bamString!=null && samtools!=null)
				{
				w.writeStartElement("div");
				w.writeAttribute("style", "text-align:center;font-size:120%;clear: both;");
				w.writeStartElement("div");
				goButton(w,bamString,"#x21DA",-30,chrom,position1);
				goButton(w,bamString,"#x21D0",-20,chrom,position1);
				goButton(w,bamString,"#x2190",-10,chrom,position1);
				goButton(w,bamString,"#x2192",1,chrom,position1);
				goButton(w,bamString,"#x21D2",10,chrom,position1);
				goButton(w,bamString,"#x21D2",20,chrom,position1);
				goButton(w,bamString,"#x21DB",30,chrom,position1);
				w.writeEndElement();//div
				w.writeEndElement();//div
				w.flush();
				bamList=bamString.split("[\n]+");
				}
			
			for(String bamFileStr:bamList)
				{
				bamFileStr=bamFileStr.trim();
				if(samtools==null || bamFileStr.isEmpty()) continue;
				File bamFile=null;
				if(bamFileStr!=null && bamFileStr.endsWith(".bam"))
					{
					bamFile=new File(bamFileStr.trim());
					if(!(bamFile.exists() && bamFile.isFile() && bamFile.canRead()))
						{
						w.writeCharacters("BAD BAM File: "+bamFileStr);
						w.writeEmptyElement("br");
						bamFile=null;
						}
					else
						{	
						File bai=new File(bamFileStr+".bai");
						if(!bai.exists())
							{
							bai=new File(bamFileStr.substring(0,bamFileStr.length()-3)+"bai");
							}
						if(!bai.exists())
							{
							w.writeCharacters("BAM File: "+bamFileStr+" is not indexed");
							w.writeEmptyElement("br");
							bamFile=null;
							}
						}
					}
				
				if(bamFile==null) continue;
				
				w.writeEmptyElement("hr");
				w.writeStartElement("h3");
				w.writeStartElement("a");
				w.writeAttribute("href", "#");
				w.writeCharacters("file://"+bamFileStr);
				w.writeEndElement();//a
				w.writeEndElement();//h3
				
				
				
				
	
				if(bamFile!=null)
					{
					
					Process proc=null;
					InputStream os=null;
					ConsummeInputStreamThread err=null;
					try
						{
						List<String> args=new ArrayList<String>();
						
						
						
						args.add(samtools.getPath());
						args.add("tview");
						args.add("-d");args.add("T");
						if(chrom!=null && position1>0)
							{
							args.add("-p");args.add(chrom+":"+position1);
							}
						
						args.add(bamFile.getPath());
						
						String ref=getPreferences().get("default.reference.path",null);
						if(ref!=null)
							{
							args.add(ref);
							}
						w.writeComment(args.toString());
						w.flush();
						
						
						proc=Runtime.getRuntime().exec(
								args.toArray(new String[args.size()]),
								new String[]{},
								null
								);
					    err=new ConsummeInputStreamThread(proc.getErrorStream());
					    err.start();
					    int c;
					    os=proc.getInputStream();
					    
						w.writeStartElement("div");
						w.writeAttribute("style", "text-align:center;");
					    w.writeStartElement("pre");
					    w.writeAttribute("class", "code");
					    while((c=os.read())!=-1)
					    	{
					    	if(c=='\n')
					    		{
					    		w.writeEmptyElement("br");
					    		}
					    	else
					    		{
					    		w.writeCharacters(String.valueOf((char)c));
					    		}
					    	if(System.out.checkError()) break;
					    	}
					    w.writeEndElement();//pre
					    w.writeEndElement();//div
					    
					    err.join();
						}
					catch(Throwable t)
						{	
						writeHTMLException(w, t);
						w.writeCharacters(String.valueOf(t.getMessage()));
						}
					finally
						{
						if(err!=null) try { err.interrupt();} catch(Exception err2){}
						CloserUtil.close(os);
						}
					}
				}
			
			printForm(w,bamString,posStr);

			
			writeHtmlFooter(w);
			w.writeEndElement();//body
			w.writeEndElement();//html
			}
		catch(Exception err)
			{
			LOG.error(err);
			}
		finally
			{
			if(w!=null)
				{
				try {w.flush();} catch(XMLStreamException err){}
				CloserUtil.close(w);
				}
			
			}
		}
	
	
	public static void main(String[] args)
		{
		new SamtoolsTviewCGI().instanceMainWithExit(args);
		}
	}
