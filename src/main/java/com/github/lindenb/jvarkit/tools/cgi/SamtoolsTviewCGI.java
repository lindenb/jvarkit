package com.github.lindenb.jvarkit.tools.cgi;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.URLEncoder;
import java.util.ArrayList;
import java.util.List;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import net.sf.samtools.util.CloserUtil;

public class SamtoolsTviewCGI extends AbstractCGICallApp
	{
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
			int shift,
			String chrom,
			int position1
			) throws XMLStreamException,IOException
		{
		if(chrom==null || bamStr==null || position1+shift<1) return;
		w.writeStartElement("a");
		w.writeAttribute("title", (position1<0?"+":"")+shift);
		w.writeAttribute("href",
				"?bam="+URLEncoder.encode(bamStr,"UTF-8")+
				"&pos="+URLEncoder.encode(chrom+":"+(position1+shift),"UTF-8")
				);
		w.writeCharacters((position1<0?"+":"")+shift);
		w.writeEndElement();
		}
		
	
	private void printForm(
			XMLStreamWriter w,
			String bamStr,
			String posStr
			) throws XMLStreamException
		{
		w.writeStartElement("form");
		w.writeAttribute("method", "GET");
		
		w.writeStartElement("label");
		w.writeAttribute("for", "bam");
		w.writeCharacters("bam File");
		w.writeEndElement();
		
		w.writeEmptyElement("input");
		w.writeAttribute("type", "text");
		w.writeAttribute("name", "bam");		
		w.writeAttribute("id", "bam");		
		w.writeAttribute("value", bamStr==null?"":bamStr.trim());		
		w.writeAttribute("placeholder", "/path/to/file.bam on server");
		w.writeAttribute("required", "true");
		
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
		w.writeAttribute("required", "true");
		
		w.writeEmptyElement("br");

		w.writeEmptyElement("input");
		w.writeAttribute("type", "submit");
		w.writeEmptyElement("br");
		
		
		w.writeEndElement();//form
		w.writeEmptyElement("hr");
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
			XMLOutputFactory xof=XMLOutputFactory.newFactory();
			w=xof.createXMLStreamWriter(System.out,"UTF-8");
			w.writeStartElement("html");
			w.writeStartElement("body");
			String bamStr=this.getString("bam");
			File bamFile=null;
			if(bamStr!=null && bamStr.trim().endsWith(".bam"))
				{
				bamFile=new File(bamStr.trim());
				if(!(bamFile.exists() && bamFile.isFile() && bamFile.canRead()))
					{
					w.writeCharacters("BAD BAM File: "+bamStr);
					w.writeEmptyElement("br");
					bamFile=null;
					}
				else
					{	
					File bai=new File(bamStr+".bai");
					if(!bai.exists())
						{
						bai=new File(bamStr.substring(0,bamStr.length()-3)+"bai");
						}
					if(!bai.exists())
						{
						w.writeCharacters("BAM File: "+bamStr+" is not indexed");
						w.writeEmptyElement("br");
						bamFile=null;
						}
					}
				}
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
				
				if(chrom==null || position1<=0)
					{
					w.writeCharacters("BAD POS: "+posStr);
					w.writeEmptyElement("br");
					}
				}
			
			
			
			printForm(w,bamStr,posStr);
			
			if(bamFile!=null)
				{
				w.writeStartElement("div");
				goButton(w,bamStr,-30,chrom,position1);
				goButton(w,bamStr,-20,chrom,position1);
				goButton(w,bamStr,-10,chrom,position1);
				goButton(w,bamStr,-1,chrom,position1);
				goButton(w,bamStr,1,chrom,position1);
				goButton(w,bamStr,10,chrom,position1);
				goButton(w,bamStr,20,chrom,position1);
				goButton(w,bamStr,30,chrom,position1);
				w.writeEndElement();//div
				
				Process proc=null;
				InputStream os=null;
				ConsummeInputStreamThread err=null;
				try
					{
					List<String> args=new ArrayList<String>();
					
					
					
					args.add(getPreferences().get("samtools.path", "samtools"));
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
					
					
					proc=Runtime.getRuntime().exec(
							args.toArray(new String[args.size()]),
							new String[]{},
							null
							);
				    err=new ConsummeInputStreamThread(proc.getErrorStream());
				    err.start();
				    int c;
				    os=proc.getInputStream();
				    w.writeStartElement("pre");
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
				    err.join();
					}
				catch(Throwable t)
					{	
					
					}
				finally
					{
					if(err!=null) try { err.interrupt();} catch(Exception err2){}
					CloserUtil.close(os);
					}
				}
			
			
			w.writeEndElement();
			w.writeEndElement();
			}
		catch(Exception err)
			{
			error(err);
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
