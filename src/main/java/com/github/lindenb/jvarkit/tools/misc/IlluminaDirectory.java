package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.PrintStream;
import java.math.BigInteger;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import net.sf.picard.PicardException;


import com.github.lindenb.jvarkit.util.illumina.FastQName;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;

public class IlluminaDirectory extends AbstractCommandLineProgram
	{
	private int ID_GENERATOR=0;

	
    
    private MessageDigest md5;

    
    private  String md5(String in)
    	{
    	if(md5==null)
	    	{
	    	  try {
	              md5 = MessageDigest.getInstance("MD5");
	          } catch (NoSuchAlgorithmException e) {
	              throw new PicardException("MD5 algorithm not found", e);
	          }
	    	}
    	 md5.reset();
         md5.update(in.getBytes());
         String s = new BigInteger(1, md5.digest()).toString(16);
         if (s.length() != 32) {
             final String zeros = "00000000000000000000000000000000";
             s = zeros.substring(0, 32 - s.length()) + s;
         }
         return s;
    	}
    
    private class Folder
    	{
    	File folder;
    	SortedMap<String, Sample> sampleMap=new TreeMap<String, IlluminaDirectory.Sample>();
    	List<Pair> undtermined=new ArrayList<Pair>();
    	void scan(File f)
    		{
    		if(f==null) return;
    		if(!f.canRead()) return;
    		IlluminaDirectory.this.info("Scanning "+f);
    		if(f.isDirectory())
    			{
    			File children[]=f.listFiles();
    			if(children==null) return;
    			for(File c:f.listFiles())
    				{
    				scan(c);
    				}
    			}
    		else
    			{
    			if(f.getName().endsWith(".fastq.gz"))
    				{
    				FastQName fq=FastQName.parse(f);
    				if(!fq.isValid())
    					{
    					warning("invalid name:"+fq);
    					return;
    					}
    				if(fq.isUndetermined())
    					{
    					for(int i=0;i< undtermined.size();++i)
	    	    			{
	    	    			Pair p=undtermined.get(i);
	    	    			if(p.complement(fq)) return;
	    	    			}
    					undtermined.add(new Pair(fq));
    					}
    				else
    					{
    					Sample sample=this.sampleMap.get(fq.getSample());
    					if(sample==null)
    						{
    						sample=new Sample();
    						sample.name=fq.getSample();
    						this.sampleMap.put(sample.name,sample);
    						}
    					sample.add(fq);
    					}
    				}
    			}
    		}
    	
    	void json(PrintStream out)
    		{
    		boolean first=true;
    		out.print("{\"directory\":\""+this.folder.getPath()+"\",\"samples\":[");
    		for(Sample S:this.sampleMap.values())
				{
    			if(!first) out.print(',');
    			first=false;
				S.json(out);
				}
    		out.print("],\"undetermined\":[");
    		first=true;
    		for(Pair p:undtermined)
				{
    			if(!first) out.print(',');
    			first=false;
				p.json(out);
				}
    		out.print("]}");
    		}
    	
    	void write(XMLStreamWriter w) throws XMLStreamException
    		{
    		w.writeStartElement("directory");
    		w.writeAttribute("path", this.folder.getPath());
    		w.writeStartElement("samples");
    		for(Sample S:this.sampleMap.values())
    			{
    			S.write(w);
    			}
    		w.writeEndElement();
    		w.writeStartElement("undetermined");
    		for(Pair p:undtermined)
    			{
    			p.write(w);
    			}
    		w.writeEndElement();
    		w.writeEndElement();
    		}
    	
    	}
    
    private class Pair
    	{
    	int id;
    	FastQName forward;
    	FastQName reverse;
    	
    	Pair(FastQName fq)
    		{
    		id=++ID_GENERATOR;
    		switch(fq.getSide())
    			{
    			case Forward:forward=fq; break;
    			case Reverse:reverse=fq; break;
    			default:throw new PicardException("bad side "+fq);
    			}
    		}
    	
    	boolean complement(FastQName other)
    		{
    		if(forward!=null && reverse!=null) return false;
    		if(forward!=null && forward.isComplementOf(other))
    			{
    			reverse=other;
    			return true;
    			}
    		else if(reverse!=null && reverse.isComplementOf(other))
    			{
    			forward=other;
    			return true;
    			}
    		return false;
    		}
    	
    	void json(PrintStream out)
    		{
    		if(forward!=null && reverse!=null)
				{
	    		
	    		out.print("{");
	    		out.print("\"id\":\"p"+this.id+"\",");
	    		out.print("\"md5pair\":\""+md5(forward.getFile().getPath()+reverse.getFile().getPath())+"\",");
	    		out.print("\"lane\":"+forward.getLane()+",");
    			if(forward.getSeqIndex()!=null)
    				{
    				out.print("\"index\":\""+forward.getSeqIndex()+"\",");
    				}
    			else
    				{
    				out.print("\"index\":null,");
    				}
	    		out.print("\"split\":"+forward.getSplit()+",");
	    		out.print("\"forward\":{");
	    		out.print("\"md5filename\":\""+md5(forward.getFile().getPath())+"\",");
	    		out.print("\"path\":\""+forward.getFile().getPath()+"\",");
	    		out.print("\"side\":"+forward.getSide().ordinal()+",");
	    		out.print("\"file-size\":"+forward.getFile().length());
	    		out.print("},\"reverse\":{");
	    		out.print("\"md5filename\":\""+md5(reverse.getFile().getPath())+"\",");
	    		out.print("\"path\":\""+reverse.getFile().getPath()+"\",");
	    		out.print("\"side\":"+reverse.getSide().ordinal()+",");
	    		out.print("\"file-size\":"+reverse.getFile().length());
	    		out.print("}}");
				}
    		else
    			{
    			FastQName F=(forward==null?reverse:forward);
    			out.print("{");
    			out.print("\"id\":\"p"+this.id+"\",");
    			out.print("\"md5filename\":\""+md5(F.getFile().getPath())+"\",");
	    		out.print("\"lane\":"+F.getLane()+",");
    			if(forward.getSeqIndex()!=null)
    				{
    				out.print("\"index\":\""+F.getSeqIndex()+"\",");
    				}
    			else
    				{
    				out.print("\"index\":null,");
    				}
	    		out.print("\"split\":"+F.getSplit()+",");
	    		out.print("\"path\":\""+F.getFile().getPath()+"\",");
	    		out.print("\"side\":"+F.getSide()+",");
    			out.print("}");
    			}
			}
    	
    	void write(XMLStreamWriter w) throws XMLStreamException
    		{
    		if(forward!=null && reverse!=null)
    			{
    			w.writeStartElement("pair");
    			w.writeAttribute("id","p"+this.id);
    			w.writeAttribute("md5",md5(forward.getFile().getPath()+reverse.getFile().getPath()));
    			w.writeAttribute("lane", String.valueOf(forward.getLane()));
    			if(forward.getSeqIndex()!=null) w.writeAttribute("index", String.valueOf(forward.getSeqIndex()));
    			w.writeAttribute("split", String.valueOf(forward.getSplit()));
    			
    			w.writeStartElement("fastq");
    			w.writeAttribute("md5filename",md5(forward.getFile().getPath()));
    			w.writeAttribute("side", String.valueOf(forward.getSide().ordinal()));
    			w.writeAttribute("path", forward.getFile().getPath());
    			w.writeAttribute("file-size",String.valueOf( forward.getFile().length()));
    			w.writeEndElement();
    			
    			
    			w.writeStartElement("fastq");
    			w.writeAttribute("md5filename",md5(reverse.getFile().getPath()));
    			w.writeAttribute("side", String.valueOf(reverse.getSide().ordinal()));
    			w.writeAttribute("path", reverse.getFile().getPath());
    			w.writeAttribute("file-size",String.valueOf( reverse.getFile().length()));
    			w.writeEndElement();
   			
    			w.writeEndElement();
    			}
    		else
    			{
    			FastQName F=(forward==null?reverse:forward);
    			
    			w.writeStartElement("single");
    			w.writeAttribute("id","p"+this.id);
    			w.writeAttribute("md5filename",md5(F.getFile().getPath()));
    			w.writeAttribute("lane", String.valueOf(F.getLane()));
    			if(forward.getSeqIndex()!=null) w.writeAttribute("index", String.valueOf(F.getSeqIndex()));
    			w.writeAttribute("split", String.valueOf(F.getSplit()));
    			
    			w.writeStartElement("fastq");
    			w.writeAttribute("side", String.valueOf(F.getSide().ordinal()));
    			w.writeAttribute("path", F.getFile().getPath());
    			w.writeAttribute("file-size",String.valueOf( F.getFile().length()));
    			w.writeEndElement();
   			
    			w.writeEndElement();
    			}
    		}
    	}
    
    private class Sample
		{
		String name;
    	List<Pair> pairs=new ArrayList<Pair>();
    	
    	private void add(FastQName fq)
    		{
    		for(int i=0;i< pairs.size();++i)
    			{
    			Pair p=pairs.get(i);
    			if(p.complement(fq)) return;
    			}
    		pairs.add(new Pair(fq));
    		}
    	
    	void write(XMLStreamWriter w) throws XMLStreamException
			{
			w.writeStartElement("sample");
			w.writeAttribute("name",this.name);
			for(Pair p:this.pairs)
				{
				p.write(w);
				}
				
			w.writeEndElement();
			}
    	
    	void json(PrintStream out)
    		{
    		out.print("{\"sample\":\""+ this.name +"\",\"files\":[");
    		for(int i=0;i< pairs.size();++i)
    			{
    			if(i>0) out.print(",");
    			pairs.get(i).json(out);
    			}
    		out.print("]}");
    		}
    	
    	
		}
    
    
    
    @Override
	public String getProgramDescription() {
		return "Scan folders and generate a structured summary of the files in JSON or XML";
		}
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -h get help (this screen)");
		out.println(" -v print version and exit.");
		out.println(" -L (level) log level. One of java.util.logging.Level . currently:"+getLogger().getLevel());
		out.println(" -J produces a JSON output");
		}
	

	@Override
	public int doWork(String[] args)
		{
		 boolean JSON=false;
		com.github.lindenb.jvarkit.util.cli.GetOpt getopt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=getopt.getopt(args, "hvL:J"))!=-1)
			{
			switch(c)
				{
				case 'h': printUsage();return 0;
				case 'v': System.out.println(getVersion());return 0;
				case 'L': getLogger().setLevel(java.util.logging.Level.parse(getopt.getOptArg()));break;
				case ':': System.err.println("Missing argument for option -"+getopt.getOptOpt());return -1;
				case 'J': JSON=true; break;
				default: System.err.println("Unknown option -"+getopt.getOptOpt());return -1;
				}
			}
		if(getopt.getOptInd()==args.length)
			{
			error("No directory given");
			return -1;
			}
				
		try
			{

			List<Folder> folders=new ArrayList<Folder>();
	    	for(int i=getopt.getOptInd();i< args.length;++i)
	    		{
	    		File f=new File(args[i]);
	    		if(!f.isDirectory()) throw new PicardException("Not a directory:"+f);
	    		Folder folder=new Folder();
	    		folder.folder=f;
	    		folder.scan(f);
	    		folders.add(folder);
	    		}
	    	if(JSON)
	    		{
	    		System.out.print("[");
	    		for(int i=0;i< folders.size();++i)
	    			{
	    			if(i>0) System.out.print(",");
	    			folders.get(i).json(System.out);
	    			}
	    		
	    		System.out.println("]");
	    		}
	    	else
	    		{
    			XMLOutputFactory xmlfactory= XMLOutputFactory.newInstance();
    			XMLStreamWriter w= xmlfactory.createXMLStreamWriter(System.out,"UTF-8");
    			w.writeStartDocument("UTF-8","1.0");
    			w.writeStartElement("illumina");
    			w.writeComment(this.getProgramCommandLine());
    			for(Folder f:folders) f.write(w);
    			w.writeEndElement();
    			w.writeEndDocument();
    			w.flush();
    			w.close();
	    		}
			
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			
			}
		return 0;
		}
   

	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new IlluminaDirectory().instanceMainWithExit(args);
		}

}
