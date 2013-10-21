package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import net.sf.picard.PicardException;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.util.Log;

import com.github.lindenb.jvarkit.util.illumina.FastQName;
import com.github.lindenb.jvarkit.util.picard.AbstractCommandLineProgram;

public class IlluminaDirectory extends AbstractCommandLineProgram
	{
	private static Log LOG=Log.getInstance(IlluminaDirectory.class);

    @Usage(programVersion="1.0")
    public String USAGE = getStandardUsagePreamble() + " Scan folders and generate a structured summary of the files. ";

	
    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME,doc="root directories",minElements=0)
    public Set<File> IN=new HashSet<File>();
    
    @Option(shortName="J",doc="json output",optional=true)
    public boolean JSON=false;
    
    
    private class Folder
    	{
    	File folder;
    	SortedMap<String, Sample> sampleMap=new TreeMap<String, IlluminaDirectory.Sample>();
    	List<Pair> undtermined=new ArrayList<Pair>();
    	void scan(File f)
    		{
    		if(f==null) return;
    		if(!f.canRead()) return;
    		LOG.debug("Scanning "+f);
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
    					LOG.warn("invalid name:"+fq);
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
				p.json(out);
				}
    		out.print("]");
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
    	FastQName forward;
    	FastQName reverse;
    	
    	Pair(FastQName fq)
    		{
    		if(fq.getSide()==1)
    			{
    			forward=fq;
    			}
    		else if(fq.getSide()==2)
    			{
    			reverse=fq;
    			}
    		else
    			{
    			throw new PicardException("bad side "+fq);
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
	    		out.print("\"path\":\""+forward.getFile().getPath()+"\",");
	    		out.print("\"side\":"+forward.getSide()+",");
	    		out.print("\"file-size\":"+forward.getFile().length());
	    		out.print("},\"reverse\":{");
	    		out.print("\"path\":\""+reverse.getFile().getPath()+"\",");
	    		out.print("\"side\":"+reverse.getSide()+",");
	    		out.print("\"file-size\":"+reverse.getFile().length());
	    		out.print("}}");
				}
    		else
    			{
    			FastQName F=(forward==null?reverse:forward);
    			out.print("{");
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
    			w.writeAttribute("lane", String.valueOf(forward.getLane()));
    			if(forward.getSeqIndex()!=null) w.writeAttribute("index", String.valueOf(forward.getSeqIndex()));
    			w.writeAttribute("split", String.valueOf(forward.getSplit()));
    			
    			w.writeStartElement("fastq");
    			w.writeAttribute("side", String.valueOf(forward.getSide()));
    			w.writeAttribute("path", forward.getFile().getPath());
    			w.writeAttribute("file-size",String.valueOf( forward.getFile().length()));
    			w.writeEndElement();
    			
    			
    			w.writeStartElement("fastq");
    			w.writeAttribute("side", String.valueOf(reverse.getSide()));
    			w.writeAttribute("path", reverse.getFile().getPath());
    			w.writeAttribute("file-size",String.valueOf( reverse.getFile().length()));
    			w.writeEndElement();
   			
    			w.writeEndElement();
    			}
    		else
    			{
    			FastQName F=(forward==null?reverse:forward);
    			
    			w.writeStartElement("single");
    			w.writeAttribute("lane", String.valueOf(F.getLane()));
    			if(forward.getSeqIndex()!=null) w.writeAttribute("index", String.valueOf(F.getSeqIndex()));
    			w.writeAttribute("split", String.valueOf(F.getSplit()));
    			
    			w.writeStartElement("fastq");
    			w.writeAttribute("side", String.valueOf(F.getSide()));
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
    protected int doWork()
    	{
    	List<Folder> folders=new ArrayList<Folder>();
    	for(File f:this.IN)
    		{
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
    		try
    			{
    			XMLOutputFactory xmlfactory= XMLOutputFactory.newInstance();
    			XMLStreamWriter w= xmlfactory.createXMLStreamWriter(System.out,"UTF-8");
    			w.writeStartDocument("UTF-8","1.0");
    			w.writeStartElement("illumina");
    			w.writeComment(getCommandLine());
    			for(Folder f:folders) f.write(w);
    			w.writeEndElement();
    			w.writeEndDocument();
    			w.flush();
    			w.close();
    			}
    		catch(Exception err)
    			{
    			err.printStackTrace();
    			LOG.error(err);
    			return -1;
    			}
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
