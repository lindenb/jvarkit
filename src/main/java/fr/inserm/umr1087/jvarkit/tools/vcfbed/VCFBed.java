package fr.inserm.umr1087.jvarkit.tools.vcfbed;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

import javax.script.Bindings;
import javax.script.Compilable;
import javax.script.CompiledScript;
import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;

import fr.inserm.umr1087.jvarkit.util.vcf.VCFUtils;


import net.sf.samtools.tabix.TabixReader;


public class VCFBed
	{
	private CompiledScript  script=null;
	private ScriptEngine engine=null;

	private PrintStream out=System.out;
	private String tabixFile;
	private TabixReader tabixReader =null;
	private Set<String> extraHeaderLine=new LinkedHashSet<String>();
	private VCFUtils vcfUtils=new VCFUtils();
	
	
	public class Context
		{
		private String chrom;
		private int pos;
		private String ref;
		private String alt;
		
		private String id=null;
		private Set<String> filters=new LinkedHashSet<String>();
		private Map<String,String> info=new LinkedHashMap<String, String>();
		
		
		public String getChrom()
			{
			return chrom;
			}
		
		
		public int getPos()
			{
			return pos;
			}
		
		
		public String getRef()
			{
			return ref;
			}
				
		public String getAlt()
			{
			return alt;
			}
		
		public String getId()
			{
			return id;
			}
		
		public void setId(String id)
			{
			this.id = id;
			}
		public Set<String> getFilterSet()
			{
			return filters;
			}
		
		public Map<String, String> getInfoMap()
			{
			return info;
			}
		}
	
	private void run(BufferedReader in) throws Exception
		{
		Bindings bindings = this.engine.createBindings();
		Pattern tab=Pattern.compile("[\t]");
		String line;
		while((line=in.readLine())!=null)
			{

			if(line.isEmpty()) continue;
			
			if(line.startsWith("#"))
				{
				if(line.startsWith("#CHROM"))
					{
					out.println("##Annotated with "+getClass()+":"+tabixFile);
					for(String s:this.extraHeaderLine)
						{
						if(!s.startsWith("##")) out.print("##");
						out.println(s);
						}
					
					out.println(line);
					continue;
					}
				out.println(line);
				continue;
				}
			
			String tokens[]=tab.split(line,9);
			if(tokens.length<8)
				{
				System.err.println("[VCFTabix] Error not enough columns in vcf: "+line);
				continue;
				}
			String chrom=tokens[0];
			Integer pos1=Integer.parseInt(tokens[1]);
			
			
			List<String> tabixrows=new ArrayList<String>();
			TabixReader.Iterator iter=tabixReader.query(chrom+":"+pos1+"-"+(pos1+1));
			String line2;
			
			while(iter!=null && (line2=iter.next())!=null)
				{
				tabixrows.add(line2);
				}
			Context ctx=new Context();
			ctx.chrom=chrom;
			ctx.pos=pos1;
			ctx.id=tokens[2];
			ctx.ref=tokens[3];
			ctx.alt=tokens[4];
			ctx.filters=this.vcfUtils.parseFilters(tokens[6]);
			ctx.info=this.vcfUtils.parseInfo(tokens[7]);
			
			bindings.put("ctx",ctx);
			bindings.put("tabix",tabixrows.toArray());
			
			script.eval(bindings);
			
			
			tokens[2]=(this.vcfUtils.isEmpty(ctx.getId())?".":ctx.getId());
			tokens[6]=this.vcfUtils.joinFilters(ctx.getFilterSet());
			tokens[7]=this.vcfUtils.joinInfo(ctx.getInfoMap());
			
			for(int i=0;i< tokens.length;++i)
				{
				if(i>0) out.print('\t');
				out.print(tokens[i]);
				}
			out.println();
			}
		
		}
	public int run(String[] args) throws Exception
		{
		String scriptStr=null;
		File scriptFile=null;

		int optind=0;
		while(optind<args.length)
			{
			if(args[optind].equals("-h"))
				{
				System.out.println("VCFBed. Author: Pierre Lindenbaum PhD. 2013.");
				System.out.println("Usage: java -jar vcftabed.jar (-e script | -F scriptfile) -f src.vcf.gz (file.vcf|stdin) " );
				System.out.println(" -f (vcf indexed with tabix) REQUIRED.");
				System.out.println(" -t ( String) add extra info header line. Can be called multiple time.");
				System.out.println(" -T ( file) file containing extra info header line. Can be called multiple time.");
				System.out.println(" -e (script).");
				System.out.println(" -E (script file ).");
				return 0;
				}
			else if(args[optind].equals("-e") && optind+1< args.length)
				{
				scriptStr=args[++optind];
				}
			else if(args[optind].equals("-E") && optind+1< args.length)
				{
				scriptFile=new File(args[++optind]);
				}
			else if(args[optind].equals("-T") && optind+1< args.length)
				{
				File extra=new File(args[++optind]);
				BufferedReader in=new BufferedReader(new FileReader(extra));
				String line;
				while((line=in.readLine())!=null)
					{
					if(line.isEmpty()) continue;
					this.extraHeaderLine.add(line);
					}
				in.close();
				}
			else if(args[optind].equals("-T") && optind+1< args.length)
				{
				this.extraHeaderLine.add(args[++optind]);
				}
			else if(args[optind].equals("-f") && optind+1< args.length)
				{
				this.tabixFile=args[++optind];
				}
			else if(args[optind].equals("--"))
				{
				optind++;
				break;
				}
			else if(args[optind].startsWith("-"))
				{
				System.err.println("Unnown option: "+args[optind]);
				return -1;
				}
			else
				{
				break;
				}
			++optind;
			}
		
		ScriptEngineManager manager = new ScriptEngineManager();
		this.engine = manager.getEngineByName("js");
		if(this.engine==null)
			{
			System.err.println("not available: javascript. Use the SUN/Oracle JDK ?");
			System.exit(-1);
			}
		
		Compilable compilingEngine = (Compilable)this.engine;
		this.script = null;
		if(scriptFile!=null)
			{
			FileReader r=new FileReader(scriptFile);
			this.script=compilingEngine.compile(r);
			r.close();
			}
		else if(scriptStr!=null)
			{
			this.script=compilingEngine.compile(scriptStr);
			}
		else
			{
			System.err.println("Script missing.");
			System.exit(-1);
			}
		
		
		if(tabixFile==null)
			{
			System.err.println("Undefined tabix File");
			return -1;
			}
		
		this.tabixReader=new TabixReader(this.tabixFile);
		
		
		
		if(optind==args.length)
			{
			this.run(new BufferedReader(new InputStreamReader(System.in)));
			}
		else if(optind+1==args.length)
			{
			String inputName=args[optind++];
			BufferedReader in=new BufferedReader(new FileReader(inputName));
			this.run(in);
			in.close();
			}
		else
			{
			System.err.println("Illegal Number of arguments");
			return -1;
			}
		return 0;
		}
	
	public static void main(String[] args) throws Exception
		{
		VCFBed app=new VCFBed();
		app.run(args);
		}
}
