package fr.inserm.umr1087.jvarkit.tools.vcftabix;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

import fr.inserm.umr1087.jvarkit.util.vcf.VCFUtils;


import net.sf.samtools.tabix.TabixReader;
import net.sf.samtools.util.BlockCompressedInputStream;


public class VCFTabix
	{
	private PrintStream out=System.out;
	private String tabixFile;
	private TabixReader tabixReader =null;
	private Set<String> infoIds=new LinkedHashSet<String>();
	private boolean replaceInfoField=true;
	private boolean refMatters=true;
	private boolean altMatters=true;
	private boolean replaceID=true;
	private String altConflictTag=null;
	private VCFUtils vcfUtils=new VCFUtils();
	
	
	
	private void run(BufferedReader in) throws IOException
		{
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
					if(!this.infoIds.isEmpty())
						{
						BlockCompressedInputStream src=new BlockCompressedInputStream(new File(tabixFile));
						BufferedReader r2=new BufferedReader(new InputStreamReader(src));
						String line2;
						while((line2=r2.readLine())!=null)
							{
							if(line2.startsWith("#CHROM")) break;
							if(!line2.startsWith("##")) break;
							if(!line2.startsWith("##INFO=")) continue;
							
							int i=line2.indexOf("ID=");
							if(i==-1) continue;
							int j=line2.indexOf(',',i+1);
							if(j==-1) j=line2.indexOf('>',i+1);
							if(j==-1) continue;
							if(this.infoIds.contains(line2.substring(i+3,j)))
								{
								this.out.println(line2);
								}
							} 
						r2.close();
						}
					if(altConflictTag!=null)
						{
						out.println("##INFO=<ID="+this.altConflictTag+",Number=1,Type=String,Description=\"Conflict ALT alleles with '"+ tabixFile+"' \">");
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
			Map<String,String> infos1=this.vcfUtils.parseInfo(tokens[7]);
		
			
			TabixReader.Iterator iter=tabixReader.query(chrom+":"+pos1+"-"+(pos1+1));
			String line2;
			
			while(iter!=null && (line2=iter.next())!=null)
				{

				String tokens2[]=tab.split(line2,9);
				
				if(tokens2.length<8)
					{
					System.err.println("[VCFTabix]. Error not enough columns in tabix line "+line2);
					return;
					}
				if(!tokens[0].equals(tokens2[0])) continue;
				if(!tokens[1].equals(tokens2[1])) continue;
				if(this.refMatters && !tokens[3].equalsIgnoreCase(tokens2[3])) continue;
				
				
				Map<String,String> infos2=this.vcfUtils.parseInfo(tokens2[7]);
				
				if(!tokens[4].equalsIgnoreCase(tokens2[4]))
					{
					if(this.altMatters) continue;
					}
				
				
				if(altConflictTag!=null && !tokens[4].equalsIgnoreCase(tokens2[4]))
					{
					infos1.put(this.altConflictTag,tokens2[4]);
					}
				
				if(this.vcfUtils.isEmpty(tokens[2]) && !this.vcfUtils.isEmpty(tokens2[2]))
					{
					tokens[2]=tokens2[2];
					}
				else if(!this.vcfUtils.isEmpty(tokens[2]) && !this.vcfUtils.isEmpty(tokens2[2]) && !tokens[2].equals(tokens2[2]))
					{
					System.err.println("[WARNING]Not same ID for:\n[WARNING]  "+line+"\n[WARNING]  "+line2);
					if(this.replaceID)
						{
						tokens[2]=tokens2[2];
						}
					}
				
				
				for(String id:this.infoIds)
					{
					if(!infos2.containsKey(id)) continue;
					if(infos1.containsKey(id) && !this.replaceInfoField)
						{
						continue;
						}
					infos1.put(id,infos2.get(id));
					}
				}
			String newinfo=this.vcfUtils.joinInfo(infos1);
			
			
			
			for(int i=0;i< tokens.length;++i)
				{
				if(i>0) out.print('\t');
				out.print(i==7?newinfo.toString():tokens[i]);
				}
			out.println();
			}
		
		}
	public int run(String[] args) throws IOException
		{
		int optind=0;
		while(optind<args.length)
			{
			if(args[optind].equals("-h"))
				{
				System.out.println("VCF Tabix. Author: Pierre Lindenbaum PhD. 2013.");
				System.out.println("Usage: java -jar vcftabix.jar -f src.vcf.gz (file.vcf|stdin) " );
				System.out.println(" -f (vcf indexed with tabix) REQUIRED.");
				System.out.println(" -T (tag String) VCF-INFO-ID optional can be used several times.");
				System.out.println(" -R doesn't use REF allele");
				System.out.println(" -A doesn't use ALT allele");
				System.out.println(" -I don't replace ID if it exists.");
				System.out.println(" -F don't replace INFO field if it exists.");
				System.out.println(" -C (TAG) use this tag in case of conflict with the ALT allele.");
				return 0;
				}
			else if(args[optind].equals("-C")&& optind+1< args.length)
				{
				this.altConflictTag=args[++optind];
				}
			else if(args[optind].equals("-F"))
				{
				this.replaceInfoField=false;
				}
			else if(args[optind].equals("-R"))
				{
				this.refMatters=false;
				}
			else if(args[optind].equals("-A"))
				{
				this.altMatters=false;
				}
			else if(args[optind].equals("-I"))
				{
				this.replaceID=false;
				}
			else if(args[optind].equals("-T") && optind+1< args.length)
				{
				this.infoIds.add(args[++optind]);
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
	
	public static void main(String[] args) throws IOException
		{
		VCFTabix app=new VCFTabix();
		app.run(args);
		}
}
