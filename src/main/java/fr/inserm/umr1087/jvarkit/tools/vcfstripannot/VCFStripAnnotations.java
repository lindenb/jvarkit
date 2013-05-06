package fr.inserm.umr1087.jvarkit.tools.vcfstripannot;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

import fr.inserm.umr1087.jvarkit.util.vcf.VCFUtils;

public class VCFStripAnnotations
	{
	private PrintStream out=System.out;
	private Set<String> infoIds=new LinkedHashSet<String>();
	private VCFUtils vcfUtils=new VCFUtils();
	private boolean reset_filters=false;
	
	
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
					out.println("##Processed with "+getClass());
					out.println(line);
					continue;
					}
				else if(line.startsWith("##INFO=<"))
					{
					String tokens[]=line.substring(8).split("[<>=,]");//UGLY ! I know...
					if(tokens.length>2 &&
						tokens[0].equals("ID") &&
						this.infoIds.contains(tokens[1])
						)
						{
						continue;
						}
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
			Map<String,String> infos1=this.vcfUtils.parseInfo(tokens[7]);
			for(String f:this.infoIds)
				{
				infos1.remove(f);
				}
			
			
			String newinfo=this.vcfUtils.joinInfo(infos1);
			if(reset_filters)
				{
				tokens[6]=".";
				}
			
			
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
				System.out.println("VCFStripAnnotations. Author: Pierre Lindenbaum PhD. 2013.");
				System.out.println("Usage: java -jar vcfstripannotations.jar -T IF1 -T F2 -F (file.vcf|stdin) " );
				System.out.println(" -T (tag String) remove this VCF-INFO-ID");
				System.out.println(" -F reset the FILTER field.");
				return 0;
				}
			else if(args[optind].equals("-F"))
				{
				this.reset_filters=true;
				}
			else if(args[optind].equals("-T") && optind+1< args.length)
				{
				this.infoIds.add(args[++optind]);
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
		VCFStripAnnotations app=new VCFStripAnnotations();
		app.run(args);
		}
}
