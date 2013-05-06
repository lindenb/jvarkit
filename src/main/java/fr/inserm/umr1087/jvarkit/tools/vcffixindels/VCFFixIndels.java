package fr.inserm.umr1087.jvarkit.tools.vcffixindels;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;


public class VCFFixIndels
	{
	private static final Logger LOG=Logger.getLogger(VCFFixIndels.class.getSimpleName());
	private PrintStream out=System.out;
	
	
	private void run(BufferedReader in) throws IOException
		{
		Pattern tab=Pattern.compile("[\t]");
		Pattern dna=Pattern.compile("[ATGCatgc]+");
		String line;
		while((line=in.readLine())!=null)
			{

			if(line.isEmpty()) continue;
			
			if(line.startsWith("#"))
				{
				if(line.startsWith("#CHROM"))
					{
					out.println("##Processed with "+getClass());
					}
				out.println(line);
				continue;
				}
			
			String tokens[]=tab.split(line,6);
			if(tokens.length<5)
				{
				System.err.println("[VCFTabix] Error not enough columns in vcf: "+line);
				continue;
				}
			if(tokens[3].length()!=tokens[4].length() &&
				dna.matcher(tokens[3]).matches() &&
				dna.matcher(tokens[4]).matches() )
				{
				int pos=Integer.parseInt(tokens[1]);
				StringBuffer ref=new StringBuffer(tokens[3].toUpperCase());
				StringBuffer alt=new StringBuffer(tokens[4].toUpperCase());
				boolean changed=false;
				//REF=AA ALT=AAT
				while(		ref.length()>1 &&
							ref.length() <alt.length() &&
							ref.charAt(0)==alt.charAt(0)
							)
					{
					changed=true;
					ref.deleteCharAt(0);
					alt.deleteCharAt(0);
					pos++;
					}
				//REF=AAT ALT=AA
				while(	alt.length()>1 &&
						alt.length() < ref.length() &&
						ref.charAt(0)==alt.charAt(0)
						)
					{
					changed=true;
					ref.deleteCharAt(0);
					alt.deleteCharAt(0);
					pos++;
					}
				
				if(changed)
					{
					LOG.info("changed "+tokens[1]+"/"+tokens[3]+"/"+tokens[4]+" to "+pos+"/"+ref+"/"+alt);
					}
				
				tokens[1]=String.valueOf(pos);
				tokens[3]=ref.toString();
				tokens[4]=alt.toString();
				}
		
			
			
			
			for(int i=0;i< tokens.length;++i)
				{
				if(i>0) out.print('\t');
				out.print(tokens[i]);
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
				System.out.println("VCFFixIndels. Author: Pierre Lindenbaum PhD. 2013.");
				System.out.println("Usage: java -jar vcffixindels.jar (options) (file.vcf|stdin) " );
				System.out.println(" -L (log-level)");
				return 0;
				}
			
			else if(args[optind].equals("-L") && optind+1< args.length)
				{
				LOG.setLevel(Level.parse(args[++optind]));
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
		LOG.setLevel(Level.OFF);
		VCFFixIndels app=new VCFFixIndels();
		app.run(args);
		}
}
