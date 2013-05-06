package fr.inserm.umr1087.jvarkit.tools.vcfbigwig;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;

public class VCFBigWig
	{
	private PrintStream out=System.out;
	private String bigWigFile;
	private BBFileReader bbFileReader=null;
	private boolean contained=true;
	private String infoId=null;
	private void run(BufferedReader in) throws IOException
		{
		Pattern tab=Pattern.compile("[\t]");
		List<Float> values=new ArrayList<Float>();

		String line;
		while((line=in.readLine())!=null)
			{
			if(line.isEmpty()) continue;
			
			if(line.startsWith("#"))
				{
				if(line.startsWith("#CHROM"))
					{
					out.println("##Annotated with "+getClass()+":"+bigWigFile);
					out.println("##INFO=<ID="+this.infoId+",Number=1,Type=Float,Description=\"Annotations from "+this.bigWigFile+"\">");

					out.println(line);
					continue;
					}
				out.println(line);
				continue;
				}
			
			String tokens[]=tab.split(line,9);
			if(tokens.length<8)
				{
				System.err.println("Error not enough columns in "+line);
				continue;
				}
			String chrom=tokens[0];
			Integer pos1=Integer.parseInt(tokens[1]);
			
			values.clear();
			
			BigWigIterator iter=this.bbFileReader.getBigWigIterator(chrom, pos1-1, chrom, pos1, this.contained);
			while(iter.hasNext())
				{
				WigItem item=iter.next();
				float v=item.getWigValue();
				values.add(v);
				
				}
			
			if(values.isEmpty())
				{
				out.println(line);
				continue;
				}
			double total=0L;
			for(Float f:values) total+=f;
			
			String newinfp=this.infoId+"="+String.format("%.2f",(float)(total/values.size()));
			String info=tokens[7];
			if(info.equals(".") || info.isEmpty())
				{
				info=newinfp;
				}
			else
				{
				info=newinfp+";"+info;
				}
			for(int i=0;i< tokens.length;++i)
				{
				if(i>0) out.print('\t');
				out.print(i==7?info:tokens[i]);
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
				System.out.println(" -f (bigwig-file) required.");
				System.out.println(" -i (String) VCF-ID optional.");
				return 0;
				}
			else if(args[optind].equals("-i") && optind+1< args.length)
				{
				this.infoId=args[++optind];
				}
			else if(args[optind].equals("-f") && optind+1< args.length)
				{
				this.bigWigFile=args[++optind];
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
		if(bigWigFile==null)
			{
			System.err.println("Undefined bigwigFile");
			return -1;
			}
		this.bbFileReader=new BBFileReader(this.bigWigFile);
		if(!this.bbFileReader.isBigWigFile())
			{
			System.err.println(this.bigWigFile+" is not a bigWIG file.");
			return -1;
			}
		
		if(infoId==null)
			{
			infoId=this.bigWigFile;
			int i=infoId.lastIndexOf(File.separator);
			if(i!=-1) infoId=infoId.substring(i+1);
			i=infoId.indexOf('.');
			infoId=infoId.substring(0,i);
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
		VCFBigWig app=new VCFBigWig();
		app.run(args);
		}
}
