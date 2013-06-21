package com.github.lindenb.jvarkit.tools.vcf2sql;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.regex.Pattern;

import org.broad.tribble.readers.AsciiLineReader;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.CommonInfo;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.GenotypesContext;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.vcf.VCFCodec;
import org.broadinstitute.variant.vcf.VCFConstants;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;
import org.broadinstitute.variant.vcf.VCFHeaderVersion;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;

@SuppressWarnings("rawtypes")
public class VcfToSql extends CommandLineProgram
	{
	private static final String COLUMN_ID="id INTEGER PRIMARY KEY AUTOINCREMENT,";
	
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+"Creates the code to insert one or more VCF into a SQL database. ";
    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="VCF files to process.",minElements=0)
	public List<File> IN=new ArrayList<File>();
    
    @Option(shortName="SFX",doc="Table suffix",optional=true)
	public String SUFFIX="";
    @Option(shortName="VEP",doc="Use  and explode VEP predictions",optional=true)
	public boolean USE_VEP=true;
    @Option(shortName="SNPEFF",doc="Use and explode SNPEFF predictions",optional=true)
	public boolean USE_SNPEFF=true;
    @Option(shortName="SQLIDX",doc="Create misc SQL Indexes.",optional=true)
	public boolean SQLINDEX=true;

    
    
    private static Log LOG=Log.getInstance(VcfToSql.class); 
    
    private PrintWriter out=new PrintWriter(System.out);
    
    @Override
    public String getVersion()
    	{
    	return "1.0";
    	}
    
    
	@Override
	protected int doWork()
		{
		try
			{
			
		    out.println( "create table if not exists FILE"+SUFFIX+"("+
		    		COLUMN_ID+
		    		"filename TEXT NOT NULL"+
		    		");");

		    out.println( "create table if not exists HEADER"+SUFFIX+"("+
		    		COLUMN_ID+
					"file_id INT NOT NULL REFERENCES FILE"+SUFFIX+"(id) ON DELETE CASCADE,"+
		    		"header TEXT"+
		    		");");

		    
		    out.println( "create table if not exists SAMPLE"+SUFFIX+"("+
			    COLUMN_ID+
			    "name TEXT NOT NULL UNIQUE"+
			   ");");
		    out.println( "create table if not exists VARIATION"+SUFFIX+"("+
				COLUMN_ID+
				"file_id INT NOT NULL REFERENCES FILE"+SUFFIX+"(id) ON DELETE CASCADE,"+
				"CHROM VARCHAR(20) NOT NULL,"+
				"POS INT NOT NULL,"+
				"START0 INT NOT NULL,"+
				"END0 INT NOT NULL,"+
				"RS_ID VARCHAR(50),"+
				"REF VARCHAR(10) NOT NULL,"+
				"QUAL FLOAT"+
			    ");");
		    
		    out.println("create table if not exists ALT"+SUFFIX+"("+
				COLUMN_ID+
				"var_id INT NOT NULL REFERENCES VARIATION"+SUFFIX+"(id) ON DELETE CASCADE,"+
				"ALT TEXT"+
				  ");");
		    out.println("create table if not exists FILTER"+SUFFIX+"("+
				COLUMN_ID+
				"var_id INT NOT NULL REFERENCES VARIATION"+SUFFIX+"(id) ON DELETE CASCADE,"+
				"FILTER varchar(50) not null"+
				 ");");

		    out.println("create table if not exists INFO"+SUFFIX+"("+
				COLUMN_ID+
				"var_id INT NOT NULL REFERENCES VARIATION"+SUFFIX+"(id) ON DELETE CASCADE,"+
				"k varchar(50) not null,"+
				"v TEXT not null"+
				 ");");
		    
		    out.println("create table if not exists EXTRAINFO"+SUFFIX+"("+
					COLUMN_ID+
					"info_id INT NOT NULL REFERENCES INFO"+SUFFIX+"(id) ON DELETE CASCADE,"+
					"type varchar(50) not null"+
					 ");");
		    
		    out.println("create table if not exists EXTRAINFOPROP"+SUFFIX+"("+
					COLUMN_ID+
					"extrainfo_id INT NOT NULL REFERENCES EXTRAINFO"+SUFFIX+"(id) ON DELETE CASCADE,"+
					"k varchar(50) not null,"+
					"v TEXT not null"+
					 ");");

		    
		    out.println(     "create table if not exists GENOTYPE"+SUFFIX+"("+
				COLUMN_ID+
				"var_id INT NOT NULL REFERENCES VARIATION"+SUFFIX+"(id) ON DELETE CASCADE,"+
				"sample_id INT NOT NULL REFERENCES SAMPLE"+SUFFIX+"(id) ON DELETE CASCADE,"+
				"A1 TEXT, A2 TEXT, dp int, ad varchar(50), gq float,pl TEXT,"+
				"is_phased SMALLINT not null,is_hom SMALLINT not null,is_homref  SMALLINT not null,is_homvar  SMALLINT not null,is_mixed  SMALLINT not null," +
				"is_nocall SMALLINT not null,is_noninformative SMALLINT not null,is_available SMALLINT not null,is_called SMALLINT not null,is_filtered  SMALLINT not null"+
				");"
				);
		    out.println("create table if not exists GTPROP"+SUFFIX+"("+
				COLUMN_ID+
				"genotype_id INT NOT NULL REFERENCES GENOTYPE"+SUFFIX+"(id) ON DELETE CASCADE,"+
				"k varchar(50) not null,"+
				"v TEXT not null"+
				 ");");

			out.println("begin transaction;");
			if(IN.isEmpty())
				{
				read(System.in,"<stdin>");
				}
			else
				{
				for(File input: IN)
					{
					InputStream in=IoUtil.openFileForReading(input);
					read(in,input.toString());
					in.close();
					}
				}
			if(SQLINDEX)
				{
				index("SAMPLE","name");
				index("EXTRAINFO","type");
				index("EXTRAINFOPROP","k");
				index("EXTRAINFOPROP","v");
				
				index("INFO","var_id");
				index("INFO","k");
				index("EXTRAINFO","info_id");
				index("EXTRAINFOPROP","extrainfo_id");
				index("GENOTYPE","var_id");
				index("GENOTYPE","sample_id");
				}
			out.println("commit;");
			out.flush();
			}
		catch(IOException err)
			{
			err.printStackTrace();
			return -1;
			}
		return 0;
		}
	
	private void index(String table,String column)
		{
		out.println("create index if not exists  "+
				(table+SUFFIX+"_"+column+"_IDX").toUpperCase() +
				" on "+table+SUFFIX+"("+column+");");
		}
	
	private void read(InputStream in,String filename)
		throws IOException
		{
		//Pattern comma=Pattern.compile("[,]");
		Pattern pipe=Pattern.compile("[\\|]");
		Pattern amp=Pattern.compile("&");

		out.println("insert into FILE"+SUFFIX+"(filename) values ("+quote(filename)+");");
		AsciiLineReader r=new AsciiLineReader(in);
		

		VCFCodec codec=new VCFCodec();
		VCFHeader header=(VCFHeader)codec.readHeader(r);
		
		String csqColumns[]=null;
		VCFInfoHeaderLine infoHeader=header.getInfoHeaderLine("CSQ");
		if(infoHeader!=null && this.USE_VEP)
			{
			LOG.info("parsing VEP "+infoHeader.getDescription());
			final String formatStr="Format: ";
			int i=infoHeader.getDescription().indexOf(formatStr);
			if(i!=-1)
				{
				csqColumns=pipe.split(infoHeader.getDescription().substring(i+formatStr.length()).trim());
				LOG.debug(Arrays.asList(csqColumns));
				}
			else
				{
				LOG.error("Cannot parse "+infoHeader.getDescription());
				}
			}
		String snpEffColumns[]=null;
		infoHeader=header.getInfoHeaderLine("EFF");
		if(infoHeader!=null && this.USE_SNPEFF)
			{
			LOG.info("parsing EFF "+infoHeader.getDescription());

			final String formatStr=".Format: '";
			final String desc=infoHeader.getDescription();
			int i=desc.indexOf(formatStr);
			if(i!=-1) i=desc.indexOf('(',i+formatStr.length());
			int j=desc.lastIndexOf(')');
			if(i!=-1 && j>i)
				{
				snpEffColumns=pipe.split(desc.substring(i+1,j).replaceAll("[ \\[\\]()\\.]", "").trim());
				LOG.info(Arrays.asList(snpEffColumns));
				}
			else
				{
				LOG.error("Cannot parse "+infoHeader.getDescription());
				}
			}
		
		String nmdColumns[]=null;
		infoHeader=header.getInfoHeaderLine("NMD");
		if(infoHeader!=null && this.USE_SNPEFF)
			{

			final String formatStr=" Format: '";
			final String desc=infoHeader.getDescription();
			int i=desc.indexOf(formatStr);
			int j=(i==-1?-1:desc.lastIndexOf('\''));

			if(i!=-1 && j>i)
				{
				nmdColumns=pipe.split(
						desc.substring(i+formatStr.length(),j).replaceAll("[ \\[\\]()\\.]", "").trim());
				}
			else
				{
				LOG.error("Cannot parse "+infoHeader.getDescription());
				}
			}
		
		String lofColumns[]=null;
		infoHeader=header.getInfoHeaderLine("LOF");
		if(infoHeader!=null && this.USE_SNPEFF)
			{

			final String formatStr=" Format: '";
			final String desc=infoHeader.getDescription();
			int i=desc.indexOf(formatStr);
			int j=(i==-1?-1:desc.lastIndexOf('\''));

			if(i!=-1 && j>i)
				{
				lofColumns=pipe.split(
						desc.substring(i+formatStr.length(),j).replaceAll("[ \\[\\]()\\.]", "").trim());
				}
			else
				{
				LOG.error("Cannot parse "+infoHeader.getDescription());
				}
			}
		

		
		for(String S:header.getSampleNamesInOrder())
			{
			out.println("insert or ignore into SAMPLE"+SUFFIX+"(name) values ("+quote(S)+");");
			}
		
		List<String> headers=new ArrayList<String>();

         for ( VCFHeaderLine line : header.getMetaDataInSortedOrder() )
         	{
             if ( VCFHeaderVersion.isFormatString(line.getKey()) )
                 continue;
             headers.add(VCFHeader.METADATA_INDICATOR+line);
             }	
         
         String chromLine=VCFHeader.HEADER_INDICATOR;
         for ( VCFHeader.HEADER_FIELDS field : header.getHeaderFields() )
         	{
         	if(!VCFHeader.HEADER_INDICATOR.equals(chromLine))chromLine+=(VCFConstants.FIELD_SEPARATOR);
         	chromLine+=(field);
         	}

         if ( header.hasGenotypingData() )
         	 {
        	 chromLine+=VCFConstants.FIELD_SEPARATOR+"FORMAT";
             for ( String sample : header.getGenotypeSamples() ) {
            	 chromLine+=VCFConstants.FIELD_SEPARATOR;
                chromLine+=sample;
             	}
         	}
         headers.add(chromLine); 

         for(String line:headers)
         	{
 			out.println(
					"insert into HEADER"+SUFFIX+
					"(file_id,header) values ("+
					"(select max(id) from FILE"+SUFFIX+"),"+
					quote(line)+");"
					);
         	}
		
		
		String line;
		
		while((line=r.readLine())!=null)
			{
			LOG.debug(line);
			
			VariantContext var=codec.decode(line);
			if(var==null)
				{
				LOG.error("Cannot parse "+line);
				continue;
				}
			//"create table if not exists FILE(id,filename text)";
			//"create table if not exists VARIATION(id,file_id,chrom,pos,start0,end0,rs_id,ref,qual)";
		
			
			out.println(
					"insert into VARIATION"+SUFFIX+
					"(file_id,chrom,pos,START0,END0,rs_id,ref,qual) values ("+
					"(select max(id) from FILE"+SUFFIX+"),"+
					quote(var.getChr())+","+
					var.getStart()+","+
					(var.getStart()-1)+","+
					var.getEnd()+","+
					(var.getID()==null || var.getID().equals(VCFConstants.EMPTY_ID_FIELD) ?"NULL":quote(var.getID()))+","+
					quote(var.getReference().getDisplayString())+","+
					(var.getPhredScaledQual()<0?"NULL":var.getPhredScaledQual())+");"
					);
			//"create table if not exists ALT(id,var_id,alt)";

			for(Allele alt: var.getAlternateAlleles())
				{
				out.println(
						"insert into ALT"+SUFFIX+"(var_id,alt) values ("+
						"(select max(id) from VARIATION"+SUFFIX+"),"+
						quote(alt.getDisplayString())+");"
						);
				}
			//"create table if not exists FILTER(id,var_id,filter)";

			for(String filter:var.getFilters())
				{
				out.println(
						"insert into FILTER"+SUFFIX+"(var_id,filter) values ("+
						"(select max(id) from VARIATION"+SUFFIX+"),"+
						quote(filter)+");"
						);
				}
			CommonInfo infos=var.getCommonInfo();
			for(String key:infos.getAttributes().keySet())
				{
				Object val=infos.getAttribute(key);
				//"create table if not exists INFO(id,var_id,k,v)";

				out.println(
						"insert into INFO"+SUFFIX+"(var_id,k,v) values ("+
						"(select max(id) from VARIATION"+SUFFIX+"),"+
						quote(key)+","+
						quote(infotoString(val))+");"
						);
					
				
				if(key.equals("CSQ") && csqColumns!=null)
					{
					List as_array=castToStringArray(val);
					
					for(Object csqs:as_array)
						{
						if(csqs.toString().isEmpty()) continue;
						String tokens[]=pipe.split(csqs.toString());
						List<String> extraInfo=new ArrayList<String>();
						for(int t=0;t<tokens.length && t<csqColumns.length;++t)
							{
							if(tokens[t].isEmpty()) continue;
							if(csqColumns[t].equals("Consequence"))
								{
								for(String pred:amp.split(tokens[t]))
									{
									if(pred.isEmpty()) continue;
									extraInfo.add(csqColumns[t]);
									extraInfo.add(pred);
									}
								
								}
							else
								{
								extraInfo.add(csqColumns[t]);
								extraInfo.add(tokens[t]);
								}
							}
						insertExtraInfos("CSQ",extraInfo);
						}
					}
				
				if(key.equals("EFF") && snpEffColumns!=null)
					{
					for(Object item:castToStringArray(val))
						{
						String snpeff=item.toString();
						if(snpeff.isEmpty()) continue;
						int opar=snpeff.indexOf('(');
						if(opar==-1) continue;
						int cpar=snpeff.lastIndexOf(')');
						if(cpar==-1) continue;
						String tokens[]=pipe.split(snpeff.substring(opar+1,cpar));
						List<String> h=new ArrayList<String>();
						h.add("Effect");
						h.add(snpeff.substring(0,opar));
						for(int t=0;t<tokens.length && t<snpEffColumns.length;++t)
							{
							if(tokens[t].isEmpty()) continue;
							h.add(snpEffColumns[t]);
							h.add(tokens[t]);
							}
						insertExtraInfos(key, h);	
						}
					}
				
				if(key.equals("NMD") && nmdColumns!=null)
					{
					
					for(Object item:castToStringArray(val))
						{
						String nmd=item.toString();
						if(nmd.isEmpty()) continue;
						String tokens[]=pipe.split(nmd);
						List<String> h=new ArrayList<String>(nmdColumns.length*2);
						for(int t=0;t<tokens.length && t<nmdColumns.length;++t)
							{
							if(tokens[t].isEmpty()) continue;
							h.add(nmdColumns[t]);
							h.add(tokens[t]);
							}
						insertExtraInfos(key, h);
						}	
					}
				
				if(key.equals("LOF") && lofColumns!=null)
					{
					
					for(Object item:castToStringArray(val))
						{
						String lof=item.toString();
						if(lof.isEmpty()) continue;
						String tokens[]=pipe.split(lof);
						List<String> h=new ArrayList<String>(lofColumns.length*2);
						for(int t=0;t<tokens.length && t<lofColumns.length;++t)
							{
							if(tokens[t].isEmpty()) continue;
							h.add(lofColumns[t]);
							h.add(tokens[t]);
							}
						insertExtraInfos(key, h);
						}	
					}
				
				
				}
			GenotypesContext genotypesCtx =var.getGenotypes();
			for(Genotype g:genotypesCtx)
				{
				//"create table if not exists GENOTYPE(id,var_id,k,v)";
				
				List<Allele> alleles=g.getAlleles();
				
				out.println(
						"insert into GENOTYPE"+SUFFIX+
						"(var_id,sample_id,A1,A2,dp,ad,gq,pl," +
						"is_phased,is_hom,is_homref,is_homvar,is_mixed," +
						"is_nocall,is_noninformative,is_available,is_called,is_filtered"+
						") values ("+
						"(select max(id) from VARIATION"+SUFFIX+"),"+
						"(select id from SAMPLE"+SUFFIX+" where name="+quote(g.getSampleName())+"),"+
						(alleles.size()==2?quote(alleles.get(0).getBaseString()):"NULL")+","+
						(alleles.size()==2?quote(alleles.get(1).getBaseString()):"NULL")+","+
						(g.hasDP()?g.getDP():"NULL")+","+
						(g.hasAD()?quote(infotoString(g.getAD())):"NULL")+","+
						(g.hasGQ()?g.getGQ():"NULL")+","+
						(g.hasPL()?quote(infotoString(g.getPL())):"NULL")+","+
						(g.isPhased()?1:0)+","+
						(g.isHom()?1:0)+","+
						(g.isHomRef()?1:0)+","+
						(g.isHomVar()?1:0)+","+
						(g.isMixed()?1:0)+","+
						(g.isNoCall()?1:0)+","+
						(g.isNonInformative()?1:0)+","+
						(g.isAvailable()?1:0)+","+
						(g.isCalled()?1:0)+","+
						(g.isFiltered()?1:0)+");"
						);
				
				for(String key:g.getExtendedAttributes().keySet())
					{
					Object val=g.getExtendedAttribute(key);
					if(val==null) continue;
					out.println(
							"insert into GTPROP"+SUFFIX+"(genotype_id,k,v) values ("+
							"(select max(id) from GENOTYPE"+SUFFIX+"),"+
							quote(key)+","+
							quote(infotoString(val))+");"
							);
					}

				}
			
			}
		
		}
	
	private String quote(String s)
		{
		if(s==null) return "NULL";
		StringBuilder b=new StringBuilder();
		b.append("\'");
		for(int i=0;i<s.length();++i)
			{
			char c=s.charAt(i);
			switch(c)
				{
				case '\'': b.append("''"); break;
				default: b.append(c); break;
				}
			}
		b.append("\'");
		return b.toString();
		}
	private void insertExtraInfos(String type,List<String> h)
		{
		boolean first=true;
		for(int i=0;i+1< h.size();i+=2)
			{			
			if(h.get(i+1).isEmpty()) continue;
			if(first)
				{
				
				out.println(
						"insert into EXTRAINFO"+SUFFIX+"(info_id,type) values ("+
						"(select max(id) from INFO"+SUFFIX+"),"+
						quote(type)+
						");"
						);
				first=false;
	
				}
			
			out.println(
					"insert into EXTRAINFOPROP"+SUFFIX+"(extrainfo_id,k,v) values ("+
					"(select max(id) from EXTRAINFO"+SUFFIX+"),"+
					quote(h.get(i))+","+
					quote(h.get(i+1))+");"
					);
			}

		}
	
	@SuppressWarnings("unchecked")
	private List castToStringArray(Object val)
		{
		if(val instanceof List)
			{
			return (List)val;
			}
		else
			{
			return new ArrayList(Collections.singleton(val.toString()));
			}
		}
	private String infotoString(Object o)
		{
		if(o instanceof int[])
			{
			int array[]=(int[])o;
			StringBuilder b=new StringBuilder();
			for(int i=0;i< array.length;++i)
				{
				if(i>0) b.append(",");
				b.append(infotoString(array[i]));
				}
			return b.toString();
			}
		if(o instanceof List)
			{
			List<?> L=List.class.cast(o);
			StringBuilder b=new StringBuilder();
			for(int i=0;i< L.size();++i)
				{
				if(i>0) b.append(",");
				b.append(infotoString(L.get(i)));
				}
			return b.toString();
			}
		return o.toString();
		}
	public static void main(String[] args)
		{
		new VcfToSql().instanceMainWithExit(args);
		}
	}
