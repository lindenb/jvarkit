package com.github.lindenb.jvarkit.tools.ws.server;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.jws.WebService;
import javax.xml.ws.Endpoint;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.util.Log;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

import org.broad.tribble.readers.TabixReader;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.vcf.VCFCodec;

@WebService(
		endpointInterface="com.github.lindenb.jvarkit.tools.ws.server.NGSService"
		)
public class NGSServiceImpl 
	extends CommandLineProgram
	implements NGSService
	{
	public Log LOG=Log.getInstance(NGSServiceImpl.class);
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+"We services for BAM and VCFs. ";
    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="file and directories.",minElements=1,optional=false)
	public List<File> IN=new ArrayList<File>();
    @Option(shortName= "HOST", doc="endpoint URL.",optional=true)
	public String HOSTNAME="http://localhost:8080/ngs";
    @Option(shortName= "MXL", doc="max number of items per output",optional=true)
	public int MAX_LINES=10000;
	
	
	private Map<String,File> bams=new HashMap<String,File>();
	private Map<String,File> vcfs=new HashMap<String,File>();
		
	
	@Override
	public List<String> getBams()
		{
		return new ArrayList<String>(bams.keySet());
		}

	@Override
	public List<String> getVCFs()
		{
		return new ArrayList<String>(vcfs.keySet());
		}

	@Override
	public Bam getBam(String filename, String chrom, int start, int end)
		{
		LOG.debug("query is "+filename+"\t"+chrom+":"+start+"-"+end);
		SAMRecordIterator iter=null;
		SAMFileReader r=null;
		long nLines=0L;
		
		try
			{
			File f=this.bams.get(filename);
			if(f==null || end<start) return null;
			Bam bam=new Bam();
			r=new  SAMFileReader(f);
			r.setValidationStringency(super.VALIDATION_STRINGENCY);
			iter=r.queryOverlapping(chrom, start, end);
			while(iter.hasNext())
				{
				++nLines;
				SAMRecord rec=iter.next();
				if(nLines<=MAX_LINES)
					{
					BamRecord br=new BamRecord();
					br.setName(rec.getReadName());
					bam.getRecord().add(br);
					}
				else
					{
					bam.setNextChrom(rec.getReferenceName());
					bam.setNextPos(rec.getAlignmentStart());
					break;
					}
				}
			return bam;
			}
		catch(Exception err)
			{
			LOG.info(err);
			return null;
			}
		finally
			{
			if(iter!=null) iter.close();
			if(r!=null) r.close();
			}
		}

	@Override
	public Vcf getVcf(String filename, String chrom, int start, int end)
		{
		LOG.debug("query is "+filename+"\t"+chrom+":"+start+"-"+end);
		TabixReader tbx=null;
		long nLines=0L;
		try
			{
			File f=this.vcfs.get(filename);
			if(f==null) return null;
			Vcf vcf=new Vcf();
			VCFCodec codec=new VCFCodec();
			tbx=new TabixReader(filename);
			TabixReader.Iterator r=tbx.query(chrom+":"+start+"-"+end);
			String line;
			while(r!=null && (line=r.next())!=null)
				{
				VariantContext ctx=codec.decode(line);
				++nLines;
				if(nLines<=this.MAX_LINES)
					{
					VcfVar var=new VcfVar();
					var.setChrom(ctx.getChr());
					var.setPos(ctx.getStart());
					var.setId(ctx.getID());
					var.setRef(ctx.getReference().getDisplayString());
					
					var.setInfos(new VcfVar.Infos());
					for(String key:ctx.getAttributes().keySet())
						{
						VcfVar.Infos.Info info=new VcfVar.Infos.Info();
						info.setKey(key);
						info.setValue(ctx.getAttributeAsString(key, ""));
						var.getInfos().getInfo().add(info);
						}
					for(Genotype g:ctx.getGenotypes())
						{
						VcfVar.Genotypes.Genotype e=new VcfVar.Genotypes.Genotype();
						VcfVar.Genotypes.Genotype.Property prop=new VcfVar.Genotypes.Genotype.Property();
						if(g.hasDP())
							{
							prop.setKey("dp");
							prop.setValue(String.valueOf(g.getDP()));
							}
						if(g.hasGQ())
							{
							prop.setKey("gq");
							prop.setValue(String.valueOf(g.getGQ()));
							}
						e.setSample(g.getSampleName());
						var.getGenotypes().getGenotype().add(e);
						}
					vcf.getVar().add(var);
					}
				else
					{
					vcf.setNextChrom(ctx.getChr());
					vcf.setNextPos(ctx.getStart());
					break;
					}
				}
			return vcf;
			}
		catch (Exception e)
			{
			LOG.error(e);
			return null;
			}
		finally
			{
			if(tbx!=null) tbx.close();
			}
		}
	private String hash(File f)
		{
		return "f"+String.valueOf(1+vcfs.size()+bams.size());
		}
	private void scan(File f)
		{
		LOG.debug("scanning "+f);
		if(f.isFile() && f.canRead())
			{
			String name=f.getName();
			if(name.endsWith(".bam"))
				{
				String name2=name.substring(0,name.length()-4);
				File bai=new File(f.getParent(),name+".bai");
				if(!bai.exists())
					{
					bai=new File(f.getParent(),name2+".bai");
					}
				if(bai.exists() && bai.lastModified()>=f.lastModified())
					{
					LOG.info("adding "+f);
					bams.put(hash(f),f);
					}
				}
			else if(f.getName().endsWith(".vcf.gz"))
				{
				File tbi=new File(f.getParent(),name+".tbi");
				if(tbi.exists() && tbi.lastModified()>=tbi.lastModified())
					{
					LOG.info("adding "+f);
					vcfs.put(hash(f),f);
					}
				}
			}
		else if(f.isDirectory())
			{
			for(File f2:f.listFiles())
				{
				scan(f2);
				}
			}
		}
	@Override
	protected int doWork()
		{
		try
			{
			for(File f:IN)
				{
				scan(f);
				}
			if(bams.isEmpty() && vcfs.isEmpty())
				{
				LOG.error("no bam or no vcf found");
				return -1;
				}
			LOG.info("publishing at "+HOSTNAME);
			Endpoint.publish(this.HOSTNAME, this);
			return 0;
			}
		catch (Exception e)
			{
			LOG.error(e);
			return -1;
			}
		finally
			{
			
			}
		}
	
	public static void main(String[] args)
		{
		NGSServiceImpl service=new NGSServiceImpl();
		service.instanceMainWithExit(args);
		}
	
	
	
	
	}
