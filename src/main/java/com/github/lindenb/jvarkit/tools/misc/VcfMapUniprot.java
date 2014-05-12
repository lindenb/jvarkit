package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.net.URLEncoder;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import htsjdk.samtools.util.CloserUtil;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter2;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.jvarkit.util.vcf.predictions.MyPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBElement;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Unmarshaller;
import javax.xml.transform.stream.StreamSource;

import org.uniprot.Entry;
import org.uniprot.FeatureType;
import org.uniprot.LocationType;
import org.uniprot.Uniprot;


public class VcfMapUniprot extends AbstractVCFFilter2
	{
	private int taxonid=9606;
	private Unmarshaller uniprotUnmarshaller=null;
	
	@SuppressWarnings("unused")
	private org.uniprot.ObjectFactory _fool_javac=null;
	
	private VcfMapUniprot()
		{
		}
	private long intervalMillisec=1000;
	private long last_fetch=0L;
	private Uniprot fetch(String gene) throws IOException,JAXBException
		{
		Exception lastErr=null;
		long diff=System.currentTimeMillis()-last_fetch;
		if(diff< intervalMillisec && intervalMillisec>0)
			{
			info("Waiting : "+diff+" ms");
			try
				{
				Thread.sleep(intervalMillisec-diff);
				}
			catch(Exception err)
				{
				
				}
			}
		last_fetch=System.currentTimeMillis();
		
		URL url=new URL("http://www.uniprot.org/uniprot/?query=organism:"+taxonid+"+"+
				URLEncoder.encode(gene,"UTF-8")
				+"&format=xml");
		
		for(int i=0;i< 20;++i)
			{
			InputStream in=null;
			try
				{
				last_fetch=System.currentTimeMillis();
				
				info(url);
				in=url.openStream();
				StreamSource is=new StreamSource(in);
				JAXBElement<Uniprot> jaxbElement=uniprotUnmarshaller.unmarshal(is, Uniprot.class);
				
				return jaxbElement.getValue();
				}
			catch(javax.xml.bind.UnmarshalException err)
				{
				return null;
				}
			catch(Exception err)
				{
				System.err.println("Failed for "+url);
				lastErr=err;
				}
			finally
				{
				CloserUtil.close(in);
				}
			}
		if(lastErr==null) lastErr=new IllegalStateException(); 
		error(lastErr);
		return null;
		}
	private CacheMap<String, List<Entry>> ensemblToUniprot=new CacheMap<String,  List<Entry>>(50);
	
	
	private List<Entry> getUniprotEntriesFromEnsembl(String ensemblId) throws IOException,JAXBException
		{
		if(ensemblId==null || ensemblId.trim().isEmpty()) return Collections.emptyList();
		if(this.ensemblToUniprot.containsKey(ensemblId))
			{
			return this.ensemblToUniprot.get(ensemblId);
			}
		Uniprot uniprot=fetch(ensemblId);
		if(uniprot!=null)
			{
			this.ensemblToUniprot.put(ensemblId,uniprot.getEntry());
			return uniprot.getEntry();
			}
		else
			{
			this.ensemblToUniprot.put(ensemblId, new ArrayList<Entry>());
			return Collections.emptyList();
			}
		
		}
	
	private String escape(String s)
		{
		if(s==null) return "";
		s=s.trim().replaceAll("[ ,=\\|;:\"\']+", "_");
		if(s.endsWith(".")) s=s.substring(0,s.length()-1).trim();
		return s;
		}
	
	private void map(List<Entry> uniprotEntries,String ens,int aapos,Set<String> annotations)
		{
		for(Entry entry:uniprotEntries)
			{
			String entryName=entry.getAccession().isEmpty()?"":entry.getAccession().get(0);
			for(FeatureType feat:entry.getFeature())
				{
				if(feat.getType()==null || feat.getType().isEmpty()) continue;
				LocationType locType=feat.getLocation();
				if(locType==null) continue;
				int pepStart,pepEnd;
				if(locType.getPosition()!=null && locType.getPosition().getPosition()!=null)
					{
					pepStart=locType.getPosition().getPosition().intValue();
					pepEnd=pepStart;
					}
				else if(locType.getBegin()!=null &&
						locType.getEnd()!=null &&
						locType.getBegin().getPosition()!=null &&
						locType.getEnd().getPosition()!=null )
					{
					pepStart=locType.getBegin().getPosition().intValue();
					pepEnd=locType.getEnd().getPosition().intValue();
					}
				else
					{
					continue;
					}
				if(pepStart<= aapos && aapos <=pepEnd)
					{
					annotations.add(
							entryName+"|"+
							pepStart+"|"+
							pepEnd+"|"+
							escape(feat.getDescription())+"|"+
							escape(feat.getType())+"|"+
							escape(ens)+"|"+
							escape(String.valueOf(aapos))
							);
					}
				}
	
			}
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/VcfMapUniprot";
		}
	
	@Override
	public String getProgramDescription() {
		return "Map uniprot features on VCF annoted with VEP or SNPEff";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println("-t (taxon-id) default:"+this.taxonid);
		out.println("-w (milliSec) wait between two uniprot queries : -1 = no-wait. default:"+this.intervalMillisec);
		super.printOptions(out);
		}

	
	@Override
	protected void doWork(VcfIterator in, VariantContextWriter out)
			throws IOException
		{
		final String TAG="UNIPROT";
		VCFHeader header=in.getHeader();
		MyPredictionParser myPredictionParser=new MyPredictionParser(header);
		SnpEffPredictionParser snpEffPredictionParser=new SnpEffPredictionParser(header);
		VepPredictionParser vepPredictionParser=new VepPredictionParser(header);
		header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
		header.addMetaDataLine(new VCFInfoHeaderLine(TAG,VCFHeaderLineCount.UNBOUNDED,VCFHeaderLineType.String,"Uniprot Features. Format: ACN|START|END|DESC|TYPE|ENSGENE|POSAA"));
		out.writeHeader(header);
		try
			{
			while(in.hasNext())
				{
				VariantContext ctx=in.next();
				VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				Set<String> annotations=new HashSet<String>();
				for(SnpEffPredictionParser.SnpEffPrediction pred:snpEffPredictionParser.getPredictions(ctx))
					{
					Integer aapos=pred.getAminoAcidPosition();
					if(aapos==null) continue;
					String ens=pred.getEnsemblTranscript();
					if(ens==null) continue;
					map(getUniprotEntriesFromEnsembl(ens),ens, aapos, annotations);
					}
				for(VepPredictionParser.VepPrediction pred:vepPredictionParser.getPredictions(ctx))
					{
					Integer aapos=pred.getAminoAcidPosition();
					if(aapos==null) continue;
					String ens=pred.getEnsemblTranscript();
					if(ens==null) continue;
					map(getUniprotEntriesFromEnsembl(ens),ens,  aapos, annotations);
					}
				for(MyPredictionParser.MyPrediction pred:myPredictionParser.getPredictions(ctx))
					{
					Integer aapos=pred.getAminoAcidPosition();
					if(aapos==null) continue;
					String ens=pred.getEnsemblTranscript();
					if(ens==null) continue;
					map(getUniprotEntriesFromEnsembl(ens),ens,  aapos, annotations);
					}
				
				
				if(!annotations.isEmpty())
					{
					vcb.attribute(TAG,new ArrayList<String>(annotations));
					}
				
				
				out.add(vcb.make());
				if(super.outCheckError()) break;
				}
			}
		catch(JAXBException err)
			{
			throw new IOException(err);
			}
		}
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"t:w:"))!=-1)
			{
			switch(c)
				{
				case 't': this.taxonid=Integer.parseInt(opt.getOptArg());break;
				case 'w': this.intervalMillisec=Long.parseLong(opt.getOptArg());break;
				default:
					{
					switch(handleOtherOptions(c, opt, null))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		
		try
			{
			JAXBContext jc = JAXBContext.newInstance("org.uniprot");
			this.uniprotUnmarshaller=jc.createUnmarshaller();
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		return super.doWork(opt.getOptInd(), args);
		}
	public static void main(String[] args) {
		new VcfMapUniprot().instanceMainWithExit(args);
		}
	
	
	@SuppressWarnings("serial")
	private static class CacheMap<K,V>
		extends LinkedHashMap<K, V>
		{
		private int capacity=10;
			
		public CacheMap(int capacity)
			{
			this.capacity=Math.max(capacity,1);
			}
		
		public int getCapacity()
			{
			return capacity;
			}
		
		@Override
		public V put(K key, V value) {
			V old= super.put(key, value);
			if(size()> getCapacity() )
				{
				K first=keySet().iterator().next();
				remove(first);
				}
			return old;
			}
		@Override
		public void putAll(Map<? extends K, ? extends V> m) {
			for(K k:m.keySet())
				{
				this.put(k,m.get(k));
				}
			}
		}

	
}
