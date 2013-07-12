package com.github.lindenb.jvarkit.tools.vcfgo;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

import javax.xml.stream.XMLStreamException;

import net.sf.picard.cmdline.Option;
import net.sf.picard.util.Log;

import org.broad.tribble.readers.LineReader;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFCodec;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLineCount;
import org.broadinstitute.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.AbstractVCFFilter;
import com.github.lindenb.jvarkit.util.go.GoTree;
import com.github.lindenb.jvarkit.util.picard.IOUtils;
import com.github.lindenb.jvarkit.util.vcf.predictions.Prediction;
import com.github.lindenb.jvarkit.util.vcf.predictions.PredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;

public class VcfGeneOntology extends AbstractVCFFilter
	{
	private GoTree goTree;
	private Map<String,Set<GoTree.Term>> name2go=new HashMap<String, Set<GoTree.Term>>();
	@Option(shortName="GOA_INPUT", doc="GOA file/URI.",optional=true)
	public String GOA="http://cvsweb.geneontology.org/cgi-bin/cvsweb.cgi/go/gene-associations/gene_association.goa_human.gz?rev=HEAD";
	@Option(shortName="GO_INPUT", doc="GOA file/URI.",optional=true)
	public String GO="http://archive.geneontology.org/latest-termdb/go_daily-termdb.rdf-xml.gz";

	private static final Log LOG=Log.getInstance(VcfGeneOntology.class);

	
	private VepPredictionParser vepPredictionParser;
	private SnpEffPredictionParser snpeffPredictionParser;
	
	private void readGO() throws IOException
		{
		LOG.info("read GO "+GO);
		try
			{
			goTree=GoTree.parse(GO);
			LOG.info("GO size:"+goTree.size());
			}
		catch(XMLStreamException err)
			{
			throw new IOException(err);
			}
		}
	
	private void readGOA() throws IOException
		{
		LOG.info("read GOA "+GOA);
		Pattern tab=Pattern.compile("[\t]");
		BufferedReader in=IOUtils.openURIForBufferedReading(GOA);
		String line;
		while((line=in.readLine())!=null)
			{
			if(line.isEmpty() || line.startsWith("!")) continue;
			String tokens[]=tab.split(line,6);
			if(tokens.length<6) continue;
			
			GoTree.Term term=goTree.getTermByAccession(tokens[4]);
			if(term==null)
				{
				
				continue;
				}
			
			
			Set<GoTree.Term> set=name2go.get(tokens[2]);
			if(set==null)
				{
				set=new HashSet<GoTree.Term>();
				name2go.put(tokens[2],set);
				}
			set.add(term);
			}
		in.close();
		LOG.info("GOA size:"+name2go.size());
		}
	private Set<String> getGeneNames(VariantContext ctx)
		{
		Set<String> names=new HashSet<String>();
		for(PredictionParser p:new PredictionParser[]{
				this.vepPredictionParser,	
				this.snpeffPredictionParser
			})
			{
			for(Prediction pred:p.getPredictions(ctx))
				{
				if(pred.getGeneName()==null || pred.getGeneName().isEmpty())
					{
					continue;
					}
				names.add(pred.getGeneName());
				}
			}
		
		return names;
		}
	@Override
	protected void doWork(LineReader in, VariantContextWriter w)
			throws IOException
		{
		readGO();
		readGOA();
		
		
		final String TAG="GOA";
		VCFCodec codeIn=new VCFCodec();		
		VCFHeader header=(VCFHeader)codeIn.readHeader(in);
		
		this.vepPredictionParser= new VepPredictionParser(header);
		this.snpeffPredictionParser= new SnpEffPredictionParser(header);
		
		
		VCFHeader h2=new VCFHeader(header.getMetaDataInInputOrder(),header.getSampleNamesInOrder());
		h2.addMetaDataLine(new VCFInfoHeaderLine(TAG,VCFHeaderLineCount.UNBOUNDED,VCFHeaderLineType.String,"GO terms from GO "+GO+" and GOA="+GOA));
		
		h2.samplesWereAlreadySorted();
		w.writeHeader(h2);
	
		String line;
		while((line=in.readLine())!=null)
			{
			VariantContext ctx=codeIn.decode(line);
			Set<String> geneNames=getGeneNames(ctx);
		
			
			List<String> atts=new ArrayList<String>();
			for(String GN:geneNames)
				{
				StringBuilder sb=new StringBuilder(GN);
				sb.append("|");
				Set<GoTree.Term> t2=this.name2go.get(GN);
				if(t2==null) continue;
				boolean first=true;
				for(GoTree.Term gt:t2)
					{
					if(!first) sb.append("&");
					sb.append(gt.getAcn());
					first=false;
					}
				atts.add(sb.toString());
				}
			if(atts.isEmpty())
				{
				w.add(ctx);
				continue;
				}
			VariantContextBuilder b=new VariantContextBuilder(ctx);

			b.attribute(TAG, atts);
			
			
			
			
			w.add(b.make());
			}
		}
	public static void main(String[] args)
		{
		new VcfGeneOntology().instanceMainWithExit(args);
		}
	}
