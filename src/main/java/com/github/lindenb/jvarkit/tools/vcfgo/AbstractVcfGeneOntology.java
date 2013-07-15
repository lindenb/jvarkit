package com.github.lindenb.jvarkit.tools.vcfgo;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

import javax.xml.stream.XMLStreamException;

import net.sf.picard.cmdline.Option;
import net.sf.picard.util.Log;

import org.broadinstitute.variant.variantcontext.VariantContext;

import com.github.lindenb.jvarkit.util.AbstractVCFFilter;
import com.github.lindenb.jvarkit.util.go.GoTree;
import com.github.lindenb.jvarkit.util.picard.IOUtils;
import com.github.lindenb.jvarkit.util.vcf.predictions.Prediction;
import com.github.lindenb.jvarkit.util.vcf.predictions.PredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;

public abstract class AbstractVcfGeneOntology extends AbstractVCFFilter
	{
	protected GoTree goTree;
	protected Map<String,Set<GoTree.Term>> name2go=new HashMap<String, Set<GoTree.Term>>();
	@Option(shortName="GOA_INPUT", doc="GOA file/URI.",optional=true)
	public String GOA="http://cvsweb.geneontology.org/cgi-bin/cvsweb.cgi/go/gene-associations/gene_association.goa_human.gz?rev=HEAD";
	@Option(shortName="GO_INPUT", doc="GOA file/URI.",optional=true)
	public String GO="http://archive.geneontology.org/latest-termdb/go_daily-termdb.rdf-xml.gz";

	protected static final Log LOG=Log.getInstance(AbstractVcfGeneOntology.class);

	
	protected VepPredictionParser vepPredictionParser;
	protected SnpEffPredictionParser snpeffPredictionParser;
	
	protected void readGO() throws IOException
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
	
	protected void readGOA() throws IOException
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
	
	protected PredictionParser[] getPredictionParsers()
		{
		return new PredictionParser[]{
				this.vepPredictionParser,	
				this.snpeffPredictionParser
			};
		}
	
	protected Set<String> getGeneNames(VariantContext ctx)
		{
		Set<String> names=new HashSet<String>();
		for(PredictionParser p:getPredictionParsers())
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
	
	}
