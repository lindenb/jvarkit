package com.github.lindenb.jvarkit.tools.vcfgo;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

import javax.xml.stream.XMLStreamException;

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalTreeMap;
import net.sf.picard.util.Log;

import org.broad.tribble.readers.LineReader;
import org.broadinstitute.variant.variantcontext.VariantContext;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.Function;
import com.github.lindenb.jvarkit.util.biomart.BiomartQuery;
import com.github.lindenb.jvarkit.util.go.GoTree;
import com.github.lindenb.jvarkit.util.picard.IntervalTreeMapFactory;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter;

public abstract class AbstractVcfGeneOntology extends AbstractVCFFilter
	{
	protected GoTree goTree;
	protected Map<String,Set<GoTree.Term>> name2go=new HashMap<String, Set<GoTree.Term>>();
	protected IntervalTreeMap<String> hgncmap;
	
	@Option(shortName="GOA_INPUT", doc="GOA file/URI.",optional=true)
	public String GOA="http://cvsweb.geneontology.org/cgi-bin/cvsweb.cgi/go/gene-associations/gene_association.goa_human.gz?rev=HEAD";
	@Option(shortName="GO_INPUT", doc="GOA file/URI.",optional=true)
	public String GO="http://archive.geneontology.org/latest-termdb/go_daily-termdb.rdf-xml.gz";
	@Option(shortName=StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc="Reference (used to reduce the number of mapped genes)",optional=true)
	public File REF=null;

	
	
	protected static final Log LOG=Log.getInstance(AbstractVcfGeneOntology.class);

	
	
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
	
	protected void loadBiomartHGNC() throws IOException
		{
		BiomartQuery q=new BiomartQuery();
		q.setDataSetName("hsapiens_gene_ensembl");
		q.setAttributes(
				"chromosome_name",
				"start_position",
				"end_position",
				"hgnc_symbol"
				);
		q.setUniqRows(true);
		
		
		IntervalTreeMapFactory<String> itf=new IntervalTreeMapFactory<String>();
		if(REF!=null)
			{
			itf.setSamSequenceDictionary(new IndexedFastaSequenceFile(REF).getSequenceDictionary());
			}
		itf.setValueFunction(new Function<String[], String>()
			{
			@Override
			public String apply(String[] param)
				{
				if(param.length<4 || param[3].isEmpty()) return null;
				if(!name2go.containsKey(param[3])) return null;
				return param[3];
				}
			});
		LOG.info("invoking biomart "+q);
		LineReader r=q.execute();
		
		this.hgncmap=itf.createIntervalMap(r);
		r.close();
		}
	
	
	protected Set<String> getGeneNames(VariantContext ctx)
		{
		return new HashSet<String>(this.hgncmap.getOverlapping(new Interval(
				ctx.getChr(),
				ctx.getStart(),
				ctx.getEnd()
				)));
		}
	
	}
