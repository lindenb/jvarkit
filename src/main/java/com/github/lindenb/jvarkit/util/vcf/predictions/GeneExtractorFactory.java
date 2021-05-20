/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/
package com.github.lindenb.jvarkit.util.vcf.predictions;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser.AnnPrediction;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParser.SnpEffPrediction;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser.VepPrediction;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;


public class GeneExtractorFactory {

public interface KeyAndGene
	extends Comparable<KeyAndGene>
	{	
	/** get the name of the entity (gene-id, transcript-id, etc... )*/
	public String getKey();
	/** name of the gene may be null , not used for comparaison */
	public String getGene();
	/** name of the extractor */
	public String getMethod();
	/** first key is method, second is gene */
	@Override
	public default int compareTo(final KeyAndGene o) {
		int i = this.getMethod().compareTo(o.getMethod());
		if(i!=0) return i;
		return this.getKey().compareTo(o.getKey());
		}
	}
	
public static class KeyAndGeneImpl implements KeyAndGene
	{
	final String key;
	final String gene;
	final String method;
	public KeyAndGeneImpl(final String key,final String gene,final String method) {
		this.key = key;
		this.gene = StringUtils.isBlank(gene)?".":gene;
		this.method = method;
		}
	@Override
	public String getKey() {
		return key;
		}
	@Override
	public String getGene() {
		return gene;
		}
	@Override
	public String getMethod() {
		return method;
		}
	@Override
	public int hashCode() {
		return this.key.hashCode()*31+this.method.hashCode();
		}
	@Override
	public boolean equals(Object obj) {
		if(obj==this) return true;
		if(obj==null || !(obj instanceof KeyAndGene)) return false;
		return this.compareTo(KeyAndGene.class.cast(obj))==0;
		}
	
	@Override
	public String toString() {
		return this.getKey();
		}
	}

	
	
public interface GeneExtractor extends Function<VariantContext, Map<KeyAndGene,Set<String>>> {
	/** get name of the tag in the INFO column */
	public String getInfoTag();
	/** get name for this extraction */
	public String getName();
	}

private abstract class  AbstractGeneExtractorImpl implements GeneExtractor {
	private String extractorName;
	AbstractGeneExtractorImpl(final String extractorName) {
		this.extractorName = extractorName;
		}
	@Override
	public String getName() {
		return this.extractorName;
		}
	@Override
	public int hashCode() {
		return this.extractorName.hashCode();
		}
	@Override
	public boolean equals(Object obj) {
		if(obj==this) return true;
		if(obj==null || !(obj instanceof GeneExtractor)) return false;
		return this.extractorName.equals(GeneExtractor.class.cast(obj).getName());
		}
	
	@Override
	public String toString() {
		return this.getName();
		}
	}

private class VepGeneExtractor   extends AbstractGeneExtractorImpl {
	private VepPredictionParser parser = null;
	private final Function<VepPrediction, String> pred2gene;
	VepGeneExtractor(final VepPredictionParser parser,final String name,final Function<VepPrediction, String> pred2gene)
		{
		super(name);
		this.parser= parser;
		this.pred2gene = pred2gene;
		}

	@Override
	public String getInfoTag() {
		return this.parser.getTag();
		}
	
	@Override
	public Map<KeyAndGene,Set<String>> apply(final VariantContext ctx ) {
		final Map<KeyAndGene,Set<String>> gene2values = new HashMap<>();
		for(final VepPrediction pred:this.parser.getPredictions(ctx)){
			final String geneName=this.pred2gene.apply(pred);
			if(StringUtils.isBlank(geneName)) continue;
			final KeyAndGene keyAndGene=new KeyAndGeneImpl(geneName,pred.getGeneName(),this.getName());
			Set<String> values = gene2values.get(keyAndGene);
			if(values==null)  {
				values = new LinkedHashSet<>();
				gene2values.put(keyAndGene,values);
				}
			values.add(pred.getOriginalAttributeAsString());
			}
		return gene2values;
		}
	}


/** bcftools csq extractor */
private class BcftoolsCsqExtractor  extends AbstractGeneExtractorImpl {
	private final BcfToolsPredictionParser parser;
	private final Function<BcfToolsPredictionParser.BcfToolsPrediction, String> pred2gene;

	BcftoolsCsqExtractor(final BcfToolsPredictionParser parser,final String name,final Function<BcfToolsPredictionParser.BcfToolsPrediction, String> pred2gene)
		{
		super(name);
		this.parser= parser;
		this.pred2gene = pred2gene;
		}
	@Override
	public String getInfoTag() {
		return parser.getTag();
		}
	
	@Override
	public Map<KeyAndGene,Set<String>> apply(final VariantContext ctx) {
		final Map<KeyAndGene,Set<String>> gene2values = new HashMap<>();
		for(final BcfToolsPredictionParser.BcfToolsPrediction pred:this.parser.getPredictions(ctx)){
			if(pred.isIntergenicRegion()) continue;
			final String geneName=this.pred2gene.apply(pred);
			if(StringUtils.isBlank(geneName)) continue;
			final KeyAndGene keyAndGene = new KeyAndGeneImpl(geneName, pred.getGeneName(),this.getName());
			Set<String> values = gene2values.get(keyAndGene);
			if(values==null)  {
				values = new LinkedHashSet<>();
				gene2values.put(keyAndGene,values);
				}
			values.add(pred.getOriginalAttributeAsString());
			}
		return gene2values;
		}
	}

/** SNPEFF/ANN extractor */
private class AnnGeneExtractor   extends AbstractGeneExtractorImpl {
	private final AnnPredictionParser parser;
	private final Function<AnnPrediction, String> pred2gene;

	AnnGeneExtractor(final AnnPredictionParser parser,final String name,final Function<AnnPrediction, String> pred2gene)
		{
		super(name);
		this.parser= parser;
		this.pred2gene = pred2gene;
		}
	@Override
	public String getInfoTag() {
		return parser.getTag();
		}
	
	@Override
	public Map<KeyAndGene,Set<String>> apply(final VariantContext ctx) {
		final Map<KeyAndGene,Set<String>> gene2values = new HashMap<>();
		for(final AnnPrediction pred:this.parser.getPredictions(ctx)){
			if(pred.isIntergenicRegion()) continue;
			final String geneName=this.pred2gene.apply(pred);
			if(StringUtils.isBlank(geneName)) continue;
			final KeyAndGene keyAndGene = new KeyAndGeneImpl(geneName, pred.getGeneName(),this.getName());
			Set<String> values = gene2values.get(keyAndGene);
			if(values==null)  {
				values = new LinkedHashSet<>();
				gene2values.put(keyAndGene,values);
				}
			values.add(pred.getOriginalAttributeAsString());
			}
		return gene2values;
		}
	}

private class SnpEffGeneExtractor   extends AbstractGeneExtractorImpl {
	private final SnpEffPredictionParser parser;
	private final Function<SnpEffPrediction, String> pred2gene;

	SnpEffGeneExtractor(final SnpEffPredictionParser parser,final String name,final Function<SnpEffPrediction, String> pred2gene)
		{
		super(name);
		this.parser= parser;
		this.pred2gene = pred2gene;
		}
	@Override
	public String getInfoTag() {
		return parser.getTag();
		}
	
	@Override
	public Map<KeyAndGene,Set<String>> apply(final VariantContext ctx) {
		final Map<KeyAndGene,Set<String>> gene2values = new HashMap<>();
		for(final SnpEffPrediction pred:this.parser.getPredictions(ctx)){
			//if(pred.isIntergenicRegion()) continue;
			final String geneName=this.pred2gene.apply(pred);
			if(StringUtils.isBlank(geneName)) continue;
			final KeyAndGene keyAndGene = new KeyAndGeneImpl(geneName, pred.getGeneName(),this.getName());
			Set<String> values = gene2values.get(keyAndGene);
			if(values==null)  {
				values = new LinkedHashSet<>();
				gene2values.put(keyAndGene,values);
				}
			values.add(pred.getOriginalAttributeAsString());
			}
		return gene2values;
		}
	}

private class SmooveExtractor extends AbstractGeneExtractorImpl {
	private final SmooveGenesParser parser;
	SmooveExtractor(final String name,final SmooveGenesParser smooveGenesParser)
		{
		super(name);
		this.parser = smooveGenesParser;
		}
	@Override
	public Map<KeyAndGene, Set<String>> apply(final VariantContext vc) {
		final Map<KeyAndGene,Set<String>> gene2values = new HashMap<>();
		for(final SmooveGenesParser.Prediction pred:this.parser.parse(vc)){
			//if(pred.isIntergenicRegion()) continue;
			final String geneName= pred.getGeneName();
			if(StringUtils.isBlank(geneName)) continue;
			final KeyAndGene keyAndGene = new KeyAndGeneImpl(geneName, geneName,this.getName());
			Set<String> values = gene2values.get(keyAndGene);
			if(values==null)  {
				values = new LinkedHashSet<>();
				gene2values.put(keyAndGene,values);
				}
			values.add(pred.getOriginalAttributeAsString());
			}
		return gene2values;
		}
	@Override
	public String getInfoTag()
		{
		return this.parser.getTag();
		}
	}

private final List<GeneExtractor> extractors =  new ArrayList<>();
/* WARNING keep that order: see constuctor */
private static List<String> AVAILABLE_EXTRACTORS_NAMES = Collections.unmodifiableList(Arrays.asList(
		"ANN/GeneId","ANN/FeatureId","ANN/GeneName",// 0 & 1 & 2
		"VEP/GeneId","VEP/Ensp","VEP/Feature",// 3 & 4 & 5
		"EFF/Gene","EFF/Transcript",// 6 & 7
		"BCSQ/gene","BCSQ/transcript",//8 & 9
		"SMOOVE" //10
		))
		;

public static final String OPT_DESC = "Gene Extractors Name. Space/semicolon/Comma separated";

public GeneExtractorFactory(final VCFHeader header) {
	
	/* WARNING keep that order: see AVAILABLE_EXTRACTORS_NAMES */

	final AnnPredictionParser annparser =new AnnPredictionParserFactory().header(header).get();
	
	extractors.add( new AnnGeneExtractor(annparser,AVAILABLE_EXTRACTORS_NAMES.get(0), P->P.getGeneId()));
	extractors.add( new AnnGeneExtractor(annparser,AVAILABLE_EXTRACTORS_NAMES.get(1), P->P.getFeatureId()));
	extractors.add( new AnnGeneExtractor(annparser,AVAILABLE_EXTRACTORS_NAMES.get(2), P->P.getGeneName()));
	
	final VepPredictionParser vepparser = new VepPredictionParserFactory().header(header).get();
	
	extractors.add( new VepGeneExtractor(vepparser,AVAILABLE_EXTRACTORS_NAMES.get(3), P->P.getEnsemblGene()));
	extractors.add( new VepGeneExtractor(vepparser,AVAILABLE_EXTRACTORS_NAMES.get(4), P->P.getENSP()));
	extractors.add( new VepGeneExtractor(vepparser,AVAILABLE_EXTRACTORS_NAMES.get(5), P->P.getFeature()));
	

	final SnpEffPredictionParser effparser = new SnpEffPredictionParser(header);
	extractors.add( new SnpEffGeneExtractor(effparser,AVAILABLE_EXTRACTORS_NAMES.get(6), P->P.getGeneName()));
	extractors.add( new SnpEffGeneExtractor(effparser,AVAILABLE_EXTRACTORS_NAMES.get(7), P->P.getEnsemblTranscript()));

	final BcfToolsPredictionParser csqParser = new BcfToolsPredictionParser(header);
	extractors.add( new BcftoolsCsqExtractor(csqParser,AVAILABLE_EXTRACTORS_NAMES.get(8), P->P.getGeneName()));
	extractors.add( new BcftoolsCsqExtractor(csqParser,AVAILABLE_EXTRACTORS_NAMES.get(9), P->P.getTranscript()));

	
	final SmooveGenesParser smooveGenesParser = new SmooveGenesParser(header);
	extractors.add( new SmooveExtractor(AVAILABLE_EXTRACTORS_NAMES.get(10), smooveGenesParser));
	}

/** return a list of all the available extractors' names */
public static List<String> getExtractorNames() {
	return AVAILABLE_EXTRACTORS_NAMES;
	}

/** get all available extractors */
public List<GeneExtractor> getAllExtractors() {
	return Collections.unmodifiableList(this.extractors);
}


public List<GeneExtractor> parse(final String arg)
	{
	if(StringUtils.isBlank(arg)) return Collections.emptyList();
	final List<GeneExtractor> L = new ArrayList<>();

	for(final String s:arg.split("[\\s,;]+")) {
		if(StringUtils.isBlank(s)) continue;
		final Optional<GeneExtractor> ex = this.getAllExtractors().
				stream().
				filter(E->E.getName().equals(s)).
				findFirst();
		if(!ex.isPresent()) {
			throw new IllegalArgumentException("Cannot find gene extractor \""+s+"\" in \""+arg+"\". Available are: " +
						this.getAllExtractors().stream().map(E->E.getName()).collect(Collectors.joining(" ")));
			}
		L.add(ex.get());
		}
	return L;
	}

}
