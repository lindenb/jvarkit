package com.github.lindenb.jvarkit.tools.plotrare;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import org.w3c.dom.DocumentFragment;
import org.w3c.dom.Element;

import com.beust.jcommander.DynamicParameter;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.css.Colors;
import com.github.lindenb.jvarkit.gatk.GATKConstants;
import com.github.lindenb.jvarkit.gtf.GtfTrack;
import com.github.lindenb.jvarkit.html.HTMLDocument;
import com.github.lindenb.jvarkit.lang.AttributeMap;
import com.github.lindenb.jvarkit.lang.HasName;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.net.Hyperlink;
import com.github.lindenb.jvarkit.pedigree.CasesControls;
import com.github.lindenb.jvarkit.samtools.util.LocatableUtils;
import com.github.lindenb.jvarkit.svg.GenomeSVGDocument;
import com.github.lindenb.jvarkit.ucsc.CytobandTrack;
import com.github.lindenb.jvarkit.util.Maps;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;

public class PlotRareVariants extends Launcher {
	private static final Logger LOG = Logger.of(PlotRareVariants.class);
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	protected Path outputFile=null;
	@Parameter(names={"-B","--bed"},description="BED file for the regions to be displayed",required=true)
	protected Path bedFile=null;
	@Parameter(names={"-gtf","--gtf"},description="GTF file for plotting genes",required=false)
	protected Path gtfFile=null;
	@DynamicParameter(names={"-D"},description="undocumented properties")
	protected Map<String,String> _properties = new HashMap<>();
	@ParametersDelegate
	protected CasesControls casesControls = new CasesControls();
	
	private static class NamedFilter implements HasName,Predicate<VariantContext> {
		private final String name;
		private final Predicate<VariantContext> pred;
		NamedFilter(final String name, final Predicate<VariantContext> pred) {
			this.name = name;
			this.pred = pred;
			}
		@Override
		public String getName() {
			return name;
			}
		@Override
		public boolean test(VariantContext t) {
			return pred.test(t);
			}
		}

	@Override
	public int doWork(List<String> args) {
		try {
			final HTMLDocument htmlDoc = new HTMLDocument();
			final Set<String> observableInfoFields = new HashSet<>(
					Arrays.asList(VCFConstants.ALLELE_COUNT_KEY,VCFConstants.ALLELE_FREQUENCY_KEY,VCFConstants.ALLELE_NUMBER_KEY,
							GATKConstants.FS_KEY,GATKConstants.MQ_KEY,GATKConstants.SOR_KEY,GATKConstants.ReadPosRankSum_KEY,GATKConstants.QD_KEY,GATKConstants.MQRankSum_KEY
							));
			casesControls.load();
			casesControls.checkHaveCasesControls();
			LOG.info("open vcf");
			try(VCFReader vcfR = VCFReaderFactory.makeDefault().open(oneAndOnlyOneFile(args))) {
				final VCFHeader header= vcfR.getHeader();
				casesControls.retain(header);
				casesControls.checkHaveCasesControls();
				
				final SAMSequenceDictionary dict= SequenceDictionaryUtils.extractRequired(header);
				final Hyperlink hyperlink = Hyperlink.compile(dict);
				List<BedLine> intervals =  new ArrayList<>();
				
				try(BedLineReader blr= new BedLineReader(this.bedFile)) {
					blr.setContigNameConverter(ContigNameConverter.fromOneDictionary(dict));
					blr.stream().forEach(R->intervals.add(R));
					}
				
				
				final List<VariantContext> variants = new ArrayList<>();

				for(Locatable loc: intervals) {
					variants.addAll(vcfR.query(loc).
							stream().
							filter(V->V.getGenotypes().stream().anyMatch(G->G.hasAltAllele() && casesControls.contains(G))).
							collect(Collectors.toList())
							);
					}

				
				Collections.sort(variants, ContigDictComparator.createLocatableComparator(dict));
				
				final List<String> all_samples = new ArrayList<>();
				all_samples.addAll(casesControls.getCases());
				all_samples.addAll(casesControls.getControls());
				Collections.sort(all_samples,(SN1,SN2)->Long.compare(
						variants.stream().map(V->V.getGenotype(SN2)).filter(G->G.hasAltAllele()).count(),
						variants.stream().map(V->V.getGenotype(SN1)).filter(G->G.hasAltAllele()).count()
						));

				{
				int i=0;
				
				while(i+1< variants.size()) {
					VariantContext ctx1 = variants.get(i);
					VariantContext ctx2 = variants.get(i+1);
					if(ctx1.contigsMatch(ctx2) &&
							ctx1.getStart()==ctx2.getStart() && 
							ctx1.getEnd()==ctx2.getEnd() && 
							ctx1.hasSameAllelesAs(ctx2)) {
						variants.remove(i+1);
						}
					else
						{
						i++;
						}
					}
				}
				LOG.info("N-variants "+variants.size());

				htmlDoc.bodyElement.appendChild(htmlDoc.h1("Rare Variants"));
				final Element divForSvg = htmlDoc.div(null,Maps.of("text-align", "center"));
				htmlDoc.bodyElement.appendChild(htmlDoc.h2("Variants Context"));
				htmlDoc.bodyElement.appendChild(divForSvg);
				
				/* TABLE OF VARIANTS **************************************************************************************/
				htmlDoc.bodyElement.appendChild(htmlDoc.h2("Variants Table"));
				HTMLDocument.Table table= htmlDoc.createTable(Arrays.asList("Variant","ID","REF","ALT","INFO",
						"CASES (N="+StringUtils.niceInt(this.casesControls.getCases().size())+")",
						"CONTROLS (N="+StringUtils.niceInt(this.casesControls.getControls().size())+")"
						));
				htmlDoc.bodyElement.appendChild(table.table);
				for(int i=0;i< variants.size();i++) {
					final VariantContext ctx = variants.get(i);
					table.set(i, 0,htmlDoc.anchor(
							ctx.getContig()+":"+StringUtils.niceInt(ctx.getStart()),
							hyperlink.apply(variants.get(i)).orElse("")
							));
					
					if(ctx.hasID()) {
						table.set(i, 1,Arrays.stream(ctx.getID().split(VCFConstants.ID_FIELD_SEPARATOR)).map(ID->{
							if(ID.startsWith("rs")) {
								return htmlDoc.anchor( ID,"https://www.ncbi.nlm.nih.gov/snp/"+ID);
								}
							else
								{
								return htmlDoc.text(ID);
						}
							}).collect(htmlDoc.joining("; ")));
						}
					table.set(i, 2,ctx.getReference().getDisplayString());
					table.set(i, 3,ctx.getAlternateAlleles().stream().map(A->A.getDisplayString()).collect(Collectors.joining(",")));
					
					table.set(i, 4,
							ctx.getAttributes().entrySet().
							stream().
							filter(KV->observableInfoFields.contains(KV.getKey())).
							sorted((A,B)->A.getKey().compareTo(B.getKey())).
							map(KV->{
								final DocumentFragment df = htmlDoc.fragment();
								df.appendChild(htmlDoc.bold(KV.getKey()));
								df.appendChild(htmlDoc.text(":"));
								df.appendChild(htmlDoc.text(KV.getValue()));
								return df;
								}).
							collect(htmlDoc.joining("; "))
							);

					
					for(int side=0;side<2;++side) {
						final List<Genotype> gL=ctx.getGenotypes(this.casesControls.get(side)).stream().filter(G->G.hasAltAllele()).collect(Collectors.toList());
						if(gL.isEmpty()) continue;
						table.set(i, 5+side, gL.stream().map(G->{
							final DocumentFragment df = htmlDoc.fragment();
							df.appendChild(htmlDoc.bold(G.getSampleName()));
							df.appendChild(htmlDoc.text(":"+G.getType().name()));
							if(G.hasDP()) {
								df.appendChild(htmlDoc.text(":"));
								df.appendChild(htmlDoc.bold("DP"));
								df.appendChild(htmlDoc.text("="));
								df.appendChild(htmlDoc.text(G.getDP()));
								}
							if(G.hasGQ()) {
								df.appendChild(htmlDoc.text(":"));
								df.appendChild(htmlDoc.bold("GQ"));
								df.appendChild(htmlDoc.text("="));
								df.appendChild(htmlDoc.text(G.getGQ()));
								}
							if(G.hasAD()) {
								df.appendChild(htmlDoc.text(":"));
								df.appendChild(htmlDoc.bold("AD"));
								df.appendChild(htmlDoc.text("="));
								df.appendChild(htmlDoc.text(Arrays.stream(G.getAD()).mapToObj(Integer::toString).collect(Collectors.joining(","))));
								}
							return df;
							}).collect(htmlDoc.joining("; ")));
						}
					}
				
				/* TABLE OF Samples **************************************************************************************/
				htmlDoc.bodyElement.appendChild(htmlDoc.h2("Samples"));
				table= htmlDoc.createTable(Arrays.asList("Sample","Type","Count","Variants"));
				htmlDoc.bodyElement.appendChild(table.table);
				for(int i=0;i< all_samples.size();++i) {
					final String sample = all_samples.get(i);
					final List<VariantContext> subL = variants.stream().filter(V->V.getGenotype(sample).hasAltAllele()).collect(Collectors.toList());
					table.set(i, 0,htmlDoc.span(sample,Maps.of("style","color:"+(casesControls.isCase(sample)?"red":"blue"))));
					table.set(i, 1,htmlDoc.span(casesControls.isCase(sample)?"CASE":"CONTROL",Maps.of("style","color:"+(casesControls.isCase(sample)?"red":"blue"))));
					table.set(i, 2,StringUtils.niceInt(subL.size()));
					table.set(i, 3,subL.stream().map(V->V.getContig()+":"+V.getStart()).collect(Collectors.joining(" ")));
					}
				/* **************************************************************************************/
				
				final List<Interval> extendedIntervals = 	intervals.stream().
					map(R->LocatableUtils.slopByFraction(R, 1.1, dict)).
					collect(LocatableUtils.mergeIntervals()).
					stream().
					map(R->new Interval(R)).
					collect(Collectors.toList())
					;
				
				
				
				final GenomeSVGDocument svgDoc = new GenomeSVGDocument(
						dict,
						extendedIntervals,
						AttributeMap.wrap(this._properties)
						);
				
				final CytobandTrack cytobandTrack = new CytobandTrack();
				cytobandTrack.setDocument(svgDoc);
				cytobandTrack.paint();
				
				final GtfTrack gtfTrack = new GtfTrack();
				gtfTrack.setGtfPath(this.gtfFile);
				gtfTrack.setDocument(svgDoc);
				gtfTrack.paint();

				final int fontSize=10;
				double y = svgDoc.lastY;
				y+=5;
				final Element g  = svgDoc.group();
				svgDoc.rootElement.appendChild(g);
				
				for(int side=0;side< 3;++side) {
					int featureHeight=20;
					final Set<String> samples = (side==0?this.casesControls.getCases():side==1?this.casesControls.getAll():this.casesControls.getControls());
					
					
					LOG.info("ploting variant side="+side);
					final Element g1 = svgDoc.group();
					g.appendChild(g1);
					g1.appendChild(svgDoc.comment("BEGIN side="+side));
					int countVariants =0;
					
					for(final VariantContext ctx:variants) {
						final Set<String> sampleWithAlt = ctx.getGenotypes().
									stream().
									filter(G->samples.contains(G.getSampleName())).
									filter(G->G.hasAltAllele()).
									map(G->G.getSampleName()).
									collect(Collectors.toSet());
									;
							
						//side == 1 must have both case and controls
						if(side==1 && !(
								ctx.getGenotypes().stream().filter(G->G.hasAltAllele()).anyMatch(G->casesControls.isCase(G)) &&
								ctx.getGenotypes().stream().filter(G->G.hasAltAllele()).anyMatch(G->casesControls.isControl(G)) 
								)) {
							sampleWithAlt.clear();//don't continue, we increase y later
							}
							
									
						if(!sampleWithAlt.isEmpty()) {
							countVariants++;
							for(GenomeSVGDocument.IntervalInfo ii: svgDoc.getIntervalInfoForInterval(ctx)) {
								if(!ctx.overlaps(ii)) continue;
								final Element rec = ii.rect(
										ii.trim(ctx),
										y,
										featureHeight,
										Maps.of(
											"stroke","none",
											"fill",(side==0?"red":(side==1?"green":"blue"))
											)
										);
								svgDoc.setTitle(rec,ctx.getContig()+":"+
										StringUtils.niceInt(ctx.getStart())+":"+
										ctx.getReference().getDisplayString()+" "+
										String.join(" ", sampleWithAlt)
										);
										
										
								g1.appendChild( svgDoc.anchor(rec, ctx));
								}// loop over ii
							}// if sampleWithAlt != empty
						}//loop over variants
					g1.appendChild(svgDoc.line(svgDoc.margin_left-5,y+featureHeight/2,svgDoc.margin_left,y+featureHeight/2));
					g1.appendChild(svgDoc.text(
							svgDoc.margin_left-5,
							y+fontSize+fontSize/2,
							(side==0?"Cases":(side==1?"Case and Controls":"Controls"))+" N="+StringUtils.niceInt(countVariants),
							Maps.of("fill",(side==0?"red":(side==2?"blue":"green")),
									"text-anchor","end",
									"font-size",fontSize+"px"
									)
							));
						
					g1.appendChild(svgDoc.comment("END "+side));
	
					y+= featureHeight;
					y+= 5;
					}// end loop for side /case+ctrl
				
				svgDoc.frame(svgDoc.lastY,y);
				svgDoc.lastY = y;
				
				svgDoc.addBanner("Sample summary", Maps.of("font-size", 30));
                y = svgDoc.lastY;
				
				{
				y+=2;
				final double y_top= y;
				final Element g1 = svgDoc.group();
				double featureHeight=12;
				final double pop_width = svgDoc.image_width*0.9;
				final double pop_x = svgDoc.margin_left + (svgDoc.image_width-pop_width)/2.0;
				final List<String> samples_with_alt = new ArrayList<>(all_samples);
				samples_with_alt.removeIf(SN->variants.stream().map(V->V.getGenotype(SN)).noneMatch(G->G.hasAltAllele()));
				
				final double sample_width= pop_width/samples_with_alt.size();
				final int maxVariantPerSamples= (int)samples_with_alt.stream().mapToLong(SN->variants.stream().map(V->V.getGenotype(SN)).filter(G->G.hasAltAllele()).count()).max().orElse(1L);
				double x=pop_x;
				y+=2;
				final List<NamedFilter> filters = Arrays.asList(
						new NamedFilter("ALL",V->true),
						new NamedFilter("SNP",V->V.isSNP()),
						new NamedFilter("Indel",V->V.isIndel()),
						new NamedFilter("Filtered",V->V.isFiltered()),
						new NamedFilter("PASS",V->!V.isFiltered()),
						new NamedFilter("Singleton",V->V.getGenotypes().stream().filter(G->G.hasAltAllele()).limit(2L).count()==1L),
						new NamedFilter("NOT Singleton",V->V.getGenotypes().stream().filter(G->G.hasAltAllele()).limit(2L).count()>1L)
						);

				for(NamedFilter filter:filters) {
					LOG.info("filter "+filter.getName());
					if( samples_with_alt.stream().
							flatMap(SN->variants.stream().filter(filter).map(V->V.getGenotype(SN))).
							filter(G->G.hasAltAllele()).
							noneMatch(G->samples_with_alt.contains(G.getSampleName()))) continue;
					
					/* frame */
					g1.appendChild(svgDoc.setTitle(svgDoc.rect(
							pop_x,
							y,
							pop_width,
							featureHeight,
							Maps.of("fill",Colors.toRGB(230),"stroke","none")
							),"Samples with Variants matching filter "+filter.getName()));
					
					/* legend */
					g1.appendChild(svgDoc.line(
							svgDoc.margin_left-5,
							y+featureHeight/2,
							svgDoc.margin_left,
							y+featureHeight/2
							));
					
					g1.appendChild(svgDoc.text(
							svgDoc.margin_left-5,
							y+fontSize,
							filter.getName()+" N="+variants.stream().filter(filter).count(),
							Maps.of(
									"font-size",fontSize+"px",
									"text-anchor", "end"
									)
							));
					x= pop_x;
					for(final String sn:samples_with_alt) {

						final int countVariantsForSample = (int)variants.stream().filter(filter).map(V->V.getGenotype(sn)).filter(G->G.hasAltAllele()).count();
						if(countVariantsForSample==0) {
							x+=sample_width;
							continue;
							}
						
						boolean is_case = this.casesControls.isCase(sn);
						double h =(countVariantsForSample/(double)maxVariantPerSamples)*featureHeight;
						final String fill_color =  Colors.shadeOf(
								(countVariantsForSample/maxVariantPerSamples),
								(is_case?"indianred":"dodgerblue"),
								(is_case?"darkred":"mediumblue")
								)
								;

						g1.appendChild(
								svgDoc.setTitle(
										svgDoc.rect(
												x, y+featureHeight-h,
												sample_width,
												h,
												Maps.of(
													"fill",fill_color,
													"stroke",(sample_width<3?"none":"orange"))
												),
										sn+ " "+ filter.getName()+" N(variants)="+countVariantsForSample
										)
								);
						x+=sample_width;
						}
					y+=featureHeight;
					y+=2;
					} /* end filter */
				y+=10;
				g1.appendChild(
						svgDoc.rect(pop_x, y_top, pop_width, y_top-y)
						);

				svgDoc.rootElement.appendChild(g1);
				svgDoc.lastY = y;
				}//end block
				
				svgDoc.finish();
				
				LOG.info("import SVG");
				divForSvg.appendChild(htmlDoc.importSVG(svgDoc));
				LOG.info("save");
				htmlDoc.saveToFileOrStdout(this.outputFile);
				LOG.info("saved");

				}
			
			return 0;
			}
		catch(Throwable err) {
			err.printStackTrace();
			LOG.error(err);
			return -1;
			}
		}
	
	public static void main(String[] args) {
		new PlotRareVariants().instanceMainWithExit(args);
	}

}
