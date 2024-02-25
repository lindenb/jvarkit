package com.github.lindenb.jvarkit.tools.plotrare;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.w3c.dom.Element;

import com.beust.jcommander.DynamicParameter;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.css.Colors;
import com.github.lindenb.jvarkit.gtf.GtfTrack;
import com.github.lindenb.jvarkit.html.HTMLDocument;
import com.github.lindenb.jvarkit.lang.AttributeMap;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.net.Hyperlink;
import com.github.lindenb.jvarkit.pedigree.CasesControls;
import com.github.lindenb.jvarkit.samtools.util.LocatableUtils;
import com.github.lindenb.jvarkit.svg.GenomeSVGDocument;
import com.github.lindenb.jvarkit.svg.SVGDocument;
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
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;

public class PlotRareVariants extends Launcher {
	private static final Logger LOG = Logger.build(PlotRareVariants.class).make();
	@Parameter(names={"-o","--output"},description="TODO")
	protected Path outputFile=null;
	@Parameter(names={"-B","--bed"},description="TODO",required=true)
	protected Path bedFile=null;
	@Parameter(names={"-gtf","--gtf"},description="TODO",required=false)
	protected Path gtfFile=null;
	@DynamicParameter(names={"-D"},description="properties")
	protected Map<String,String> _properties = new HashMap<>();
	@ParametersDelegate
	protected CasesControls casesControls = new CasesControls();
	
	@Override
	public int doWork(List<String> args) {
		try {
			final HTMLDocument htmlDoc = new HTMLDocument();
			
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
				final Set<String> regionsIds0 = intervals.stream().map(L->L.getOrDefault(3,"")).collect(Collectors.toSet());
				
				
				final List<VariantContext> variants = new ArrayList<>();
					
				for(Locatable loc: intervals) {
					variants.addAll(vcfR.query(loc).toList());
					}
				Collections.sort(variants, ContigDictComparator.createLocatableComparator(dict));
				
				

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
				
				final HTMLDocument.Table table= htmlDoc.createTable(Arrays.asList("Variant","REF","ALT","INFO",
						"CASES (N="+StringUtils.niceInt(this.casesControls.getCases().size())+")",
						"CONTROLS (N="+StringUtils.niceInt(this.casesControls.getControls().size())+")"
						));
				htmlDoc.bodyElement.appendChild(table.table);
				for(int i=0;i< variants.size();i++) {
					final VariantContext ctx = variants.get(i);
					table.set(i, 0,htmlDoc.anchor(
							ctx.getContig()+":"+StringUtils.niceInt(i),
							hyperlink.apply(variants.get(i)).orElse("")
							));
					table.set(i, 1,ctx.getReference().getDisplayString());
					table.set(i, 2,ctx.getAlternateAlleles().stream().map(A->A.getDisplayString()).collect(Collectors.joining(",")));
					
					for(int side=0;side<2;++side) {
						List<Genotype> gL=ctx.getGenotypes(this.casesControls.get(side)).stream().filter(G->G.hasAltAllele()).collect(Collectors.toList());
						if(gL.isEmpty()) continue;
						table.set(i, 4+side, gL.stream().map(A->A.getSampleName()).collect(Collectors.joining(" "))+" N="+gL.size());
						}
					}
				
				
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

				
				double y = svgDoc.lastY;
				final Element g  = svgDoc.group();
				svgDoc.rootElement.appendChild(g);
				
				for(int side=0;side< 2;++side) {
					int featureHeight=10;
					final Set<String> samples = this.casesControls.get(side);
					final Element g1 = svgDoc.group();
					g.appendChild(g1);
					g1.setAttribute("transform", "translate(0,"+(y)+")");
					g1.appendChild(svgDoc.comment("BEGIN "+(side==0?"CASES":"CONTROLS")));
					int countVariants =0;
					
					for(final VariantContext ctx:variants) {
						final List<String> sampleWithAlt = ctx.getGenotypes().
								stream().
								filter(G->samples.contains(G.getSampleName())).
								filter(G->G.hasAltAllele()).
								map(G->G.getSampleName()).
								collect(Collectors.toList());
								;
						if(sampleWithAlt.isEmpty()) continue;
						countVariants++;
						final double height2 = featureHeight/(double)sampleWithAlt.size();
						
						for(GenomeSVGDocument.IntervalInfo ii: svgDoc.getIntervalInfoForInterval(ctx)) {
							final Element g2 = svgDoc.group(Maps.of(
									"stroke", "none",
									"opacity","0.7"
									));
							g1.appendChild(g2);
							for(int i=0;i< sampleWithAlt.size();i++) {
									final Element rec = ii.rect(
											ii.trim(ctx),
											i*height2,
											height2,
											Maps.of(
												"fill",(side==0?
														Colors.shadeOf(i/(float)sampleWithAlt.size(),"crimson","carmine"):
														Colors.shadeOf(i/(float)sampleWithAlt.size(),"royalblue","blue"))
												)
											);
									svgDoc.setTitle(rec,ctx.getContig()+":"+
											StringUtils.niceInt(ctx.getStart())+":"+
											ctx.getReference().getDisplayString()+" "+
											sampleWithAlt.get(i)+
											" N="+sampleWithAlt.size()
											);
									g2.appendChild(rec);
									}
							g1.appendChild( svgDoc.anchor(g2, ctx));
							}
						}
				g1.appendChild(svgDoc.text(
						svgDoc.margin_left-5,
						0,
						(side==0?"Cases":"Controls")+" N="+StringUtils.niceInt(countVariants),
						Maps.of("fill",(side==0?"red":"blue"),"text-anchor","end","font-size","10px")
						));
					
				g1.appendChild(svgDoc.comment("END "+(side==0?"CASES":"CONTROLS")));

				y+= featureHeight;
				y+= 5;
				}// end loop for side /case+ctrl
				
				svgDoc.frame(svgDoc.lastY,y);
				svgDoc.lastY = y;
				
				{
				y+=10;
				final Element g1 = svgDoc.group();
				double featureHeight=12;
				final double pop_width = svgDoc.image_width*0.9;
				final double pop_x = svgDoc.margin_left + svgDoc.image_width-pop_width;
				final double sample_width= pop_width/(casesControls.getAll().size());
				final int maxVariantPerSamples= (int)casesControls.getAll().stream().mapToLong(SN->variants.stream().map(V->V.getGenotype(SN)).filter(G->G.hasAltAllele()).count()).max().orElse(1L);
				double x=pop_x;
				for(int side=0;side<2;++side) {
					for(final String sn: casesControls.get(side)) {
						g1.appendChild(
								svgDoc.setTitle(
										svgDoc.rect(x, y, sample_width, featureHeight,Maps.of("fill", (side==0?"red":"blue"),"stroke","orange")),
										sn
										)
								);
						int countVariantsForSample = (int)variants.stream().map(V->V.getGenotype(sn)).filter(G->G.hasAltAllele()).count();
						if(countVariantsForSample>0) {
							double h =(countVariantsForSample/(double)maxVariantPerSamples)*featureHeight;
							String color = Colors.shadeOf((countVariantsForSample/(float)maxVariantPerSamples), "black","darkgray");
							g1.appendChild(
									svgDoc.setTitle(
											svgDoc.rect(x,  y+featureHeight+featureHeight-h, sample_width, h,Maps.of("fill",color,"stroke","none")),
											sn+" "+countVariantsForSample+" Variant(s)."
											)
									);
							}
						x+=sample_width;
						}
					}
				y+=10;
				g1.appendChild(
						svgDoc.rect(pop_x, y, pop_width, featureHeight*2+10,Maps.of("fill","none","stroke","gray"))
						);
				y+=featureHeight*3;
				
				svgDoc.rootElement.appendChild(g1);
				svgDoc.frame(svgDoc.lastY,y);
				svgDoc.lastY = y;
				}
				
				svgDoc.finish();
				
				
				htmlDoc.bodyElement.appendChild(htmlDoc.importSVG(svgDoc));
				LOG.info("save");
				svgDoc.saveToFileOrStdout(this.outputFile);
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
