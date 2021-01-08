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


History:

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;
import java.util.function.IntFunction;
import java.util.stream.Collectors;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.pedigree.Pedigree;
import com.github.lindenb.jvarkit.pedigree.PedigreeParser;
import com.github.lindenb.jvarkit.pedigree.Sample;
import com.github.lindenb.jvarkit.samtools.util.IntervalListProvider;
import com.github.lindenb.jvarkit.stream.HtsCollectors;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.structure.Exon;
import com.github.lindenb.jvarkit.util.bio.structure.Gene;
import com.github.lindenb.jvarkit.util.bio.structure.GtfReader;
import com.github.lindenb.jvarkit.util.bio.structure.Transcript;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.ns.XLINK;
import com.github.lindenb.jvarkit.util.svg.SVG;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;
/**
BEGIN_DOC

## Example

```

# run vcf2svg 
$ java -jar dist/vcf2svg.jar  -r "chr22:100-2000"  --gtf gtf.txt.gz input.vcf -o out.zip
```

## Screenshot

https://twitter.com/yokofakun/status/851875435948462080

![screenshot](https://pbs.twimg.com/media/C9J4LeoXkAEqvIN.jpg)




END_DOC
 */
@Program(name="vcf2svg",
	description="write a vcf to svg , with gene context",
	keywords={"vcf","svg","xlm","visualization"},
	modificationDate="20200715",
	creationDate="20170411"
	)
public class VcfToSvg extends Launcher {
private static final Logger LOG=Logger.build(VcfToSvg.class).make();

@Parameter(names={"-o","--out"},description=ArchiveFactory.OPT_DESC,required=true)
private Path outputPath=null;
@Parameter(names={"-r","--interval"},description=IntervalListProvider.OPT_DESC,converter=IntervalListProvider.StringConverter.class,splitter=NoSplitter.class,required=true)
private IntervalListProvider intervalListProvider = IntervalListProvider.unspecified();

@Parameter(names={"-g","--gtf"},description=GtfReader.OPT_DESC)
private Path gtfPath = null;

@Parameter(names={"-m","--manifest"},description="Manifest bed file containing the names of the files.")
private Path manifestFile=null;
@Parameter(names={"-gw","--genotypeWidth"},description="Genotype square width")
private int genotype_width=10;
@Parameter(names={"--nonCoding"},description="Ignore Non-coding transcripts")
private boolean removeNonCoding = false;
@Parameter(names={"--exon","--exons"},description="Only keep variants in exons")
private boolean variantsInExonOnly = false;
@Parameter(names={"--alphaFILTER"},description="Variant having FILTER!=PASS opacity (0== hide FILTERED variants)")
private double variantFILTEREDOpacity = 1.0;
@Parameter(names={"--alphaINDEL"},description="Variant INDEL opacity (0== hide INDEL variants) ")
private double variantIndelOpacity = 1.0;
@Parameter(names={"--pedigree"},description="Optional pedigree. "+PedigreeParser.OPT_DESC)
private Path pedPath=null;




private void title(final XMLStreamWriter w,final String title) throws XMLStreamException
	{
	w.writeStartElement("title");
	w.writeCharacters(title);
	w.writeEndElement();
	}

@Override
public int doWork(final List<String> args) {
	VCFReader r=null;
	OutputStream outputStream=null;
	XMLStreamWriter w=null;
	PrintWriter manifestW=null;
	ArchiveFactory archiveFactory = null;
	try {
		r= VCFReaderFactory.makeDefault().open(Paths.get(oneAndOnlyOneFile(args)), true);
		final VCFHeader header=r.getHeader();
		final List<String> samples= new ArrayList<>(header.getSampleNamesInOrder());
	
		
		
		final SAMSequenceDictionary dict= SequenceDictionaryUtils.extractRequired(header);
		intervalListProvider.dictionary(dict);
		
		/* read gtf if any */
		final IntervalTreeMap<Gene> geneMap = new IntervalTreeMap<>();
		if(this.gtfPath!=null) {
			try(GtfReader gtfReader=new GtfReader(this.gtfPath)) {
				gtfReader.setContigNameConverter(ContigNameConverter.fromOneDictionary(dict));
				 gtfReader.getAllGenes().
					stream().
					filter(G->!this.removeNonCoding || G.getTranscripts().stream().anyMatch(T->T.isCoding())).
					forEach(G->geneMap.put(new Interval(G),G));
				}
			}
		archiveFactory = ArchiveFactory.open(this.outputPath);
		if(manifestFile!=null)
			{
			manifestW= IOUtils.openPathForPrintWriter(this.manifestFile);
			}
		else
			{
			manifestW= new PrintWriter(new NullOuputStream());
			}
		
		final Pedigree pedigree;
		if(this.pedPath==null) {
			pedigree = PedigreeParser.empty();
		} else {
			pedigree = new PedigreeParser().parse(this.pedPath);
		}

		final Path tmpSvg = Files.createTempFile("vcf.", ".svg");
		final XMLOutputFactory xof=XMLOutputFactory.newInstance();
		
		for(final Locatable interval : intervalListProvider.dictionary(dict).stream().collect(Collectors.toList()))
			{
			final List<VariantContext> variants =  r.query(interval).stream().
					filter(V->this.variantFILTEREDOpacity>0 ||!V.isFiltered()).					
					filter(V->this.variantIndelOpacity>0 ||!V.isIndel()).					
					collect(Collectors.toCollection(ArrayList::new));
			
			if(variants.isEmpty()) continue;
			

			final List<Transcript> transcripts =  geneMap.
					getOverlapping(interval).
					stream().
					flatMap(G->G.getTranscripts().stream()).
					filter(T->!this.removeNonCoding || T.isCoding()).
					collect(Collectors.toList());

			variants.removeIf(V->this.gtfPath!=null && this.variantsInExonOnly && transcripts.stream().flatMap(T->T.getExons().stream()).noneMatch(EX->EX.overlaps(V)));

			if(variants.isEmpty()) continue;

			
			final String geneId =  transcripts.stream().
					map(T->T.getGene().getId()).
					collect(Collectors.toSet()).stream().
					collect(HtsCollectors.oneAndOnlyOne()).
					orElse(null);
			final String geneName =  transcripts.stream().
					map(T->T.getGene().getGeneName()).
					collect(Collectors.toSet()).stream().
					collect(HtsCollectors.oneAndOnlyOne()).
					orElse(null);

			
			
			outputStream=IOUtils.openPathForWriting(tmpSvg);
			w=xof.createXMLStreamWriter(outputStream);

			
			
			
			double featureHeight=10;
			double TRANSCRIPT_HEIGHT=featureHeight; 

			final int all_genotypes_width = variants.size()*this.genotype_width;
           
	        final int drawinAreaWidth=Math.max(all_genotypes_width,1000);
	
			final int interline_weight=6;
			final int margin_top=10;
			final int margin_bottom=10;
			final int margin_right=100;
			final int margin_left=100;
			
			w.writeStartDocument("UTF-8", "1.0");
			
            w.writeStartElement("svg");
            w.writeDefaultNamespace(SVG.NS);
            w.writeNamespace("xlink", XLINK.NS);
            w.writeAttribute("version", "1.1");
            w.writeAttribute("width",String.valueOf(margin_right+margin_right+drawinAreaWidth));
            w.writeAttribute("height",String.valueOf(
            		margin_top+margin_bottom+
            		transcripts.size()*TRANSCRIPT_HEIGHT+
            		interline_weight*featureHeight+
            		samples.size()*this.genotype_width
            		));
            title(w,interval.getContig()+":"+interval.getStart()+"-"+interval.getEnd());
            
            w.writeStartElement("desc");
            w.writeCharacters("generated with "+getProgramName()+"\n"+
            		"Author: Pierre Lindenbaum PhD. @yokofakun .");
            w.writeEndElement();
            
            //defs
            w.writeStartElement("defs");
			
            
    		//genotypes
    		
    		w.writeStartElement("g");
				w.writeAttribute("id","g_"+GenotypeType.HOM_REF); //
					w.writeEmptyElement("rect");
						w.writeAttribute("style","fill:lime;stroke;none;");
						w.writeAttribute("x","0");
						w.writeAttribute("y","0" );
						w.writeAttribute("width",String.valueOf(this.genotype_width));
						w.writeAttribute("height",String.valueOf(this.genotype_width));
    		w.writeEndElement();
    		
    		w.writeStartElement("g");
				w.writeAttribute("id","g_"+GenotypeType.NO_CALL); //
					w.writeEmptyElement("rect");
						w.writeAttribute("style","fill:silver;stroke;gray;");
						w.writeAttribute("x","0");
						w.writeAttribute("y","0" );
						w.writeAttribute("width",String.valueOf(this.genotype_width));
						w.writeAttribute("height",String.valueOf(this.genotype_width));
    		w.writeEndElement();
    		
			w.writeStartElement("g");
				w.writeAttribute("id","g_"+GenotypeType.HOM_VAR); //
					w.writeEmptyElement("rect");
						w.writeAttribute("style","fill:crimson;stroke;none;");
						w.writeAttribute("x","0");
						w.writeAttribute("y","0" );
						w.writeAttribute("width",String.valueOf(this.genotype_width));
						w.writeAttribute("height",String.valueOf(this.genotype_width));
    		w.writeEndElement();	
    		
    		w.writeStartElement("g");
			w.writeAttribute("id","g_"+GenotypeType.MIXED); //
				w.writeEmptyElement("rect");
					w.writeAttribute("style","fill:pink;stroke;none;");
					w.writeAttribute("x","0");
					w.writeAttribute("y","0" );
					w.writeAttribute("width",String.valueOf(this.genotype_width));
					w.writeAttribute("height",String.valueOf(this.genotype_width));
		    w.writeEndElement();	
    		
    		w.writeStartElement("g");
			w.writeAttribute("id","g_"+GenotypeType.UNAVAILABLE); //
				w.writeEmptyElement("rect");
					w.writeAttribute("style","fill:gray;stroke;none;");
					w.writeAttribute("x","0");
					w.writeAttribute("y","0" );
					w.writeAttribute("width",String.valueOf(this.genotype_width));
					w.writeAttribute("height",String.valueOf(this.genotype_width));
		    w.writeEndElement();	

		    
    		w.writeStartElement("g");
    		w.writeAttribute("id","g_"+GenotypeType.HET); //
    			w.writeEmptyElement("rect");
    			w.writeAttribute("style","fill:lime;stroke;black;");
    			w.writeAttribute("x", "0" );
    			w.writeAttribute("y", "0" );
        		w.writeAttribute("width",String.valueOf(genotype_width));
        		w.writeAttribute("height",String.valueOf(genotype_width));
    			w.writeEmptyElement("polygon");
    			w.writeAttribute("style","fill:crimson;stroke;black;");
    			w.writeAttribute("points","0,0 "+genotype_width+",0 0,"+genotype_width+" 0,0");
    		w.writeEndElement();

    		
    		//strand
    		w.writeEmptyElement("polyline");
    			w.writeAttribute("id","strandF");
    			w.writeAttribute("points", "-5,-5 0,0 -5,5" );
    		
    		w.writeEmptyElement("polyline");
    			w.writeAttribute("id","strandR");
    			w.writeAttribute("points", "5,-5 0,0 5,5" );
    		
    		//gradients
    			w.writeStartElement("linearGradient");
    				w.writeAttribute("id","grad01");
    				w.writeAttribute("x1","50%");
    				w.writeAttribute("x2","50%");
    				w.writeAttribute("y1","0%");
    				w.writeAttribute("y2","100%");
    				w.writeEmptyElement("stop");
    					w.writeAttribute("offset","0%");
    					w.writeAttribute("style","stop-color:black;stop-opacity:1;");
    				w.writeEmptyElement("stop");
    					w.writeAttribute("offset","50%");
    					w.writeAttribute("style","stop-color:white;stop-opacity:1;");
    				w.writeEmptyElement("stop");
    					w.writeAttribute("offset","100%");
    					w.writeAttribute("style","stop-color:black;stop-opacity:1;");
    			w.writeEndElement();
    		
    		w.writeEndElement();//defs

    		w.writeStartElement("style");
    		w.writeCharacters(
    				"svg {fill:none; stroke:black;}\n"+
    				"text {fill:black;stroke:none;font-size:"+ (featureHeight/1.5) +"px;}\n"+
    				".ruler-label { stroke:red;}\n"+
    				".frame { stroke:black;fill:none;}\n"+
    				".kgexon {fill:url(#grad01);stroke:black;}\n"+
    				".gcpercent {fill:url(#grad02);stroke:black;}"+
    				".coverage {fill:url(#grad03);stroke:black;}"+
    				".kgcds {fill:yellow;stroke:black;opacity:0.7;}\n"+
    				".variant{stroke:none;fill:red;opacity:0.2;}\n"+
    				".xaxis{stroke:gray;fill:none;opacity:0.2;}\n"+
    				".postick{font-size:9px;stroke:black;stroke-width:1;}"
    				);
    		w.writeEndElement();//style

            
            	
			final IntFunction<Integer> trim= t-> Math.max(interval.getStart(), Math.min(interval.getEnd(), t));
			
			final IntFunction<Double> baseToPixel=  t->margin_left + drawinAreaWidth*(t-(double)interval.getStart())/((double)interval.getLengthOnReference());
				
			final IntFunction<Double> variantIndexToPixel= idx-> {
					final double variant_width= drawinAreaWidth/(double)variants.size();
					final double midx=variant_width*idx+variant_width/2.0;
					return margin_left+ midx-genotype_width/2.0;
					};
			
			
			final Function<VariantContext,String> variantTitle= V-> (V.getContig().startsWith("chr")?V.getContig().substring(3):V.getContig())+":"+V.getStart()+" "+V.getReference().getDisplayString();
			
			/** title */
			double y=0;
			w.writeStartElement("text");
			w.writeAttribute("x","0");
			w.writeAttribute("y",String.valueOf(featureHeight));
			w.writeCharacters(interval.toString());
			w.writeEndElement();
			y+= featureHeight;
			
			for(final Transcript g:transcripts)
				{
				int cdsHeigh= 5;
				double exonHeight=TRANSCRIPT_HEIGHT-5;
				double midY=TRANSCRIPT_HEIGHT/2;
		
				w.writeStartElement("g");
				
				
				
				w.writeAttribute("transform", "translate(0,"+y+")");
				
				title(w, g.getId());
				
				w.writeStartElement("text");
				w.writeAttribute("x",String.valueOf(margin_left-10));
				w.writeAttribute("y",String.valueOf(featureHeight));
				w.writeAttribute("style","text-anchor:end;");
				w.writeCharacters(g.getId());
				w.writeEndElement();
				
				/* transcript line */
				w.writeEmptyElement("line");
				w.writeAttribute("class","kgtr");
				w.writeAttribute("x1",String.valueOf(baseToPixel.apply(trim.apply(g.getTxStart()))));
				w.writeAttribute("y1",String.valueOf(midY));
				w.writeAttribute("x2",String.valueOf(baseToPixel.apply(trim.apply(g.getTxEnd()))));
				w.writeAttribute("y2",String.valueOf(midY));
				
				
				
				
				/* strand symbols */
				for(double pixX=0;
					pixX< drawinAreaWidth;
					pixX+=30)
					{
					double pos0= interval.getStart()+(pixX/(double)drawinAreaWidth)*interval.getLengthOnReference();
						if(pos0+1< g.getTxStart()) continue;
						if(pos0> g.getTxEnd()) break;
						w.writeEmptyElement("use");
						w.writeAttribute("class","kgstrand");
						w.writeAttribute("xlink", XLINK.NS, "href", "#strand"+(g.isPositiveStrand()?"F":"R"));
						w.writeAttribute("x",String.valueOf(margin_left + pixX));
						w.writeAttribute("y",String.valueOf(midY));
					}
			
				/* exons */
				for(final Exon exon:g.getExons())
						{
						if(exon.getStart()+1>=  interval.getEnd()) continue;
						if(exon.getEnd()<= interval.getStart()) continue;
						w.writeStartElement("rect");
						w.writeAttribute("class","kgexon");
						
						w.writeAttribute("x",String.valueOf(baseToPixel.apply(trim.apply(exon.getStart()))));
						w.writeAttribute("y",String.valueOf(midY-exonHeight/2));
						w.writeAttribute("width",String.valueOf(baseToPixel.apply(trim.apply(exon.getEnd()))-baseToPixel.apply((trim.apply(exon.getStart())))));
						w.writeAttribute("height",String.valueOf(exonHeight));
						title(w,exon.getName());
						w.writeEndElement();
						}
			
				/* coding line */
				if(!g.isNonCoding() && g.hasCodonStartDefined() && g.hasCodonStopDefined())
					{
					final double codonx1 = baseToPixel.apply(trim.apply(g.getLeftmostCodon().get().getStart()));
					final double codonx2 = baseToPixel.apply(trim.apply(g.getRightmostCodon().get().getEnd()));
					w.writeEmptyElement("rect");
					w.writeAttribute("class","kgcds");
					w.writeAttribute("x",String.valueOf(codonx1));
					w.writeAttribute("y",String.valueOf(midY-cdsHeigh/4.0));
					w.writeAttribute("width",String.valueOf(baseToPixel.apply((int)(codonx2 - codonx1))));
					w.writeAttribute("height",String.valueOf(cdsHeigh/2.0));
					}
				
				//String label=String.format("%15s", g.getName());
				//w.writeEmptyElement("path");
				//double fontHeight=Math.min(10,0.8*TRANSCRIPT_HEIGHT);
				//w.writeAttribute("d",this.hershey.svgPath(label,-insets.left,midY-fontHeight/2,insets.left*0.9,fontHeight));

				
				w.writeEndElement();
				w.writeCharacters("\n");
				y+=featureHeight;
				}
			
			/* draw lines to variants */
			for(int vidx=0;vidx < variants.size();++vidx)
				{
				final VariantContext vc=variants.get(vidx);
				double x1 =  baseToPixel.apply(vc.getStart());
				double x2 =  baseToPixel.apply(vc.getEnd());
				final double y2= y+featureHeight*interline_weight;
				w.writeStartElement("polygon");
				w.writeAttribute("style","fill:"+(vidx%2==0?"ghostwhite":"lavender")+";stroke:black;opacity:0.6;stroke-width:0.5;");
				w.writeAttribute("points",
						 ""+x1+","+(y-featureHeight/2.0)+
						" "+x2+","+(y-featureHeight/2.0)+
						" "+variantIndexToPixel.apply(vidx)+","+y2+
						" "+(variantIndexToPixel.apply(vidx)+this.genotype_width)+","+y2
						);
				title(w,variantTitle.apply(vc));
				w.writeEndElement();
				}
			for(int vidx=0;vidx < variants.size();++vidx)
				{
				final VariantContext vc=variants.get(vidx);
				final double y2= y+featureHeight*interline_weight;
				w.writeStartElement("text");
				w.writeAttribute("transform",
						"translate("+(String.valueOf(variantIndexToPixel.apply(vidx)+genotype_width/2.0))+","+String.valueOf(y2-5)+") "+
						"rotate(-45)");
				w.writeAttribute("x", "0");
				w.writeAttribute("y", "0");
				w.writeAttribute("class", "postick");
				w.writeCharacters(variantTitle.apply(vc));
				w.writeEndElement();
				w.writeCharacters("\n");
				}
			y+=featureHeight*interline_weight;
			
			w.writeStartElement("g");
			/* step 0: affected, 1: unaffected, 2: others */
			for(int step=0;step<3;++step) {
				for(final String sample: samples )
					{
					final Sample individual = pedigree.getSampleById(sample);
					if(step==0 && (individual==null || !individual.isAffected())) continue;
					if(step==1 && (individual==null || !individual.isUnaffected())) continue;
					if(step==2 && individual!=null && individual.isStatusSet()) continue;
					
					w.writeStartElement("g");//sample
					switch(step) {
						case 0: w.writeAttribute("style", "hue-rotate(195deg);");break;
						case 1: w.writeAttribute("style", "hue-rotate(45deg);");break;
						default:break;
						}
					
					for(int vidx=0;vidx < variants.size();++vidx)
						{
						final VariantContext vc=variants.get(vidx);
						final Genotype g=vc.getGenotype(sample);
						double opacity= 1.0;
						
						if(vc.isIndel())  opacity*=this.variantIndelOpacity;
						if(vc.isFiltered())  opacity*=this.variantFILTEREDOpacity;
						if(opacity>1) opacity=1;
						if(opacity<=0) continue;
						if(opacity<1)
							{
							w.writeStartElement("g");
							w.writeAttribute("style","opacity:"+opacity+";");
							}
						
						w.writeEmptyElement("use");
						w.writeAttribute("x",String.valueOf(variantIndexToPixel.apply(vidx)));
						w.writeAttribute("y",String.valueOf(y));
						w.writeAttribute("xlink", XLINK.NS, "href", "#g_"+g.getType());
						
						if(opacity<1)
							{
							w.writeEndElement();
							}
						
						}
					w.writeCharacters("\n");
					w.writeStartElement("text");
					w.writeAttribute("x",String.valueOf(margin_left-10));
					w.writeAttribute("y",String.valueOf(y + this.genotype_width/2.0));
					w.writeAttribute("style","text-anchor:end;");
					w.writeCharacters(sample);
					w.writeEndElement();//text
					
					w.writeEndElement();//g for sample
					
					y+=this.genotype_width;
					}
			}
			w.writeCharacters("\n");
			
			
			
			w.writeEndDocument();
			w.writeCharacters("\n");
			w.flush();
			w.close();
			
			
			final String md5 = StringUtils.md5(interval.getContig()+":"+interval.getStart()+"-"+interval.getEnd());
			final String filename =  md5.substring(0,2) + File.separatorChar + md5.substring(2) + File.separator+
					interval.getContig().replaceAll("[/\\:]", "_") +
					"_"+
					interval.getStart() +
					"_"+
					interval.getEnd() +
					(StringUtils.isBlank(geneName)?"":"."+geneName.replaceAll("[/\\:]", ""))+
					(StringUtils.isBlank(geneId)?"":"."+geneId.replaceAll("[/\\:]", ""))+
					".svg"
					;

			
			OutputStream os = archiveFactory.openOuputStream(filename);
			IOUtils.copyTo(tmpSvg, os);
			os.flush();
			os.close();
			
			Files.delete(tmpSvg);
			
			
			manifestW.print(interval.getContig());
			manifestW.print('\t');
			manifestW.print(interval.getStart()-1);
			manifestW.print('\t');
			manifestW.print(interval.getEnd());
			manifestW.print('\t');
			manifestW.print(transcripts.stream().map(G->G.getGene().getId()).collect(Collectors.toSet()).stream().collect(Collectors.joining(";")));
			manifestW.print('\t');
			manifestW.print(transcripts.stream().map(G->G.getGene().getGeneName()).collect(Collectors.toSet()).stream().collect(Collectors.joining(";")));
			manifestW.print('\t');
			manifestW.print(transcripts.stream().map(G->G.getId()).collect(Collectors.toSet()).stream().collect(Collectors.joining(";")));
			manifestW.print('\t');
			manifestW.print((archiveFactory.isTarOrZipArchive()?"":this.outputPath.toString()+File.separator)+filename);
			manifestW.print('\t');
			manifestW.println(variants.size());
			}
		r.close();
		manifestW.flush();
		manifestW.close();
		manifestW=null;
		
		archiveFactory.close();
		archiveFactory=null;

		return 0;
		}
	catch(final Throwable err)
		{
		LOG.error(err);
		return -1;
		}
	finally
		{
		CloserUtil.close(archiveFactory);
		CloserUtil.close(r);
		CloserUtil.close(outputStream);
		CloserUtil.close(manifestW);
		}
	}
	
	
public static void main(final String[] args) {
	new VcfToSvg().instanceMainWithExit(args);
	}
}
