package com.github.lindenb.jvarkit.tools.genome2svg.beans;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;
import java.util.function.IntToDoubleFunction;

import org.w3c.dom.Element;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.tools.genome2svg.SVGContext;
import com.github.lindenb.jvarkit.ucsc.Cytoband;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

public class CytobandTrack extends Track {
	private Function<Locatable,String> locatable2url = null;

	@Override
	public void paint(final SVGContext ctx) {
		Function<Locatable,String> loc2url = this.locatable2url;
		final List<Cytoband> entries= new ArrayList<>();
		String ucsc_name = "";
		BufferedReader br=null;
		try {
			
			if(getReference()!=null) {
				final SAMSequenceDictionary dict = getReference().getSAMSequenceDictionary();
				if(dict!=null) {
					if(SequenceDictionaryUtils.isGRCh37(dict)) {
						ucsc_name = "hg19";
						}
					else if(SequenceDictionaryUtils.isGRCh38(dict)) {
						ucsc_name = "hg38";
						}
					}
				}
			if(!StringUtils.isBlank(ucsc_name) && loc2url==null) {
				final String final_ucsc_name= ucsc_name;
				loc2url = L->"http://genome.ucsc.edu/cgi-bin/hgTracks?org=Human&db="+final_ucsc_name+"&position=";
				}
			
			
			if(getPath()==null) {
				if(StringUtils.isBlank(ucsc_name)) return;
				final String url = "https://hgdownload.cse.ucsc.edu/goldenpath/"+ucsc_name+"/database/cytoBand.txt.gz";
				br = IOUtils.openURIForBufferedReading(url);
				}
			else
				{
				br = IOUtils.openPathForBufferedReading(getPath().asPath());
				}
			for(;;) {
				final String line = br.readLine();
				if(line==null) break;
				final Cytoband e = Cytoband.of(line);
				if(!e.contigsMatch(ctx.loc)) continue;
				entries.add(e);
				}
			}
		catch(final IOException err) {
			throw new RuntimeIOException(err);
			}
		finally
			{
			if(br!=null) try {br.close();} catch(Throwable err) {}
			}

		
		
		
		if(entries.isEmpty()) return;
		if(loc2url==null) {
			loc2url= L->"";
			}
		
		ctx.appendStyle(
			".ctgborderp {fill:url(#grad01);stroke:green;}" +
			".ctgborderq {fill:url(#grad01);stroke:green;}" +
			".ctglabel {text-anchor:middle;stroke:none;fill:darkgrey;font: bold 10px Verdana, Helvetica, Arial, sans-serif;}" +
			".cytoband {fill:silver;stroke:none;}" +
			".bedlabel {stroke:red;fill:none;text-anchor:start;font: normal 7px Verdana, Helvetica, Arial, sans-serif;}" +
			".maintitle {stroke:none;fill:darkgrey;text-anchor:middle;font: normal 12px Verdana, Helvetica, Arial, sans-serif;}" +
			".ctgback {fill:gainsboro;stroke:none;filter:url(#filter01);}"
			);
		
		final double contig_len = entries.stream().mapToLong(T->T.getEnd()).max().orElse(1L);
		
		final IntToDoubleFunction loc2pix = P->(P/contig_len)*ctx.image_width;
		final double featureHeight=20;	
		final Element g = ctx.element("g");
		ctx.tracksNode.appendChild(g);
		final List<String> names= new ArrayList<>();
		names.add(ctx.loc.getContig()+":");
		for(Cytoband e:entries) {
			final Element r = ctx.element("rect");
			double x1 = loc2pix.applyAsDouble(e.getStart());
			double x2 = loc2pix.applyAsDouble(e.getEnd()+1);
			r.setAttribute("style", e.getCssColor());
			r.setAttribute("x", format(x1));
			r.setAttribute("y", format(ctx.y));
			r.setAttribute("width", format(x2-x1));
			r.setAttribute("height", format(featureHeight));
			g.appendChild(ctx.anchor(r, loc2url.apply(e)));
			r.appendChild(ctx.title(e.getName()));
			
			if(CoordMath.overlaps(e.getStart(), e.getEnd(), ctx.loc.getStart(), ctx.loc.getEnd())) {
				names.add(e.getName());
				}
			}
		Element frame = ctx.element("rect");
		frame.setAttribute("style","fill:none;stroke:red;");
		frame.setAttribute("x", format(loc2pix.applyAsDouble(ctx.loc.getStart())));
		frame.setAttribute("y", format(ctx.y));
		frame.setAttribute("width", format(loc2pix.applyAsDouble(ctx.loc.getEnd())-loc2pix.applyAsDouble(ctx.loc.getStart())));
		frame.setAttribute("height", format(featureHeight));
		g.appendChild(frame);
		
		final Element legend = ctx.element("text",String.join(" ", names));
		legend.setAttribute("x", format(-5));
		legend.setAttribute("y", format(ctx.y));
		g.appendChild(legend);
		
		ctx.y += featureHeight+2;
	}


}
