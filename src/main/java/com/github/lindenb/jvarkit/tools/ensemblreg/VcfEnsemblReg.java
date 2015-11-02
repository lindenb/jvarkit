/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.ensemblreg;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.seekablestream.SeekableStreamFactory;

import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BedFeature;
import org.broad.igv.bbfile.BigBedIterator;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;

import com.github.lindenb.jvarkit.util.igv.SeekableStreamAdaptor;
import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class VcfEnsemblReg extends AbstractVcfEnsemblReg
	{
	private static class Track
		{
		String id;
		String type=null;
		String parent=null;
		String shortLabel=null;
		String longLabel=null;
		URL url=null;
		}
	
	
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(VcfEnsemblReg.class);

	@Override
	public Command createCommand() {
		return new MyCommand();
		}

	static private class MyCommand extends AbstractVcfEnsemblReg.AbstractVcfEnsemblRegCommand
		{    

	private List<Track> tracks=new ArrayList<Track>();

	private static final String BrightRed="187,0,0";
	private static final String LightRed="255,105,105";
	private static final String Orange="250,202,0";
	private static final String Blue="10,190,254";
	private static final String Gold="209,157,0";
	private static final String Yellow="255,252,4";
	private static final String LightGray="225,225,225";
	private static final String Gray="127,127,127";
	private static final String DarkGreen="0,176,80";
	

	/** http://ngs.sanger.ac.uk/production/ensembl/regulation/hg19/segmentation_summaries/ */
	private String segway_17SegmentationSummaries(String color)
	{
	if(color==null) return null;
		switch(color)
		{
		case BrightRed: return  "ActivePromoters";
		case LightRed :return "ProximalEnhancer";
		case Orange : return "DistalEnhancer";
		case Blue : return "DistalCTF";
		case DarkGreen: return " TranscriptionAssociated";
		case Gray: return "PolycombRepressed";
		case LightGray: return null;//  Weak signal  | Heterochromatin/Repetitive/Copy Number Variation  
		default: System.err.println("UNKNOWN Seg COLOR "+color);return null;
		}
	}
	/** http://ngs.sanger.ac.uk/production/ensembl/regulation/hg19/projected_segmentations/ */
	private String projectedSegments(String color)
	{
	if(color==null) return null;
		switch(color)
		{
		case BrightRed: return  "PredictedActivePromoters";
		case LightRed :return "PredicteActivedPromoterFlankingRegions";
		case Orange : return "PredictedActiveEnhancers";
		case Blue : return "ActiveCTCFBindingSite";
		case Gold : return "UnannotatedActiveTFBS";
		case Yellow: return "UnannotatedActiveOpenChromatinRegions" ;
		case LightGray: return "InactiveRegions"; 
		default: LOG.warn("projectedSegments: undefined color: "+color); return null;
		}
	}
	/** http://ngs.sanger.ac.uk/production/ensembl/regulation/hg19/overview/ */
	private String regBuildOverview(String color)
	{
	if(color==null) return null;
		switch(color)
		{
		case BrightRed: return  "PredictedPromoters";
		case LightRed :return "PredictedPromoterFlankingRegions";
		case Orange : return "PredictedEnhancers";
		case Blue : return "CTCFBindingSite";
		case Gold : return "UnannotatedTFBS";
		case Yellow: return "UnannotatedOpenChromatinRegions" ;
		default: LOG.warn("regBuildOverview: undefined color: "+color);return null;
		}
	}
	
	/** http://ngs.sanger.ac.uk/production/ensembl/regulation/hg19/segmentations/ */
	private String segway_17CellSegments(String color)
		{
		if(color==null) return null;
			switch(color)
			{
			case BrightRed: return  "ActivePromoters";
			case LightRed :return "ProximalEnhancer";
			case Orange : return "DistalEnhancer";
			case Blue : return "DistalCTF";
			case DarkGreen: return "TranscriptionAssociated";
			case Gray: return "PolycombRepressed";
			case LightGray: return null;//  Weak signal  | Heterochromatin/Repetitive/Copy Number Variation  
			default: System.err.println("UNKNOWN Seg COLOR "+color);return null;
			}
		}

	private void annotate(Track track,File inf,File outf) throws IOException
		{
		boolean contained=false;
		LOG.warn("Processing "+track.id+" ("+track.shortLabel+") "+track.url);
		VcfIterator in=VCFUtils.createVcfIteratorFromFile(inf);
		VCFHeader header=in.getHeader();
		VCFInfoHeaderLine info=null;
		
		
		SeekableStream sstream=SeekableStreamFactory.getInstance().getStreamFor(track.url);
		BBFileReader bigFile = new BBFileReader(track.url.toString(), new SeekableStreamAdaptor(sstream));
		
		VariantContextWriter w1=VCFUtils.createVariantContextWriter(outf);
		
		if(bigFile.isBigWigFile())
			{
			info = new VCFInfoHeaderLine(track.id, 1, VCFHeaderLineType.Float,String.valueOf(track.longLabel)+" "+track.url);
			}
		else
			{
			info = new VCFInfoHeaderLine(track.id, 1, VCFHeaderLineType.String,String.valueOf(track.longLabel)+" "+track.url);
			}
		
		header.addMetaDataLine(info);
		w1.writeHeader(in.getHeader());

		
		
		while(in.hasNext())
			{
			VariantContext ctx=in.next();
			String chrom=ctx.getContig();
			if(!chrom.startsWith("chr")) chrom="chr"+chrom;
			if(!chrom.matches("(chrX|chrY|chr[0-9]|chr1[0-9]|chr2[12])"))
				{
				w1.add(ctx);
				}
			else if(bigFile.isBigWigFile())
				{
				BigWigIterator iter=bigFile.getBigWigIterator(
						chrom,
						ctx.getStart()-1,
						chrom,
						ctx.getStart(),
						contained
						);
				Float wigValue=null;
				while(iter!=null && iter.hasNext() && wigValue==null)
					{
					WigItem item=iter.next();
					wigValue=item.getWigValue();
					}
				if(wigValue==null)
					{
					w1.add(ctx);
					continue;
					}
				VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				vcb.attribute(track.id, wigValue);
				w1.add(vcb.make());
				}
			else
				{
				BigBedIterator iter = bigFile.getBigBedIterator(
						chrom,
						ctx.getStart()-1,
						chrom,
						ctx.getStart(),
						contained
						);
				Set<String> bedValues=new HashSet<String>();
				while(iter!=null && iter.hasNext())
					{
					BedFeature item=iter.next();
					String rest[]=item.getRestOfFields();
					if(rest==null || rest.length!=6)
						{
						System.err.println(track.id+" "+Arrays.toString(item.getRestOfFields()));
						continue;
						}
					String color=null;
					if(track.parent!=null)
						{
						if(track.parent.startsWith("Segway_17SegmentationSummaries"))
							{
							color = segway_17SegmentationSummaries(rest[5]);
							}
						else if(track.parent.startsWith("ProjectedSegments") )
							{
							color = projectedSegments(rest[5]);
							}
						else if(track.parent.startsWith("RegBuildOverview") )
							{
							color = regBuildOverview(rest[5]);
							}
						else if(track.parent.startsWith("Segway_17CellSegments") )
							{
							color = segway_17CellSegments(rest[5]);
							}
						else
							{
							System.err.println("Unknown parent:"+track.parent);
							}
						}
					
					
					if(color==null) continue;
					bedValues.add(rest[0]+"|"+color);
					} 
				if(bedValues.isEmpty())
					{
					w1.add(ctx);
					continue;
					}
				StringBuilder sb=new StringBuilder();
				for(String s:bedValues) { if(sb.length()!=0) sb.append(",");sb.append(s);}
				VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				vcb.attribute(track.id, sb.toString());
				w1.add(vcb.make());
				}
			}
		
		sstream.close();
		in.close();
		w1.close();
		}
	@Override
		protected Collection<Throwable> doVcfToVcf(String inputName,
				VcfIterator in, VariantContextWriter out) throws IOException
		{
		File tmpDir= this.getTmpDirectories().get(0);
		File tmp1=File.createTempFile("tmp_",".vcf", tmpDir);
		File tmp2=File.createTempFile("tmp_",".vcf", tmpDir);
		tmp1.deleteOnExit();
		tmp2.deleteOnExit();
		
		VariantContextWriter w1=VCFUtils.createVariantContextWriter(tmp1);
		VCFHeader header=in.getHeader();
		addMetaData(header);
		w1.writeHeader(header);
		while(in.hasNext())
			{
			w1.add(in.next());
			}
		in.close();
		w1.close();
		
	
		for(Track track:this.tracks)
			{
			if(track.url==null) continue;
			annotate(track, tmp1,tmp2);
			tmp1.delete();

			tmp2.renameTo(tmp1);
			}
		in=VCFUtils.createVcfIteratorFromFile(tmp1);
		out.writeHeader(in.getHeader());
		while(in.hasNext())
			{
			out.add(in.next());
			}
		in.close();
		tmp1.delete();
		tmp2.delete();
		return RETURN_OK;
		}
	
	private void parseTrackDB(URL url) throws IOException
		{
		LOG.info("Parsing "+url);
		BufferedReader r=new BufferedReader(new InputStreamReader(url.openStream(),"UTF-8"));
		String line;
		Track track=null;
		for(;;)
			{
			line=r.readLine();
			
			while(line!=null &&
					!line.isEmpty() &&
					Character.isWhitespace(line.charAt(0)))
					line=line.substring(1);
				
			if(line!=null && line.isEmpty()) continue;
			
			if(	line==null ||
				line.trim().startsWith("track "))
				{
				if(track!=null &&
						track.url!=null &&
						track.type!=null
						)
					{
					this.tracks.add(track);
					}
				if(line==null) break;
				track=null;
				}
			int w=line.indexOf(' ');
			if(w==-1)
				{
				LOG.warn("No whitespace in "+line+" ?");
				continue;
				}
			if(track==null) track=new Track();
			String key=line.substring(0,w).trim();
			String val=line.substring(w+1).trim();
			switch(key)
				{
				case "track":track.id=val;break;
				case "shortLabel":track.shortLabel=val;break;
				case "longLabel":track.longLabel=val;break;
				case "bigDataUrl":track.url=new URL(url,val);break;
				case "type":track.type=val;break;
				case "parent":track.parent=val.split("[ \t]")[0];break;
				default:break;
				}

			}
		
		CloserUtil.close(r);
		}

	@Override
		protected Collection<Throwable> call(String inputName) throws Exception
			{
			try
				{
				parseTrackDB(new URL(trackDBUrl));
				
				return doVcfToVcf(inputName);
				}
			catch(Exception err)
				{
				return wrapException(err);
				}
			}
		}
	
	public static void main(String[] args)
		{
		new VcfEnsemblReg().instanceMainWithExit(args);
		}
}
