package com.github.lindenb.jvarkit.tools.bam2graphics;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.Shape;
import java.awt.geom.GeneralPath;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Function;
import java.util.function.Predicate;

import javax.imageio.ImageIO;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.Hershey;
import com.github.lindenb.jvarkit.util.bio.IntervalParser;
import com.github.lindenb.jvarkit.util.bio.samfilter.SamFilterParser;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;
import com.github.lindenb.jvarkit.util.swing.ColorUtils;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public abstract class AbstractBam2Raster extends Launcher{
	protected static final Color ALMOST_BLACK = new Color(20,20,20);
	protected static final Color ALMOST_WHITE = new Color(240,240,240);

	private static final Logger LOG = Logger.build(AbstractBam2Raster.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	protected File outputFile = null;
	@Parameter(names={"--groupby"},description="Group Reads by")
	protected SAMRecordPartition groupBy=SAMRecordPartition.sample;
	@Parameter(names={"-r","--region"},description="Restrict to that region. "+IntervalParser.OPT_DESC,required=true)
	protected String regionStr = null;
	@Parameter(names={"-w","--width"},description="Image width")
	protected int WIDTH = 1000 ;
	@Parameter(names={"-clip","--clip"},description="Show clipping")
	protected boolean showClip=false;
	@Parameter(names={"--limit","--maxrows"},description="Limit number of rows to 'N' lines. negative: no limit.")
	protected int maxRows=-1;
	@Parameter(names={"-depth","--depth"},description="Depth track height.")
	protected int depthSize=100;
	@Parameter(names={"-srf","--samRecordFilter"},description=SamFilterParser.FILTER_DESCRIPTION,converter=SamFilterParser.StringConverter.class)
	protected SamRecordFilter samRecordFilter = SamFilterParser.buildDefault();
	@Parameter(names={"-minh","--minh"},description="Min. distance between two reads.")
	protected int minDistanceBetweenPairs=2;
	@Parameter(names={"--spaceyfeature"},description="number of pixels between features")
	protected int spaceYbetweenFeatures=4;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION)
	protected File referenceFile = null;
	@Parameter(names={"--highlight"},description="hightligth those positions.",converter=com.beust.jcommander.converters.IntegerConverter.class)
	protected Set<Integer> highlightPositions = new HashSet<>() ;
	@Parameter(names={"-V","--variants","--vcf"},description="VCF files used to fill the position to hightlight with POS")
	protected List<String> variants=new ArrayList<>();
	private enum AlphaHandler {all_opaque,handler1}
	@Parameter(names={"--mapqopacity"},description="How to handle the MAPQ/ opacity of the reads. all_opaque: no opacity, handler 1: transparency under MAPQ=60")
	private AlphaHandler alpha_handler=AlphaHandler.handler1;

	protected final Function<SAMRecord,Color> samRecord2color = new ColorUtils.SAMRecordColorExtractor();

	
	protected final Function<SAMRecord,Float> samRecord2alpha = R->{
		final int mapq = R.getMappingQuality();
		switch(alpha_handler)
			{
			case handler1:
				if( mapq ==  SAMRecord.UNKNOWN_MAPPING_QUALITY) return 0.5f;
				if( mapq >= 60 ) return 1f;
				if(mapq >50 && mapq <= 60) return 0.9f;
				if(mapq >40 && mapq <= 50) return 0.8f;
				if(mapq >30 && mapq <= 40) return 0.7f;
				if(mapq >20 && mapq <= 30) return 0.6f;
				if(mapq >10 && mapq <= 20) return 0.5f;
				return 0.4f;
			case all_opaque: return 1f;
			default: return 1f;
			}
		
	};
	
	
	protected Interval interval=null;
	protected IndexedFastaSequenceFile indexedFastaSequenceFile=null;
	protected final Hershey hersheyFont=new Hershey();

	
	protected final Function<SAMRecord, Integer> readLeft = R->showClip?R.getUnclippedStart():R.getAlignmentStart();
	protected final Function<SAMRecord, Integer> readRight = R->showClip?R.getUnclippedEnd():R.getAlignmentEnd();
	protected final Function<Integer, Double> pos2pixel = POS->((double)(POS - this.interval.getStart())/(double)this.interval.length()) * this.WIDTH;
	protected final Function<Integer, Integer> pixel2pos = PIX-> interval.getStart()+(int)(((double)(PIX)/(double)WIDTH) * (double)this.interval.length());
	protected final Function<SAMRecord, Double> left2pixel = R->pos2pixel.apply(readLeft.apply(R));
	protected final Function<SAMRecord, Double> right2pixel = R->pos2pixel.apply(readRight.apply(R));
	protected final Predicate<Locatable> isInInterval = LOC->{
		if(!LOC.getContig().equals(this.interval.getContig())) return false;
		if(LOC.getEnd()<this.interval.getStart()) return false;
		if(LOC.getStart()>this.interval.getEnd()) return false;
		return true;
		};

	protected void saveImages(final Collection<BufferedImage> imgs) throws IOException 
		{
		int image_width= imgs.stream().mapToInt(P->P.getWidth()).max().getAsInt();
		int image_height= imgs.stream().mapToInt(P->P.getHeight()).sum();

		final BufferedImage img= new BufferedImage(image_width, image_height, BufferedImage.TYPE_INT_RGB);
		final Graphics2D g=img.createGraphics();
		g.setRenderingHint(
				RenderingHints.KEY_RENDERING,
				RenderingHints.VALUE_RENDER_QUALITY
				);

		int y=0;
		for(final BufferedImage subImg : imgs) {
			g.drawImage(subImg,0,y,null);
			y+=subImg.getHeight();
			}
		g.dispose();

		if(this.outputFile==null)
			{
			ImageIO.write(img, "PNG", stdout());
			}
		else
			{
			LOG.info("saving to "+this.outputFile);
			final String format=(this.outputFile.getName().toLowerCase().endsWith(".png")?"PNG":"JPG");
			ImageIO.write(img,format, this.outputFile);
			}
		}
	protected void loadVCFs() {
		for(final String vcfFile: IOUtils.unrollFiles(variants)) {
			final VCFFileReader vcfFileReader = new VCFFileReader(new File(vcfFile),true);
			final CloseableIterator<VariantContext> r=vcfFileReader.query(this.interval.getContig(), this.interval.getStart(), this.interval.getEnd());
			while(r.hasNext())
				{
				this.highlightPositions.add(r.next().getStart());
				}
			r.close();
			vcfFileReader.close();
			}
		}
	protected  Shape createTriange(double cx,double cy,double r,double angle)
		{
		GeneralPath gp = new GeneralPath();
		for(int i=0;i< 3;++i)
			{
			double x= cx + Math.cos(angle)*r;
			double y= cy + Math.sin(angle)*r;
			if(i==0) 
				{
				gp.moveTo(x, y);
				}
			else
				{
				gp.lineTo(x, y);
				}
			angle+=(Math.PI*2.0)/3.0;
			}
		gp.closePath();
		return gp;
		}
}
