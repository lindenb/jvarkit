package com.github.lindenb.jvarkit.tools.genscan;

import java.awt.AlphaComposite;
import java.awt.Color;
import java.awt.Composite;
import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.Insets;
import java.awt.Toolkit;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;

import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JScrollPane;


import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.Hershey;
import com.github.lindenb.jvarkit.util.cli.GetOpt;

public abstract class AbstractGeneScan extends AbstractCommandLineProgram
	{
	protected Dimension screenSize=new Dimension(1000,300);
    //@Option(shortName="DC",doc="Do not display unseen chromosome",optional=false)	
    //@Option(shortName="DC",doc="Do not display unseen chromosome",optional=false)	
	protected boolean DISCARD_UNSEEN_CHROM=false;
    //@Option(shortName="DB",doc="ignore chromosomes left/right side if there's no data there.",optional=false)	
	protected boolean DISCARD_BOUNDS=false;
    private Hershey hershey=new Hershey();
    protected Insets insets=new Insets(20,150, 30, 30);
    protected MinMaxDouble minMaxY=new MinMaxDouble();
    /**all samples in the same chart */
    private boolean merge_sample=false;
    /** user min Value */
    private Double user_min_value=null;
    /** user max Value */
    private Double user_max_value=null;
    
	protected List<Sample> samples=new ArrayList<Sample>();
	protected List<ChromInfo> chromInfos=new ArrayList<ChromInfo>();
	
	protected class MinMaxDouble
		{
		double min_value=Double.MAX_VALUE;
		double max_value=-Double.MAX_VALUE;
		MinMaxDouble()
			{
			this(Double.MAX_VALUE,-Double.MAX_VALUE);
			}
		MinMaxDouble(double m,double M)
			{
			this.min_value= m;
			this.max_value= M;
			}
		boolean isValid()
			{
			return this.min_value < this.max_value &&
					!Double.isNaN(min_value) &&
					!Double.isNaN(max_value)
					;
 			}
		public void setMin(double  min_value)
			{
			this.min_value= min_value;
			}
		public void setMax(double  max_value)
			{
			this.max_value= max_value;
			}
		public double getMax()
			{
			return max_value;
			}
		public double getMin()
			{
			return min_value;
			}
		public void visit(double v)
			{
			if(Double.isNaN(v)) return;
			this.min_value= Math.min(v,min_value);
			this.max_value= Math.max(v,max_value);
			}
		public double getAmplitude()
			{
			return getMax()-getMin();
			}
		public boolean contains(double v)
			{
			return getMin()<=v && v <= getMax();
			}
		public double getFraction(double v)
			{
			return (v-getMin())/(getMax()-getMin());
			}
		public void merge(MinMaxDouble mM)
			{
			this.visit(mM.getMin());
			this.visit(mM.getMax());
			}
		@Override
		public String toString()
			{
			return "["+getMin()+" - "+ getMax()+"]";
			}
		}
	
	protected class Sample
		{
		int sample_id=-1;
		String name="";
		MinMaxDouble minmax=new MinMaxDouble();
		double y;
		double height;
		boolean visible=true;

		}

    
    
	protected class ChromInfo 
		{
		int tid=-1;
		String sequenceName=null;
		Double dictSequenceLength=null;//sequence set via SamSequenceDict
		double x;
		double width;
		MinMaxDouble minmaxBase=new MinMaxDouble();
		boolean visible=true;
		
		private Rectangle2D.Double getFrame()
			{
			return new Rectangle2D.Double(
					x,insets.top,
					width,
					screenSize.getHeight()-(insets.top+insets.bottom)
					);
			}
		
		void draw(Graphics2D g)
			{
			Rectangle2D.Double r=getFrame();
			g.setColor(Color.GRAY);
			g.draw(r);
			
			double optW=Math.min(this.getSequenceName().length()*10,this.width);
			Rectangle2D.Double titleRec=new Rectangle2D.Double(
					(x+width/2)-optW/2.0,
					insets.top-22,
					optW,
					20);
			
			hershey.paint(g, this.getSequenceName(), titleRec);
			}
		
		String getSequenceName()
			{
			return this.sequenceName;
			}
		
		int getSequenceLength()
			{
			if(DISCARD_BOUNDS)
				{
				return (int)minmaxBase.getAmplitude();
				}
			else
				{
				return this.getSequenceLength();
				}
			}
		}
	
	
	
	protected AbstractGeneScan()
		{
		
		}
	protected abstract List<ChromInfo> getChromInfos();

	protected List<Sample> getSamples()
		{
		return this.samples;
		}
	
	protected abstract void drawPoints(Graphics2D g);
	
	
	
	protected void prepareData()
		{
		info("prepare sample data");
		info(this.minMaxY);
		int i=0;
		int num_visible_samples=samples.size();
		while(i<samples.size())
			{
			Sample sample=samples.get(i);
			if(sample.minmax.isValid())
				{
				this.minMaxY.merge(sample.minmax);
				info("merge "+sample.minmax+" "+this.minMaxY);
				++i;
				}
			else
				{
				warning("No valid data for "+sample.name);
				sample.visible=false;
				num_visible_samples--;
				}
			}
		
		info(this.minMaxY);
		
		if(this.user_min_value!=null)
			{
			this.minMaxY.setMin(this.user_min_value);
			}
		if(this.user_max_value!=null)
			{
			this.minMaxY.setMax(this.user_max_value);
			}
		
		if(!this.minMaxY.isValid())
			{
			warning("invalid y values");
			this.minMaxY.setMax(1.0);
			this.minMaxY.setMin(0.0);
			}

		if(this.minMaxY.getMin()==this.minMaxY.getMax())
			{
			warning(" y min= y max");
			this.minMaxY.setMax(this.minMaxY.getMin()+1.0);
			}
		
		info(this.minMaxY);
		
		if(this.merge_sample)
			{
			for(Sample sample:this.samples)
				{
				if(!sample.visible) continue;
				sample.y=this.insets.top;
				sample.height=this.screenSize.height-(this.insets.top+this.insets.bottom);
				}
			}
		else
			{
			double y=this.insets.top;
			double h=(this.screenSize.height-(this.insets.top+this.insets.bottom))/(double)num_visible_samples;
			
			for(Sample sample:this.samples)
				{
				if(!sample.visible) continue;
				sample.y=y;
				sample.height=h-2;
				y+=h;
				}
				
			}
		info("prepare chrom data");
		long genome_length=0L;
		i=0;
		while(i< this.chromInfos.size())
			{
			ChromInfo ci=this.chromInfos.get(i);
			if(DISCARD_UNSEEN_CHROM && !ci.minmaxBase.isValid())
				{
				warning("no data for chrom "+ci.sequenceName);
				ci.visible=false;
				continue;
				}
			if(!DISCARD_BOUNDS)
				{
				ci.minmaxBase.setMin(0);
				}
			
			if(!DISCARD_BOUNDS && ci.dictSequenceLength!=null)
				{
				ci.minmaxBase.setMax(ci.dictSequenceLength);
				}
		
			if(ci.minmaxBase.getMin()==ci.minmaxBase.getMax())
				{
				ci.minmaxBase.setMax(1+ci.minmaxBase.getMin());
				}
			
			if(!DISCARD_BOUNDS)
				{
				genome_length+=(long)ci.minmaxBase.getMax();
				}
			else
				{
				genome_length+=(long)ci.minmaxBase.getAmplitude();
				}
			++i;
			}
		double x=this.insets.left;
		final double W=this.screenSize.getWidth()-(this.insets.left+this.insets.right);
		for(i=0;i< this.chromInfos.size();++i)
			{
			ChromInfo ci=this.chromInfos.get(i);
			if(!ci.visible) continue;
			ci.x=x;
			if(!DISCARD_BOUNDS)
				{
				ci.width+=W*(ci.minmaxBase.getMax()/genome_length);
				}
			else
				{
				ci.width+=W*(ci.minmaxBase.getAmplitude()/genome_length);
				}
			x+=ci.width;
			}
		
		}

	
	protected BufferedImage makeImage()
		{
		prepareData();
		NumberFormat fmt=NumberFormat.getInstance();
		info("make image");
		BufferedImage img=new BufferedImage(
				this.screenSize.width,
				this.screenSize.height,
				BufferedImage.TYPE_INT_RGB
				);
		
		Graphics2D g=img.createGraphics();
		g.setColor(Color.WHITE);
		g.fillRect(0, 0,this.screenSize.width, this.screenSize.height);
		
		
		for(int i=0;i< getChromInfos().size();++i)
			{
			if(!getChromInfos().get(i).visible) continue;
			g.setColor(i%2==0?Color.WHITE:Color.LIGHT_GRAY);
			g.fill(getChromInfos().get(i).getFrame());
			}
			
		
		//draw axis y for each seample
		for(Sample sample:this.samples)
			{
			if(!sample.visible) continue;
			g.setColor(Color.BLACK);
			
			for(int i=0;i<= 10;++i)
				{
				double v=minMaxY.getMin()+ ((minMaxY.getAmplitude())/10.0)*i;
				double y=sample.y+sample.height-sample.minmax.getFraction(v)*sample.height;
				double font_size=Math.min(12,sample.height/(10.0+2));
				
				g.setColor(Color.LIGHT_GRAY);
				g.draw(new Line2D.Double(
						insets.left-5, y,
						this.screenSize.width-insets.right,
						y));
				
				Rectangle2D.Double r=new Rectangle2D.Double(
						0,y-font_size,
						insets.left,
						font_size
						);
				g.setColor(Color.BLACK);
				this.hershey.paint(g, fmt.format(v), r);
				}
			
			if(merge_sample) break;
			}
		/* print sample label */
		if(!merge_sample)
			{
			for(Sample sample:this.samples)
				{
				if(!sample.visible) continue;
				Rectangle2D.Double r=new Rectangle2D.Double(
						this.screenSize.width-(double)this.insets.right,
						sample.y,
						(double)insets.right,
						Math.min(12.0, sample.height)
						);
				this.hershey.paint(g,sample.name, r);
				}
			}
		float OPACITY=0.9f;
		Composite oldComposite=g.getComposite();
		g.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER,(float)OPACITY));
		drawPoints(g);
		g.setComposite(oldComposite);
		
		for(ChromInfo ci:getChromInfos())
			{
			if(!ci.visible) continue;
			ci.draw(g);
			}
		
		
		g.dispose();	
		return img;
		}
	
	@Override
	protected GetOptStatus handleOtherOptions(int c, GetOpt opt,String args[])
		{
		switch(c)
			{
			case GetOpt.LONG_OPT:
					{
					String lo=opt.getLongOpt();
					if("min-y".equals(lo))
						{
						this.user_min_value=Double.parseDouble(opt.increaseOptind(args));
						return GetOptStatus.OK;
						}
					else if("max-y".equals(lo))
						{
						this.user_max_value=Double.parseDouble(opt.increaseOptind(args));
						return GetOptStatus.OK;
						}
					break;
					}
			}
		return super.handleOtherOptions(c, opt, args);
		}
	
	protected void showGui(BufferedImage img)
		{
		JLabel label=new JLabel(new ImageIcon(img));
		Dimension d=Toolkit.getDefaultToolkit().getScreenSize();
		JScrollPane scroll=new JScrollPane(label);
		scroll.setPreferredSize(new Dimension(
				Math.min(img.getWidth(),d.width-50),
				Math.min(img.getHeight(),d.height-50))
				);
		JOptionPane.showMessageDialog(null, scroll);
		}
	}