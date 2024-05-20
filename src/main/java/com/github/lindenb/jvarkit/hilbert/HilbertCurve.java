/*
The MIT License (MIT)

Copyright (c) 2020 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.hilbert;

import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.List;

public class HilbertCurve {
	private final double squareWidth;
	private final int recursionLevel;
	private final long maximumLength;
	private double genomicSizePerCurveUnit=1.0;
	
	public HilbertCurve(
			final double squareWidth,
			final long maximumLength,
			int recursionLevel
			)
		{
		this.squareWidth = squareWidth;
		this.recursionLevel = recursionLevel;
		this.maximumLength = maximumLength;
		
		// get maximum number of segments 
		class Count implements Handler {
			//double curveLength=0.0;
			int count=0;
			@Override
			public void segment(Point2D.Double beg, Point2D.Double end, long chromStart, long chromEnd) {
				//curveLength+= beg.distance(end);
				count++;
				}
			}
		final Count evalCurve=new Count();
		this.visit(evalCurve);
		
		
		this.genomicSizePerCurveUnit = maximumLength/evalCurve.count;		
		}
	
	/** return maximum length */
	public long size() {
		return maximumLength;
		}
	
	public static interface Handler
		{
		public default void initialize() {
			}
		public void segment(
				Point2D.Double beg,
	    		Point2D.Double end,
	    		long chromStart,
	    		long chromEnd
	    		);
		public default void finish() {
			}
		}
	
	public List<Point2D.Double> getPoints(long start,long end) {
		return visit(new PointsCollectors(start, end)).points;
		}
	
	
	public <T extends Handler> T visit(final T handler)
		{
		final Visitor visitor=new Visitor();
		handler.initialize();
		visitor.visit(handler);
		handler.finish();
		return handler;
		}
		
	
    private class Visitor
		{
	    /** last index of the position in the genome */
		double prev_base=0;
	    /** size (in pb) of an edge */
	    private final double d_base=HilbertCurve.this.genomicSizePerCurveUnit;
		/** size of an edge */
	    private double dist=-1;
	    /** last time we plot a segment, point at end */
	    private Point2D.Double prevPoint;
	    /** curve length so far */
	   // private  double curveLength=0.0;
	    
	    
	    Visitor()
	    	{	    	
	    	this.dist = HilbertCurve.this.squareWidth;
	    	for (int i= HilbertCurve.this.recursionLevel; i>0; i--)
	            {
	            this.dist /= 2.0;
	            }
	    	this.prevPoint=new Point2D.Double(
	    			this.dist/2.0,
	    			this.dist/2.0
	    			);
	        this.prev_base=0L;
	    	}
	    
	   private void segment(
	    		final Point2D.Double beg,
	    		final Point2D.Double end,
	    		long chromStart,
	    		long chromEnd,
	    		Handler h
	    		)
		   	{
		   	//this.curveLength += beg.distance(end);
		    if(h!=null) h.segment(beg, end, chromStart, chromEnd);
		   	}
	    
	   private Handler visit(final Handler handler)
	    	{
	    	this.HilbertU(recursionLevel,handler);
	    	return handler;
	    	}	
	    
	    private void lineRel(double deltaX, double deltaY,final Handler handler)
	    	{
	    	final Point2D.Double point=new Point2D.Double(
	    		this.prevPoint.x + deltaX,
	    		this.prevPoint.y + deltaY
	    		);
	        final long chromStart=(long)prev_base;
	        final long chromEnd=(long)(chromStart+ d_base);
	        segment(prevPoint, point, chromStart, chromEnd,handler);
	        this.prevPoint= point;
	        this.prev_base= chromEnd;;
	    	}
		//make U shaped curve       
	    private void HilbertU(int level,final Handler handler)
	        {
	        if (level <= 0) return;
	        HilbertD(level-1,handler);    this.lineRel(0, dist,handler);
	        HilbertU(level-1,handler);    this.lineRel(dist, 0,handler);
	        HilbertU(level-1,handler);    this.lineRel(0, -dist,handler);
	        HilbertC(level-1,handler);
	        }
	 
	    //make D shaped rule
	    private void HilbertD(int level,final Handler handler)
	    	{
	    	if (level <= 0) return;
	        HilbertU(level-1,handler);    this.lineRel(dist, 0,handler);
	        HilbertD(level-1,handler);    this.lineRel(0, dist,handler);
	        HilbertD(level-1,handler);    this.lineRel(-dist, 0,handler);
	        HilbertA(level-1,handler);
			}
	 
	    // make C shaped
	    private void HilbertC(int level,final Handler handler)
	    	{
	    	if (level <= 0) return;
	        HilbertA(level-1,handler);    this.lineRel(-dist, 0,handler);
	        HilbertC(level-1,handler);    this.lineRel(0, -dist,handler);
	        HilbertC(level-1,handler);    this.lineRel(dist, 0,handler);
	        HilbertU(level-1,handler);
	    	}
	 	//make A shaped
	    private void HilbertA(int level,final Handler handler) {
	    	if (level <= 0) return;
	        HilbertC(level-1,handler);    this.lineRel(0, -dist,handler);
	        HilbertA(level-1,handler);    this.lineRel(-dist, 0,handler);
	        HilbertA(level-1,handler);    this.lineRel(0, dist,handler);
	        HilbertD(level-1,handler);
	    	}
		}
    
    private static class PointsCollectors implements HilbertCurve.Handler {
    	final long genomicStart;
    	final long genomicEnd;
    	final List<Point2D.Double> points = new ArrayList<>();
    	PointsCollectors(long genomicStart, long genomicEnd){
    		this.genomicStart = genomicStart;
    		this.genomicEnd = genomicEnd;
    		}
    	    	
    	@Override
    	public void segment(
         		final Point2D.Double beg,
         		final Point2D.Double end,
         		final long chromStart,
         		final long chromEnd
         		)
    		{
    		if(this.genomicStart > chromEnd) return;
    		if(this.genomicEnd < chromStart) return;
    		
    		long p0= Math.max(this.genomicStart,chromStart);
    		long p1= Math.min(this.genomicEnd,chromEnd);
    		double segsize=chromEnd-chromStart;
    		if(segsize==0) segsize=1.0;
    		
    		double r0= (p0-(double)chromStart)/segsize;
    		double r1= (p1-(double)chromStart)/segsize;
    		
    		double x0= beg.x + (end.x - beg.x)*r0;
    		double y0= beg.y + (end.y - beg.y)*r0;
    		double x1= beg.x + (end.x - beg.x)*r1;
    		double y1= beg.y + (end.y - beg.y)*r1;
    		
    		if(points.isEmpty()) {
    			points.add(new Point2D.Double(x0, y0));
    			}
    		points.add(new Point2D.Double(x1, y1));
    		}
    	}

	}
