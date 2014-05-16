package com.github.lindenb.jvarkit.util.bio.circular;

import java.awt.geom.Arc2D;
import java.awt.geom.Point2D;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

public class CircularContext
	{
	private SAMSequenceDictionary dictionary;
	private Point2D center=new Point2D.Double(0, 0);
	
	public CircularContext(SAMSequenceDictionary dictionary)
		{
		this.dictionary=dictionary;
		}
	
	public SAMSequenceDictionary getDictionary()
		{
		return dictionary;
		}
	
	public Point2D getCenter()
		{
		return center;
		}
	
	public double getCenterX()
		{
		return getCenter().getX();
		}	
	public double getCenterY()
		{
		return getCenter().getY();
		}	
	
	public void setCenter(double x,double y)
		{
		this.center = new Point2D.Double(x,y);
		}
	public double convertPositionToRadian(int tid,int pos0)
		{
		return convertPositionToRadian(getDictionary().getSequence(tid),pos0);
		}
	public double convertPositionToRadian(String chrom,int pos0)
		{
		return convertPositionToRadian(getDictionary().getSequence(chrom),pos0);
		}
	public double convertPositionToRadian(SAMSequenceRecord ssr,int pos0)
		{
		if(ssr==null) throw new IllegalArgumentException("ssr==null");
		if(pos0<0)  throw new IllegalArgumentException("pos<0");
		if(pos0>ssr.getSequenceLength())  throw new IllegalArgumentException("pos>seq.length");
		double radperpb=(2*Math.PI)/this.getDictionary().getReferenceLength();
		double a=0.0;
		for(int i=0;i< ssr.getSequenceIndex();++i)
			{
			a+=radperpb*getDictionary().getSequence(i).getSequenceLength();
			}
		a+=radperpb*pos0;
		//convert to clock direction		
		return  a-Math.PI/2.0;
		}
	public Point2D convertPositionToPoint(int tid,int pos0,double radius)
		{
		return convertPositionToPoint(getDictionary().getSequence(tid), pos0, radius);
		}

	public Point2D convertPositionToPoint(String chrom,int pos0,double radius)
		{
		return convertPositionToPoint(getDictionary().getSequence(chrom), pos0, radius);
		}
	
	public Point2D convertPositionToPoint(SAMSequenceRecord ssr,int pos0,double radius)
		{
		Point2D c=getCenter();
		double angle=convertPositionToRadian(ssr, pos0);
		return new Point2D.Double(
				c.getX()+ radius* Math.cos(angle),
				c.getY()+ radius* Math.sin(angle)
				);
		}
	
	public Arc2D getArc(String chrom,int start0,int end0,double radius,int arcType)
		{
		return getArc(getDictionary().getSequence(chrom), start0, end0, radius, arcType);
		}

	
	public Arc2D getArc(int tid,int start0,int end0,double radius,int arcType)
		{
		return getArc(getDictionary().getSequence(tid), start0, end0, radius, arcType);
		}
	
	
	public Arc2D getArc(SAMSequenceRecord ssr,int start0,int end0,double radius,int arcType)
		{
		Point2D c=getCenter();
		double a1= convertPositionToRadian(ssr, start0);
		double a2= convertPositionToRadian(ssr, end0);
		Arc2D.Double arc=new Arc2D.Double(
				c.getX()-radius,
				c.getY()-radius,
				radius*2.0,
				radius*2.0,
				-rad2deg(a1),
				-rad2deg(a2-a1),
				arcType
				);
		
		return arc;
		}
	private static final double rad2deg(double a)
		{
		return (a/Math.PI)*180.0;
		}
	
	
	}
