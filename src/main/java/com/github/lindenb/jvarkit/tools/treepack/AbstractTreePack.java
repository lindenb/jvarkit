package com.github.lindenb.jvarkit.tools.treepack;

import java.awt.geom.Rectangle2D;
import java.util.Collection;



public abstract class AbstractTreePack implements TreePack
	{
	private AbstractTreePack parent=null;
	private Rectangle2D bounds=new Rectangle2D.Double();
	
	protected AbstractTreePack(AbstractTreePack parent)
		{
		this.parent=parent;
		}
	
	public AbstractTreePack getParent()
		{
		return parent;
		}
	
	
	public int getDepth()
		{
		int d=0;
		AbstractTreePack p=this;
		while(p.parent!=null)
			{
			p=p.parent;
			++d;
			}
		return d;
		}
	
	@Override
	public void setBounds(Rectangle2D bounds)
		{
		this.bounds=bounds;
		}

	@Override
	public Rectangle2D getBounds()
		{
		return this.bounds;
		}

	@Override
	public abstract double getWeight();
	
	
	
	
	}
