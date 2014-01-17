package com.github.lindenb.jvarkit.tools.treepack;

import java.awt.geom.Rectangle2D;



public abstract class AbstractTreePackCommandLine implements TreePack
	{
	private AbstractTreePackCommandLine parent=null;
	private Rectangle2D bounds=new Rectangle2D.Double();
	
	protected AbstractTreePackCommandLine(AbstractTreePackCommandLine parent)
		{
		this.parent=parent;
		}
	
	public AbstractTreePackCommandLine getParent()
		{
		return parent;
		}
	
	
	public int getDepth()
		{
		int d=0;
		AbstractTreePackCommandLine p=this;
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
