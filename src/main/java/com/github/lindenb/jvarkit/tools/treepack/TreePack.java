package com.github.lindenb.jvarkit.tools.treepack;
import java.awt.geom.Rectangle2D;

public interface TreePack
	{
	public void setBounds(Rectangle2D bounds);
	public Rectangle2D getBounds();
	public double getWeight();
	}
