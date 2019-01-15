/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.util.swing;

import java.awt.Component;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.Insets;
import java.awt.LayoutManager;

import javax.swing.JSeparator;

public class FormLayout implements LayoutManager
{
private int marginLeft=200;
private int spacingx=5;
private int spacingy=spacingx*2;
@Override
public void addLayoutComponent(String name, Component comp) {
	//ignore
	}

@Override
public void layoutContainer(Container parent)
	{
	synchronized (parent.getTreeLock())
	    {
		 Insets insets = parent.getInsets();
		 int y=insets.top;
		 final int n= parent.getComponentCount();
		 int i=0;
		 while(i<n)
		 	{
			Component c= parent.getComponent(i); 
			if(c instanceof JSeparator)
				{
				y += spacingy;
				c.setBounds(
						insets.left,
						y,
						parent.getWidth()-(insets.left+insets.right),
						spacingy/2
						);
				y += spacingy;
				++i;
				continue;
				}
			
			
			Dimension d = c.getPreferredSize();
			int rowHeight=  d.height;
			c.setBounds(
					insets.left,
					y,
					marginLeft,
					d.height
					);
			
			if(i+1<n)
				{
				i++;
				c= parent.getComponent(i); 
				d = c.getPreferredSize();
				int x= insets.left+marginLeft+spacingx;
				int width = parent.getWidth()-(x+insets.right);
				if(width<=marginLeft) width=marginLeft;
				rowHeight= Math.max(rowHeight, d.height);
				c.setBounds(
						x,
						y,
						width,
						d.height
						);
				
				if(i+1<n) y += spacingy;
				}
			y+= rowHeight;
			++i;
		 	}
		 y+=insets.bottom;
	    }
	}
@Override
public Dimension minimumLayoutSize(Container parent)
	{
	synchronized (parent.getTreeLock())
    {
	 int width=marginLeft;
	 Insets insets = parent.getInsets();
	 int y=insets.top;
	 final int n= parent.getComponentCount();
	 int i=0;
	 while(i<n)
	 	{
		Component c= parent.getComponent(i); 
		if(c instanceof JSeparator)
			{
			y += spacingy*2;
			++i;
			continue;
			}
		
		
		Dimension d = c.getPreferredSize();
		int rowHeight=  d.height;
		
		if(i+1<n)
			{
			i++;
			c= parent.getComponent(i); 
			d = c.getPreferredSize();
			rowHeight= Math.max(rowHeight, d.height);
			width = Math.max(width, marginLeft+spacingx+d.width);
			if(i+1<n) y += spacingy;
			}
		y+= rowHeight;
		++i;
	 	}
	 y+=insets.bottom;
	 return new Dimension(width+insets.left+insets.right, y);
    }
	}
@Override
public Dimension preferredLayoutSize(Container parent) {
	return minimumLayoutSize(parent);
	}
@Override
public void removeLayoutComponent(Component parent) {
	//ignore
	}
}