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
package com.github.lindenb.jvarkit.util.hershey;

import java.util.List;


import javafx.scene.canvas.GraphicsContext;

/** moved JFX to its own package to allow tools to compile without oracle */
public class JfxHershey extends Hershey {
	public void paint(
			final GraphicsContext g,
			final String s,
			final double x, final double y,
			final double width, 
			final double height
			)
		{
		if(s.isEmpty() || width==0 || height==0) return;
		
		double dx=width/s.length();
		for(int i=0;i < s.length();++i)
			{
			final List<PathOp> array=charToPathOp(s.charAt(i));
			
			for(int n=0;n< array.size();++n)
				{
				PathOp p2=transform(array.get(n));
				if(p2.operator==Op.MOVETO) continue;
				PathOp p1=transform(array.get(n-1));
				
				double x1=(p1.x/this.getScaleX())*dx + x+dx*i +dx/2.0;
				double y1=(p1.y/this.getScaleY())*height +y +height/2.0 ;
				
				double x2=(p2.x/this.getScaleX())*dx + x+dx*i +dx/2.0;
				double y2=(p2.y/this.getScaleY())*height +y +height/2.0 ;
				
				g.strokeLine(x1, y1, x2, y2);
				}
			}
		}
}
