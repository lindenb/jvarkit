/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.math;

/**
 * Nice Scale for X/Y axis 
 * from https://stackoverflow.com/questions/8506881/
 * 
 * @author lindenb
 *
 */
public class NiceScale {
	private static final int DEFAULT_NTICKS=10;
    private final int nticks;
    private final double tickSpacing;
    private final double niceMin;
    private final double niceMax;

    public NiceScale(final MinMaxDouble minmax,int nticks) {
    	this(
    		minmax.isEmpty()?0:minmax.getMinAsDouble(),
    		minmax.isEmpty()?1:minmax.getMaxAsDouble(),
    		nticks
    		);
    	}
    public NiceScale(final MinMaxDouble minmax) {
    	this(minmax,DEFAULT_NTICKS);
    	}
    
    public NiceScale(final double min_value, final double max_value) {
    	this(min_value,max_value,DEFAULT_NTICKS);
    	}
    
    public NiceScale(double min,double max,int nticks) {
        this.nticks=nticks;
        if(nticks<1) throw new IllegalArgumentException("nticks<1");
        if(min>max) throw new IllegalArgumentException("min>max");
        if(min==max) {
        	min-=0.5;
        	max+=0.5;
        	}
        
        final double range       = niceNum(max - min, false);
        tickSpacing = niceNum(range / ((double)this.nticks - 1.0), true);
        niceMin     = Math.floor(min / tickSpacing) * tickSpacing;
        niceMax     = Math.ceil(max / tickSpacing) * tickSpacing;
    	}

    private static double niceNum(final double RANGE, final boolean ROUND) {
        double exponent;     // exponent of RANGE
        double fraction;     // fractional part of RANGE
        double niceFraction; // nice, rounded fraction

        exponent = Math.floor(Math.log10(RANGE));
        fraction = RANGE / Math.pow(10, exponent);

        if (ROUND) {
            if (fraction < 1.5)
                niceFraction = 1;
            else if (fraction < 3)
                niceFraction = 2;
            else if (fraction < 7)
                niceFraction = 5;
            else
                niceFraction = 10;
        } else {
            if (fraction <= 1)
                niceFraction = 1;
            else if (fraction <= 2)
                niceFraction = 2;
            else if (fraction <= 5)
                niceFraction = 5;
            else
                niceFraction = 10;
        }
        return niceFraction * Math.pow(10, exponent);
    	}


    public double getTickSpacing() { return tickSpacing; }

    public double getMin() { return niceMin; }

    public double getMax() { return niceMax; }
    
    /** return number of ticks */
    public int size()  {
    	return this.nticks;
    }
    
    /** return value of the n-th tick */
    public double get(int i) {
    	if(i<0 || i>=size()) throw new IndexOutOfBoundsException("i="+i+"/"+size());
    	return getMin()+getTickSpacing()*i;
    	}

    
    public double[] getTicks() {
    	final double[] array = new double[size()];
    	for(int i=0;i< array.length;i++) {
    		array[i]=get(i);
    		}
    	return array;
    	}
    
    @Override
    public String toString() {
    	return "NiceScale("+getMin()+","+getMax()+")";
    	}
}
