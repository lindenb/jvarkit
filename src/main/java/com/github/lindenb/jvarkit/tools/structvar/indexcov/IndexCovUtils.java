/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.structvar.indexcov;

/** utilies for indexcov */
public class IndexCovUtils {
public static final String TRESHOLD_OPT_DESC = "DUP if 1.5-x<=depth<=1.5+x . HET_DEL if 0.5-x<=depth<=0.5+x HOM_DEL if 0.0-x<=depth<=0.0+x.";
public static final double DEFAULT_TRESHOLD=0.05;
private final double treshold;
public enum SvType {HOM_DEL(0),HET_DEL(0.5),REF(1.0),HET_DUP(1.5),HOM_DUP(2.0),AMBIGOUS(-1);
	final double theoritical;
	SvType(final double theoritical) {
		this.theoritical = theoritical;
		}
	public double getTheoriticalDepth() {
		return theoritical;
		}
	public int getGenotypeQuality(final double normDepth) {
		if(this.equals(AMBIGOUS)) return 0;
		double gq = Math.abs(this.getTheoriticalDepth()-normDepth);
		gq = Math.min(0.5, gq);
		gq = gq * gq;
		gq = gq / 0.25;
		gq = 99 * (1.0 - gq);
		return (int)gq;
		}
	
	public boolean isAmbigous() {
		return this.equals(AMBIGOUS);
		}
	public boolean isVariant() {
		switch(this) {
			case REF: return false;
			case AMBIGOUS: return false;
			default: return true;
			}
		}
	public boolean isDeletion() {
		switch(this) {
			case HOM_DEL: return true;
			case HET_DEL: return true;
			default: return false;
			}
		}
	public boolean isDuplication() {
		switch(this) {
			case HOM_DUP: return true;
			case HET_DUP: return true;
			default: return false;
			}
		}
	}
public IndexCovUtils(final double treshold) {
	this.treshold = treshold;
	if(this.treshold<0f || this.treshold>=0.25f) {
		throw new IllegalArgumentException("Bad treshold 0 < "+this.treshold+" >=0.25 ");
		}
	}

public double getTreshold() {
	return treshold;
	}
private boolean testLimit(double normDepth,double limit) {
	return Math.abs(normDepth-limit)<= getTreshold();
	}
public SvType getType(double normDepth) {
	if(isHomDel(normDepth)) return SvType.HOM_DEL;
	if(isHetDel(normDepth)) return SvType.HET_DEL;
	if(isHomDup(normDepth)) return SvType.HOM_DUP;
	if(isHetDup(normDepth)) return SvType.HET_DUP;
	if(isRef(normDepth)) return SvType.REF;
	return SvType.AMBIGOUS;
	}
public boolean isHomDel(double normDepth) { return testLimit(normDepth,0.0);}
public boolean isHetDel(double normDepth) { return testLimit(normDepth,0.5) || (!isHomDel(normDepth) && normDepth<=0.5);}
public boolean isHomDup(double normDepth) { return testLimit(normDepth,2.0) || normDepth>2.0 ;}
public boolean isHetDup(double normDepth) { return testLimit(normDepth,1.5) || (!isHomDup(normDepth) && normDepth>=1.5) ;}
public boolean isRef(double normDepth) { return testLimit(normDepth,1.0);}

public boolean isVariant(double normDepth) {
	return getType(normDepth).isVariant();
	}
}
