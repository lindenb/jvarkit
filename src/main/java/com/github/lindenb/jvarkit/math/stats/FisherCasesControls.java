/*
The MIT License (MIT)

Copyright (c) 2024 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.math.stats;

import java.util.Collection;
import java.util.HashSet;
import java.util.Objects;
import java.util.Set;
import java.util.function.Consumer;
import java.util.function.DoubleSupplier;

import com.github.lindenb.jvarkit.pedigree.CasesControls;

import htsjdk.variant.variantcontext.Genotype;


/* http://lh3lh3.users.sourceforge.net/fisher.shtml 
 * https://github.com/molgenis/systemsgenetics/blob/master/genetica-libraries/src/main/java/umcg/genetica/math/stats/FisherExactTest.java
 * 
 * 
 * */

public class FisherCasesControls implements DoubleSupplier,Consumer<String> {
	private final CasesControls casesControls;
	private final Set<String> seen = new HashSet<>();
    public FisherCasesControls(final CasesControls casesControls) {
    	this.casesControls=Objects.requireNonNull(casesControls);
    	}
    public FisherCasesControls reset() {
    	this.seen.clear();
    	return this;
    	}
    
    public boolean isSeen(final String sn) {
    	return this.seen.contains(sn);
    	}
    public void acceptAll(Collection<String> samplesNames) {
    	for(final String sn:samplesNames) accept(sn);
    	}
    @Override
    public void accept(String sampleName) {
    	this.seen.add(sampleName);
    	}

    public void accept(Genotype gt) {
    	if(gt.hasAltAllele()) this.accept(gt.getSampleName());
    	}
    
    public int getCasesAltCount() {
    	return (int)this.seen.stream().filter(sn->this.casesControls.isCase(sn)).count();
    	}
    public int getControlsAltCount() {
    	return (int)this.seen.stream().filter(sn->this.casesControls.isControl(sn)).count();
    	}
    public int getCasesRefCount() {
    	return this.casesControls.getCasesCount()-getCasesAltCount();
    	}
    public int getControlsRefCount() {
    	return this.casesControls.getControlsCount()- getControlsAltCount();
    	}
    
    public FisherExactTest getFisherExactTest() {
    	final FisherExactTest fisher = FisherExactTest.compute(
    			getCasesAltCount(),getCasesRefCount(),
    			getControlsAltCount(),getControlsRefCount()
				);
    	return fisher;
    	}
    
    @Override
    public double getAsDouble() {
    	return getFisherExactTest().getAsDouble();
    	}
	}
