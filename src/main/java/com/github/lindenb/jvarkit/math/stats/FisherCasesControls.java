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
import java.util.OptionalDouble;
import java.util.Set;
import java.util.function.Consumer;
import java.util.function.DoubleSupplier;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.pedigree.CasesControls;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;


/**
 *  Wrapper for Fisher and Case+Controls
 * */

public class FisherCasesControls implements DoubleSupplier,Consumer<String> {
	private final CasesControls casesControls;
	private final Set<String> seen = new HashSet<>();
    public FisherCasesControls(final CasesControls casesControls) {
    	this.casesControls=Objects.requireNonNull(casesControls);
    	}
    /** reset count */
    public FisherCasesControls reset() {
    	this.seen.clear();
    	return this;
    	}
    /** did we see sample sn ? */
    public boolean isSeen(final String sn) {
    	return this.seen.contains(sn);
    	}
    public void acceptAll(Collection<String> samplesNames) {
    	for(final String sn:samplesNames) accept(sn);
    	}
    @Override
    public void accept(final String sampleName) {
    	this.seen.add(sampleName);	
    	}
    
    public void accept(final VariantContext vc) {
    	for(final String sn:this.casesControls.getAll()) {
    		if(isSeen(sn)) continue;
    		final Genotype gt= vc.getGenotype(sn);
    		if(gt==null) continue;
    		this.accept(gt);
    		}
    	}
    public void accept(final Genotype gt) {
    	if(gt.hasAltAllele()) this.accept(gt.getSampleName());
    	}
    
    /** return odd ratio https://en.wikipedia.org/wiki/Odds_ratio */
    public OptionalDouble getOddRatio() {
    	final int nCases = this.casesControls.getCasesCount();
    	if(nCases==0) return OptionalDouble.empty();
    	final int nCtrls = this.casesControls.getControlsCount();
    	if(nCtrls==0) return OptionalDouble.empty();
    	final double f1 = getCasesAltCount()/(double)nCases;
    	final double f2 = getControlsAltCount()/(double)nCtrls;
    	if(f2==0) return OptionalDouble.empty();
    	return OptionalDouble.of(f1/f2);
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
    /** return name of case samples with alt */
    public Set<String> getCasesAlt() {
    	return this.seen.stream().filter(sn->this.casesControls.isCase(sn)).collect(Collectors.toSet());
    	}
    /** return name of ctrl samples with alt */
    public Set<String> getControlsAlt() {
    	return this.seen.stream().filter(sn->this.casesControls.isControl(sn)).collect(Collectors.toSet());
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
