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
package com.github.lindenb.jvarkit.variant.swing;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Vector;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.pedigree.Pedigree;
import com.github.lindenb.jvarkit.pedigree.Sample;
import com.github.lindenb.jvarkit.pedigree.Status;
import com.github.lindenb.jvarkit.util.swing.AbstractGenericTable;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;

@SuppressWarnings("serial")
public class SwingVCFGenotypesTableModel extends AbstractGenericTable<Genotype>{
private final Vector<ColumnInfo> columns = new Vector<>();
private final Pedigree pedigree;
private class ColumnInfo {
	String name;
	Class<?> clazz;
	Function<Genotype, Object> extractor;
	}

public SwingVCFGenotypesTableModel(final Pedigree pedigree) {
	super();
	this.pedigree=pedigree;
	}

public void setVariant(final VariantContext ctx) {
final List<Genotype> genotypes;
if(ctx==null || !ctx.hasGenotypes()) {
	genotypes = Collections.emptyList();
	}
else {
	genotypes = ctx.getGenotypes();
	}
setRows(genotypes);
}

@Override
synchronized public void  setRows(final List<Genotype> genotypes) {
	if(genotypes==null) {
	setRows(Collections.emptyList());
		}
    this.columns.clear();
	    
	    final Vector<ColumnInfo> columns = new Vector<>();
	    ColumnInfo ci = new ColumnInfo();
	    ci.clazz = String.class;
	    ci.name= "Sample";
	    ci.extractor = GT->GT.getSampleName();
	    columns.add(ci);
	    if(this.pedigree!=null) {
		    ci = new ColumnInfo();
	        ci.clazz = String.class;
	        ci.name= "Status";
	        ci.extractor = GT->{
	        		final Sample sn= this.pedigree.getSampleById(GT.getSampleName());
	        		if( sn==null || !sn.isStatusSet()) return null;
	        		if(sn.getStatus().equals(Status.affected)) return "[*]";
	        		if(sn.getStatus().equals(Status.unaffected)) return "[ ]";
	        		return null;
	        	};
	        columns.add(ci);
	    	}
	    
	    
	    if(genotypes.stream().anyMatch(G->G.isAvailable())) {
	    	ci = new ColumnInfo();
	        ci.clazz = String.class;
	        ci.name= VCFConstants.GENOTYPE_KEY;
	        ci.extractor = GT->GT.getAlleles().stream().
	        		map(A->A.getDisplayString()).
	        		collect(Collectors.joining(GT.isPhased()?Genotype.PHASED_ALLELE_SEPARATOR:Genotype.UNPHASED_ALLELE_SEPARATOR));
	        columns.add(ci);
	    }
	    
	    ci = new ColumnInfo();
        ci.clazz = String.class;
        ci.name= "Type";
        ci.extractor = GT->GT.getType().name();
        columns.add(ci);
	    
	    if(genotypes.stream().anyMatch(G->G.hasDP())) {
	    	ci = new ColumnInfo();
	        ci.clazz = String.class;
	        ci.name= VCFConstants.DEPTH_KEY;
	        ci.extractor = GT->GT.hasDP()?GT.getDP():null;
	        columns.add(ci);
	    }
	    if(genotypes.stream().anyMatch(G->G.hasGQ())) {
	    	ci = new ColumnInfo();
	        ci.clazz = String.class;
	        ci.name= VCFConstants.GENOTYPE_QUALITY_KEY;
	        ci.extractor = GT->GT.hasGQ()?GT.getGQ():null;
	        columns.add(ci);
	    }
	    if(genotypes.stream().anyMatch(G->G.hasAD())) {
	    	ci = new ColumnInfo();
	        ci.clazz = String.class;
	        ci.name= VCFConstants.GENOTYPE_ALLELE_DEPTHS;
	        ci.extractor = GT->GT.hasAD()?Arrays.stream(GT.getAD()).mapToObj(P->String.valueOf(P)).collect(Collectors.joining(",")):null;
	        columns.add(ci);
	    }
	    if(genotypes.stream().anyMatch(G->G.hasPL())) {
	    	ci = new ColumnInfo();
	        ci.clazz = String.class;
	        ci.name= VCFConstants.GENOTYPE_PL_KEY;
	        ci.extractor = GT->GT.hasPL()?Arrays.stream(GT.getPL()).mapToObj(P->String.valueOf(P)).collect(Collectors.joining(",")):null;
	        columns.add(ci);
	    }
	    for(final String att: genotypes.stream().
	    		flatMap(G->G.getExtendedAttributes().keySet().stream()).
	    		collect(Collectors.toSet())) {
	    	ci = new ColumnInfo();
	        ci.clazz = Object.class;
	        ci.name= att;
	        ci.extractor = GT->GT.hasExtendedAttribute(att)?GT.getExtendedAttribute(att):null;
	        columns.add(ci);
	    	}
	    this.columns.addAll(columns);
		
    
    super.rows.clear();
    super.rows.addAll(genotypes);
	fireTableStructureChanged();
	}

@Override
public int getColumnCount() {
	return this.columns.size();
	}

@Override
public Class<?> getColumnClass(int column) {
	return this.columns.get(column).clazz;
	}
@Override
public String getColumnName(int column) {
	return this.columns.get(column).name;
	}

@Override
public Object getValueOf(final Genotype F, int column) {
	return this.columns.get(column).extractor.apply(F); 
	}
}
