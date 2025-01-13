/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.pedigree.Pedigree;
import com.github.lindenb.jvarkit.pedigree.Sample;
import com.github.lindenb.jvarkit.pedigree.Status;
import com.github.lindenb.jvarkit.swing.AbstractGenericTableModel;
import com.github.lindenb.jvarkit.swing.ColumnDef;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;

@SuppressWarnings("serial")
public class SwingVCFGenotypesTableModel extends AbstractGenericTableModel<Genotype>{
private final Vector<ColumnDef<Genotype>> columns = new Vector<>();
private final Pedigree pedigree;

public SwingVCFGenotypesTableModel() {
	this(null);
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
setGenotypes(genotypes);
}

public void setGenotypes(final List<Genotype> genotypes) {
	setRows(genotypes==null?Collections.emptyList():genotypes);
}

@Override
synchronized public void  setRows(final List<Genotype> genotypes) {
	if(genotypes==null) {
	setRows(Collections.emptyList());
		}
    this.columns.clear();
	    
	    final Vector<ColumnDef<Genotype>> columns = new Vector<>();
	    ColumnDef<Genotype> ci = new ColumnDef<Genotype>("Sample",String.class,GT->GT.getSampleName());
	    columns.add(ci);
	    if(this.pedigree!=null) {
		    ci = new ColumnDef<Genotype>(
	    		"Status",
	    		String.class,
		    	GT->{
	        		final Sample sn= this.pedigree.getSampleById(GT.getSampleName());
	        		if( sn==null || !sn.isStatusSet()) return null;
	        		if(sn.getStatus().equals(Status.affected)) return "[*]";
	        		if(sn.getStatus().equals(Status.unaffected)) return "[ ]";
	        		return null;
		    	});
	        columns.add(ci);
	    	}
	    
	    
	    if(genotypes.stream().anyMatch(G->G.isAvailable())) {
	    	ci = new ColumnDef<Genotype>(
	    			VCFConstants.GENOTYPE_KEY,
	    			String.class,
	    			GT->GT.getAlleles().stream().
	        		map(A->A.getDisplayString()).
	        		collect(Collectors.joining(GT.isPhased()?Genotype.PHASED_ALLELE_SEPARATOR:Genotype.UNPHASED_ALLELE_SEPARATOR))
	    			);
	    	columns.add(ci);
	    }
	    
	    ci = new ColumnDef<Genotype>("Type",String.class,GT->GT.getType().name());
        columns.add(ci);
	    
	    if(genotypes.stream().anyMatch(G->G.hasDP())) {
	    	ci = new ColumnDef<Genotype>(VCFConstants.DEPTH_KEY, String.class, GT->GT.hasDP()?GT.getDP():null);
	        columns.add(ci);
	    }
	    if(genotypes.stream().anyMatch(G->G.hasGQ())) {
	    	ci = new ColumnDef<Genotype>(VCFConstants.GENOTYPE_QUALITY_KEY, String.class, GT->GT.hasGQ()?GT.getGQ():null);
	        columns.add(ci);
	    }
	    if(genotypes.stream().anyMatch(G->G.hasAD())) {
	    	ci = new ColumnDef<Genotype>(VCFConstants.GENOTYPE_ALLELE_DEPTHS, String.class, GT->GT.hasAD()?Arrays.stream(GT.getAD()).mapToObj(P->String.valueOf(P)).collect(Collectors.joining(",")):null);
	        columns.add(ci);
	    }
	    if(genotypes.stream().anyMatch(G->G.hasPL())) {
	    	ci = new ColumnDef<Genotype>(VCFConstants.GENOTYPE_PL_KEY,  String.class, GT->GT.hasPL()?Arrays.stream(GT.getPL()).mapToObj(P->String.valueOf(P)).collect(Collectors.joining(",")):null);
	        columns.add(ci);
	    }
	    for(final String att: genotypes.stream().
	    		flatMap(G->G.getExtendedAttributes().keySet().stream()).
	    		collect(Collectors.toSet())) {
	    	ci = new ColumnDef<Genotype>(att, Object.class, GT->GT.hasExtendedAttribute(att)?GT.getExtendedAttribute(att):null);
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
	return this.columns.get(column).getColumnClass();
	}
@Override
public String getColumnName(int column) {
	return this.columns.get(column).getColumnName();
	}

@Override
public Object getValueOf(final Genotype F, int column) {
	return this.columns.get(column).getExtractor().apply(F); 
	}
}
