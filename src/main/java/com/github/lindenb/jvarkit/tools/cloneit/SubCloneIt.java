package com.github.lindenb.jvarkit.tools.cloneit;

import java.util.ArrayList;
import java.util.List;

public class SubCloneIt extends AbstractCloneIt {
private class HemiStrategy
	{
	SiteWithTreatment left;
	SiteWithTreatment right;
	}

private class Strategy
	{
	HemiStrategy left;
	HemiStrategy riht;
	}
private VectorPlasmid vector;
private InsertPlasmid insert;

private List<Strategy> run() {
	List<Enzyme> enzymes =new ArrayList<>();
	for(final Enzyme e: enzymes)
		{
		long count_out = insert.getSites().stream().
				filter(S->S.getEnzyme().equals(e)).
				filter(S->S.getMaxPosition()< insert.getPolylinker5().getStart() || S.getPosition() > insert.getPolylinker3()).
				count();
		}
	
	
	vector.getSites().stream().filter(S->S.getPosition()>0).forEach(SITE_V_5->{
		
		vector.getSites().stream().filter(S->S.getPosition()>SITE_V_5.getMaxPosition() && S.getPosition()<0).forEach(SITE_V_3->{
			
			
			
			});
		
		});
	
	}

}
