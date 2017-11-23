package com.github.lindenb.jvarkit.tools.cloneit;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.log.Logger;

public class SubCloneIt extends AbstractCloneIt {
	private static final Logger LOG = Logger.build(SubCloneIt.class).make();
	
	@ParametersDelegate
	private InsertPlasmidLoader insertLoader=new InsertPlasmidLoader();
	@ParametersDelegate
	private VectorPlasmidLoader vectorLoader=new VectorPlasmidLoader();
	
	private float enzymeWeight=5f;
	private InsertPlasmid insert=null;
	private VectorPlasmid vector=null;
	
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

private class SearchV3
	{
	private int index=0;
	private Enzyme enz;
	private SearchV5 searchV5;
	
	List<Enzyme> getEnzymes() { return this.searchV5.enzymes;}
	
	void run(final SearchV5 searchV5) {
		this.searchV5 = searchV5;
		for(this.index=0;this.index<this.getEnzymes().size();++this.index)
			{
			this.enz = this.getEnzymes().get(this.index);
			
			vector.getSites().stream().
				filter(S->S.getEnzyme().equals(this.enz)).
				filter(S->S.getMaxPosition()<=vector.getPolylinker().getEnd()).forEach(S->{
					});

			}
		}
	}

private class SearchV5
	{
	private List<Enzyme> enzymes;
	private int index=0;
	private Enzyme enz;
	void run() {
		final SearchV3 searchv3=new SearchV3();
		for(this.index=0;this.index<this.enzymes.size();++this.index)
			{
			this.enz = this.enzymes.get(this.index);
			vector.getSites().stream().
				filter(S->S.getEnzyme().equals(enz)).
				filter(S->S.getPosition()>=vector.getPolylinker().getStart()).
				filter(S->S.getMaxPosition()<=vector.getPolylinker().getEnd()).forEach(S->{
					searchv3.run(this);
					});
			}
		}
	}

private List<Strategy> run() {
	List<Enzyme> enzymes =new ArrayList<>();
	for(final Enzyme e: enzymes)
		{
		long count_out = insert.getSites().stream().
				filter(S->S.getEnzyme().equals(e)).
				filter(S->S.getMaxPosition()< insert.getPolylinker5().getStart() /* todo || S.getPosition() > insert.getPolylinker3()*/).
				count();
		}
	
	
	vector.getSites().stream().filter(S->S.getPosition()>0).forEach(SITE_V_5->{
		
		vector.getSites().stream().filter(S->S.getPosition()>SITE_V_5.getMaxPosition() && S.getPosition()<0).forEach(SITE_V_3->{
			
			
			
			});
		
		});
	return null;//TODO
	}

@Override
public int doWork(final List<String> args) {
	if(!args.isEmpty()) {
		LOG.error("too many arguments");
		return -1;
	}
	try 
		{
		final List<Enzyme> enzymes = getRebase().stream().
				filter(E->E.isPalindromic()).
				filter(E->E.getWeight()>this.enzymeWeight).
				map(E->new EnzymeImpl(E)).collect(Collectors.toList());
		if(enzymes.isEmpty()) {
			LOG.error("enzyme list is empty");
			return -1;
		}
		this.insert  = this.insertLoader.create();
		this.vector  = this.vectorLoader.create();
		this.insert.digest(enzymes);
		this.vector.digest(enzymes);

		
		
		return 0;
		}
	catch(final Exception err) {
		return -1;
		}
	
	}

public static void main(String[] args) {
	new SubCloneIt().instanceMainWithExit(args);
	}
}
