package com.github.lindenb.jvarkit.util.jcommander;

import java.lang.reflect.Field;
import java.lang.reflect.Type;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import com.beust.jcommander.JCommander;

public interface JCommanderVisit {
	public JCommanderVisit visit(final Object...array) ;

	public static JCommanderVisit of(final JCommander cmd) {
		final JCommanderVisit v = new JCommanderVisit() {	
		private final Collection<Object> visited = new HashSet<>();
		@Override
		public JCommanderVisit visit(final Object...array) {
			if(array==null) return this;
			for(final Object o:array) {
				if(o==null) return this;
				if(this.visited.stream().anyMatch(O->O==o)) return this;
				cmd.addObject(o);
				try {
					for(final Field fld:o.getClass().getFields()) {
						Object o2=fld.get(o);
						if(o2==null) continue;
						if(o2 instanceof JCommanderVisit) {
							this.visit(o2);
						}
					}
					
				} catch(Throwable err) {
					
				}
				
				
				}
			return this;
			}
		};
	return v;
	}
}
