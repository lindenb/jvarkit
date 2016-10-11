package com.github.lindenb.jvarkit.util.command;

import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class Arguments {
@Argument(opt="")
private final List<Object> objects=new ArrayList<>();

public Arguments register(final Object...objs) {
	if(objs==null || objs.length==0) return this;
	for(final Object obj:objs)
		{
		if(obj==null) continue;
		boolean f=false;
		for(final Object o: this.objects) if(o==obj){ f=true; break;}
		if(!f) this.objects.add(obj);
		}
	return this;
	}

List<String> parseOptions(final String argcargv[]) {
	final List<String> args = new ArrayList<>(Arrays.asList(argcargv));
	final Set<String> visited=new HashSet<>();
	scanAll(new ArgScanner() {
		@Override
		void process(Object o, Field f, Argument arg) {
			boolean multiple=false;
			boolean has_arg =false;
			
			List<String> search=new ArrayList<>();
			if(!(arg.opt()==null || arg.opt().isEmpty())) {
				search.add("-"+arg.opt());
				}
			if(!(arg.longopt()==null || arg.longopt().isEmpty())) {
				search.add("--"+arg.longopt());
				}
			for(String opt: search) {
				int optind = -1;
				for(int i=0;i< args.size();++i) {
					String x = args.get(i);
					if(!x.startsWith(opt)) continue;
					int colon = x.indexOf(":");
					String left=(colon==-1?opt:opt.substring(colon+1));
					String right=(colon==-1?opt:opt.substring(colon+1));
					optind=i;
					break;
				}
				if(optind==-1) continue;
				if(!multiple && visited.contains(opt)) {
					throw new RuntimeException("Option used twice");
					}
				if(has_arg) {
					Constructor<?> ctor = f.getType().getConstructor(String.class);
					Object v = ctor.newInstance(args.get(optind+1));
					f.set(o, v);
					}
				visited.addAll(search);
				}
			}
		});
	
	return args;
	}

private void scanAll(final ArgScanner scanner) {
	for(final Object o: this.objects) {
		scanObject(o,o.getClass(),scanner);
		}
	}
private void scanObject(final Object o,final Class<?> clazz,final ArgScanner scanner)
	{
	for(final Field field :clazz.getDeclaredFields()) {
		final Argument argument = field.getAnnotation(Argument.class);
		if(argument==null) continue;
		scanner.process(o, field, argument);
		}
	final Class<?> parent = clazz.getSuperclass();
	if(parent==null || parent == Object.class) return;
	scanObject(o,parent,scanner);
	}

private abstract class ArgScanner {
	abstract void process(Object o,Field f,final Argument argument);
	}

}