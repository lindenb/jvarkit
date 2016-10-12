package com.github.lindenb.jvarkit.util.command;

import java.io.File;
import java.io.PrintWriter;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.function.Predicate;


public class Arguments {
@Argument
private static boolean print_hidden_options=false;

private final List<ArgumentDef> argumentdefs = new ArrayList<>();
private Predicate<ArgumentDef> argumentFilter=null;
public Arguments() {
	register(Arguments.class);
	}

public ArgumentDef findArgumentDefByOpt(final String opt) {
	if(opt==null) return null;
	for(final ArgumentDef ad: this.argumentdefs) {
		if(ad.hasOpt() && opt.equals(ad.getOpt())) return ad;
	}
	return null;
}

public ArgumentDef findArgumentDefByLongOpt(final String longopt) {
	if(longopt==null) return null;
	for(final ArgumentDef ad: this.argumentdefs) {
		if(ad.hasLongOpt() && longopt.equals(ad.getLongOpt())) return ad;
	}
	return null;
}

public Arguments register(final Object obj) {
	if(obj!=null) {
		register(obj,obj.getClass());
		}
	return this;
	}

public Arguments register(final Class<?> clazz) {
	register(null,clazz);
	return this;
	}

public Arguments filter(Predicate<ArgumentDef> filter) {
	this.argumentFilter  = filter;
	int i=0;
	while(this.argumentFilter!=null && i< this.argumentdefs.size()) {
		if( !this.argumentFilter.test(this.argumentdefs.get(i))) {
			this.argumentdefs.remove(i);
			}
		else
			{
			++i;
			}
		}
	return this;
	}


private void register(final Object obj,final Class<?> clazz) {
	/** loop over the field of the class */
	for(final Field field: clazz.getDeclaredFields()) {
		//if no object, we're only looking a static fields
		if(obj==null && !java.lang.reflect.Modifier.isStatic(field.getModifiers())) continue;
		//if real object, we're not looking a static fields
		if(obj!=null && java.lang.reflect.Modifier.isStatic(field.getModifiers())) continue;
		//retrieve argument
		final Argument argument  = field.getAnnotation(Argument.class);
		//no argument
		if(argument==null) continue;
		//create new argument def
		final ArgumentDef def = new  ArgumentDef(obj, field,argument);
		//check we can modify it
		if(!def.isMultiple() && java.lang.reflect.Modifier.isFinal(field.getModifiers())) {
			throw new RuntimeException("not a collection and flagged final : "+def);
			}
		//discard if a argument filter is set
		if( argumentFilter!=null && !argumentFilter.test(def)) continue;
		//discard if already registerd
		if( argumentdefs.contains(def) ) {
			System.err.println("Argument already registered: "+def);
			continue;
			}
		//add to list of arguments
		this.argumentdefs.add(def);
		}
	//observe parent class
	final Class<?> parentClass = clazz.getSuperclass();
	if(parentClass==null || parentClass == Object.class) return;
	register(obj,parentClass);
	
	}

List<String> parseOptions(final String argcargv[]) {
	final List<String> args = new ArrayList<>(Arrays.asList(argcargv));
	for(final ArgumentDef def: this.argumentdefs) {
		
		}
	return args;
	}

public void printHelp(final PrintWriter w) {
w.print(this.argumentdefs);
w.flush();
}


private static class Test {
@Argument
boolean f1;
@Argument
Boolean f2;
@Argument
Collection<File> f3 = new ArrayList<>();
}

public static void main(String[] args)
	{
	args=new String[]{};
	Arguments arguments = new Arguments().register(new Test());
	arguments.printHelp(new PrintWriter(System.out));
	for(ArgumentDef tt:arguments.argumentdefs)
		System.err.println(tt);
	}
}