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
	for(final Field field: clazz.getDeclaredFields()) {
		if(obj==null && !java.lang.reflect.Modifier.isStatic(field.getModifiers())) continue;
		if(obj!=null && java.lang.reflect.Modifier.isStatic(field.getModifiers())) continue;
		if(java.lang.reflect.Modifier.isFinal(field.getModifiers())) continue;
		final Argument argument  = field.getAnnotation(Argument.class);
		if(argument==null) continue;
		final ArgumentDef def = new  ArgumentDef(obj, field,argument);
		if( argumentFilter!=null && !argumentFilter.test(def)) continue;
		if( argumentdefs.contains(def) ) {
			System.err.println("Argument already registered: "+def);
			continue;
			}
		this.argumentdefs.add(def);
		}
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