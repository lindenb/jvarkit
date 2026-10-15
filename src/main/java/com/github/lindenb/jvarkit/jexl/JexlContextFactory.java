package com.github.lindenb.jvarkit.jexl;

import java.util.HashMap;
import java.util.Map;
import java.util.Objects;
import java.util.function.Function;

import org.apache.commons.jexl2.JexlContext;

/**
 * Factory for JEXL contexts
 */
public class JexlContextFactory<T>  {
	private final Map<String, Function<T,Object>> extractors;
	public JexlContextFactory() {
		this.extractors = new HashMap<>();
		}
	/** copy constructor */
	public JexlContextFactory(final JexlContextFactory<T> cp) {
		this.extractors = new HashMap<>(cp.extractors);
		}

	/** invoke copy constructor */
	@Override
	public JexlContextFactory<T>  clone() {
		return new JexlContextFactory<T>(this);
		}
	
	/** 
	 * add extractor to extract property 'key' from oject T
	 * @return this
	 *
	 */
	public JexlContextFactory<T> add(final String key,Function<T,Object> extractor ) {
		this.extractors.put(key, extractor);
		return this;
		}
	
	/**  implementation of a JEXL context */
	private static class MyContext<T> implements JexlContext {
		private final T object;
		private final Map<String, Function<T,Object>> extractors;
		MyContext(final T object,Map<String, Function<T,Object>> extractors) {
			this.object = object;
			this.extractors = new HashMap<>(extractors);
			}
		@Override
		public boolean has(final String key) {
			return extractors.containsKey(key);
			}
		@Override
		public Object get(final String key) {
			if(this.extractors.containsKey(key)) {
				return this.extractors.get(key).apply(this.object);
				}
			throw new IllegalArgumentException("cannot invoke jexl context for \""+key+"\". Available keys are :"+String.join(", ", this.extractors.keySet()));
			}
		@Override
		public void set(final String key, final Object arg1) {
			throw new IllegalArgumentException("cannot set jexl context for \""+key+"\". This is a read-only context");
			}
		@Override
		public int hashCode() {
			return Objects.hash(this.object,this.extractors);
			}
		@Override
		public String toString() {
			return "JexlContextFactory:MyContext("+String.valueOf(this.object)+")";
			}
		}

	/** create a new JexlContext from the object */
	public JexlContext make(final T object) {
		return new  MyContext<>(object,this.extractors);
		}
	
	@Override
	public String toString() {
		return "JexlContextFactory["+String.join("| ", this.extractors.keySet())+"]";
		}
	}
