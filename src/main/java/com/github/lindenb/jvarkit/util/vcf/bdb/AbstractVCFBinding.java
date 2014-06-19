package com.github.lindenb.jvarkit.util.vcf.bdb;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


import com.sleepycat.bind.tuple.BigDecimalBinding;
import com.sleepycat.bind.tuple.BigIntegerBinding;
import com.sleepycat.bind.tuple.BooleanBinding;
import com.sleepycat.bind.tuple.ByteBinding;
import com.sleepycat.bind.tuple.CharacterBinding;
import com.sleepycat.bind.tuple.DoubleBinding;
import com.sleepycat.bind.tuple.FloatBinding;
import com.sleepycat.bind.tuple.IntegerBinding;
import com.sleepycat.bind.tuple.LongBinding;
import com.sleepycat.bind.tuple.ShortBinding;
import com.sleepycat.bind.tuple.StringBinding;
import com.sleepycat.bind.tuple.TupleBinding;
import com.sleepycat.bind.tuple.TupleInput;
import com.sleepycat.bind.tuple.TupleOutput;

/**
 * Base class for serializing some java.lang.Object in a VCF (in INFO and Genotype)
 * @author lindenb
 *
 * @param <T>
 */
@SuppressWarnings({"rawtypes","unchecked"})
public abstract class AbstractVCFBinding<T> extends TupleBinding<T>
	{
	/** array of tuple com.sleepycat.bind.tuple.*Binding to serialize some simple objects */
	private static final List<TupleBinding> INDEX_TO_BINDING=new ArrayList<TupleBinding>();

	/** map java.langClass to an index in INDEX_TO_BINDING */
	private static final Map<Class,Integer> CLASS_TO_BINDING_INDEX=new HashMap<Class, Integer>();
	static 
		{
		addClassBinding(Byte.class,Byte.TYPE,new ByteBinding());
		addClassBinding(Short.class,Short.TYPE,new ShortBinding());
		addClassBinding(Integer.class,Integer.TYPE,new IntegerBinding());
		addClassBinding(Long.class,Long.TYPE,new LongBinding());
		addClassBinding(Float.class,Float.TYPE,new FloatBinding());
		addClassBinding(Double.class,Double.TYPE,new DoubleBinding());
		addClassBinding(Character.class,Character.TYPE,new CharacterBinding());
		addClassBinding(Boolean.class,Boolean.TYPE,new BooleanBinding());
		addClassBinding(String.class,null,new StringBinding());
		addClassBinding(BigDecimal.class,null,new BigDecimalBinding());
		addClassBinding(BigInteger.class,null,new BigIntegerBinding());
		}
	/** add this opcode to tell that and object is an 'array-of' */
	private static final byte OP_ARRAY=(byte)100;
	
	/**put binding for this primitive */
	private static <T> void  addClassBinding(Class<T> c1,Class<T> c2,TupleBinding<T> bind)
		{
		CLASS_TO_BINDING_INDEX.put(c1, INDEX_TO_BINDING.size());
		INDEX_TO_BINDING.add(bind);
		if(c2!=null) addClassBinding(c2,null,bind);
		}
	
	protected AbstractVCFBinding()
		{
		}
	
	
	/** unserialize any object using opcodes */
	protected Object readAttribute(TupleInput in)
		{
		byte opcode=in.readByte();
		/* this is an array */
		if(opcode>=OP_ARRAY)
			{
			int n_items=in.readInt();
			TupleBinding bind=INDEX_TO_BINDING.get(opcode-OP_ARRAY);
			List array=new ArrayList(n_items);
			for(int i=0;i< n_items;++i)
				{
				array.add(bind.entryToObject(in));
				}
			return array;
			}
		/* this is a simple type */
		TupleBinding bind=INDEX_TO_BINDING.get(opcode);
		return bind.entryToObject(in);
		}
	
	/** serialize any object using opcodes */
	protected void writeAttribute(TupleOutput out,Object o)
		{
		/* this is an array */
		if(o.getClass().isArray())
			{
			Class clazz=o.getClass().getComponentType();
			int idx=CLASS_TO_BINDING_INDEX.get(clazz);
			TupleBinding bind=INDEX_TO_BINDING.get(idx);
			out.writeByte((byte)(idx+OP_ARRAY));/* add opcode for array */
			Object array[]=(Object[])o;
			out.writeInt(array.length);
			for(int i=0;i< array.length;++i)
				{
				bind.objectToEntry(array[i], out);
				}
			}
		/* this is a list */
		else if(o instanceof List)
			{
			List array=(List)o;
			Class clazz=array.isEmpty()?Byte.class:array.get(0).getClass();
			int idx=CLASS_TO_BINDING_INDEX.get(clazz);
			TupleBinding bind=INDEX_TO_BINDING.get(idx);
			out.writeByte((byte)(idx+OP_ARRAY));/* add opcode for array */
			out.writeInt(array.size());
			for(int i=0;i< array.size();++i)
				{
				bind.objectToEntry(array.get(i), out);
				}
			}
		/* this is a simple object */
		else
			{
			int idx= CLASS_TO_BINDING_INDEX.get(o.getClass());
			TupleBinding bind=INDEX_TO_BINDING.get(idx);
			out.writeByte((byte)idx);
			bind.objectToEntry(o, out);
			}
		}
	}
