package com.github.lindenb.jvarkit.tools.jmx;


import java.text.DecimalFormat;
import java.text.NumberFormat;

import htsjdk.samtools.util.Locatable;

import javax.management.Attribute;
import javax.management.AttributeChangeNotification;
import javax.management.AttributeList;
import javax.management.AttributeNotFoundException;
import javax.management.DynamicMBean;
import javax.management.InvalidAttributeValueException;
import javax.management.MBeanAttributeInfo;
import javax.management.MBeanConstructorInfo;
import javax.management.MBeanException;
import javax.management.MBeanInfo;
import javax.management.MBeanOperationInfo;
import javax.management.MBeanParameterInfo;
import javax.management.NotificationBroadcasterSupport;
import javax.management.ReflectionException;
import javax.management.RuntimeOperationsException;

/* http://docs.oracle.com/cd/E19698-01/816-7609/6mdjrf83d/index.html */
public class LocatableStreamInfo
	extends NotificationBroadcasterSupport
	implements DynamicMBean
	{
	public enum Status{IDLE,RUNNING,BREAK,ABORT};
	private String projectName="";
	private String lastPosition="";
	private long startMillisec =-1L;
	private long countRecords=0L;
	private Status control=Status.IDLE;
	
	public LocatableStreamInfo(String projectName)
		{
		this.projectName=projectName;
		}
	
	public LocatableStreamInfo()
		{
		this("undefined");
		}
	
	public String getProjectName() {
		return projectName;
		}
	
	public void setProjectName(String projectName) {
		this.projectName = projectName;
		}
	
	public void setLastPosition(String lastPosition) {
		this.lastPosition = lastPosition;
		}
	public String getLastPosition() {
		return lastPosition;
		}
    
    
	public String getDuration()
		{
		if(this.startMillisec<0L) return "";
		NumberFormat timeFmt = new DecimalFormat("00");
		long seconds =  (System.currentTimeMillis() - this.startMillisec) / 1000;
        final long s = seconds % 60;
        final long allMinutes = seconds / 60;
        final long m = allMinutes % 60;
        final long h = allMinutes / 60;

        return timeFmt.format(h) + ":" + timeFmt.format(m) + ":" + timeFmt.format(s);
		}
	
	public long getCount() {
		return countRecords;
	}
	
	public String getState() {
		return this.control.name();
	}
	
	public Status watch(Locatable locatable)
		{
		if(this.startMillisec<0L)
			{
			synchronized (LocatableStreamInfo.class) 
				{
				this.startMillisec = System.currentTimeMillis();
				if(this.control.equals(Status.IDLE)) this.control= Status.RUNNING;		
				}
			
			}
		if(locatable.getContig()!=null)
			{
			this.lastPosition = locatable.getContig()+":"+locatable.getStart();
			}
		else
			{
			this.lastPosition =  "unmapped";
			}
		this.countRecords++;
		/** no need to send messages 
		
			sendNotification(new AttributeChangeNotification(this,
	                0,
	                System.currentTimeMillis(),
	                "Chromosome",
	                "Chromosome",
	                "java.lang.String",
	                locatable.getContig(),
	                prevContig==null?".":prevContig
		            ));
			prevContig=locatable.getContig();
		*/
		
		return this.control;
		}
	
	
	
	/*  metadata information about attributes and operations. */
	private MBeanInfo dMBeanInfo = buildDynamicMBeanInfo();	
	
	private MBeanInfo buildDynamicMBeanInfo()
	 	{		
		return new MBeanInfo(
				LocatableStreamInfo.class.getName(),
				"Mon premier MBean Dynamic",
				new MBeanAttributeInfo[]{
					new MBeanAttributeInfo("ProjectName","java.lang.String","Project Name", true,false,false),
					new MBeanAttributeInfo("LastPosition","java.lang.String","Last Genomic Position", true,false,false),
					new MBeanAttributeInfo("Duration","java.lang.String","Duration", true,false,false),
					new MBeanAttributeInfo("Count","java.lang.Long","Number of records seen so far", true,false,false),
					new MBeanAttributeInfo("State","java.lang.String","Streaming State (IDLE/RUNNING/BREAK/ABORT)", true,false,false),
				},
				new MBeanConstructorInfo[]{
					new MBeanConstructorInfo(
							"SimpleDynamic",
							"SimpleDynamic",
							new MBeanParameterInfo[0]
							)
					
				},
				new MBeanOperationInfo[]{
					new MBeanOperationInfo("doAbort", "Abort the loop with error", new MBeanParameterInfo[0],Void.TYPE.getName(),MBeanOperationInfo.ACTION),
					new MBeanOperationInfo("doBreak", "Break the loop and continue", new MBeanParameterInfo[0],Void.TYPE.getName(),MBeanOperationInfo.ACTION),
					
				},
				null);
	 	}
	
	
	
	@Override
	public synchronized Object getAttribute(String attribute_name)
			throws AttributeNotFoundException, MBeanException,
			ReflectionException
		{
	    // Check attribute_name to avoid NullPointerException later on
	    if (attribute_name == null) 
	        throw new RuntimeOperationsException(new IllegalArgumentException("attName cannot be null"));
	    // Call the corresponding getter for a recognized attribute_name
	    if (attribute_name.equals("ProjectName")) return getProjectName();
	    if (attribute_name.equals("LastPosition")) return getLastPosition();
	    if (attribute_name.equals("Duration")) return getDuration();
	    if (attribute_name.equals("Count")) return getCount();
	    if (attribute_name.equals("State")) return getState();
	    // If attribute_name has not been recognized
	    throw(new AttributeNotFoundException(
	        "Cannot find " + attribute_name + " attribute in " + LocatableStreamInfo.class.getName()));
		}

	@Override
	public synchronized void setAttribute(Attribute attribute)
			throws AttributeNotFoundException, InvalidAttributeValueException,
			MBeanException, ReflectionException {
	    // Check attribute to avoid NullPointerException later on
	    if (attribute == null) {
	        throw new RuntimeOperationsException(
	            new IllegalArgumentException("Attribute cannot be null"),
	            "Cannot invoke a setter of " +  LocatableStreamInfo.class.getName() +
	                " with null attribute");
	    	}
	    
        throw new AttributeNotFoundException(
            "Attribute " + attribute.getName() + " read only or not found in " +
                this.getClass().getName());
   	
		}

	
	@Override
	public synchronized AttributeList getAttributes(String[] attributeNames) {
	    if (attributeNames == null) {
	        throw new RuntimeOperationsException(
	            new IllegalArgumentException(
	                "attributeNames[] cannot be null"),
	            "Cannot invoke a getter of " + LocatableStreamInfo.class.getName());
	    }
	    AttributeList resultList = new AttributeList();
	    // build the result attribute list
	    for (String attName:attributeNames)
	    	{
	        try {        
	            Object value = getAttribute(attName);     
	            resultList.add(new Attribute(attName,value));
	        	}
	        catch (Exception e)
	        	{
	            e.printStackTrace();
	        	}
	    	}
	    return(resultList);
	    }

	@Override
	public synchronized AttributeList setAttributes(AttributeList attributes) {
	    // Check attributes to avoid NullPointerException later on
	    if (attributes == null) {
	        throw new RuntimeOperationsException(
	            new IllegalArgumentException(
	                "AttributeList attributes cannot be null"),
	            "Cannot invoke a setter of " + LocatableStreamInfo.class.getName());
	    }
	    AttributeList resultList = new AttributeList();
	    // try to set each attribute and add to result list if successful
	    for (Object i: attributes) {
	        Attribute attr = (Attribute) i;
	        try {
	            setAttribute(attr);
	            String name = attr.getName();
	            Object value = getAttribute(name); 
	            resultList.add(new Attribute(name,value));
	        } catch(Exception e) {
	            e.printStackTrace();
	        }
	    }
	    return(resultList);	}

	@Override
	public synchronized Object invoke(String operationName, Object[] params, String[] signature)
			throws MBeanException, ReflectionException {
	    // Check operationName to avoid NullPointerException later on
	    if (operationName == null) {
	        throw new RuntimeOperationsException(
	            new IllegalArgumentException(
	                "Operation name cannot be null"),
	            "Cannot invoke a null operation in " + LocatableStreamInfo.class.getName());
	    }

	    // Call the corresponding operation for a recognized name
	    if (operationName.equals("doBreak")){
	        // this code is specific to the internal "reset" method:
	        doBreak();     // no parameters to check
	        return null; // and no return value
	    }else if (operationName.equals("doAbort")){
	        // this code is specific to the internal "reset" method:
	        doAbort();     // no parameters to check
	        return null; // and no return value
	    } else { 
	        // unrecognized operation name:
	        throw new ReflectionException(
	            new NoSuchMethodException(operationName), 
	            "Cannot find the operation " + operationName +
	                " in " + LocatableStreamInfo.class.getName());
	    }
	    
	}

	@Override
	public MBeanInfo getMBeanInfo() {
		return dMBeanInfo;
	}

	private synchronized void setControl(Status ctrl)
		{
		if(this.control.equals(Status.ABORT)) return;
		if(this.control.equals(Status.BREAK) && !ctrl.equals(Status.ABORT)) return;
		
		AttributeChangeNotification acn =
	            new AttributeChangeNotification(this,
                0,
                System.currentTimeMillis(),
                "State changed",
                "State",
                "java.lang.String",
                this.control.name(),
                ctrl.name()
	            );
		this.control=ctrl;
		sendNotification(acn); 
		}
	
	public synchronized void doBreak()
		{
		setControl(Status.BREAK);
		}
	
	public synchronized void doAbort()
		{
		setControl(Status.ABORT);
		}
	}
