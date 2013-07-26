package com.github.lindenb.jvarkit.tools.ws.rmi;

import java.rmi.Remote;
import java.rmi.RemoteException;
import java.util.List;

import com.github.lindenb.jvarkit.tools.ws.WSProject;

public interface NGSService extends Remote
	{
	public static final String SERVICE_NAME="ngsservice";
	public  List<? extends WSProject> getProjects() throws RemoteException;
	public  WSProject getProjectById(String id) throws RemoteException;
	}
