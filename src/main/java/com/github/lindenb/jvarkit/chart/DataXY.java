package com.github.lindenb.jvarkit.chart;

import java.awt.geom.Point2D;
import java.io.IOException;

import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.google.gson.stream.JsonWriter;

public class DataXY extends AbstractDataY {
	double x;
	public DataXY(final double x, final double y) {
		super(y);
		this.x = x;
		}
	public DataXY(Point2D pt) {
		this(pt.getX(),pt.getY());
		}
	
	
	public double getX() {
		return x;
		}
	
	void saveMultiQC(final JsonWriter w) throws IOException {
		w.beginArray();
		w.value(getX());
		w.value(getY());
		w.endArray();
		}
	
	void saveXml(XMLStreamWriter w) throws IOException,XMLStreamException {
		w.writeEmptyElement("point");
		w.writeAttribute("x", String.valueOf(getX()));
		w.writeAttribute("y", String.valueOf(getY()));
		}
	}