//
// This file was generated by the Eclipse Implementation of JAXB, v4.0.5 
// See https://eclipse-ee4j.github.io/jaxb-ri 
// Any modifications to this file will be lost upon recompilation of the source schema. 
//


package gov.nih.nlm.ncbi.gb;

import jakarta.xml.bind.annotation.XmlAccessType;
import jakarta.xml.bind.annotation.XmlAccessorType;
import jakarta.xml.bind.annotation.XmlElement;
import jakarta.xml.bind.annotation.XmlRootElement;
import jakarta.xml.bind.annotation.XmlType;


/**
 * 
 */
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "", propOrder = {
    "gbInterval"
})
@XmlRootElement(name = "GBAltSeqItem_interval")
public class GBAltSeqItemInterval {

    @XmlElement(name = "GBInterval", required = true)
    protected GBInterval gbInterval;

    /**
     * Gets the value of the gbInterval property.
     * 
     * @return
     *     possible object is
     *     {@link GBInterval }
     *     
     */
    public GBInterval getGBInterval() {
        return gbInterval;
    }

    /**
     * Sets the value of the gbInterval property.
     * 
     * @param value
     *     allowed object is
     *     {@link GBInterval }
     *     
     */
    public void setGBInterval(GBInterval value) {
        this.gbInterval = value;
    }

}
