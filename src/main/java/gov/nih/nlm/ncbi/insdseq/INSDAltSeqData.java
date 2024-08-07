//
// This file was generated by the Eclipse Implementation of JAXB, v4.0.5 
// See https://eclipse-ee4j.github.io/jaxb-ri 
// Any modifications to this file will be lost upon recompilation of the source schema. 
//


package gov.nih.nlm.ncbi.insdseq;

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
    "insdAltSeqDataName",
    "insdAltSeqDataItems"
})
@XmlRootElement(name = "INSDAltSeqData")
public class INSDAltSeqData {

    @XmlElement(name = "INSDAltSeqData_name", required = true)
    protected String insdAltSeqDataName;
    @XmlElement(name = "INSDAltSeqData_items")
    protected INSDAltSeqDataItems insdAltSeqDataItems;

    /**
     * Gets the value of the insdAltSeqDataName property.
     * 
     * @return
     *     possible object is
     *     {@link String }
     *     
     */
    public String getINSDAltSeqDataName() {
        return insdAltSeqDataName;
    }

    /**
     * Sets the value of the insdAltSeqDataName property.
     * 
     * @param value
     *     allowed object is
     *     {@link String }
     *     
     */
    public void setINSDAltSeqDataName(String value) {
        this.insdAltSeqDataName = value;
    }

    /**
     * Gets the value of the insdAltSeqDataItems property.
     * 
     * @return
     *     possible object is
     *     {@link INSDAltSeqDataItems }
     *     
     */
    public INSDAltSeqDataItems getINSDAltSeqDataItems() {
        return insdAltSeqDataItems;
    }

    /**
     * Sets the value of the insdAltSeqDataItems property.
     * 
     * @param value
     *     allowed object is
     *     {@link INSDAltSeqDataItems }
     *     
     */
    public void setINSDAltSeqDataItems(INSDAltSeqDataItems value) {
        this.insdAltSeqDataItems = value;
    }

}
