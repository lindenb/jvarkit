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
    "gbFeatureKey",
    "gbFeatureLocation",
    "gbFeatureIntervals",
    "gbFeatureOperator",
    "gbFeaturePartial5",
    "gbFeaturePartial3",
    "gbFeatureQuals",
    "gbFeatureXrefs"
})
@XmlRootElement(name = "GBFeature")
public class GBFeature {

    @XmlElement(name = "GBFeature_key", required = true)
    protected String gbFeatureKey;
    @XmlElement(name = "GBFeature_location", required = true)
    protected String gbFeatureLocation;
    @XmlElement(name = "GBFeature_intervals")
    protected GBFeatureIntervals gbFeatureIntervals;
    @XmlElement(name = "GBFeature_operator")
    protected String gbFeatureOperator;
    @XmlElement(name = "GBFeature_partial5")
    protected GBFeaturePartial5 gbFeaturePartial5;
    @XmlElement(name = "GBFeature_partial3")
    protected GBFeaturePartial3 gbFeaturePartial3;
    @XmlElement(name = "GBFeature_quals")
    protected GBFeatureQuals gbFeatureQuals;
    @XmlElement(name = "GBFeature_xrefs")
    protected GBFeatureXrefs gbFeatureXrefs;

    /**
     * Gets the value of the gbFeatureKey property.
     * 
     * @return
     *     possible object is
     *     {@link String }
     *     
     */
    public String getGBFeatureKey() {
        return gbFeatureKey;
    }

    /**
     * Sets the value of the gbFeatureKey property.
     * 
     * @param value
     *     allowed object is
     *     {@link String }
     *     
     */
    public void setGBFeatureKey(String value) {
        this.gbFeatureKey = value;
    }

    /**
     * Gets the value of the gbFeatureLocation property.
     * 
     * @return
     *     possible object is
     *     {@link String }
     *     
     */
    public String getGBFeatureLocation() {
        return gbFeatureLocation;
    }

    /**
     * Sets the value of the gbFeatureLocation property.
     * 
     * @param value
     *     allowed object is
     *     {@link String }
     *     
     */
    public void setGBFeatureLocation(String value) {
        this.gbFeatureLocation = value;
    }

    /**
     * Gets the value of the gbFeatureIntervals property.
     * 
     * @return
     *     possible object is
     *     {@link GBFeatureIntervals }
     *     
     */
    public GBFeatureIntervals getGBFeatureIntervals() {
        return gbFeatureIntervals;
    }

    /**
     * Sets the value of the gbFeatureIntervals property.
     * 
     * @param value
     *     allowed object is
     *     {@link GBFeatureIntervals }
     *     
     */
    public void setGBFeatureIntervals(GBFeatureIntervals value) {
        this.gbFeatureIntervals = value;
    }

    /**
     * Gets the value of the gbFeatureOperator property.
     * 
     * @return
     *     possible object is
     *     {@link String }
     *     
     */
    public String getGBFeatureOperator() {
        return gbFeatureOperator;
    }

    /**
     * Sets the value of the gbFeatureOperator property.
     * 
     * @param value
     *     allowed object is
     *     {@link String }
     *     
     */
    public void setGBFeatureOperator(String value) {
        this.gbFeatureOperator = value;
    }

    /**
     * Gets the value of the gbFeaturePartial5 property.
     * 
     * @return
     *     possible object is
     *     {@link GBFeaturePartial5 }
     *     
     */
    public GBFeaturePartial5 getGBFeaturePartial5() {
        return gbFeaturePartial5;
    }

    /**
     * Sets the value of the gbFeaturePartial5 property.
     * 
     * @param value
     *     allowed object is
     *     {@link GBFeaturePartial5 }
     *     
     */
    public void setGBFeaturePartial5(GBFeaturePartial5 value) {
        this.gbFeaturePartial5 = value;
    }

    /**
     * Gets the value of the gbFeaturePartial3 property.
     * 
     * @return
     *     possible object is
     *     {@link GBFeaturePartial3 }
     *     
     */
    public GBFeaturePartial3 getGBFeaturePartial3() {
        return gbFeaturePartial3;
    }

    /**
     * Sets the value of the gbFeaturePartial3 property.
     * 
     * @param value
     *     allowed object is
     *     {@link GBFeaturePartial3 }
     *     
     */
    public void setGBFeaturePartial3(GBFeaturePartial3 value) {
        this.gbFeaturePartial3 = value;
    }

    /**
     * Gets the value of the gbFeatureQuals property.
     * 
     * @return
     *     possible object is
     *     {@link GBFeatureQuals }
     *     
     */
    public GBFeatureQuals getGBFeatureQuals() {
        return gbFeatureQuals;
    }

    /**
     * Sets the value of the gbFeatureQuals property.
     * 
     * @param value
     *     allowed object is
     *     {@link GBFeatureQuals }
     *     
     */
    public void setGBFeatureQuals(GBFeatureQuals value) {
        this.gbFeatureQuals = value;
    }

    /**
     * Gets the value of the gbFeatureXrefs property.
     * 
     * @return
     *     possible object is
     *     {@link GBFeatureXrefs }
     *     
     */
    public GBFeatureXrefs getGBFeatureXrefs() {
        return gbFeatureXrefs;
    }

    /**
     * Sets the value of the gbFeatureXrefs property.
     * 
     * @param value
     *     allowed object is
     *     {@link GBFeatureXrefs }
     *     
     */
    public void setGBFeatureXrefs(GBFeatureXrefs value) {
        this.gbFeatureXrefs = value;
    }

}