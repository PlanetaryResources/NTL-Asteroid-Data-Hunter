/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package gov.nasa.asteroid.hunter.models;


/**
 * Represents a NEO.
 *
 * <p>
 * <strong>Thread-Safety:</strong> This class is not thread safe since it is
 * mutable.
 * </p>
 *
 * @author albertwang, TCSASSEMBLER
 * @version 1.0
 */
public class NEO {
    /**
     * Represents object designation.
     *
     * Required.
     */
    private String objectDesignation;

    /**
     * Represents the right ascension.
     */
    private String rightAscension;
    
    private double rightAscensionValue;

    /**
     * Represents the declination.
     */
    private String declination;
    
    private double declinationValue;

    /**
     * Represents the V value.
     */
    private double v;

    /**
     * Represents the offset right ascension.
     */
    private String offsetRightAscension;

    /**
     * Represents the offset declination.
     */
    private String offsetDeclination;

    /**
     * Represents the motion right ascension.
     */
    private String motionRightAscension;

    /**
     * Represents the motion declination.
     */
    private String motionDeclination;

    /**
     * Represents the orbit.
     */
    private String orbit;

    /**
     * <p>
     * The default constructor for class <code>NEO</code>.
     * </p>
     */
    public NEO() {
        // does nothing
    }

    /**
     * <p>
     * Gets object designation.
     * </p>
     * @return object designation.
     */
    public String getObjectDesignation() {
        return objectDesignation;
    }

    /**
     * <p>
     * Sets object designation.
     * </p>
     * @param objectDesignation object designation.
     */
    public void setObjectDesignation(String objectDesignation) {
        this.objectDesignation = objectDesignation;
    }

    /**
     * <p>
     * Gets the right ascension.
     * </p>
     * @return the right ascension.
     */
    public String getRightAscension() {
        return rightAscension;
    }

    /**
     * <p>
     * Sets the right ascension.
     * </p>
     * @param rightAscension the right ascension.
     */
    public void setRightAscension(String rightAscension) {
        this.rightAscension = rightAscension;
    }

    /**
     * <p>
     * Gets the right ascension.
     * </p>
     * @return the right ascension.
     */
    public String getDeclination() {
        return declination;
    }

    /**
     * <p>
     * Sets the right ascension.
     * </p>
     * @param declination the right ascension.
     */
    public void setDeclination(String declination) {
        this.declination = declination;
    }

    /**
     * <p>
     * Gets the V value.
     * </p>
     * @return the V value.
     */
    public double getV() {
        return v;
    }

    /**
     * <p>
     * Sets the V value.
     * </p>
     * @param v the V value.
     */
    public void setV(double v) {
        this.v = v;
    }

    /**
     * <p>
     * Gets the offset right ascension.
     * </p>
     * @return the offset right ascension.
     */
    public String getOffsetRightAscension() {
        return offsetRightAscension;
    }

    /**
     * <p>
     * Sets the offset right ascension.
     * </p>
     * @param offsetRightAscension the offset right ascension.
     */
    public void setOffsetRightAscension(String offsetRightAscension) {
        this.offsetRightAscension = offsetRightAscension;
    }

    /**
     * <p>
     * Gets the offset declination.
     * </p>
     * @return the offset declination.
     */
    public String getOffsetDeclination() {
        return offsetDeclination;
    }

    /**
     * <p>
     * Sets the offset declination.
     * </p>
     * @param offsetDeclination the offset declination.
     */
    public void setOffsetDeclination(String offsetDeclination) {
        this.offsetDeclination = offsetDeclination;
    }

    /**
     * <p>
     * Gets the motion right ascension.
     * </p>
     * @return the motion right ascension.
     */
    public String getMotionRightAscension() {
        return motionRightAscension;
    }

    /**
     * <p>
     * Sets the motion right ascension.
     * </p>
     * @param motionRightAscension the motion right ascension.
     */
    public void setMotionRightAscension(String motionRightAscension) {
        this.motionRightAscension = motionRightAscension;
    }

    /**
     * <p>
     * Gets the motion declination.
     * </p>
     * @return the motion declination.
     */
    public String getMotionDeclination() {
        return motionDeclination;
    }

    /**
     * <p>
     * Sets the motion declination.
     * </p>
     * @param motionDeclination the motion declination.
     */
    public void setMotionDeclination(String motionDeclination) {
        this.motionDeclination = motionDeclination;
    }

    /**
     * <p>
     * Gets the orbit.
     * </p>
     * @return the orbit.
     */
    public String getOrbit() {
        return orbit;
    }

    /**
     * <p>
     * Sets the orbit.
     * </p>
     * @param orbit the orbit.
     */
    public void setOrbit(String orbit) {
        this.orbit = orbit;
    }

    public double getRightAscensionValue() {
        return rightAscensionValue;
    }

    public void setRightAscensionValue(double rightAscensionValue) {
        this.rightAscensionValue = rightAscensionValue;
    }

    public double getDeclinationValue() {
        return declinationValue;
    }

    public void setDeclinationValue(double declinationValue) {
        this.declinationValue = declinationValue;
    }
}

