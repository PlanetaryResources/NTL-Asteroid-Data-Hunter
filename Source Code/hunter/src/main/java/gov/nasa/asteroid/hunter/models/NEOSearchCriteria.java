/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package gov.nasa.asteroid.hunter.models;

import java.util.Date;

import javax.validation.constraints.NotNull;
import javax.validation.constraints.Pattern;

import org.hibernate.validator.constraints.NotBlank;
import org.hibernate.validator.constraints.Range;

/**
 * Represents NEO search criteria.
 *
 * JSR303 validation is needed for this class.
 *
 * <p>
 * <strong>Thread-Safety:</strong> This class is not thread safe since it is
 * mutable.
 * </p>
 *
 * @author albertwang, TCSASSEMBLER
 * @version 1.0
 */
public class NEOSearchCriteria extends BaseSearchParameters {
    /**
     * Represents the date time.
     *
     * Required.
     */
    @NotNull
    private Date date;

    /**
     * Represents the right ascension.
     *
     * Required. In format of "01 36 41.7", "01" here is hour, 36 here is
     * minute, 41.7 here is second.
     */
    @Pattern(regexp = "^[0-5][0-9] [0-5][0-9] [0-5][0-9](.\\d+)?$")
    private String rightAscension;

    /**
     * Represents the declination.
     *
     * Required. In format of "+01 48.3".
     * 1 here is second.
     */
    @Pattern(regexp = "^.* [0-9]*(.\\d+)?$")
    private String declination;

    /**
     * Represents the observations.
     *
     * Required. Observation 80-column record format
     * (http://www.minorplanetcenter.net/iau/info/ObsFormat.html)
     */
    @Pattern(regexp = "^.* .*$")
    private String observations;

    /**
     * Represents the radius.
     *
     * Required. 5 ~ 300.
     */
    @Range(min = 5, max = 300)
    private double radius;

    /**
     * Represents the V value.
     *
     * Required. >= 1.0.
     */
    @Range(min = 1)
    private double v;

    /**
     * Represents the observatory code.
     *
     * Required. A valid observatory code
     * (http://en.wikipedia.org/wiki/List_of_observatory_codes)
     */
    @NotBlank
    private String observatoryCode;

    /**
     * <p>
     * The default constructor for class <code>NEOSearchCriteria</code>.
     * </p>
     */
    public NEOSearchCriteria() {
        // does nothing
    }

    /**
     * <p>
     * Gets the date time.
     * </p>
     * @return the date time.
     */
    public Date getDate() {
        return date;
    }

    /**
     * <p>
     * Sets the date time.
     * </p>
     * @param date the date time.
     */
    public void setDate(Date date) {
        this.date = date;
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
     * Gets the declination.
     * </p>
     * @return the declination.
     */
    public String getDeclination() {
        return declination;
    }

    /**
     * <p>
     * Sets the declination.
     * </p>
     * @param declination the declination.
     */
    public void setDeclination(String declination) {
        this.declination = declination;
    }

    /**
     * <p>
     * Gets the observations.
     * </p>
     * @return the observations.
     */
    public String getObservations() {
        return observations;
    }

    /**
     * <p>
     * Sets the observations.
     * </p>
     * @param observations the observations.
     */
    public void setObservations(String observations) {
        this.observations = observations;
    }

    /**
     * <p>
     * Gets the radius.
     * </p>
     * @return the radius.
     */
    public double getRadius() {
        return radius;
    }

    /**
     * <p>
     * Sets the radius.
     * </p>
     * @param radius the radius.
     */
    public void setRadius(double radius) {
        this.radius = radius;
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
     * Gets the observatory code.
     * </p>
     * @return the observatory code.
     */
    public String getObservatoryCode() {
        return observatoryCode;
    }

    /**
     * <p>
     * Sets the observatory code.
     * </p>
     * @param observatoryCode the observatory code.
     */
    public void setObservatoryCode(String observatoryCode) {
        this.observatoryCode = observatoryCode;
    }
}
