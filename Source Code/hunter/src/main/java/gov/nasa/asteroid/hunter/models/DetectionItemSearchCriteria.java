/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package gov.nasa.asteroid.hunter.models;

import java.util.Date;

import org.hibernate.validator.constraints.Range;

/**
 * Represents detection item search criteria.
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
public class DetectionItemSearchCriteria extends BaseSearchParameters {
    /**
     * Represents the date range start (time part will be ignored).
     *
     * Optional.
     */
    private Date dateStart;

    /**
     * Represents the date range end (time part will be ignored).
     *
     * Optional.
     */
    private Date dateEnd;

    /**
     * Represents the hour range start.
     *
     * Optional. [0, 23]
     */
    @Range(min = 0, max = 23)
    private Integer hourStart;

    /**
     * Represents the hour range end.
     *
     * Optional. [0, 23]
     */
    @Range(min = 0, max = 23)
    private Integer hourEnd;

    /**
     * Represents the minute range start.
     *
     * Optional. [0, 59]
     */
    @Range(min = 0, max = 59)
    private Integer minuteStart;

    /**
     * Represents the minute range end.
     *
     * Optional. [0, 59]
     */
    @Range(min = 0, max = 59)
    private Integer minuteEnd;

    /**
     * Indicates whether the detection item is submitted.
     *
     * Optional.
     */
    private Boolean submitted;

    /**
     * <p>
     * The default constructor for class <code>DetectionItemSearchCriteria</code>.
     * </p>
     */
    public DetectionItemSearchCriteria() {
        // does nothing
    }

    /**
     * <p>
     * Gets the date range start (time part will be ignored).
     * </p>
     * @return the date range start (time part will be ignored).
     */
    public Date getDateStart() {
        return dateStart;
    }

    /**
     * <p>
     * Sets the date range start (time part will be ignored).
     * </p>
     * @param dateStart the date range start (time part will be ignored).
     */
    public void setDateStart(Date dateStart) {
        this.dateStart = dateStart;
    }

    /**
     * <p>
     * Gets the date range end (time part will be ignored).
     * </p>
     * @return  the date range end (time part will be ignored).
     */
    public Date getDateEnd() {
        return dateEnd;
    }

    /**
     * <p>
     * Sets the date range end (time part will be ignored).
     * </p>
     * @param dateEnd the date range end (time part will be ignored).
     */
    public void setDateEnd(Date dateEnd) {
        this.dateEnd = dateEnd;
    }

    /**
     * <p>
     * Gets the hour range start.
     * </p>
     * @return the hour range start.
     */
    public Integer getHourStart() {
        return hourStart;
    }

    /**
     * <p>
     * Sets the hour range start.
     * </p>
     * @param hourStart the hour range start.
     */
    public void setHourStart(Integer hourStart) {
        this.hourStart = hourStart;
    }

    /**
     * <p>
     * Gets the hour range end.
     * </p>
     * @return the hour range end.
     */
    public Integer getHourEnd() {
        return hourEnd;
    }

    /**
     * <p>
     * Sets the hour range end.
     * </p>
     * @param hourEnd the hour range end.
     */
    public void setHourEnd(Integer hourEnd) {
        this.hourEnd = hourEnd;
    }

    /**
     * <p>
     * Gets the minute range start.
     * </p>
     * @return the minute range start.
     */
    public Integer getMinuteStart() {
        return minuteStart;
    }

    /**
     * <p>
     * Sets the minute range start.
     * </p>
     * @param minuteStart the minute range start.
     */
    public void setMinuteStart(Integer minuteStart) {
        this.minuteStart = minuteStart;
    }

    /**
     * <p>
     * Gets the minute range end.
     * </p>
     * @return the minute range end.
     */
    public Integer getMinuteEnd() {
        return minuteEnd;
    }

    /**
     * <p>
     * Sets the minute range end.
     * </p>
     * @param minuteEnd the minute range end.
     */
    public void setMinuteEnd(Integer minuteEnd) {
        this.minuteEnd = minuteEnd;
    }

    /**
     * <p>
     * Gets whether the detection item is submitted.
     * </p>
     * @return whether the detection item is submitted.
     */
    public Boolean isSubmitted() {
        return submitted;
    }

    /**
     * <p>
     * Sets whether the detection item is submitted.
     * </p>
     * @param submitted whether the detection item is submitted.
     */
    public void setSubmitted(Boolean submitted) {
        this.submitted = submitted;
    }
}

