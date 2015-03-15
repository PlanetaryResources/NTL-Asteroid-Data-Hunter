/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package gov.nasa.asteroid.hunter.models;

import javax.persistence.Column;
import javax.persistence.EmbeddedId;
import javax.persistence.Entity;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.persistence.MapsId;
import javax.persistence.Table;
import javax.validation.constraints.Max;
import javax.validation.constraints.Min;

import com.fasterxml.jackson.annotation.JsonIgnore;


/**
 * Represents a detection item frame.
 *
 * <p>
 * <strong>Thread-Safety:</strong> This class is not thread safe since it is
 * mutable.
 * </p>
 *
 * @author albertwang, TCSASSEMBLER
 * @version 1.0
 */
@Entity
@Table(name = "detection_item_frame")
public class DetectionItemFrame {

    /**
     * <p>
     * Represents the composite primary key.
     * </p>
     * <p>
     * It is required.
     * </p>
     */
    @EmbeddedId
    private DetectionItemFramePk detectionItemFramePk;

    /**
     * <p>
     * Represents detection item of this frame.
     * </p>
     * <p>
     * It is required.
     * </p>
     */
    @MapsId("detectionItemId")
    @JoinColumn(name = "detection_item_id", referencedColumnName = "id")
    @ManyToOne
    @JsonIgnore
    private DetectionItem detectionItem;

    /**
     * Represents the right ascension.
     *
     * Required. [0, 23] hours.
     */
    @Column(name = "right_ascension", nullable = false)
    @Min(0)
    @Max(23)
    private double rightAscension;

    /**
     * Represents the declination.
     *
     * Required. -90 ~ 90 degrees
     */
    @Column(name = "declination", nullable = false)
    @Min(-90)
    @Max(90)
    private double declination;

    /**
     * Represents the rough magnitude.
     *
     * Required. Float number.
     */
    @Column(name = "rough_magnitude", nullable = false)
    private double roughMagnitude;

    /**
     * Represents the observation latitude.
     *
     * Required. -90 ~ 90 degrees
     */
    @Column(name = "observation_latitude", nullable = false)
    @Min(-90)
    @Max(90)
    private double observationLatitude;

    /**
     * Represents the observation longitude.
     *
     * Required. -180 ~ 180 degrees
     */
    @Column(name = "observation_longitude", nullable = false)
    @Min(-180)
    @Max(180)
    private double observationLongitude;

    /**
     * Represents the file path of the visualization image.
     *
     * Required.
     */
    @Column(name = "visualization_image", nullable = false)
    private String visualizationImage;
    
    @Column(name = "x", nullable = false)
    private int x;
    
    @Column(name = "y", nullable = false)
    private int y;

    /**
     * <p>
     * The default constructor for class <code>DetectionItemFrame</code>.
     * </p>
     */
    public DetectionItemFrame() {
        // does nothing
    }

    /**
     * <p>
     * Gets detection item of this frame.
     * </p>
     * @return detection item of this frame.
     */
    public DetectionItem getDetectionItem() {
        return detectionItem;
    }

    /**
     * <p>
     * Sets detection item of this frame.
     * </p>
     * @param detectionItem detection item of this frame.
     */
    public void setDetectionItem(DetectionItem detectionItem) {
        this.detectionItem = detectionItem;
    }

    /**
     * <p>
     * Gets the composite primary key.
     * </p>
     * @return the composite primary key.
     */
    public DetectionItemFramePk getDetectionItemFramePk() {
        return detectionItemFramePk;
    }

    /**
     * <p>
     * Sets the composite primary key.
     * </p>
     * @param detectionItemFramePk the composite primary key.
     */
    public void setDetectionItemFramePk(DetectionItemFramePk detectionItemFramePk) {
        this.detectionItemFramePk = detectionItemFramePk;
    }

    /**
     * <p>
     * Gets the right ascension.
     * </p>
     * @return the right ascension.
     */
    public double getRightAscension() {
        return rightAscension;
    }

    /**
     * <p>
     * Sets the right ascension.
     * </p>
     * @param rightAscension the right ascension.
     */
    public void setRightAscension(double rightAscension) {
        this.rightAscension = rightAscension;
    }

    /**
     * <p>
     * Gets the declination.
     * </p>
     * @return the declination.
     */
    public double getDeclination() {
        return declination;
    }

    /**
     * <p>
     * Sets the declination.
     * </p>
     * @param declination the declination.
     */
    public void setDeclination(double declination) {
        this.declination = declination;
    }

    /**
     * <p>
     * Gets the rough magnitude.
     * </p>
     * @return the rough magnitude.
     */
    public double getRoughMagnitude() {
        return roughMagnitude;
    }

    /**
     * <p>
     * Sets the rough magnitude.
     * </p>
     * @param roughMagnitude the rough magnitude.
     */
    public void setRoughMagnitude(double roughMagnitude) {
        this.roughMagnitude = roughMagnitude;
    }

    /**
     * <p>
     * Gets the observation latitude.
     * </p>
     * @return the observation latitude.
     */
    public double getObservationLatitude() {
        return observationLatitude;
    }

    /**
     * <p>
     * Sets the observation latitude.
     * </p>
     * @param observationLatitude the observation latitude.
     */
    public void setObservationLatitude(double observationLatitude) {
        this.observationLatitude = observationLatitude;
    }

    /**
     * <p>
     * Gets the observation longitude.
     * </p>
     * @return  the observation longitude.
     */
    public double getObservationLongitude() {
        return observationLongitude;
    }

    /**
     * <p>
     * Sets  the observation longitude.
     * </p>
     * @param observationLongitude  the observation longitude.
     */
    public void setObservationLongitude(double observationLongitude) {
        this.observationLongitude = observationLongitude;
    }

    /**
     * <p>
     * Gets the file path of the visualization image.
     * </p>
     * @return the file path of the visualization image.
     */
    public String getVisualizationImage() {
        return visualizationImage;
    }

    /**
     * <p>
     * Sets the file path of the visualization image.
     * </p>
     * @param visualizationImage the file path of the visualization image.
     */
    public void setVisualizationImage(String visualizationImage) {
        this.visualizationImage = visualizationImage;
    }

    public int getX() {
        return x;
    }

    public void setX(int x) {
        this.x = x;
    }

    public int getY() {
        return y;
    }

    public void setY(int y) {
        this.y = y;
    }
}

