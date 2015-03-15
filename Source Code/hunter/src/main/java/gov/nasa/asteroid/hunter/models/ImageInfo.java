package gov.nasa.asteroid.hunter.models;

import java.util.Date;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.Table;

@Entity
@Table(name = "image_info")
public class ImageInfo extends IdentifiableEntity {
    @Column(name = "time_obs", nullable = false)
    private String timeOBS;
    
    @Column(name = "date_obs", nullable = false)
    private String dateOBS;
    
    @Column(name = "image_create_date", nullable = false)
    private String imageCreateDate;
    
    @Column(name = "ra", nullable = false)
    private String rightAscension;
    
    @Column(name = "dec", nullable = false)
    private String declination;
    
    @Column(name = "timestamp", nullable = false)
    private Date timestamp;
    
    @Column(name = "observatory_code", nullable = false)
    private String observatoryCode;
    
    public String getDateOBS() {
        return dateOBS;
    }
    public void setDateOBS(String dateOBS) {
        this.dateOBS = dateOBS;
    }
    public String getTimeOBS() {
        return timeOBS;
    }
    public void setTimeOBS(String timeOBS) {
        this.timeOBS = timeOBS;
    }
    public String getImageCreateDate() {
        return imageCreateDate;
    }
    public void setImageCreateDate(String imageCreateDate) {
        this.imageCreateDate = imageCreateDate;
    }

    public Date getTimestamp() {
        return timestamp;
    }
    public void setTimestamp(Date timestamp) {
        this.timestamp = timestamp;
    }
    public String getRightAscension() {
        return rightAscension;
    }
    public void setRightAscension(String rightAscension) {
        this.rightAscension = rightAscension;
    }
    public String getDeclination() {
        return declination;
    }
    public void setDeclination(String declination) {
        this.declination = declination;
    }
    public String getObservatoryCode() {
        return observatoryCode;
    }
    public void setObservatoryCode(String observatoryCode) {
        this.observatoryCode = observatoryCode;
    }
}
