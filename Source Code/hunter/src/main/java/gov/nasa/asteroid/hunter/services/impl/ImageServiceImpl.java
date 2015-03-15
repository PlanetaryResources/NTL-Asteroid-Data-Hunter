package gov.nasa.asteroid.hunter.services.impl;

import gov.nasa.asteroid.hunter.models.ImageInfo;
import gov.nasa.asteroid.hunter.services.ServiceException;
import gov.nasa.asteroid.tester.AsteroidDetectorTester.FITSImage;
import gov.nasa.asteroid.tester.ImageService;

import java.util.Date;
import java.util.List;

import javax.persistence.EntityManager;
import javax.persistence.TypedQuery;

import org.apache.log4j.Logger;
import org.springframework.transaction.annotation.Transactional;

public class ImageServiceImpl extends BasePersistenceService implements ImageService {

    /**
     * <P>
     * Represents the name of the class for logging.
     * </p>
     */
    private static final String CLASS_NAME = ImageServiceImpl.class.getName();
    
    public ImageServiceImpl() {
        
    }
    
    @Transactional
    @Override
    public long addImage(FITSImage image) {
        final Logger logger = super.getLogger();
        final String signature = CLASS_NAME + ".addImage(FITSImage)";
        
        ImageInfo info = new ImageInfo();
        info.setDateOBS(image.getDateOBS());
        info.setTimeOBS(image.getTimeOBS());
        info.setDeclination(image.getDEC());
        info.setRightAscension(image.getRA());
        info.setImageCreateDate(image.getImageCreateTime());
        info.setTimestamp(new Date());
        info.setObservatoryCode(image.getObservatoryCode());
        
        try {
            persist(logger, signature, info);
            return info.getId();
        } catch (ServiceException e) {
            // ignore the exception
        }
        return 0;
    }

    @Override
    public boolean imageExists(FITSImage image) {
        EntityManager em = super.getEntityManager();
        TypedQuery<ImageInfo> query = em.createQuery("SELECT e FROM ImageInfo e " +
        		"WHERE e.dateOBS = :dateOBS AND e.timeOBS = :timeOBS AND e.declination = :declination AND e.rightAscension = :rightAscension", ImageInfo.class);
        query.setParameter("dateOBS", image.getDateOBS());
        query.setParameter("timeOBS", image.getTimeOBS());
        query.setParameter("declination", image.getDEC());
        query.setParameter("rightAscension", image.getRA());
        return query.getResultList().size() > 0;
    }
    
    @Override
    public <T> T getLatestImage(Class<T> clazz) {
        EntityManager em = super.getEntityManager();
        TypedQuery<T> query = em.createQuery("SELECT e FROM ImageInfo e ORDER BY e.timestamp DESC", clazz);
        query.setMaxResults(1);
        List<T> result = query.getResultList();
        if (result.size() == 0) {
            return null;
        }
        return result.get(0);
    }

}
