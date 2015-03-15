/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package gov.nasa.asteroid.hunter.services.impl;

import gov.nasa.asteroid.hunter.ConfigurationException;
import gov.nasa.asteroid.hunter.Helper;
import gov.nasa.asteroid.hunter.LoggingHelper;
import gov.nasa.asteroid.hunter.models.DetectionItem;
import gov.nasa.asteroid.hunter.models.DetectionItemFrame;
import gov.nasa.asteroid.hunter.models.DetectionItemFramePk;
import gov.nasa.asteroid.hunter.models.DetectionItemSearchCriteria;
import gov.nasa.asteroid.hunter.models.NEO;
import gov.nasa.asteroid.hunter.models.NEOSearchCriteria;
import gov.nasa.asteroid.hunter.models.SearchResult;
import gov.nasa.asteroid.hunter.models.SortType;
import gov.nasa.asteroid.hunter.services.AsteroidDetectionService;
import gov.nasa.asteroid.hunter.services.DetectionSessionService;
import gov.nasa.asteroid.hunter.services.EntityNotFoundException;
import gov.nasa.asteroid.hunter.services.NEOService;
import gov.nasa.asteroid.hunter.services.ServiceException;
import gov.nasa.asteroid.tester.AsteroidDetectorTester;
import gov.nasa.asteroid.tester.AsteroidDetectorTester.ImageSet;
import gov.nasa.asteroid.tester.DetectionProgressListener;
import gov.nasa.asteroid.tester.ImageAlreadyExistException;
import gov.nasa.asteroid.tester.ImageService;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.annotation.PostConstruct;
import javax.persistence.Query;

import org.apache.log4j.Logger;
import org.springframework.transaction.annotation.Transactional;

/**
 * <p>
 * This class is the implementation of AsteroidDetectionService.
 * </p>
 * <p>
 * This service provides methods to detect asteroids and access asteroid
 * detection results.
 * </p>
 *
 * <p>
 * <strong>Thread Safety:</strong> This class is effectively thread safe
 * (injected configurations are not considered as thread safety factor).
 * </p>
 *
 * @author albertwang, TCSASSEMBLER
 * @version 1.0
 */
public class AsteroidDetectionServiceImpl extends BasePersistenceService
    implements AsteroidDetectionService {
    

    /**
     * <P>
     * Represents the name of the class for logging.
     * </p>
     */
    private static final String CLASS_NAME = AsteroidDetectionServiceImpl.class.getName();

    /**
     * <p>
     * Represents the minutes in an hour.
     * </p>
     */
    private static final int MINUTES_PER_HOUR = 60;

    /**
     * <p>
     * Represents the time end parameter.
     * </p>
     */
    private static final String PARAM_TIME_END = "timeEnd";

    /**
     * <p>
     * Represents the time start parameter.
     * </p>
     */
    private static final String PARAM_TIME_START = "timeStart";

    /**
     * <p>
     * Represents the time end query where cause.
     * </p>
     */
    private static final String WHERE_CAUSE_TIME_END = 
            " AND HOUR(e.timestamp) * 60 + MINUTE(e.timestamp) <= :timeEnd";

    /**
     * <p>
     * Represents the time start query where cause.
     * </p>
     */
    private static final String WHERE_CAUS_TIME_START = 
            " AND HOUR(e.timestamp) * 60 +  MINUTE(e.timestamp) >= :timeStart";

    /**
     * <p>
     * Represents the maximum minute.
     * </p>
     */
    private static final int MAX_MINUTE = 59;

    /**
     * <p>
     * Represents the maximum hour.
     * </p>
     */
    private static final int MAX_HOUR = 23;
    
    /**
     * <p>
     * Represents the suffix of the image name.
     * </p>
     */
    private static final String IMAGE_NAME_SUFFIX = ".arch.H";

    /**
     * <p>
     * Represents the middle of the name of the image.
     * </p>
     */
    private static final String IMAGE_NAME_MIDDLE = "_000";

    /**
     * <p>
     * Represents the folder name time format.
     * </p>
     */
    private static final String FOLDER_DATE_FORMAT = "yyyyMMddhhmmssS";

    /**
     * <p>
     * Represents the resource path of the trained data.
     * </p>
     */
    private static final String TRAINED_DATA_RESOURCE_PATH = "/trained-data";

    /**
     * <p>
     * Represents the prefix of the image name.
     * </p>
     */
    private static final String IMAGE_NAME_PREFIX = "INPUT";

    /**
     * <p>
     * Represents the test data file name.
     * </p>
     */
    private static final String TESTDATA_FILE_NAME = "testdata.txt";

    /**
     * <p>
     * Represents the number of the frames.
     * </p>
     */
    private static final int NUMBER_OF_FRAMES = 4;

    /**
     * <p>
     * Represents the parameter name of submitted.
     * </p>
     */
    private static final String PARAM_SUBMITTED = "submitted";

    /**
     * <p>
     * Represents the parameter name of date end.
     * </p>
     */
    private static final String PARAM_DATE_END = "timestampEnd";

    /**
     * <p>
     * Represents the parameter name of date start.
     * </p>
     */
    private static final String PARAM_DATE_START = "timestampStart";


    /**
     * <p>
     * Represents the where cause for date end.
     * </p>
     */
    private static final String WHERE_CAUSE_DATE_END = " AND e.timestamp < :timestampEnd";

    /**
     * <p>
     * Represents the where cause for date start.
     * </p>
     */
    private static final String WHERE_CAUSE_DATE_START = " AND e.timestamp >= :timestampStart";

    /**
     * <p>
     * Represents the prefix of the where cause
     * </p>
     */
    private static final String WHERE_CAUSE_PREFIX = "1 = 1";

    /**
     * <p>
     * Represents max hour value.
     * </p>
     */
    private static final int MAX_HOUR_VALUE = 23;

    /**
     * <p>
     * Represents the max minute value
     * </p>
     */
    private static final int MAX_MINUTE_VALUE = 59;

    /**
     * <p>
     * Represents the options to run the asteroid command.
     * </p>
     */
    private static final String ASTEROIRD_COMMAND_OPTIONS = " --mode test --folder %1$s";



    /**
     * <p>
     * Represents the full path of the asteroids detection program.
     * It should be something like "/home/test/detector".
     * </p>
     *
     * Required. Not null/empty.
     */
    private String asteroidsDetectionExe;

    /**
     * <p>
     * Represents the base directory used to save image files and other asteroid
     * detection data files.
     * </p>
     *
     * Required. Not null/empty.
     */
    private String baseDirectory;

    /**
     * <p>
     * Represents the max number of answers returned by the detector.
     * </p>
     */
    private int maxNumberOfAnswers = 100000;


    /**
     * <p>
     * Represents the injected detection session service.
     * </p>
     */
    private DetectionSessionService detectionSessionService;

    private ImageService imageService;
    
    private NEOService neoService;
    
    /**
     * <p>
     * This is the default constructor of <code>AsteroidDetectionServiceImpl</code>.
     * </p>
     */
    public AsteroidDetectionServiceImpl() {
        // does nothing
    }

    /**
     * <p>
     * Check if all required fields are initialized properly.
     * </p>
     * Note, in this class, asteroidsDetectionExe and baseDirectory should be not null or empty.
     * @throws ConfigurationException if any required field is not initialized properly.
     */
    @PostConstruct
    @Override
    protected void checkConfiguration() {
        // prepare for logging
        Logger logger = getLogger();
        final String signature = CLASS_NAME + ".checkConfiguration()";

        // log the entrance
        LoggingHelper.logEntrance(logger, signature, null, null);

        // call the super checkConfiguration
        super.checkConfiguration();

        // validate the injections
        Helper.checkConfigurationNullOrEmpty(logger, signature,
                "asteroidsDetectionExe", asteroidsDetectionExe);
        Helper.checkConfigurationNullOrEmpty(logger, signature,
                "baseDirectory", baseDirectory);
        Helper.checkNull(logger, signature, "detectionSessionService", detectionSessionService);

        Helper.checkPositive(logger, signature, "maxNumberOfAnswers", maxNumberOfAnswers);

        // log the exit
        LoggingHelper.logExit(logger, signature, null);

    }

    /**
     * <p>
     * This method is used to search asteroid detection items.
     * </p>
     *
     * @param criteria the search criteria.
     *
     * @return the search result.
     *
     * @throws IllegalArgumentException if criteria is null or invalid
     * (see <code>DetectionItemSearchCriteria</code> document for details).
     * @throws ServiceException if any other error occurred during the operation
     *
     */
    @Override
    public SearchResult<DetectionItem> searchDetectionItems(DetectionItemSearchCriteria criteria)
            throws ServiceException {

        // prepare for logging
        Logger logger = getLogger();
        final String signature = CLASS_NAME + ".searchDetectionItems(DetectionItemSearchCriteria)";

        // log the entrance
        LoggingHelper.logEntrance(logger, signature, new String[] { "criteria" }, new Object[] { criteria });

        // validate the parameters
        Helper.checkBaseSearchParameters(logger, signature, "criteria", criteria,
                Arrays.asList("id", "time", "timestamp", "date", "neo", 
                        "frame.rightAscension", 
                        "frame.declination", 
                        "frame.observationLatitude", 
                        "frame.observationLongitude", PARAM_SUBMITTED, "knownByMPC", "significance"));

        checkRange(logger, signature, "creteria.minuteStart", criteria.getMinuteStart(), 0, MAX_MINUTE_VALUE);
        checkRange(logger, signature, "creteria.minuteEnd", criteria.getMinuteEnd(), 0, MAX_MINUTE_VALUE);
        checkRange(logger, signature, "creteria.hourStart", criteria.getHourStart(), 0, MAX_HOUR_VALUE);
        checkRange(logger, signature, "creteria.hourEnd", criteria.getHourEnd(), 0, MAX_HOUR_VALUE);

        StringBuffer sb = new StringBuffer(WHERE_CAUSE_PREFIX);
        Map<String, Object> queryParameters = new HashMap<String, Object>();

        // Append criteria
        Calendar calendar = Calendar.getInstance();
        if (criteria.getDateStart() != null) {
            calendar.setTime(criteria.getDateStart());
            calendar.set(Calendar.HOUR_OF_DAY, 0);
            calendar.set(Calendar.MINUTE, 0);
            calendar.set(Calendar.SECOND, 0);
            calendar.set(Calendar.MILLISECOND, 0);

            sb.append(WHERE_CAUSE_DATE_START);
            queryParameters.put(PARAM_DATE_START, calendar.getTime());
        }
        if (criteria.getDateEnd() != null) {
            calendar.setTime(criteria.getDateEnd());
            calendar.add(Calendar.DATE, 1);
            calendar.set(Calendar.HOUR_OF_DAY, 0);
            calendar.set(Calendar.MINUTE, 0);
            calendar.set(Calendar.SECOND, 0);
            calendar.set(Calendar.MILLISECOND, 0);
            sb.append(WHERE_CAUSE_DATE_END);
            queryParameters.put(PARAM_DATE_END, calendar.getTime());
        }
        
        int hourStart = 0;
        int miniteStart = 0;
        int hourEnd = MAX_HOUR;
        int miniteEnd = MAX_MINUTE;
        
        if (criteria.getHourStart() != null || criteria.getMinuteStart() != null) {
            if (criteria.getHourStart() != null) {
                hourStart = criteria.getHourStart();
            }
            if (criteria.getMinuteStart() != null) {
                miniteStart = criteria.getMinuteStart();
            }
            sb.append(WHERE_CAUS_TIME_START);
            queryParameters.put(PARAM_TIME_START, hourStart * MINUTES_PER_HOUR + miniteStart);
        }
        
        if (criteria.getHourEnd() != null || criteria.getMinuteEnd() != null) {
            if (criteria.getHourEnd() != null) {
                hourEnd = criteria.getHourEnd();
            }
            if (criteria.getMinuteEnd() != null) {
                miniteEnd = criteria.getMinuteEnd();
            }
            sb.append(WHERE_CAUSE_TIME_END);
            queryParameters.put(PARAM_TIME_END, hourEnd * MINUTES_PER_HOUR + miniteEnd);
        }
        
        if (criteria.isSubmitted() != null) {
            sb.append(" AND e.submitted = :submitted");
            queryParameters.put("submitted", criteria.isSubmitted());
        }
        
        SearchResult<DetectionItem> result = null;
        String queryString = null;
        String countQueryString = null;
        
        // modify the order by
        if (criteria.getSortType() == null) {
            criteria.setSortType(SortType.ASC);
        }
        if ("date".equals(criteria.getSortBy())) {
            queryString = "SELECT e FROM DetectionItem e WHERE %1$s ORDER BY TRUNC(timestamp) %2$s, id ASC";
            queryString = String.format(queryString,  sb.toString(), criteria.getSortType().name());
            countQueryString = "SELECT COUNT(e) FROM DetectionItem e WHERE %1$s";
            countQueryString = String.format(countQueryString,  sb.toString());
        } else if ("time".equals(criteria.getSortBy())) {
            queryString = "SELECT e FROM DetectionItem e WHERE %1$s ORDER BY (HOUR(e.timestamp) * 60 + MINUTE(e.timestamp)) %2$s, id ASC";
            queryString = String.format(queryString,  sb.toString(), criteria.getSortType().name());
            countQueryString = "SELECT COUNT(e) FROM DetectionItem e WHERE %1$s";
            countQueryString = String.format(countQueryString,  sb.toString());
        } else if (criteria.getSortBy() != null && criteria.getSortBy().startsWith("frame.")) {
            String sortBy = criteria.getSortBy().substring(6);
            queryString = "SELECT e FROM DetectionItemFrame f INNER JOIN f.detectionItem e WHERE %1$s AND f.detectionItemFramePk.frame = 0 ORDER BY f.%2$s %3$s, e.id ASC";
            queryString = String.format(queryString,  sb.toString(), sortBy, criteria.getSortType().name());
            countQueryString = "SELECT COUNT(e) FROM DetectionItemFrame f INNER JOIN f.detectionItem e WHERE %1$s AND f.detectionItemFramePk.frame = 0";
            countQueryString = String.format(countQueryString,  sb.toString());
        }
        
        try {
            if (queryString != null && countQueryString != null) {
                result =  searchWithQuery(criteria, queryString, countQueryString, queryParameters, DetectionItem.class);
            } else {
                result = search(criteria, sb.toString(), queryParameters, DetectionItem.class);
            }
        } catch (ServiceException e) {
            // log and re-throw
            throw LoggingHelper.logException(logger, signature, e);
        }
        
        // log the exit
        LoggingHelper.logExit(logger, signature, new Object[] { result });

        return result;
    }

    /**
     * <p>
     * This method is used to retrieve a detection item.
     * </p>
     *
     * @param id the ID
     * @return  the detection item, null will be returned if there is no such entity.
     * @throws ServiceException if there are any errors.
     * @throws IllegalArgumentException if the id is not positive.
     */
    public DetectionItem getDetectionItem(long id) throws ServiceException {
        // prepare for logging
        Logger logger = getLogger();
        final String signature = CLASS_NAME + ".getDetectionItem(long)";

        // log the entrance
        LoggingHelper.logEntrance(logger, signature, new String[] { "id" }, new Object[] { id });

        // check the parameters
        Helper.checkPositive(logger, signature, "id", id);

        try {
            DetectionItem result = getEntityManager().find(DetectionItem.class, id);

            // log the exit
            LoggingHelper.logExit(logger, signature, new Object[] { result });

            return result;
        } catch (IllegalArgumentException e) {
            throw LoggingHelper.logException(logger, signature, new ServiceException("Failed to get the result.", e));
        }
    }

    /**
     * <p>
     * This method is used to submit a detection item (mark the detection item as submitted).
     * </p>
     *
     * @param id the ID
     *
     * @throws IllegalArgumentException if id is not positive
     * @throws EntityNotFoundException if there is no such detection item
     * @throws ServiceException if any other error occurred during the operation
     *
     */
    @Transactional
    @Override
    public void submitDetectionItem(long id) throws ServiceException {
        // prepare for logging
        Logger logger = getLogger();
        final String signature = CLASS_NAME + ".submitDetectionItem(long)";

        // log the entrance
        LoggingHelper.logEntrance(logger, signature, new String[] { "id" }, new Object[] { id });

        // check the parameters
        Helper.checkPositive(logger, signature, "id", id);

        try {
            DetectionItem item = getDetectionItem(id);
            if (item == null) {
                throw LoggingHelper.logException(logger, signature, new EntityNotFoundException(
                        "Detection Item with id " + id + " does not exist."));
            }
            // mark it submitted
            item.setSubmitted(true);

            // update the item.
            merge(logger, signature, item);

            // log the exit
            LoggingHelper.logExit(logger, signature, null);

        } catch (ServiceException e) {
            // log and re-throw
            throw LoggingHelper.logException(logger, signature, e);
        }
    }

    @Transactional
    @Override
    public void deleteAllDetectionItems() throws ServiceException {
        String hql = String.format("delete from %s", "DetectionItemFrame");
        Query query = getEntityManager().createQuery(hql);
        query.executeUpdate();
        
        hql = String.format("delete from %s", "DetectionItem");
        query = getEntityManager().createQuery(hql);
        query.executeUpdate();
        
        hql = String.format("delete from %s", "ImageInfo");
        query = getEntityManager().createQuery(hql);
        query.executeUpdate();
    }
    
    /**
     * <p>
     * This method is used to detect asteroids in FITS images.
     * </p>
     *
     * <p>
     * This method should use READ_UNCOMMITTED transaction isolation level so
     * that progress information can be read while this method is executing.
     * </p>
     *
     * @param fitsImages
     *            the FITS image files, must be an array of exactly 4 items
     * @param detectionSessionId
     *            the ID of the detection session.
     * @return the detection result
     *
     * @throws IllegalArgumentException
     *             if fitsImages is null, or does not contain 4 elements, or any
     *             item is null or detectionSessionId is not positive.
     * @throws ServiceException
     *             if any other error occurred during the operation
     * @throws ImageAlreadyExistException 
     *
     */
    @Transactional
    @Override
    public List<DetectionItem> detectAsteroids(List<File> fitsImages, long detectionSessionId, boolean forced, String observatoryCode)
        throws ServiceException, ImageAlreadyExistException {

        // prepare for logging
        Logger logger = getLogger();
        final String signature = CLASS_NAME + ".detectAsteroids(List<File>, long)";

        // log the entrance
        LoggingHelper.logEntrance(logger, signature,
                new String[] {"fitsImages", "detectionSessionId"},
                new Object[] {fitsImages, detectionSessionId});

        // validate the parameters
        Helper.checkNull(logger, signature, "fitsImages", fitsImages);
        if (fitsImages.size() != NUMBER_OF_FRAMES) {
            throw LoggingHelper.logException(logger, signature,
                    new IllegalArgumentException("fitsImages should contains 4 files."));
        }
        String fileExtension = null;
        for (File file : fitsImages) {
            Helper.checkNull(logger, signature, "Element in fitsImages", file);
            String extension = getExtension(file.getName());
            if (extension == null) {
                throw LoggingHelper.logException(logger, signature,
                        new IllegalArgumentException("Unsupported extension.")); 
            } else if (fileExtension == null) {
                fileExtension = extension;
            } else if (!fileExtension.equals(extension)) {
                throw LoggingHelper.logException(logger, signature,
                        new IllegalArgumentException("There are multiple file types in uploaded files.")); 
            }
        }
        Helper.checkPositive(logger, signature, "detectionSessionId", detectionSessionId);

        // Copy input FITS images to the working directory
        File directory = prepareDataDirectory(logger, signature, fitsImages, new Date(), fileExtension);

        // create the tester
        AsteroidDetectorTester tester = createAsteroidDetectorTester(detectionSessionId, directory, fileExtension);
        tester.setForced(forced);
        tester.setObservatoryCode(observatoryCode);
        
        // Execute AsteroidDetectorTester to get the entities
        List<AsteroidDetectorTester.Entity> entities;
        try {
            entities = tester.doExec();
        } catch (ImageAlreadyExistException e) {
            // the image is already exist
            throw LoggingHelper.logException(logger, signature, e);
        } catch (Exception e) {
            throw LoggingHelper.logException(logger, signature, new ServiceException("Error occurs on execute.", e));
        }

        ImageSet imageSet = tester.getCurrentImageSet();
        double imageWidth = imageSet.img[0].W;
        double imageHeight = imageSet.img[0].H;
        HashMap<String, List<NEO>> mpcItems = new HashMap<String, List<NEO>>();
        // get the neos of the 4 frame
        for (int i = 0; i < 4; i++) {
            // read the image create time
            DateFormat sdf = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss.SSS");
            Date date;
            try {
                date = sdf.parse(imageSet.img[i].getDateOBS() + " " + imageSet.img[i].getTimeOBS());
            } catch (ParseException e) {
                throw LoggingHelper.logException(logger, signature, new ServiceException("Cannot parse the date.", e));
            }
            
            NEOSearchCriteria criteria = new NEOSearchCriteria();
            criteria.setDate(date);
            criteria.setRightAscension(imageSet.img[i].getRA().replaceAll(":", " "));
            criteria.setDeclination(imageSet.img[i].getDEC().replaceAll(":", " "));
            criteria.setObservatoryCode(observatoryCode);
            criteria.setRadius(129);
            criteria.setV(22);
            
            SearchResult<NEO> neos = neoService.search(criteria);
            for (NEO neo : neos.getValues()) {
                List<NEO> list = mpcItems.get(neo.getObjectDesignation());
                if (list == null) {
                    list = new ArrayList<NEO>();
                    mpcItems.put(neo.getObjectDesignation(), list);
                }
                list.add(neo);
            }
        }
        Set<String> keySet = new HashSet<String> ();
        keySet.addAll(mpcItems.keySet());
        for (String key : keySet) {
            List<NEO> list = mpcItems.get(key);
            if (list.size() != 4) {
                mpcItems.remove(key);
            }
        }
        System.out.println("Number of objects in MPC: " + mpcItems.keySet().size());
        
        // persist the entities
        List<DetectionItem> items = persistEntities(logger, signature, tester.getImageCreateTime(), entities, mpcItems, 
                observatoryCode, imageWidth, imageHeight, tester.getImageId());
        // mark finished.
        detectionSessionService.updateDetectionSessionProgress(detectionSessionId, 1.0f);
        
        // log the exit
        LoggingHelper.logExit(logger, signature, new Object[] {items});

        return items;
    }

    private String getExtension(String fileName) {
        if (fileName.endsWith(".arch.H")) {
            return ".arch.H";
        } else if (fileName.endsWith(".fit")) {
            return ".fit";
        } else if (fileName.endsWith(".fits")) {
            return ".fits";
        } else {
            return null;
        }
    }
    /**
     * <p>
     * Persist the detected items to Database.
     * </p>
     * @param logger the logger for logging.
     * @param signature the signature of the caller for logging.
     * @param timestamp the timestamp.
     * @param entities the detected entities.
     * @param mpcItems 
     * @param observatoryCode 
     * @param imageHeight 
     * @param imageWidth 
     * @param imageId 
     * @return the detection items list.
     * @throws ServiceException if there are any errors.
     */
    private List<DetectionItem> persistEntities(Logger logger, String signature, Date timestamp,
            List<AsteroidDetectorTester.Entity> entities, Map<String, List<NEO>> mpcItems, String observatoryCode, double imageWidth, 
            double imageHeight, long imageId) throws ServiceException {
        // check the entity
        if (entities == null) {
            throw LoggingHelper.logException(logger, signature,
                    new ServiceException("Entities cannot be null."));
        }
        // Create results
        List<DetectionItem> items = new ArrayList<DetectionItem>();

        for (AsteroidDetectorTester.Entity entity : entities) {
            DetectionItem item = new DetectionItem();
            item.setTimestamp(timestamp);
            item.setNeo(entity.isNeo);
            item.setSignificance(entity.significance);
            item.setObservatoryCode(observatoryCode);
            item.setImageHeight(imageWidth);
            item.setImageWidth(imageHeight);
            item.setImageId(imageId);
            List<DetectionItemFrame> frames = new ArrayList<DetectionItemFrame>();

            for (int i = 0; i < NUMBER_OF_FRAMES; i++) {
                DetectionItemFrame frame = new DetectionItemFrame();

                DetectionItemFramePk detectionItemFramePk = new DetectionItemFramePk();
                detectionItemFramePk.setFrame(i);
                frame.setDetectionItemFramePk(detectionItemFramePk);
                frame.setRightAscension(entity.RA[i]);
                frame.setDeclination(entity.DEC[i]);
                frame.setRoughMagnitude(entity.mag[i]);
                frame.setX(entity.x[i]);
                frame.setY(entity.y[i]);
                frame.setVisualizationImage(entity.getVisualizationImages()[i]);
                frame.setObservationLatitude(entity.getLatitudes()[i]);
                frame.setObservationLongitude(entity.getLongitudes()[i]);
                frame.setDetectionItem(item);
                
                frames.add(frame);
            }
            item.setFrames(frames);

            // check if it is existing in MPC
            for (String key : mpcItems.keySet()) {
                List<NEO> list = mpcItems.get(key);
                double sum = 0;
                for (int i = 0; i < NUMBER_OF_FRAMES; i++) {
                    double RA = entity.RA[i];
                    double DEC = entity.DEC[i];
                    NEO neo = list.get(i);
                    double RA_neo = neo.getRightAscensionValue();
                    double DEC_neo = neo.getDeclinationValue();
                    sum += (RA - RA_neo) * (RA - RA_neo) + (DEC - DEC_neo) * (DEC - DEC_neo);
                    if (sum > 0.01) {
                        break;
                    }
                }
                if (sum < 0.01) {
                    item.setSubmitted(true);
                }
            }
            items.add(item);
            // Persist in database
            persist(logger, signature, item);
        }

        // return the result
        return items;
    }

    /**
     * <p>
     * Creates the asteroid detector tester and sets the options.
     * </p>
     * @param detectionSessionId the session id.
     * @param directory the working directory.
     *
     * @return the created asteroid detector tester.
     */
    private AsteroidDetectorTester createAsteroidDetectorTester(
        final long detectionSessionId, File directory, String fileExtension) {
        
        String trainedDataFolder = new File(getClass().getResource(
                TRAINED_DATA_RESOURCE_PATH).getFile()).getAbsolutePath() + File.separator;
        String cmd = "\"" + asteroidsDetectionExe + "\"" + String.format(ASTEROIRD_COMMAND_OPTIONS, trainedDataFolder);

        // Create AsteroidDetectorTester
        AsteroidDetectorTester tester = new AsteroidDetectorTester();
        tester.setVisualize(true);
        tester.setVisnr(-1);
        tester.setFolder(directory.getAbsolutePath() + File.separator);
        tester.setTestFile(new File(directory, TESTDATA_FILE_NAME).getAbsolutePath());
        tester.setExecCommand(cmd);
        tester.setTrainFolder(trainedDataFolder);
        tester.setMaxNumberOfAnswer(maxNumberOfAnswers);
        tester.setImageService(imageService);
        tester.setFileExtension(fileExtension);
        final DetectionSessionService service = detectionSessionService;

        // Set progress listener
        tester.setDetectionProgressListener(new DetectionProgressListener() {
            public void updateProgress(float progress) {
                // update the progress
                try {
                    service.updateDetectionSessionProgress(detectionSessionId, progress);
                } catch (ServiceException e) {
                    // ignore
                }
            }
        });

        // return the tester
        return tester;
    }

    /**
     * <p>
     * Creates the working directory to store the necessary data.
     * </p>
     *
     * @param logger
     *            the logger for logging.
     * @param signature
     *            the name of the caller for logging.
     * @param fitsImages
     *            the fitsImages files.
     * @param timestamp
     *            the timestamp.
     * @param fileExtension 
     * @return the path of the directory.
     * @throws ServiceException
     *             if there are any errors.
     */
    private File prepareDataDirectory(Logger logger, String signature, List<File> fitsImages, Date timestamp, String fileExtension)
            throws ServiceException {

        // create the working directory
        File directory = new File(baseDirectory, new SimpleDateFormat(FOLDER_DATE_FORMAT).format(timestamp));
        if (!directory.mkdirs()) {
            throw LoggingHelper.logException(logger, signature,
                    new ServiceException("Cannot create the working directory:" + directory.getAbsolutePath()));
        }

        FileWriter fw = null;
        try {
            // copy the fits images
            int index = 0;
            for (File fitsImage : fitsImages) {
                org.springframework.util.FileCopyUtils.copy(
                        fitsImage, new File(directory,
                                IMAGE_NAME_PREFIX + IMAGE_NAME_MIDDLE + (++index) + fileExtension));
            }

            // Generate test file, only one image set
            fw = new FileWriter(new File(directory, TESTDATA_FILE_NAME));
            fw.write(IMAGE_NAME_PREFIX);
        } catch (IOException e) {
            throw LoggingHelper.logException(logger, signature,
                    new ServiceException("Error occurs while copying the files or writing test data file.", e));
        } finally {
            if (fw != null) {
                try {
                    fw.close();
                } catch (IOException e) {
                    // ignore
                }
            }
        }
        // return the created directory
        return directory;
    }



    /**
     * <p>
     * This method is used to mark a detection item as known by MPC.
     * </p>
     *
     * @param id the ID
     *
     * @throws IllegalArgumentException if id is not positive
     * @throws EntityNotFoundException if there is no such detection item
     * @throws ServiceException if any other error occurred during the operation
     *
     */
    @Transactional
    @Override
    public void markDetectionItemKnownByMPC(long id) throws ServiceException {

        // prepare for logging
        Logger logger = getLogger();
        final String signature = CLASS_NAME + ".markDetectionItemKnownByMPC(long)";

        // log the entrance
        LoggingHelper.logEntrance(logger, signature, new String[] { "id" }, new Object[] { id });

        // check the parameters
        Helper.checkPositive(logger, signature, "id", id);

        try {
            DetectionItem item = getDetectionItem(id);
            if (item == null) {
                throw LoggingHelper.logException(logger, signature, new EntityNotFoundException(
                        "Detection Item with id " + id + " does not exist."));
            }
            // mark it
            item.setKnownByMPC(true);

            // update the item.
            merge(logger, signature, item);

            // log the exit
            LoggingHelper.logExit(logger, signature, null);

        } catch (ServiceException e) {
            // log and re-throw
            throw LoggingHelper.logException(logger, signature, e);
        }

    }

    /**
     * <p>
     * Sets the the full path of the asteroids detection program.
     * </p>
     * @param asteroidsDetectionExe the full path of the asteroids detection program.
     */
    public void setAsteroidsDetectionExe(String asteroidsDetectionExe) {
        this.asteroidsDetectionExe = asteroidsDetectionExe;
    }

    /**
     * <p>
     * Sets the base directory used to save image files and other asteroid
     * detection data files.
     * </p>
     * @param baseDirectory the base directory.
     */
    public void setBaseDirectory(String baseDirectory) {
        this.baseDirectory = baseDirectory;
    }

    /**
     * <p>
     * Sets the max number of answers.
     * </p>
     * @param maxNumberOfAnswers max number of answers.
     */
    public void setMaxNumberOfAnswers(int maxNumberOfAnswers) {
        this.maxNumberOfAnswers = maxNumberOfAnswers;
    }

    /**
     * <p>
     * Checks if a value is within the range.
     * </p>
     * @param logger the logger for logging.
     * @param signature the signature of the caller method for logging.
     * @param name the name of the value.
     * @param value the value.
     * @param minValue the minimum value.
     * @param maxValue the maximum value.
     *
     * @throws IllegalArgumentException if the value is not in the range.
     */
    private static void checkRange(Logger logger, String signature,
            String name, Integer value, int minValue, int maxValue) {
        if (value != null && (value < minValue || value > maxValue)) {
            throw LoggingHelper.logException(logger, signature, new IllegalArgumentException(
                    name + " is not in the correct range:[" + minValue + "," + maxValue + "]."));
        }
    }

    /**
     * <p>
     * Sets the detection session service.
     * </p>
     * @param detectionSessionService the detection session service to set.
     */
    public void setDetectionSessionService(DetectionSessionService detectionSessionService) {
        this.detectionSessionService = detectionSessionService;
    }

    public ImageService getImageService() {
        return imageService;
    }

    public void setImageService(ImageService imageService) {
        this.imageService = imageService;
    }

    public NEOService getNeoService() {
        return neoService;
    }

    public void setNeoService(NEOService neoService) {
        this.neoService = neoService;
    }

}

