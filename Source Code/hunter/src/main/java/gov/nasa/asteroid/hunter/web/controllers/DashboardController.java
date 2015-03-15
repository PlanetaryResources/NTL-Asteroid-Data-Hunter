/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package gov.nasa.asteroid.hunter.web.controllers;

import gov.nasa.asteroid.hunter.ConfigurationException;
import gov.nasa.asteroid.hunter.Helper;
import gov.nasa.asteroid.hunter.LoggingHelper;
import gov.nasa.asteroid.hunter.models.DetectionItem;
import gov.nasa.asteroid.hunter.models.DetectionItemSearchCriteria;
import gov.nasa.asteroid.hunter.models.DetectionSession;
import gov.nasa.asteroid.hunter.models.SearchResult;
import gov.nasa.asteroid.hunter.models.SortType;
import gov.nasa.asteroid.hunter.services.AsteroidDetectionService;
import gov.nasa.asteroid.hunter.services.DetectionSessionService;
import gov.nasa.asteroid.hunter.services.ServiceException;
import gov.nasa.asteroid.tester.ImageAlreadyExistException;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import javax.annotation.PostConstruct;
import javax.servlet.http.HttpServletResponse;

import org.apache.log4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.http.HttpStatus;
import org.springframework.stereotype.Controller;
import org.springframework.ui.Model;
import org.springframework.web.bind.annotation.ExceptionHandler;
import org.springframework.web.bind.annotation.PathVariable;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RequestMethod;
import org.springframework.web.bind.annotation.RequestParam;
import org.springframework.web.bind.annotation.ResponseBody;
import org.springframework.web.bind.annotation.ResponseStatus;
import org.springframework.web.multipart.MultipartFile;

/**
 * <p>
 * Spring MVC Controller for dashboard page.
 * </p>
 * 
 * <p>
 * 
 * <strong>Thread Safety:</strong> This class is effectively thread safe
 * (injected configurations are not considered as thread safety factor).
 * 
 * </p>
 * 
 * @author albertwang, TCSASSEMBLER
 * @version 1.0
 */
@Controller
public class DashboardController extends BaseController {
    
    /**
     * <p>
     * Represents the name of this class for logging.
     * </p>
     */
    private static final String CLASS_NAME = DashboardController.class.getName();
    
    /**
     * <p>
     * Represents the prefix for temp files.
     * </p>
     */
    private static final String TEMP_FILE_PREFIX = "temp_";
    
    /**
     * <p>
     * Represents the number of frames.
     * </p>
     */
    private static final int NUMBERS_OF_FRAMES = 4;

    /**
     * <p>
     * Represents the dashboard view.
     * </p>
     */
    private static final String VIEW_DASHBOARD_VIEW_NAME = "viewDashboard";

    /**
     * <p>
     * Represents the latest detection item .
     * </p>
     */
    private static final String LATEST_DETECTION_ITEM_ATTRIBUTE_NAME = "latestDetectionItem";

    /**
     * <p>
     * Represents attribute name of the search result.
     * </p>
     */
    private static final String SEARCH_RESULT_ATTRIBUTE_NAME = "searchResult";

    /**
     * <p>
     * Represents the time stamp column to sort by.
     * </p>
     */
    private static final String TIMESTAMP_COLUMN = "timestamp";
    
    /**
     * <p>
     * Represents the AsteroidDetectionService used for detecting asteroids.
     * </p>
     * <p>
     * It is Required.
     * </p>
     */
    @Autowired
    private AsteroidDetectionService asteroidDetectionService;

    /**
     * <p>
     * Represents the detection section service.
     * </p>
     * <p>
     * It is required.
     * </p>
     */
    @Autowired
    private DetectionSessionService detectionSectionService;
    
    /**
     * <p>
     * Represents the number of detection items to display on dashboard.
     * </p>
     * <p>
     * It is required, and should be positive integer.
     * </p>
     */
    @Value("${controller.detection_item_count}")
    private int detectionItemCount;

    /**
     * Empty constructor.
     */
    public DashboardController() {
        // does nothing
    }

    /**
     * <p>
     * Check if all required fields are initialized properly.
     * </p>
     * 
     * Note, asteroidDetectionService, detectionSectionService should be injected, 
     * and detectionItemCount should be positive.
     * 
     * @throws ConfigurationException if any required field is not initialized properly.
     */
    @PostConstruct
    protected void checkConfiguration() {
        // prepare for logging
        final Logger logger = getLogger();
        final String signature = CLASS_NAME + ".checkConfiguration()";
        
        LoggingHelper.logEntrance(logger, signature, null, null);
        
        super.checkConfiguration();
        
        Helper.checkConfigurationNull(logger, signature, "asteroidDetectionService", asteroidDetectionService);
        Helper.checkConfigurationNull(logger, signature, "detectionSectionService", detectionSectionService);
        
        if (detectionItemCount <= 0) {
            throw LoggingHelper.logException(logger, signature, 
                    new ConfigurationException("Detection Item Count should be positive."));
        }
        
        // log the exit
        LoggingHelper.logExit(logger, signature, null);
    }

    /**
     * <p>
     * Handles the root path.
     * </p>
     * @return the direct path
     * @throws ServiceException if there are any error
     */
    @RequestMapping(value = "/", method = RequestMethod.GET)
    @ResponseStatus(HttpStatus.OK)
    public void index(HttpServletResponse response) throws ServiceException {
     // prepare for logging
        final Logger logger = getLogger();
        final String signature = CLASS_NAME + ".index(HttpServletResponse)";
        // log the entrance
        LoggingHelper.logEntrance(logger, signature, null, null);
        try {
            response.sendRedirect("/dashboard");
        } catch (IOException e) {
            throw LoggingHelper.logException(logger, signature, 
                    new ServiceException("Failed to do the redirection.", e));
        }
        // log the exit
        LoggingHelper.logExit(logger, signature, null);
    }
    
    /**
     * <p>
     * This method processes the request to view dashboard.
     * </p>
     * 
     * @throws IllegalArgumentException if the model is null. 
     * @throws ServiceException if there are any error.
     * 
     */
    @RequestMapping(value = "dashboard", method = RequestMethod.GET)
    @ResponseStatus(HttpStatus.OK)
    public String view(Model model) throws ServiceException {
        // prepare for logging
        final Logger logger = getLogger();
        final String signature = CLASS_NAME + ".view(Model)";
        // log the entrance
        LoggingHelper.logEntrance(logger, signature, new String[] { "model" }, new Object[] { model });

        Helper.checkNull(logger, signature, "model", model);

        /// get the search result
        DetectionItemSearchCriteria criteria = new DetectionItemSearchCriteria();
        criteria.setSortBy(TIMESTAMP_COLUMN);
        criteria.setSortType(SortType.DESC);
        criteria.setPageNumber(1);
        criteria.setPageSize(detectionItemCount);
        
        SearchResult<DetectionItem> searchResult;
        try {
            searchResult = asteroidDetectionService.searchDetectionItems(criteria);
        } catch (ServiceException e) {
            // log and re-throw
            throw LoggingHelper.logException(logger, signature, e);
        }

        model.addAttribute(SEARCH_RESULT_ATTRIBUTE_NAME, searchResult);

        DetectionItem latestDetectionItem = null;
        
        if (searchResult.getValues().size() > 0) {
            long lastestItemId = searchResult.getValues().get(0).getId();
            try {
                latestDetectionItem = asteroidDetectionService.getDetectionItem(lastestItemId);
            } catch (ServiceException e) {
                // log and re-throw
                throw LoggingHelper.logException(logger, signature, e);
            }
        }
        
        model.addAttribute(LATEST_DETECTION_ITEM_ATTRIBUTE_NAME, latestDetectionItem);

        String viewName = VIEW_DASHBOARD_VIEW_NAME;

        // log the exit
        LoggingHelper.logExit(logger, signature, new Object[] { viewName });

        return viewName;
    }

    /**
     * <p>
     * This method processes the request to detect asteroids.
     * </p>
     * 
     * @param files
     *            the upload files for fits images
     * @param id
     *            the session id.
     * 
     * @return the detected asteroids (Spring Framework can be configured to
     *         serialize it to JSON)
     * @throws ImageAlreadyExistException 
     * 
     * @throws IllegalArgumentException
     *             if the files is null or not 4 files, or the id is not
     *             positive.
     * 
     */
    @RequestMapping(value="detectAsteroids", method=RequestMethod.POST, produces = "application/json")
    @ResponseBody
    public List<DetectionItem> detectAsteroids(@RequestParam("files") List<MultipartFile> files, 
        @RequestParam("id") long id, @RequestParam(required = false, defaultValue = "false", value = "forced") boolean forced, @RequestParam("observatoryCode") String observatoryCode) 
        throws ServiceException, ImageAlreadyExistException {
        
        final Logger logger = getLogger();
        final String signature = CLASS_NAME + ".detectAsteroids(List<MultipartFile>, long)";
        
        // log the entrance
        LoggingHelper.logEntrance(logger, signature, new String[] {"files", "id"}, new Object[] {files, id});
        
        // validate the parameters
        Helper.checkNull(logger, signature, "files", files);
        Helper.checkPositive(logger, signature, "id", id);
        if (files.size() != NUMBERS_OF_FRAMES) {
            throw LoggingHelper.logException(logger, signature, 
                    new IllegalArgumentException("files should have exactly 4 files."));
        }
        
        File[] imageFiles = new File[NUMBERS_OF_FRAMES];
        try {
            for (int i = 0; i < NUMBERS_OF_FRAMES; i++) {
                imageFiles[i] = File.createTempFile(TEMP_FILE_PREFIX, files.get(i).getOriginalFilename());
                files.get(i).transferTo(imageFiles[i]);
            }
            
            List<DetectionItem> result = asteroidDetectionService.detectAsteroids(Arrays.asList(imageFiles), id, forced, observatoryCode);
            
            // log the exit
            LoggingHelper.logExit(logger, signature, new Object[] { result });

            return result;
            
        } catch (IOException e) {
            // wrap log and re-throw
            throw LoggingHelper.logException(logger, signature, 
                    new ServiceException("IO error occurs while processing the uploaded files", e));
        } catch (ServiceException e) {
            // log and re-throw
            throw LoggingHelper.logException(logger, signature, e);
        } catch (ImageAlreadyExistException e) {
            throw LoggingHelper.logException(logger, signature, e);
        } finally {
            // remove the temp files
            for (File file : imageFiles) {
                if (file == null) {
                    continue;
                }
                file.delete();
            }
        }
    }

    @ResponseStatus(value=HttpStatus.CONFLICT, reason="Data violation")  // 409
    @ExceptionHandler(ImageAlreadyExistException.class)
    public void imageConflict() {
      // Nothing to do
      System.out.println("Image already exists.");
    }
    
    /**
     * <p>
     * This method processes the request to create a new detection session.
     * </p>
     * 
     * @return the detection session. (Spring Framework can be configured to
     *         serialize it to JSON)
     * 
     * @throws ServiceException
     *             if there are any error.
     */
    @RequestMapping(value = "detectionSessions", method = RequestMethod.POST, produces = "application/json")
    @ResponseBody
    public DetectionSession newDetectionSession() throws ServiceException {
        // prepare for logging
        final Logger logger = getLogger();
        final String signature = CLASS_NAME + ".newDetectionSession()";
        // log the entrance
        LoggingHelper.logEntrance(logger, signature, null, null);
        
        DetectionSession result;
        try {
            result = detectionSectionService.newDetectionSession();
        } catch (ServiceException e) {
            // log and re-throw
            throw LoggingHelper.logException(logger, signature, e);
        }
        
        // log the exit
        LoggingHelper.logExit(logger, signature, new Object[] {result});
        
        return result;
    }

    /**
     * <p>
     * This method processes the request to get detection session.
     * </p>
     * 
     * @param id
     *            the id of the session
     * @return the detection session. (Spring Framework can be configured to
     *         serialize it to JSON)
     * @throws ServiceException
     *             if there are any error.
     */
    @RequestMapping(value = "detectionSessions/{id}", method = RequestMethod.GET, produces = "application/json")
    @ResponseBody
    public DetectionSession getDetectionSession(@PathVariable long id) throws ServiceException {
        // prepare for logging
        final Logger logger = getLogger();
        final String signature = CLASS_NAME + ".getDetectionSession(long)";
        
        // log the entrance
        LoggingHelper.logEntrance(logger, signature, new String[] { "id" }, new Object[] { id });
        Helper.checkPositive(logger, signature, "id", id);
        
        DetectionSession result;
        try {
            result = detectionSectionService.getDetectionSession(id);
        } catch (ServiceException e) {
            throw LoggingHelper.logException(logger, signature, e);
        }
        
        // log the exit
        LoggingHelper.logExit(logger, signature, new Object[] {result});
        
        return result;
    }
}

