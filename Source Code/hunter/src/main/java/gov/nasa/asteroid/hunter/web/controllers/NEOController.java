/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package gov.nasa.asteroid.hunter.web.controllers;

import gov.nasa.asteroid.hunter.ConfigurationException;
import gov.nasa.asteroid.hunter.Helper;
import gov.nasa.asteroid.hunter.LoggingHelper;
import gov.nasa.asteroid.hunter.models.ImageInfo;
import gov.nasa.asteroid.hunter.models.NEO;
import gov.nasa.asteroid.hunter.models.NEOSearchCriteria;
import gov.nasa.asteroid.hunter.models.SearchResult;
import gov.nasa.asteroid.hunter.services.NEOService;
import gov.nasa.asteroid.hunter.services.ServiceException;
import gov.nasa.asteroid.tester.ImageService;

import javax.annotation.PostConstruct;
import javax.validation.Valid;

import org.apache.log4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.format.annotation.DateTimeFormat;
import org.springframework.http.HttpStatus;
import org.springframework.stereotype.Controller;
import org.springframework.ui.Model;
import org.springframework.web.bind.annotation.ModelAttribute;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RequestMethod;
import org.springframework.web.bind.annotation.ResponseBody;
import org.springframework.web.bind.annotation.ResponseStatus;

/**
 * <p>
 * Spring MVC Controller for NEOs.
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
@Controller
public class NEOController extends BaseController {
    
    /**
     * <p>
     * Represents the name of this class for logging.
     * </p>
     */
    private final static String CLASS_NAME = NEOController.class.getName();
    
    /**
     * <p>
     * Represents the NEOService used for searching NEOs.
     * </p>
     * <p>
     * It is required.
     * </p>
     */
    @Autowired
    private NEOService neoService;

    @Autowired
    private ImageService imageService;
    
    /**
     * <p>
     * The default constructor of this class.
     * </p>
     */
    public NEOController() {
        // does nothing
    }

    /**
     * <p>
     * Check if all required fields are initialized properly.
     * </p>
     * 
     * @throws ConfigurationException if any required field is not initialized properly.
     */
    @PostConstruct
    protected void checkConfiguration()  {
        final Logger logger = getLogger();
        final String signature = CLASS_NAME + ".checkConfiguration()";
        // log the entrance
        LoggingHelper.logEntrance(logger, signature, null, null);
        super.checkConfiguration();
        
        Helper.checkConfigurationNull(logger, signature, "neoService", neoService);
        // log the exit
        LoggingHelper.logExit(logger, signature, null);
    }

    /**
     * <p>
     * This method processes the request to list NEOs.
     * </p>
     * @return the view name
     */
    @RequestMapping(value = "neos", method = RequestMethod.GET)
    public String list(Model model) {
        final Logger logger = getLogger();
        final String signature = CLASS_NAME + ".list()";
        // log the entrance
        LoggingHelper.logEntrance(logger, signature, null, null);

        String viewName = "listNEOs";
        
        ImageInfo imageInfo = imageService.getLatestImage(ImageInfo.class);
        if (imageInfo != null) {
            model.addAttribute("lastImageInfo", imageInfo);
        }
        // log and return
        LoggingHelper.logExit(logger, signature, new Object[] { viewName });

        return viewName;
    }

    /**
     * <p>
     * This method processes the AJAX request to search NEOs.
     * </p>
     * 
     * @param criteria
     *            the search criteria.
     * 
     * @return the search result (Spring Framework can be configured to
     *         serialize it to JSON)
     * 
     * @throws IllegalArgumentException
     *             if any argument is null.
     * @throws ServiceException
     *             if there are any error.
     * 
     */
    @RequestMapping(value = "search/neos", method = RequestMethod.GET, produces = "application/json")
    @ResponseStatus(HttpStatus.OK)
    @ResponseBody
    public SearchResult<NEO> search(
            @Valid @ModelAttribute("criteria") @DateTimeFormat(pattern = "MM/dd/yyyy") NEOSearchCriteria criteria)
            throws ServiceException {
        final Logger logger = getLogger();
        final String signature = CLASS_NAME + ".search(NEOSearchCriteria)";
        LoggingHelper.logEntrance(logger, signature, new String[] { "criteria" }, new Object[] { criteria });
        // validate the parameters
        Helper.checkNull(logger, signature, "criteria", criteria);
        SearchResult<NEO> result;
        try {
            result = neoService.search(criteria);
        } catch (ServiceException e) {
            // log and re-throw
            throw LoggingHelper.logException(logger, signature, e);
        }
        LoggingHelper.logExit(logger, signature, new Object[] { result });
        return result;
    }
}

