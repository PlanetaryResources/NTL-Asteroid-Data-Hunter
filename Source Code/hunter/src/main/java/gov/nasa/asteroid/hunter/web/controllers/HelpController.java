/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package gov.nasa.asteroid.hunter.web.controllers;

import gov.nasa.asteroid.hunter.ConfigurationException;
import gov.nasa.asteroid.hunter.Helper;
import gov.nasa.asteroid.hunter.LoggingHelper;
import gov.nasa.asteroid.hunter.models.HelpItem;
import gov.nasa.asteroid.hunter.models.HelpItemSearchCriteria;
import gov.nasa.asteroid.hunter.models.SearchResult;
import gov.nasa.asteroid.hunter.services.HelpService;
import gov.nasa.asteroid.hunter.services.ServiceException;

import javax.annotation.PostConstruct;
import javax.validation.Valid;

import org.apache.log4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.HttpStatus;
import org.springframework.stereotype.Controller;
import org.springframework.ui.Model;
import org.springframework.web.bind.annotation.ModelAttribute;
import org.springframework.web.bind.annotation.PathVariable;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RequestMethod;
import org.springframework.web.bind.annotation.ResponseBody;
import org.springframework.web.bind.annotation.ResponseStatus;

/**
 * <p>
 * Spring MVC Controller for help contents.
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
public class HelpController extends BaseController {
    
    /**
     * <p>
     * Represents the name of this class for logging.
     * </p>
     */
    private static final String CLASS_NAME = HelpController.class.getName();
    
    /**
     * <p>
     * Represents the HelpService used for accessing help content.
     * </p>
     * <p>
     * It is required.
     * </p>
     */
    @Autowired
    private HelpService helpService;

    /**
     * <p>
     * The default constructor for HelpController.
     * </p>
     */
    public HelpController() {
        // does nothing
    }

    /**
     * <p>
     * Check if all required fields are initialized properly.
     * </p>
     * Note, in this class, helpService must be injected.
     * @throws ConfigurationException if any required field is not initialized properly.
     */
    @PostConstruct
    protected void checkConfiguration() {
        // prepare for logging
        final Logger logger = getLogger();
        final String signature = CLASS_NAME + ".checkConfiguration()";
        // log the entrance
        LoggingHelper.logEntrance(logger, signature, null, null);
        
        super.checkConfiguration();
        Helper.checkConfigurationNull(logger, signature, "helpService", helpService);
        
        // log the exit
        LoggingHelper.logExit(logger, signature, null);
    }

    /**
     * <p>
     * This method processes the request to list help items.
     * </p>
     * 
     * @param model
     *            the Spring MVC model
     * 
     * @return the view name
     * 
     * @throws ServiceException
     *             if there are any errors.
     */
    @RequestMapping(value = "helpItems", method = RequestMethod.GET)
    public String list(Model model) throws ServiceException {
        final Logger logger = getLogger();
        final String signature = CLASS_NAME + ".list(Model)";
        LoggingHelper.logEntrance(logger, signature, new String[] { "model" }, new Object[] { model });

        model.addAttribute("helpTopics", helpService.getHelpTopics());
        String viewName = "listHelpItems";
        LoggingHelper.logExit(logger, signature, new Object[] { viewName });
        return viewName;
    }

    /**
     * <p>
     * This method processes the request to view help item.
     * </p>
     * 
     * @param id
     *            the id of the help item
     * @param model
     *            the model
     * 
     * @return the view name
     * 
     * @throws IllegalArgumentException
     *             if any argument is null or id is not positive
     * @throws ServiceException
     *             if there are any error.
     * 
     */
    @RequestMapping(value = "helpItems/{id}", method = RequestMethod.GET)
    public String view(@PathVariable long id, Model model) throws ServiceException {
        final Logger logger = getLogger();
        final String signature = CLASS_NAME + ".view(long, Model)";
        LoggingHelper.logEntrance(logger, signature, new String[] { "id", "model" }, new Object[] { id, model });
        // validate the parameters
        Helper.checkPositive(logger, signature, "id", id);

        try {
            model.addAttribute("helpItem", helpService.getHelpItem(id));
            model.addAttribute("helpTopics", helpService.getHelpTopics());
        } catch (ServiceException e) {
            // log and re-throw
            throw LoggingHelper.logException(logger, signature, e);
        }

        String viewName = "viewHelpItem";
        // log the exit
        LoggingHelper.logExit(logger, signature, new Object[] { viewName });

        return viewName;
    }

    /**
     * <p>
     * This method processes the AJAX request to search help items.
     * </p>
     * 
     * @param criteria
     *            the search criteria
     * 
     * @return the search result (Spring Framework can be configured to
     *         serialize it to JSON)
     * 
     * @throws IllegalArgumentException
     *             if any argument is null.
     * @throws ServiceException
     *             if there are any other errors.
     */
    @RequestMapping(value = "search/helpItems", method = RequestMethod.GET, produces = "application/json")
    @ResponseStatus(HttpStatus.OK)
    @ResponseBody
    public SearchResult<HelpItem> search(@Valid @ModelAttribute HelpItemSearchCriteria criteria) throws ServiceException {
        final Logger logger = getLogger();
        final String signature = CLASS_NAME + ".search(HelpItemSearchCriteria)";
        LoggingHelper.logEntrance(logger, signature, new String[] { "criteria" }, new Object[] { criteria });

        try {
            SearchResult<HelpItem> result = helpService.searchHelpItems(criteria);
            LoggingHelper.logExit(logger, signature, new Object[] { result });
            return result;
        } catch (ServiceException e) {
            // log and re-throw
            throw LoggingHelper.logException(logger, signature, e);
        }
    }
}

