/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package gov.nasa.asteroid.hunter.web.controllers;

import gov.nasa.asteroid.hunter.ConfigurationException;
import gov.nasa.asteroid.hunter.Helper;
import gov.nasa.asteroid.hunter.LoggingHelper;
import gov.nasa.asteroid.hunter.models.DetectionItem;
import gov.nasa.asteroid.hunter.models.DetectionItemFrame;
import gov.nasa.asteroid.hunter.models.DetectionItemSearchCriteria;
import gov.nasa.asteroid.hunter.models.SearchResult;
import gov.nasa.asteroid.hunter.services.AsteroidDetectionService;
import gov.nasa.asteroid.hunter.services.EntityNotFoundException;
import gov.nasa.asteroid.hunter.services.ServiceException;

import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;
import java.net.URLEncoder;
import java.text.SimpleDateFormat;

import javax.annotation.PostConstruct;
import javax.servlet.http.HttpServletResponse;
import javax.validation.Valid;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVPrinter;
import org.apache.commons.io.IOUtils;
import org.apache.log4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.format.annotation.DateTimeFormat;
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
 * Spring MVC Controller for detection items.
 * </p>
 * <p>
 * <strong>Thread Safety:</strong> This class is effectively thread safe
 * (injected configurations are not considered as thread safety factor).
 * </p>
 * 
 * @author albertwang, TCSASSEMBLER
 * @version 1.0
 */
@Controller
public class DetectionItemController extends BaseController {
    
    /**
     * <p>
     * Represents the name of the class for logging.
     * </p>
     */
    private static final String CLASS_NAME = DetectionItemController.class.getName();

    /**
     * <p>
     * Represents the label for No.
     * </p>
     */
    private static final String NO_LABEL = "No";

    /**
     * <p>
     * Represents the label for Yes.
     * </p>
     */
    private static final String YES_LABEL = "Yes";

    /**
     * <p>
     * Represents the date format.
     * </p>
     */
    private static final String DATE_OUTPUT_FORMAT = "MM/dd/yyyy";

    /**
     * <p>
     * Represents the view detection item view name.
     * </p>
     */
    private static final String VIEW_DETECTION_ITEM_VIEW_NAME = "viewDetectionItem";

    /**
     * <p>
     * Represents the new comet reports email address attribute.
     * </p>
     */
    private static final String NEW_COMET_REPORTS_EMAIL_ADDRESS_ATTRIBUTE = "newCometReportsEmailAddress";

    /**
     * <p>
     * Represents the observation email address attribute.
     * </p>
     */
    private static final String OBSERVATIONS_EMAIL_ADDRESS_ATTRIBUTE = "observationsEmailAddress";

    /**
     * <p>
     * Represents the new comet report mail to attribute.
     * </p>
     */
    private static final String NEW_COMET_REPORT_MAIL_TO_ATTRIBUTE = "newCometReportMailTo";

    /**
     * <p>
     * Represents the observations mail to attribute.
     * </p>
     */
    private static final String OBSERVATIONS_MAIL_TO_ATTRIBUTE = "observationsMailTo";

    /**
     * <p>
     * Represents the detection item attribute.
     * </p>
     */
    private static final String DETECTION_ITEM_ATTRIBUTE = "detectionItem";

    /**
     * <p>
     * Represents list items view name.
     * </p>
     */
    private static final String LIST_ITEMS_VIEW_NAME = "listDetectionItems";
    
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
     * Represents the email address to which observations email will be sent.
     * </p>
     * <p>
     * Required, non-null/empty.
     * </p>
     */
    @Value("${controller.observationsEmailAddress}")
    private String observationsEmailAddress;

    /**
     * <p>
     * Represents the email address to which new comet reports email will be sent.
     * </p>
     * <p>
     * Required, non-null/empty.
     * </p>
     */
    @Value("${controller.newCometReportsEmailAddress}")
    private String newCometReportsEmailAddress;

    /**
     * <p>
     * Represents the observations email subject.
     * </p>
     * <p>
     * Required, non-null/empty.
     * </p>
     */
    @Value("${controller.observationsEmailSubject}")
    private String observationsEmailSubject;

    /**
     * <p>
     * Represents the new comet reports email subject.
     * </p>
     * <p>
     * Required, non-null/empty.
     * </p>
     */
    @Value("${controller.newCometReportsEmailSubject}")
    private String newCometReportsEmailSubject;

    /**
     * <p>
     * Creates the instance of this controller.
     * </p>
     */
    public DetectionItemController() {
        // does nothing
    }

    /**
     * <p>
     * Check if all required fields are initialized properly.
     * </p>
     * 
     * Note, in this controller, asteroidDetectionService should be injected,
     * newCometReportsEmailAddress should be not null or empty,
     * newCometReportsEmailSubject should be not null or empty,
     * observationsEmailSubject should be not null or empty.
     * 
     * @throws ConfigurationException
     *             if any required field is not initialized properly.
     */
    @PostConstruct
    protected void checkConfiguration() {
        // prepare for logging
        final Logger logger = getLogger();
        final String signature = CLASS_NAME + ".checkConfiguration()";
        // log the entrance
        LoggingHelper.logEntrance(logger, signature, null, null);
        super.checkConfiguration();

        Helper.checkConfigurationNull(logger, signature, "asteroidDetectionService", asteroidDetectionService);
        Helper.checkConfigurationNullOrEmpty(logger, signature, NEW_COMET_REPORTS_EMAIL_ADDRESS_ATTRIBUTE,
                newCometReportsEmailAddress);
        Helper.checkConfigurationNullOrEmpty(logger, signature, "newCometReportsEmailSubject",
                newCometReportsEmailSubject);
        Helper.checkConfigurationNullOrEmpty(logger, signature, 
                "observationsEmailSubject", observationsEmailSubject);
        Helper.checkConfigurationNullOrEmpty(logger, signature, 
                OBSERVATIONS_EMAIL_ADDRESS_ATTRIBUTE, observationsEmailAddress);
        
        // log the exit
        LoggingHelper.logExit(logger, signature, null);
    }

    /**
     * <p>
     * This method processes the request to list detection items.
     * </p>
     * 
     * @return the detection items view name.
     */
    @RequestMapping(value = "detectionItems", method = RequestMethod.GET)
    public String list() {
        // prepare for logging
        final Logger logger = getLogger();
        final String signature = CLASS_NAME + ".list()";

        LoggingHelper.logEntrance(logger, signature, null, null);

        String viewName = LIST_ITEMS_VIEW_NAME;

        // log the exit
        LoggingHelper.logExit(logger, signature, new Object[] { viewName });

        return viewName;
    }

    /**
     * <p>
     * This method processes the request to view detection item.
     * </p>
     * 
     * @param id the ID
     * @param model the Spring MVC model
     * 
     * @return the view name
     * 
     * @throws IllegalArgumentException if any argument is null or id is not positive
     * @throws ServiceException if there are any other exception.
     */
    @RequestMapping(value = "detectionItems/{id}", method = RequestMethod.GET)
    public String view(@PathVariable long id, Model model)  throws  ServiceException {
        
        // prepare for logging
        final Logger logger = getLogger();
        final String signature = CLASS_NAME + ".view(long, Model)";
        
        // log the entrance
        LoggingHelper.logEntrance(logger, signature, new String[] {"id", "model"}, new Object[] {id, model});
        
        // get the detected item
        DetectionItem item;
        try {
            item = asteroidDetectionService.getDetectionItem(id);
        } catch (ServiceException e) {
            // log and re-throw
            throw LoggingHelper.logException(logger, signature, e);
        }
        
        // set the detection item
        model.addAttribute(DETECTION_ITEM_ATTRIBUTE, item);
        
        String emailBody = buildDetectionItemBody(item);
        
        try {
            String mailToString = buildMailTo(observationsEmailAddress, observationsEmailSubject, emailBody);
            model.addAttribute(OBSERVATIONS_MAIL_TO_ATTRIBUTE, mailToString);
            model.addAttribute("observationsMailAddress", observationsEmailAddress);
            
            model.addAttribute("observationsMailSubject", URLEncoder.encode(observationsEmailSubject, "UTF-8"));
            model.addAttribute("observationsMailBody", URLEncoder.encode(emailBody, "UTF-8"));
            
            mailToString = buildMailTo(newCometReportsEmailAddress, newCometReportsEmailSubject, emailBody);
            model.addAttribute(NEW_COMET_REPORT_MAIL_TO_ATTRIBUTE, mailToString);
        } catch (ServiceException e) {
            // log and re-throw
            throw LoggingHelper.logException(logger, signature, e);
        } catch (UnsupportedEncodingException e) {
            // log and re-throw
            throw LoggingHelper.logException(logger, signature, new ServiceException("error", e));
        }
        
        model.addAttribute(OBSERVATIONS_EMAIL_ADDRESS_ATTRIBUTE, observationsEmailAddress);
        model.addAttribute(NEW_COMET_REPORTS_EMAIL_ADDRESS_ATTRIBUTE, newCometReportsEmailAddress);
        
        String viewName = VIEW_DETECTION_ITEM_VIEW_NAME;
        // log the exit
        LoggingHelper.logExit(logger, signature, new Object[] {viewName});
        
        return viewName;
    }

    /**
     * <p>
     * Builds the detection item body.
     * </p>
     * 
     * @param item
     *            the item.
     * @return the body string.
     */
    private String buildDetectionItemBody(DetectionItem item) {
        String[] names = {"DATE", "NEO?", "SUBMITTED?", "KNOWN BY MPC?", "RIGHT ASCENSION", "DECLINATION",
                "ROUGH MAGNITUDE", "OBSERVATORY CODE"};

        DetectionItemFrame frame = item.getFrames().get(0);

        String date = new SimpleDateFormat(DATE_OUTPUT_FORMAT).format(item.getTimestamp());
        String isNeo = item.isNeo() ? YES_LABEL : NO_LABEL;
        String isSubmitted = item.isSubmitted() ? YES_LABEL : NO_LABEL;
        String isKnownByMPC = item.isKnownByMPC() ? YES_LABEL : NO_LABEL;
        String rightAscension = double2String(frame.getRightAscension());
        String declination = double2String(frame.getDeclination());
        String roughMagnitude = double2String(frame.getRoughMagnitude());
        String observatoryCode = item.getObservatoryCode();

        String[] values = { date, isNeo, isSubmitted, isKnownByMPC, rightAscension, declination, roughMagnitude,
                observatoryCode };
        StringBuffer sb = new StringBuffer();
        for (int i = 0; i < names.length; i++) {
            String name = names[i];
            String value = values[i];
            sb.append(name).append(":").append(value).append('\n');
        }
        return sb.toString();
    }
    
    /**
     * <p>
     * Formats a double value to string.
     * </p>
     * @param v the value to format
     * @return the formatted string.
     */
    private String double2String(double v) {
        return String.format("%.2f", v);
    }
    
    /**
     * <p>
     * Builds the mail to string.
     * </p>
     * 
     * @param emailAddress the email address.
     * @param subject the email subject.
     * @param body the email body.
     * 
     * @return the mail to string.
     * 
     * @throws ServiceException if there are any error.
     */
    private String buildMailTo(String emailAddress, String subject, String body) 
        throws ServiceException {
        
        StringBuffer sb = new StringBuffer("mailto:").append(emailAddress).append("?");
        try {
            sb.append("subject=").append(URLEncoder.encode(subject, "UTF-8"));
            sb.append("&body=").append(URLEncoder.encode(body, "UTF-8"));
        } catch (UnsupportedEncodingException e) {
            throw new ServiceException("Failed to encode the subject or body.", e);
        }
        
        return sb.toString();
    }
    
    
    /**
     * <p>
     * This method processes the AJAX request to search detection items.
     * </p>
     * @param criteria the search criteria
     * @return the search result (Spring Framework can be configured to serialize it to JSON)
     * @throws IllegalArgumentException if the criteria is null.
     * @throws ServiceException if there are any other erros.
     */
    @RequestMapping(value = "search/detectionItems", method = RequestMethod.GET, produces = "application/json")
    @ResponseStatus(HttpStatus.OK)
    @ResponseBody
    public SearchResult<DetectionItem> search(
            @ModelAttribute("criteria") @DateTimeFormat(pattern = "MM/dd/yyyy hh:mm") 
            @Valid DetectionItemSearchCriteria criteria)
            throws ServiceException {
        // prepare for logging
        final Logger logger = getLogger();
        final String signature = CLASS_NAME + ".search(DetectionItemSearchCriteria)";
        LoggingHelper.logEntrance(logger, signature, new String[] { "criteria" }, new Object[] { criteria });
        // check the parameters
        Helper.checkNull(logger, signature, "criteria", criteria);
        try {
            SearchResult<DetectionItem> result = asteroidDetectionService.searchDetectionItems(criteria);

            // log the result and exit
            LoggingHelper.logExit(logger, signature, new Object[] { result });

            return result;
        } catch (ServiceException e) {
            // log and re-throw
            throw LoggingHelper.logException(logger, signature, e);
        }
    }

    /**
     * <p>
     * This method processes the AJAX request to submit a detection item.
     * </p>
     * 
     * @param id the id of the detection item.
     * 
     * @throws IllegalArgumentException : if id is not positive or userIds is empty/null
     * @throws ServiceException if there are any error.
     */
    @RequestMapping(value = "detectionItems/{id}/submit", method = RequestMethod.POST)
    @ResponseStatus(HttpStatus.OK)
    public void submitDetectionItem(@PathVariable long id) throws ServiceException {
        final Logger logger = getLogger();
        final String signature = CLASS_NAME + ".submitDetectionItem(long)";
        LoggingHelper.logEntrance(logger, signature, new String[] {"id"}, new Object[] {id});
        // validate the parameters
        Helper.checkPositive(logger, signature, "id", id);
        
        try {
            asteroidDetectionService.submitDetectionItem(id);
            // log and exit
            LoggingHelper.logExit(logger, signature, null);
            
        } catch (ServiceException e) {
            // log and re-throw the exception
            throw LoggingHelper.logException(logger, signature, e);
        }
        
    }
    
    /**
     * <p>
     * This method processes the AJAX request to mark a detection item as known by MPC.
     * </p>
     * 
     * @param id the id of the detection item.
     * 
     * @throws IllegalArgumentException if id is not positive or userIds is empty/null
     * @throws ServiceException if there are any error.
     */
    @RequestMapping(value = "detectionItems/{id}/known", method = RequestMethod.POST)
    @ResponseStatus(HttpStatus.OK)
    public void markDetectionItemKnownByMPC(@PathVariable long id) throws ServiceException {
        
        final Logger logger = getLogger();
        final String signature = CLASS_NAME + ".markDetectionItemKnownByMPC(long)";
        LoggingHelper.logEntrance(logger, signature, new String[] {"id"}, new Object[] {id});
        // validate the parameters
        Helper.checkPositive(logger, signature, "id", id);
        
        try {
            asteroidDetectionService.markDetectionItemKnownByMPC(id);
            // log and exit
            LoggingHelper.logExit(logger, signature, null);
            
        } catch (ServiceException e) {
            // log and re-throw the exception
            throw LoggingHelper.logException(logger, signature, e);
        }
        
    }
    
    /**
     * <p>
     * Gets the frame image.
     * </p>
     * 
     * @param id the id of the detection item.
     * @param frameIndex the frame index of the image.
     * 
     * @return the image content.
     * 
     * @throws IllegalArgumentException if the id is not positive or the frameIndex is not positive.
     * @throws ServiceException if there are any error.
     */
    @RequestMapping(value = "/frameimage/{id}/frame/{frame}/image", 
            method = RequestMethod.GET, produces = "image/png")
    @ResponseBody
    public byte[] getDetectionItemFrameImage(@PathVariable("id") long id, @PathVariable("frame") int frameIndex) 
        throws ServiceException {
        // prepare for logging
        final Logger logger = getLogger();
        final String signature = CLASS_NAME + ".getDetectionItemFrameImage(long, int)";
        // log the entrance
        LoggingHelper.logEntrance(logger, signature, 
                new String[] {"id", "frameIndex"}, new Object[] {id, frameIndex});
        
        // validate the parameters
        Helper.checkPositive(logger, signature, "id", id);
        if (frameIndex < 0) {
            throw LoggingHelper.logException(logger, signature, 
                    new IllegalArgumentException("frameIndex cannot be negative."));
        }
        
      
        // get the detect item by id
        DetectionItem item;
        try {
            item = asteroidDetectionService.getDetectionItem(id);
            if (item == null) {
                throw LoggingHelper.logException(logger, signature, 
                        new EntityNotFoundException("Detection item for id: " + id + " not found"));
            }
        } catch (ServiceException e) {
            // log and re-throw
            throw LoggingHelper.logException(logger, signature, e);
        }
        
        String imagePath = null;
        DetectionItemFrame theFrame = null;
        for (DetectionItemFrame frame : item.getFrames()) {
            if (frame.getDetectionItemFramePk().getFrame() == frameIndex) {
                // get the frame
                theFrame = frame;
            }
            if (frame.getDetectionItemFramePk().getFrame() == 0) {
                // get the frame (always use the first frame as background image)
                imagePath = frame.getVisualizationImage();
            }
        }
        if (imagePath == null || theFrame == null) {
            throw LoggingHelper.logException(logger, signature, 
                    new EntityNotFoundException("Image not found for id: " + id + " and frame: " + frameIndex));
        }
        
        InputStream is = null;
        try {
            // Map<String, Object> params = new HashMap<String, Object>();

            // set optional parameters if you like
            /*params.put(ImagingConstants.BUFFERED_IMAGE_FACTORY,
                    new ManagedImageBufferedImageFactory());*/

            // params.put(ImagingConstants.PARAM_KEY_VERBOSE, Boolean.TRUE);

            // read image
            // BufferedImage nbi = Imaging.getBufferedImage(new File(imagePath));
            // Image image = Toolkit.getDefaultToolkit().getImage(imagePath);
            /*BufferedImage nbi = ImageIO.read(new File(imagePath));
            // draw the circle
            int circleSize = 160;
            Graphics2D ng = (Graphics2D)nbi.getGraphics();
            ng.setStroke(new BasicStroke(circleSize / 10));
            ng.setColor(Color.RED);
            ng.drawOval((int)(theFrame.getX()-circleSize/2), (int)(theFrame.getY()-circleSize/2), circleSize, circleSize);
            ByteArrayOutputStream baos = new ByteArrayOutputStream();
            ImageIO.write(nbi, "jpeg", baos);
            ng.dispose();
            byte[] bytes = baos.toByteArray();*/
            is = new FileInputStream(new File(imagePath));
            byte[] bytes = IOUtils.toByteArray(is);
            LoggingHelper.logExit(logger, signature, null);
            return bytes;
        } catch (IOException e) {
            throw LoggingHelper.logException(logger, signature, 
                    new ServiceException("Error occur while reading the image.", e));
        } /*catch (ImageReadException e) {
            throw LoggingHelper.logException(logger, signature, 
                    new ServiceException("Error occur while reading the image.", e));
        }*/ finally {
            if (is != null) {
                try {
                    is.close();
                } catch (IOException e) {
                    // ignore
                }
            }
        }
    }
    
    /**
     * Converts a given Image into a BufferedImage
     *
     * @param img The Image to be converted
     * @return The converted BufferedImage
     */
    public static BufferedImage toBufferedImage(Image img)
    {
        if (img instanceof BufferedImage)
        {
            return (BufferedImage) img;
        }

        // Create a buffered image with transparency
        BufferedImage bimage = new BufferedImage(img.getWidth(null), img.getHeight(null), BufferedImage.TYPE_INT_ARGB);

        // Draw the image on to the buffered image
        Graphics2D bGr = bimage.createGraphics();
        bGr.drawImage(img, 0, 0, null);
        bGr.dispose();

        // Return the buffered image
        return bimage;
    }
    
    @RequestMapping(value = "/detectionItems/download", 
            method = RequestMethod.GET, produces = "text/csv")
    public void generateDetectionItemsCSV(HttpServletResponse response) throws ServiceException {
        // prepare for logging
        final Logger logger = getLogger();
        final String signature = CLASS_NAME + ".generateDetectionItemsCSV(HttpServletResponse)";
        
        response.setHeader("Content-Disposition", "attachment;filename=download.csv");
        response.setHeader("Content-Type", "text/csv");
        
        // get all items
        SearchResult<DetectionItem> result = asteroidDetectionService.searchDetectionItems(new DetectionItemSearchCriteria());
        
        // print the CSV header
        CSVPrinter csvPrinter = null;
        try {
            csvPrinter = new CSVPrinter(new OutputStreamWriter(response.getOutputStream()), CSVFormat.DEFAULT);
            csvPrinter.printRecord("DATE", "TIME", "RIGHT ASCENSION", "DECLINATION", "OBSERVATORY CODE", "SIGNIFICANCE", "EXISTS AT MPC?");
            for (DetectionItem item : result.getValues()) {
                csvPrinter.print(new SimpleDateFormat("MM/dd/yyyy").format(item.getTimestamp()));
                csvPrinter.print(new SimpleDateFormat("hh:mm a").format(item.getTimestamp()));
                csvPrinter.print(item.getFrames().get(0).getRightAscension());
                csvPrinter.print(item.getFrames().get(0).getDeclination());
                csvPrinter.print(item.getObservatoryCode());
                csvPrinter.print(item.getSignificance());
                csvPrinter.print(item.isSubmitted());
                csvPrinter.println();
            }
        } catch (IOException e) {
            throw LoggingHelper.logException(logger, signature, 
                    new ServiceException("Error occur while exporting the CSV.", e));
        } finally {
            if (csvPrinter != null) {
                try {
                    csvPrinter.close();
                } catch (IOException e) {
                    // ignore
                }
            }
        }
        
    }
    
    @RequestMapping(value = "/detectionItems/deleteAll", 
            method = RequestMethod.DELETE)
    @ResponseStatus(HttpStatus.OK)
    public void deleteAllDetectionItems() throws ServiceException {
        asteroidDetectionService.deleteAllDetectionItems();
    }
    
}

