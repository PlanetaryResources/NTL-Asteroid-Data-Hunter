/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package gov.nasa.asteroid.hunter.services.impl;

import gov.nasa.asteroid.hunter.ConfigurationException;
import gov.nasa.asteroid.hunter.Helper;
import gov.nasa.asteroid.hunter.LoggingHelper;
import gov.nasa.asteroid.hunter.models.NEO;
import gov.nasa.asteroid.hunter.models.NEOSearchCriteria;
import gov.nasa.asteroid.hunter.models.SearchResult;
import gov.nasa.asteroid.hunter.models.SortType;
import gov.nasa.asteroid.hunter.services.NEOService;
import gov.nasa.asteroid.hunter.services.ServiceException;

import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.UnsupportedEncodingException;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLEncoder;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import javax.annotation.PostConstruct;

import org.apache.log4j.Logger;

/**
 * <p>
 * This class is the implementation of NEOService.
 * </p>
 *
 * <p>
 * This service provides a method to search known NEOs.
 * </p>
 * <p>
 * <strong>Thread Safety:</strong> This class is effectively thread safe
 * (injected configurations are not considered as thread safety factor).
 * </p>
 *
 * @author albertwang, TCSASSEMBLER
 * @version 1.0
 */
public class NEOServiceImpl extends BaseService implements NEOService {


    /**
     * <p>
     * Represents the name of the class for logging.
     * </p>
     */
    private static final String CLASS_NAME = NEOServiceImpl.class.getName();

    /**
     * <p>
     * Represents the empty string.
     * </p>
     */
    private static final String EMPTY_STRING = "";

    /**
     * <p>
     * Represents the seconds in a minute.
     * </p>
     */
    private static final double SECONDS_PER_MINUTE = 60;

    /**
     * <p>
     * Represents the seconds in an hour.
     * </p>
     */
    private static final double SECONDS_PER_HOUR = 60 * 60;

    /**
     * <p>
     * Represents the seconds in a day.
     * </p>
     */
    private static final double SECONDS_PER_DAY = 60 * 60 * 24;


    /**
     * <p>
     * Represents the parameters for others.
     * </p>
     */
    private static final String PARAM_OTHERS = "&sort=d&mot=h&tmot=s&pdes=u&needed=f&ps=n&type=p";

    /**
     * <p>
     * Represents the parameter for oc.
     * </p>
     */
    private static final String PARAM_OC = "&oc=";

    /**
     * <p>
     * Represents the parameter for limit.
     * </p>
     */
    private static final String PARAM_LIMIT = "&limit=";

    /**
     * <p>
     * Represents the parameter for radius.
     * </p>
     */
    private static final String PARAM_RADIUS = "&radius=";

    /**
     * <p>
     * Represents the parameter for TextArea.
     * </p>
     */
    private static final String PARAM_TEXT_AREA = "&TextArea=";

    /**
     * <p>
     * Represents the parameter for decl.
     * </p>
     */
    private static final String PARAM_DECL = "&decl=";

    /**
     * <p>
     * Represents the parameter for ra.
     * </p>
     */
    private static final String PARAM_RA = "&ra=";


    /**
     * <p>
     * Represents the which parameters.
     * </p>
     */
    private static final String PARAM_WHICH = "&which=";

    /**
     * <p>
     * Represents the pos for which.
     * </p>
     */
    private static final String WHICH_VALUE_POS = "pos";

    /**
     * <p>
     * Represents the obs for which.
     * </p>
     */
    private static final String WHICH_VALUE_OBS = "obs";

    /**
     * <p>
     * Represents the parameter for day.
     * </p>
     */
    private static final String PARAM_DAY = "&day=";

    /**
     * <p>
     * Represents the parameter for month.
     * </p>
     */
    private static final String PARAM_MONTH = "&month=";

    /**
     * <p>
     * Represents the parameter for year.
     * </p>
     */
    private static final String PARAM_YEAR = "year=";

    /**
     * <p>
     * Represents the degree character for second.
     * </p>
     */
    private static final String DEGREE_SECOND_CHAR = "\"";

    /**
     * <p>
     * Represents the degree character for minute.
     * </p>
     */
    private static final String DEGREE_MINUTE_CHART = "'";

    /**
     * <p>
     * Represents the degree character.
     * </p>
     */
    private static final String DEGREE_CHAR = "&#176;";

    /**
     * <p>
     * Represents letter for second.
     * </p>
     */
    private static final String SECOND_LETTER = "s";

    /**
     * <p>
     * Represents letter for minute.
     * </p>
     */
    private static final String MINUTE_LETTER = "m";

    /**
     * <p>
     * Represents letter for hour.
     * </p>
     */
    private static final String HOUR_LETTER = "h";

    /**
     * <p>
     * Represents the default value of V if V is invalid.
     * </p>
     */
    private static final int EMPTY_DEFAULT_V = 99;

    /**
     * <p>
     * Represents the end index of the orbit column.
     * </p>
     */
    private static final int END_INDEX_ORBIT = 87;

    /**
     * <p>
     * Represents the start index of the orbit column.
     * </p>
     */
    private static final int START_INDEX_ORBIT = 80;

    /**
     * <p>
     * Represents the end index of the montion decl column.
     * </p>
     */
    private static final int START_INDEX_MOTION_DECL = 74;

    /**
     * <p>
     * Represents the end index of the motion ra column.
     * </p>
     */
    private static final int START_INDEX_MOTION_RA = 67;

    /**
     * <p>
     * Represents the end index of the decl column.
     * </p>
     */
    private static final int START_INDEX_OFFSET_DECL = 59;

    /**
     * <p>
     * Represents the end index of the ra column.
     * </p>
     */
    private static final int START_INDEX_OFFSET_RA = 52;

    /**
     * <p>
     * Represents the end index of the v column.
     * </p>
     */
    private static final int START_INDEX_V = 47;

    /**
     * <p>
     * Represents the end index of the decl column.
     * </p>
     */
    private static final int START_INDEX_DECL = 36;

    /**
     * <p>
     * Represents the end index of the ra column.
     * </p>
     */
    private static final int START_INDEX_RA = 25;

    /**
     * <p>
     * Represents the end index of the object designation column.
     * </p>
     */
    private static final int START_INDEX_OBJECT_DESIGNATION = 0;

    /**
     * <p>
     * Represents the end of the pre html tag.
     * </p>
     */
    private static final String PRE_TAG_END = "</pre>";

    /**
     * <p>
     * Represents the start of the pre html tag.
     * </p>
     */
    private static final String PRE_TAG_START = "<pre>";

    /**
     * <p>
     * Represents the US lang.
     * </p>
     */
    private static final String EN_US_LANG = "en-US";

    /**
     * <p>
     * Represents the UTF-8 encoding.
     * </p>
     */
    private static final String UTF_8_ENCODING = "UTF-8";

    /**
     * <p>
     * Represents the content type.
     * </p>
     */
    private static final String CONTENT_TYPE = "application/x-www-form-urlencoded";

    /**
     * <p>
     * Represents content language
     * </p>
     */
    private static final String HEADER_CONTENT_LANGUAGE = "Content-Language";

    /**
     * <p>
     * Represents the content length.
     * </p>
     */
    private static final String HEADER_CONTENT_LENGTH = "Content-Length";

    /**
     * <p>
     * Represents the content type.
     * </p>
     */
    private static final String HEADER_CONTENT_TYPE = "Content-Type";

    /**
     * <p>
     * Represents the post method
     * </p>
     */
    private static final String METHOD_POST = "POST";

    /**
     * <p>
     * Represents the minimum value of v.
     * </p>
     */
    private static final int MIN_V = 1;

    /**
     * <p>
     * Represents the max value of the radius.
     * </p>
     */
    private static final int MAX_RADIUS = 300;

    /**
     * <p>
     * Represents the min value of the radius.
     * </p>
     */
    private static final int MIN_RADIUS = 5;

    /**
     * <p>
     * Represents the sort by column orbit.
     * </p>
     */
    private static final String SORT_BY_ORBIT = "orbit";

    /**
     * <p>
     * Represents the sort by column motion declination.
     * </p>
     */
    private static final String SORT_BY_MOTION_DECLINATION = "motionDeclination";

    /**
     * <p>
     * Represents the sort by column montion right ascension.
     * </p>
     */
    private static final String SORT_BY_MOTION_RIGHT_ASCENSION = "motionRightAscension";

    /**
     * <p>
     * Represents the sort by column offsest declination.
     * </p>
     */
    private static final String SORT_BY_OFFSET_DECLINATION = "offsetDeclination";

    /**
     * <p>
     * Represents the sort by column right asension.
     * </p>
     */
    private static final String SORT_BY_OFFSET_RIGHT_ASCENSION = "offsetRightAscension";

    /**
     * <p>
     * Represents the sort by column v.
     * </p>
     */
    private static final String SORT_BY_V = "v";

    /**
     * <p>
     * Represents the sort by column declination.
     * </p>
     */
    private static final String SORT_BY_DECLINATION = "declination";

    /**
     * <p>
     * Represents the sort by ascension.
     * </p>
     */
    private static final String SORT_BY_RIGHT_ASCENSION = "rightAscension";

    /**
     * <p>
     * Represents the sort by column object designation.
     * </p>
     */
    private static final String SORT_BY_OBJECT_DESIGNATION = "objectDesignation";


    /**
     * <p>
     * Represents the MP Checker CGI URL.
     * </p>
     * <p>
     * Optional. Default to "http://scully.cfa.harvard.edu/cgi-bin/mpcheck.cgi".
     * </p>
     * <p>
     * Not null/empty.
     * </p>
     */
    private String mpCheckerCGIURL = "http://scully.cfa.harvard.edu/cgi-bin/mpcheck.cgi";

    /**
     * <p>
     * Represents the max entries in the cached.
     * </p>
     */
    private static final int MAX_ENTRIES = 200;
    
    /**
     * <p>
     * Represents the cache.
     * </p>
     */
    private static Map <String, List<NEO>> cache = new LinkedHashMap<String, List<NEO>> (MAX_ENTRIES + 1, .75F, true) {
        public boolean removeEldestEntry(Map.Entry<String, List<NEO>> eldest) {
            return size() > MAX_ENTRIES;
        }
    };
    
    /**
     * <p>
     * This is the default class for <code>NEOServiceImpl</code>.
     * </p>
     */
    public NEOServiceImpl() {
        // does nothing
    }

    /**
     * <p>
     * Check if all required fields are initialized properly.
     * </p>
     * <p>
     * Note, in this class, the field is mpCheckerCGIURL. It is optional, but it
     * required to be not null or empty.
     *
     * </p>
     *
     * @throws ConfigurationException  if any required field is not initialized properly.
     */
    @PostConstruct
    protected void checkConfiguration() {
        // prepare for logging
        final String signature = CLASS_NAME + ".checkConfiguration()";
        Logger logger = getLogger();

        // log the entrance
        LoggingHelper.logEntrance(logger, signature, null, null);

        // check the injection in the base class
        super.checkConfiguration();

        // validate the parameters
        Helper.checkConfigurationNullOrEmpty(logger, signature, "mpCheckerCGIURL", mpCheckerCGIURL);

        // log the exit
        LoggingHelper.logExit(logger, signature, null);
    }

    /**
     * <p>
     * Sets the MP checker CGUI URL.
     * </p>
     * @param mpCheckerCGIURL the MP Checker CGI URL.
     */
    public void setMpCheckerCGIURL(String mpCheckerCGIURL) {
        this.mpCheckerCGIURL = mpCheckerCGIURL;
    }

    /**
     * <p>
     * Searches the known NEOs.
     * </p>
     * @param criteria the search criteria
     *
     * @return the search result
     *
     * @throws IllegalArgumentException if criteria is null or invalid
     * (see <code>NEOSearchCriteria</code> document for detail)
     * @throws ServiceException if any other error occurred during the operation
     *
     */
    @Override
    public SearchResult<NEO> search(NEOSearchCriteria criteria) throws ServiceException {
        // prepare for logging
        final String signature = CLASS_NAME + ".search(NEOSearchCriteria)";
        Logger logger = getLogger();

        // log the entrance
        LoggingHelper.logEntrance(logger, signature, new String[] { "criteria" }, new Object[] { criteria });

        // validate the parameter
        Helper.checkBaseSearchParameters(logger, signature, "criteria", criteria, Arrays.asList(
                SORT_BY_OBJECT_DESIGNATION, SORT_BY_RIGHT_ASCENSION, SORT_BY_DECLINATION, SORT_BY_V,
                SORT_BY_OFFSET_RIGHT_ASCENSION, SORT_BY_OFFSET_DECLINATION, SORT_BY_MOTION_RIGHT_ASCENSION,
                SORT_BY_MOTION_DECLINATION, SORT_BY_ORBIT));

        // validate the required fields
        Helper.checkNull(logger, signature, "criteria.date", criteria.getDate());
        if (criteria.getObservations() == null) {
            Helper.checkNull(logger, signature, "criteria.rightAscension", criteria.getRightAscension());
            Helper.checkNull(logger, signature, "criteria.declination", criteria.getDeclination());
        }
        if (criteria.getRadius() < MIN_RADIUS || criteria.getRadius() > MAX_RADIUS) {
            throw LoggingHelper.logException(logger, signature,
                    new IllegalArgumentException("criteria.radius should in range [5,300]"));
        }
        if (criteria.getV() < MIN_V) {
            throw LoggingHelper.logException(logger, signature,
                    new IllegalArgumentException("criteria.v should >= 1."));
        }
        Helper.checkNullOrEmpty(logger, signature,
                "criteria.observatoryCode", criteria.getObservatoryCode());


        // Prepare POST parameters
        String requestParameters = buildRequestParameters(logger, signature, criteria);

        // Call MPChecker CGI service
        List<NEO> neos = null;
        synchronized(cache) {
            // check if in the cache
            if (cache.containsKey(requestParameters)) {
                neos = cache.get(requestParameters);
            } 
        }
        
        if (neos == null) {
            // no cache
            neos = requestNEOs(logger, signature, requestParameters);
            synchronized(cache) {
                cache.put(requestParameters, neos);
            }
        }
        

        // Sort results
        if (criteria.getSortBy() != null) {
            sortNEOs(neos, criteria);
        }

        // Construct SearchResult
        SearchResult<NEO> result = new SearchResult<NEO>();
        result.setPageSize(criteria.getPageSize());
        result.setPageNumber(criteria.getPageNumber());
        result.setSortBy(criteria.getSortBy());
        result.setSortType(criteria.getSortType());
        result.setTotal(neos.size());

        if (criteria.getPageNumber() > 0) {
            int totalPageCount = (neos.size() + criteria.getPageSize() - 1) / criteria.getPageSize();
            result.setTotalPages(totalPageCount);

            int startIndex = criteria.getPageSize() * (criteria.getPageNumber() - 1);
            
            int endIndex = criteria.getPageSize() * criteria.getPageNumber();
            if (endIndex > neos.size()) {
                endIndex = neos.size();
            }
            if (startIndex >= neos.size()) {
                result.setValues(new ArrayList<NEO>());
            } else {
                result.setValues(neos.subList(startIndex, endIndex));
            }
        } else {
            result.setTotalPages(neos.size() > 0 ? 1 : 0);
            result.setValues(neos);
        }

        // log the exit
        LoggingHelper.logExit(logger, signature, new Object[] { result });

        return result;
    }

    /**
     * <p>
     * Sorts the neos array.
     * </p>
     * @param neos the neos array for sorting.
     * @param criteria the criteria to sort.
     */
    private void sortNEOs(List<NEO> neos, NEOSearchCriteria criteria) {
        final String sortBy = criteria.getSortBy();
        final boolean asc = criteria.getSortType() != null
                ? (criteria.getSortType() == SortType.ASC) : true;
        Collections.sort(neos, new Comparator<NEO>() {
            /**
             * <p>
             * Compares the 2 NEO objects.
             * </p>
             * @param o1 the first object
             * @param o2 the second object
             * @return < 0 if they are in ascending order,
             * or 0 if they are equal or > 0 if they are in descending order.
             */
            public int compare(NEO o1, NEO o2) {
                int result = 0;
                if (SORT_BY_OBJECT_DESIGNATION.equals(sortBy)) {
                    result = o1.getObjectDesignation().compareTo(o2.getObjectDesignation());
                } else if (SORT_BY_RIGHT_ASCENSION.equals(sortBy)) {
                    result = o1.getRightAscension().compareTo(o2.getRightAscension());
                } else if (SORT_BY_DECLINATION.equals(sortBy)) {
                    result = o1.getDeclination().compareTo(o2.getDeclination());
                } else if (SORT_BY_V.equals(sortBy)) {
                    result = Double.compare(o1.getV(), o2.getV());
                } else if (SORT_BY_OFFSET_RIGHT_ASCENSION.equals(sortBy)) {
                    String offsetRightAscension1 = o1.getOffsetRightAscension();
                    String offsetRightAscension2 = o2.getOffsetRightAscension();
                    if ((offsetRightAscension1.endsWith("E") || offsetRightAscension1.endsWith("W"))
                            && (offsetRightAscension2.endsWith("E") || offsetRightAscension2.endsWith("W"))) {
                        char l1 = offsetRightAscension1.charAt(offsetRightAscension1.length() - 1);
                        char l2 = offsetRightAscension2.charAt(offsetRightAscension2.length() - 1);
                        if (l1 == l2) {
                            double v1 = new Double(offsetRightAscension1.substring(0, offsetRightAscension1.length() - 1));
                            double v2 = new Double(offsetRightAscension2.substring(0, offsetRightAscension2.length() - 1));
                            result = Double.compare(v1, v2);
                        } else if (l1 == 'E') {
                            result = -1;
                        } else {
                            result = 1;
                        }
                    } else {
                        result = o1.getOffsetRightAscension().compareTo(o2.getOffsetRightAscension());
                    }
                    
                } else if (SORT_BY_OFFSET_DECLINATION.equals(sortBy)) {
                    String offsetDeclination1 = o1.getOffsetDeclination();
                    String offsetDeclination2 = o2.getOffsetDeclination();
                    if ((offsetDeclination1.endsWith("S") || offsetDeclination1.endsWith("N"))
                            && (offsetDeclination2.endsWith("S") || offsetDeclination2.endsWith("N"))) {
                        char l1 = offsetDeclination1.charAt(offsetDeclination1.length() - 1);
                        char l2 = offsetDeclination2.charAt(offsetDeclination2.length() - 1);
                        if (l1 == l2) {
                            double v1 = new Double(offsetDeclination1.substring(0, offsetDeclination1.length() - 1));
                            double v2 = new Double(offsetDeclination2.substring(0, offsetDeclination2.length() - 1));
                            result = Double.compare(v1, v2);
                        } else if (l1 == 'S') {
                            result = -1;
                        } else {
                            result = 1;
                        }
                    } else {
                        result = o1.getOffsetDeclination().compareTo(o2.getOffsetDeclination());
                    }
                } else if (SORT_BY_MOTION_RIGHT_ASCENSION.equals(sortBy)) {
                    result = Double.compare(getMotionValue(o1.getMotionRightAscension()), getMotionValue(o2.getMotionRightAscension()));
                    
                } else if (SORT_BY_MOTION_DECLINATION.equals(sortBy)) {
                    result = Double.compare(getMotionValue(o1.getMotionDeclination()), getMotionValue(o2.getMotionDeclination()));
                } else if (SORT_BY_ORBIT.equals(sortBy)) {
                    result = o1.getOrbit().compareTo(o2.getOrbit());
                }
                return asc ? result : -result;
            }
        });
    }

    /**
     * <p>
     * Gets the motion value.
     * </p>
     * @return the motion value.
     */
    private double getMotionValue(String str) {
        if (str.endsWith("+") || str.endsWith("-")) {
            return new Double(str.charAt(str.length() - 1) + str.substring(0, str.length() - 1));
        } else {
            return new Double(str);
        }
    }
    
    /**
     * <p>
     * Requests NEOs via the HTTP connection.
     * </p>
     * @param logger the logger for logging.
     * @param signature the name of the caller for logging.
     * @param requestParameters the request parameters.
     * @return the list of NEO instances.
     * @throws ServiceException if there are any error.
     */
    private List<NEO> requestNEOs(Logger logger, String signature, String requestParameters)
        throws ServiceException {
        System.out.println("Sending request to MPC:" + requestParameters);
        List<NEO> neos = new ArrayList<NEO>();

        HttpURLConnection connection = null;
        DataOutputStream wr = null;
        BufferedReader rd = null;

        try {
            
            connection = (HttpURLConnection) new URL(mpCheckerCGIURL).openConnection();
            connection.setRequestMethod(METHOD_POST);
            connection.setRequestProperty(HEADER_CONTENT_TYPE, CONTENT_TYPE);
            connection.setRequestProperty(HEADER_CONTENT_LENGTH,
                    Integer.toString(requestParameters.getBytes(UTF_8_ENCODING).length));
            connection.setRequestProperty(HEADER_CONTENT_LANGUAGE, EN_US_LANG);

            connection.setUseCaches(false);
            connection.setDoInput(true);
            connection.setDoOutput(true);

            wr = new DataOutputStream(connection.getOutputStream ());
            wr.writeBytes(requestParameters);
            wr.flush ();

            // close it immediately to show server that all data sent
            wr.close ();
            wr = null;

            InputStream is = connection.getInputStream();
            rd = new BufferedReader(new InputStreamReader(is));
            String line;
            boolean resultDetected = false;

            // NEOs are shown between the all pairs of <pre> and </pre>
            boolean inPreSection = false;
            while((line = rd.readLine()) != null) {
                
                if (line.contains(PRE_TAG_START)) {
                    inPreSection = true;
                } else if (line.contains(PRE_TAG_END)) {
                    inPreSection = false;
                    resultDetected = false;
                } else if (line.trim().length() == 0 && inPreSection) {
                    resultDetected = true;
                } else if (resultDetected) {
                    // parse the line content to NEO instance
                    try {
                        neos.add(parseNeo(line));
                    } catch (ServiceException e) {
                        // log and re-throw
                        throw LoggingHelper.logException(logger, signature, e);
                    }
                }
            }
            
        } catch (MalformedURLException e) {
            throw LoggingHelper.logException(
                    logger, signature, new ServiceException("The URL is invalid.", e));
        } catch (IOException e) {
            throw LoggingHelper.logException(
                    logger, signature, new ServiceException("IO errors occur when getting the response.", e));
        } finally {
            if (connection != null) {
                connection.disconnect();
            }

            if (wr != null) {
                try {
                    wr.close();
                } catch (IOException e) {
                    // ignore
                }
            }

            if (rd != null) {
                try {
                    rd.close();
                } catch (IOException e) {
                    // ignore
                }
            }
        }
        return neos;
    }

    /**
     * <p>
     * Parses a string to NEO instance.
     * </p>
     * @param line the line of the string to parse
     * @return the NEO instance.
     * @throws ServiceException if there are any format error
     */
    private NEO parseNeo(String line) throws ServiceException {
        NEO neo = new NEO();

        String objectDesignation = subString(line, START_INDEX_OBJECT_DESIGNATION, START_INDEX_RA).trim();
        String[] raParts = subString(line, START_INDEX_RA, START_INDEX_DECL).split(" ");
        String[] declParts = subString(line, START_INDEX_DECL, START_INDEX_V).split(" ");
        String vString = subString(line, START_INDEX_V, START_INDEX_OFFSET_RA).trim();
        String offsetRA = subString(line, START_INDEX_OFFSET_RA, START_INDEX_OFFSET_DECL).trim();
        String offsetDecl = subString(line, START_INDEX_OFFSET_DECL, START_INDEX_MOTION_RA).trim();
        String motionRA = subString(line, START_INDEX_MOTION_RA, START_INDEX_MOTION_DECL).trim();
        String motionDecl = subString(line, START_INDEX_MOTION_DECL, START_INDEX_ORBIT).trim();
        String orbit = subString(line, START_INDEX_ORBIT, END_INDEX_ORBIT).trim();

        // object designation
        neo.setObjectDesignation(objectDesignation);

        // RA
        neo.setRightAscension(raParts[0] + HOUR_LETTER + raParts[1] + MINUTE_LETTER + raParts[2] + SECOND_LETTER);
        double part1 = Double.parseDouble(raParts[0]);
        double part2 = Double.parseDouble(raParts[1]);
        double part3 = Double.parseDouble(raParts[2]);
        double v = Math.abs(part1) + part2 / 60 + part3 / 3600;
        if (part1 < 0) {
            v = -v;
        }
        v *= 15;
        neo.setRightAscensionValue(v);
        
        // declination
        neo.setDeclination(declParts[0] + DEGREE_CHAR + declParts[1] + DEGREE_MINUTE_CHART + declParts[2]
                + DEGREE_SECOND_CHAR);
        part1 = Double.parseDouble(declParts[0]);
        part2 = Double.parseDouble(declParts[1]);
        part3 = Double.parseDouble(declParts[2]);
        v = Math.abs(part1) + part2 / 60 + part3 / 3600;
        if (part1 < 0) {
            v = -v;
        }
        neo.setDeclinationValue(v);
        
        
        // V
        if (vString.length() == 0) {
            neo.setV(EMPTY_DEFAULT_V);
        } else {
            try {
                neo.setV(Double.parseDouble(vString));
            } catch (NumberFormatException e) {
                throw new ServiceException("Invalid Format to parse v.", e);
            }
        }

        // offset RA
        neo.setOffsetRightAscension(offsetRA);
        // offset declination
        neo.setOffsetDeclination(offsetDecl);

        // motion RA
        neo.setMotionRightAscension(motionRA);
        neo.setMotionDeclination(motionDecl);

        // orbit
        neo.setOrbit(orbit);

        // return the result
        return neo;
    }

    /**
     * <p>
     * Builds the request parameters.
     * </p>
     * @param logger the logger for logging.
     * @param signature the signature for logging.
     * @param criteria the criteria.
     * @return the built string.
     * @throws ServiceException if there are any error.
     */
    private String buildRequestParameters(Logger logger, String signature, NEOSearchCriteria criteria)
        throws ServiceException {

        StringBuffer sb = new StringBuffer();

        Calendar calendar = Calendar.getInstance();
        calendar.setTime(criteria.getDate());
        // year, month, day
        sb.append(PARAM_YEAR).append(calendar.get(Calendar.YEAR));
        sb.append(PARAM_MONTH).append(calendar.get(Calendar.MONTH) + 1);
        double day = calendar.get(Calendar.DATE)
                + (calendar.get(Calendar.HOUR_OF_DAY) * SECONDS_PER_HOUR + calendar.get(Calendar.MINUTE)
                        * SECONDS_PER_MINUTE + calendar.get(Calendar.SECOND)) / SECONDS_PER_DAY;

        sb.append(PARAM_DAY).append(String.format("%.3f", day));

        try {

            if (criteria.getObservations() == null) {
                // search by position
                sb.append(PARAM_WHICH).append(WHICH_VALUE_POS);
                sb.append(PARAM_RA).append(URLEncoder.encode(criteria.getRightAscension(), UTF_8_ENCODING));
                sb.append(PARAM_DECL).append(URLEncoder.encode(criteria.getDeclination(), UTF_8_ENCODING));
                sb.append(PARAM_TEXT_AREA);
            } else {
                // search by observations
                sb.append(PARAM_WHICH).append(WHICH_VALUE_OBS);
                sb.append(PARAM_TEXT_AREA).append(URLEncoder.encode(criteria.getObservations(), UTF_8_ENCODING));
                sb.append(PARAM_RA);
                sb.append(PARAM_DECL);
            }
        } catch (UnsupportedEncodingException e) {
            throw LoggingHelper.logException(logger, signature,
                    new ServiceException("Not supporting the encoding", e));
        }

        // radius
        sb.append(PARAM_RADIUS).append(String.format("%.2f", criteria.getRadius()));
        // v limit
        sb.append(PARAM_LIMIT).append(String.format("%.2f", criteria.getV()));
        // observatory code
        sb.append(PARAM_OC).append(criteria.getObservatoryCode());
        // other parameters are left to default values
        sb.append(PARAM_OTHERS);

        // return the parameters
        return sb.toString();
    }

    /**
     * <p>
     * Helper method for extracting the substring.
     * </p>
     * @param string the string to extract
     * @param startIndex the start index of the substring.
     * @param endIndex the end index of the substring
     * @return the substring.
     */
    private String subString(String string, int startIndex, int endIndex) {
        if (startIndex >= string.length()) {
            return EMPTY_STRING;
        }
        if (endIndex > string.length()) {
            endIndex = string.length();
        }
        return string.substring(startIndex, endIndex);
    }
}

