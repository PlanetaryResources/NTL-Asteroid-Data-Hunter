/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package gov.nasa.asteroid.hunter.services.impl;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import gov.nasa.asteroid.hunter.models.NEO;
import gov.nasa.asteroid.hunter.models.NEOSearchCriteria;
import gov.nasa.asteroid.hunter.models.SearchResult;
import gov.nasa.asteroid.hunter.services.NEOService;
import gov.nasa.asteroid.hunter.services.ServiceException;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Date;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.junit4.SpringJUnit4ClassRunner;

import com.fasterxml.jackson.databind.ObjectMapper;

/**
 * <p>
 * This class is the unit tests for class <code>NEOService</code>.
 * </p>
 * @author TCSASSEMBLER
 * @version 1.0
 */
@RunWith(SpringJUnit4ClassRunner.class)
@ContextConfiguration(locations = { "classpath:applicationContext.xml" })
public class NEOServiceImplTest {

    /**
     * <p>
     * Represents the test file for observations.
     * </p>
     */
    private static final String TEST_OBSERVATIONS = "/observations.txt";

    /**
     * <p>
     * Represents the neo service for tests.
     * </p>
     */
    @Autowired
    private NEOService neoService;

    /**
     * <p>
     * Searches the observations.
     * </p>
     * @throws Exception to JUnit.
     */
    @Test
    public void testSearchByObservations() throws Exception {
        NEOSearchCriteria criteria = new NEOSearchCriteria();
        criteria.setDate(new Date());
        criteria.setObservations(readFile(TEST_OBSERVATIONS));
        criteria.setRadius(200);
        criteria.setV(20);
        criteria.setObservatoryCode("500");
        SearchResult<NEO> result = neoService.search(criteria);
        assertTrue("Should have value.", result.getValues().size() > 0);
        for (NEO neo : result.getValues()) {
            assertTrue("Should have correct V.",
                    Double.compare(neo.getV(), 20) <= 0 || Double.compare(neo.getV(), 99) == 0);
        }

        System.out.println(new ObjectMapper().writeValueAsString(result));
    }

    /**
     * <p>
     * Searches the position.
     * </p>
     * @throws Exception to JUnit.
     */
    @Test
    public void testSearchByPosition() throws Exception {
        NEOSearchCriteria criteria = new NEOSearchCriteria();
        criteria.setDate(new Date());
        criteria.setRadius(200);
        criteria.setV(20);
        criteria.setObservatoryCode("500");
        criteria.setRightAscension("14 03 08");
        criteria.setDeclination("+01 48.3");
        SearchResult<NEO> result = neoService.search(criteria);
        assertTrue("Should have value.", result.getValues().size() > 0);
        for (NEO neo : result.getValues()) {
            assertTrue("Should have correct V.",
                    Double.compare(neo.getV(), 20) <= 0 || Double.compare(neo.getV(), 99) == 0);
        }
        System.out.println(new ObjectMapper().writeValueAsString(result));
    }

    /**
     * <p>
     * Searches if there are no result.
     * </p>
     * @throws Exception to JUnit.
     */
    @Test
    public void testSearchNoResult() throws Exception {
        NEOSearchCriteria criteria = new NEOSearchCriteria();
        criteria.setDate(new Date());
        criteria.setRadius(5);
        criteria.setV(1);
        criteria.setObservatoryCode("500");
        criteria.setRightAscension("14 03 08");
        criteria.setDeclination("+01 48.3");
        SearchResult<NEO> result = neoService.search(criteria);
        assertTrue("Should have no value.", result.getValues().size() == 0);
        assertTrue("Page number should be zero for no result.", result.getTotalPages() == 0);
        assertTrue("Total Number should be zero for no result.", result.getTotal() == 0);
    }

    /**
     * <p>
     * Tests the search paging is correct or not.
     * </p>
     * @throws Exception to JUnit.
     */
    @Test
    public void testSearchPaging() throws Exception {
        // search all first
        Date now = new Date();
        NEOSearchCriteria criteria = new NEOSearchCriteria();
        criteria.setDate(now);
        criteria.setRadius(200);
        criteria.setV(20);
        criteria.setObservatoryCode("500");
        criteria.setRightAscension("14 03 08");
        criteria.setDeclination("+01 48.3");
        SearchResult<NEO> allResult = neoService.search(criteria);

        // get the second page, 10 per page
        criteria.setPageNumber(2);
        criteria.setPageSize(10);
        SearchResult<NEO> pageResult = neoService.search(criteria);
        assertEquals("Should have 10 results.", 10, pageResult.getValues().size());
        assertEquals("Same for total number.", allResult.getTotal(), pageResult.getTotal());
        // check the items
        for (int i = 0; i < 10; i++) {
            NEO item1 = allResult.getValues().get(i + 10);
            NEO item2 = pageResult.getValues().get(i);
            assertEquals("They should have equal name.", item1.getObjectDesignation(), item2.getObjectDesignation());
        }
    }

    /**
     * <p>
     * Tests the search for sort by objectDesignation
     * </p>
     * @throws ServiceException
     */
    @Test
    public void testSearchSortByObjectDesignation() throws ServiceException {
        // search all first
        Date now = new Date();
        NEOSearchCriteria criteria = new NEOSearchCriteria();
        criteria.setDate(now);
        criteria.setRadius(200);
        criteria.setV(20);
        criteria.setObservatoryCode("500");
        criteria.setRightAscension("14 03 08");
        criteria.setDeclination("+01 48.3");
        criteria.setSortBy("objectDesignation");

        SearchResult<NEO> allResult = neoService.search(criteria);
        String prev = null;
        for (NEO neo : allResult.getValues()) {
            if (prev != null) {
                assertTrue("Should in asending order.", prev.compareTo(neo.getObjectDesignation()) <= 0);
            }
            prev = neo.getObjectDesignation();
        }
    }

    /**
     * <p>
     * Tests the search for sort by declination
     * </p>
     * @throws ServiceException
     */
    @Test
    public void testSearchSortBydeclination() throws ServiceException {
        // search all first
        Date now = new Date();
        NEOSearchCriteria criteria = new NEOSearchCriteria();
        criteria.setDate(now);
        criteria.setRadius(200);
        criteria.setV(20);
        criteria.setObservatoryCode("500");
        criteria.setRightAscension("14 03 08");
        criteria.setDeclination("+01 48.3");
        criteria.setSortBy("declination");

        SearchResult<NEO> allResult = neoService.search(criteria);
        String prev = null;
        for (NEO neo : allResult.getValues()) {
            if (prev != null) {
                assertTrue("Should in asending order.", prev.compareTo(neo.getDeclination()) <= 0);
            }
            prev = neo.getDeclination();
        }
    }

    /**
     * <p>
     * Tests the search for sort by V
     * </p>
     * @throws ServiceException
     */
    @Test
    public void testSearchSortByV() throws ServiceException {
        // search all first
        Date now = new Date();
        NEOSearchCriteria criteria = new NEOSearchCriteria();
        criteria.setDate(now);
        criteria.setRadius(200);
        criteria.setV(20);
        criteria.setObservatoryCode("500");
        criteria.setRightAscension("14 03 08");
        criteria.setDeclination("+01 48.3");

        criteria.setSortBy("v");

        SearchResult<NEO> allResult = neoService.search(criteria);
        double prev = -1;
        for (NEO neo : allResult.getValues()) {
            if (prev >= 0) {
                assertTrue("Should in asending order.", Double.compare(prev, neo.getV()) <= 0);
            }
            prev = neo.getV();
        }
    }


    /**
     * <p>
     * Helper method to read the file content.
     * </p>
     * @throws Exception to JUnit.
     */
    private String readFile(String fileName) throws IOException {
        // transfer to class path file
        fileName = getClass().getResource(fileName).getFile();
        BufferedReader br = new BufferedReader(new FileReader(fileName));
        try {
            StringBuilder sb = new StringBuilder();
            String line = br.readLine();
            while (line != null) {
                if (sb.length() > 0) {
                    sb.append("\n");
                }
                sb.append(line);
                line = br.readLine();
            }
            return sb.toString();
        } finally {
            br.close();
        }
    }
}
