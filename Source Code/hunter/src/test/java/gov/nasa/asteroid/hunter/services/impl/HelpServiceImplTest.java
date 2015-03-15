/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package gov.nasa.asteroid.hunter.services.impl;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import gov.nasa.asteroid.hunter.models.HelpItem;
import gov.nasa.asteroid.hunter.models.HelpItemSearchCriteria;
import gov.nasa.asteroid.hunter.models.HelpTopic;
import gov.nasa.asteroid.hunter.models.SearchResult;
import gov.nasa.asteroid.hunter.models.SortType;
import gov.nasa.asteroid.hunter.services.HelpService;

import java.util.List;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.junit4.SpringJUnit4ClassRunner;

/**
 * <p>
 * This class is the unit tests for class <code>HelpServiceImpl</code>.
 * </p>
 * @author TCSASSEMBLER
 * @version 1.0
 */
@RunWith(SpringJUnit4ClassRunner.class)
@ContextConfiguration(locations = { "classpath:applicationContext.xml" })
public class HelpServiceImplTest extends PersistenceBaseTest {

    /**
     * <p>
     * Represents the help service for unit tests.
     * </p>
     */
    @Autowired
    private HelpService helpService;

    /**
     * <p>
     * Gets the help topics.
     * </p>
     * @throws Exception to JUnit.
     */
    @Test
    public void testGetHelpTopics() throws Exception {
        for (int i = 0; i < 10; i++) {
            HelpTopic topic = new HelpTopic();
            topic.setName("topic name" + i);
            persist(topic);
        }

        List<HelpTopic> helpTopics = helpService.getHelpTopics();
        assertEquals("Should have 10 topics.", 10, helpTopics.size());
        for (HelpTopic topic : helpTopics) {
            assertTrue("topic name is invalid.", topic.getName().contains("topic name"));
        }

    }

    /**
     * <p>
     * Tests to search help items.
     * </p>
     * @throws Exception to JUnit.
     */
    @Test
    public void testSearchHelpItems() throws Exception {
        HelpTopic topic = new HelpTopic();
        topic.setName("test");
        persist(topic);
        for (int i = 0; i < 10; i++) {
            HelpItem item = new HelpItem();
            item.setContent("a" + i);
            item.setTitle("b" + i);
            item.setTopic(topic);
            persist(item);
        }
        HelpItemSearchCriteria criteria = new HelpItemSearchCriteria();
        SearchResult<HelpItem> items = helpService.searchHelpItems(criteria);
        assertEquals("items invalid.", 10, items.getValues().size());
    }

    /**
     * <p>
     * Tests to search help items.
     * </p>
     * @throws Exception to JUnit.
     */
    @Test
    public void testSearchHelpItemsSearchForKeywoard() throws Exception {
        HelpTopic topic = new HelpTopic();
        topic.setName("test");
        persist(topic);
        for (int i = 0; i < 10; i++) {
            HelpItem item = new HelpItem();
            item.setContent("a " + i);
            item.setTitle("b " + i);
            item.setTopic(topic);
            persist(item);
        }
        HelpItemSearchCriteria criteria = new HelpItemSearchCriteria();
        criteria.setKeyword("5");
        SearchResult<HelpItem> items = helpService.searchHelpItems(criteria);
        assertTrue("should have value.", items.getValues().size() > 0);
        for (HelpItem item : items.getValues()) {
            assertTrue("should contains key word.", item.getContent().contains("5")
                    || item.getTitle().contains("5"));
        }
    }

    /**
     * <p>
     * Tests to search help items.
     * </p>
     * @throws Exception to JUnit.
     */
    @Test
    public void testSearchHelpItemsSearchForTopic() throws Exception {
        HelpTopic topic = new HelpTopic();
        topic.setName("test");
        persist(topic);
        for (int i = 0; i < 10; i++) {
            HelpItem item = new HelpItem();
            item.setContent("a " + i);
            item.setTitle("b " + i);
            item.setTopic(topic);
            persist(item);
        }
        HelpItemSearchCriteria criteria = new HelpItemSearchCriteria();
        criteria.setTopicId(topic.getId());
        SearchResult<HelpItem> items = helpService.searchHelpItems(criteria);
        assertTrue("should have value.", items.getValues().size() > 0);
        for (HelpItem item : items.getValues()) {
            assertTrue("should contains key word.", item.getTopic().getId() == topic.getId());
        }
    }

    /**
     * <p>
     * Tests the pagination of searching the help items.
     * </p>
     * @throws Exception to JUnit.
     */
    @Test
    public void testSearchHelpItemsPage() throws Exception {
        HelpTopic topic = new HelpTopic();
        topic.setName("test");
        persist(topic);
        for (int i = 0; i < 100; i++) {
            HelpItem item = new HelpItem();
            item.setContent("a");
            item.setTitle("b");
            item.setTopic(topic);
            persist(item);
        }
        HelpItemSearchCriteria criteria = new HelpItemSearchCriteria();
        criteria.setPageNumber(2);
        criteria.setPageSize(5);
        SearchResult<HelpItem> items = helpService.searchHelpItems(criteria);
        assertEquals("items invalid.", 5, items.getValues().size());
    }

    /**
     * <p>
     * Tests order by for search help items.
     * </p>
     * @throws Exception to JUnit.
     */
    @Test
    public void testSearchHelpItemsOrderBy() throws Exception {
        HelpTopic topic = new HelpTopic();
        topic.setName("test");
        persist(topic);
        for (int i = 0; i < 100; i++) {
            HelpItem item = new HelpItem();
            item.setContent("a");
            item.setTitle("b" + i);
            item.setTopic(topic);
            persist(item);
        }
        HelpItemSearchCriteria criteria = new HelpItemSearchCriteria();
        criteria.setSortBy("title");
        criteria.setSortType(SortType.DESC);
        SearchResult<HelpItem> items = helpService.searchHelpItems(criteria);
        assertEquals("items invalid.", 100, items.getValues().size());
        HelpItem prev = null;
        for (HelpItem item : items.getValues()) {
            if (prev != null) {
                assertTrue("Should desc", item.getTitle().compareTo(prev.getTitle()) <= 0);
            }
            prev = item;
        }
    }

    /**
     * <p>
     * Tests order by for search help items (in asc order).
     * </p>
     * @throws Exception to JUnit.
     */
    @Test
    public void testSearchHelpItemsASC() throws Exception {
        HelpTopic topic = new HelpTopic();
        topic.setName("test");
        persist(topic);
        for (int i = 0; i < 100; i++) {
            HelpItem item = new HelpItem();
            item.setContent("a");
            item.setTitle("b" + i);
            item.setTopic(topic);
            persist(item);
        }
        HelpItemSearchCriteria criteria = new HelpItemSearchCriteria();
        criteria.setSortBy("title");
        criteria.setSortType(SortType.ASC);
        SearchResult<HelpItem> items = helpService.searchHelpItems(criteria);
        assertEquals("items invalid.", 100, items.getValues().size());
        HelpItem prev = null;
        for (HelpItem item : items.getValues()) {
            if (prev != null) {
                assertTrue("Should desc", item.getTitle().compareTo(prev.getTitle()) >= 0);
            }
            prev = item;
        }
    }


    /**
     * <p>
     * Tests get help item.
     * </p>
     * @throws Exception to JUnit.
     */
    @Test
    public void testGetHelpItem() throws Exception {
        HelpTopic topic = new HelpTopic();
        topic.setName("test");
        persist(topic);
        long itemId = 0;
        for (int i = 0; i < 10; i++) {
            HelpItem item = new HelpItem();
            item.setContent("a" + i);
            item.setTitle("b" + i);
            item.setTopic(topic);
            persist(item);
            itemId = item.getId();
        }
        HelpItem item = helpService.getHelpItem(itemId);
        assertEquals("item's topic is not correct.", "test", item.getTopic().getName());
        assertEquals("item's content is not correct.", "a" + 9, item.getContent());
    }

}
