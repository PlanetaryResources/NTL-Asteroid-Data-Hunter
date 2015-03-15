<%--
  - Author: TCSASSEMBLER
  - Version: 1.0
  - Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
  -
  - Description: The search detection items pages.
--%>
<%@taglib uri="http://www.springframework.org/tags/form" prefix="form" %>

<%@ page import="com.fasterxml.jackson.databind.ObjectMapper" %>

<%@include file="header.jsp" %>

<script>

</script>
        <div class="page-content search-page">
            <div class="round-box search-page-content">
                <div class="title-panel">
                    <img src='../i/search-icon.png' /><span class="title-label">Object Identification History</span>
    <div class="clear-history-btn">
    </div>
                    <a class="download-btn" href="/detectionItems/download">
                    </a>

    <div class="search-filter">
    </div>
                </div>
                <div class="seperator"></div>
                <div class="filter-panel hidden">
                    <form:form method="get" modelAttribute="criteria" action="/search/detectionItems">
                        <span>DATE:</span><input name="dateStart" id="dateStart" type="text" class="blue-input width120 dateselection"/>-<input name="dateEnd" id="dateEnd" type="text" class="blue-input width120 dateselection"/>
                        <span>TIME:</span><input name="timeStart" id="timeStart" type="text" class="blue-input width80"/>-<input name="timeEnd" id = "timeEnd" type="text" class="blue-input width80"/>
                        <span>SUBMITTED:</span><select id="submittedOptions" class="width80"><option>All</option><option>Yes</option><option>No</option></select>
                        <a href="javascript:void(0)" class="do-filter-btn">FILTER</a>
                    </form:form>
                     <div class="seperator"></div>
                </div>
                <div class="asteroid-table-container search-list-table">
                    <div class="scrollwrapper">
                    <table>
                    	<colgroup>
                    		<col style="width:12%">
                    		<col style="width:12%">
                    		<col style="width:16%">
                    		<col style="width:14%">
                    		<col style="width:17%">
                    		<col style="width:17%">
                    		<col style="width:12%">
                    	</colgroup>
                        <thead>
                            <tr>
                                <th class="sortable" data-sort-by="date">IMAGE CREATE DATE</th>
                                <th class="sortable" data-sort-by="time">IMAGE CREATE TIME</th>
                                <th class="sortable" data-sort-by="frame.rightAscension">RIGHT ASCENSION</th>
                                <th class="sortable" data-sort-by="frame.declination">DECLINATION</th>
                                <th class="sortable" data-sort-by="frame.observationLatitude">OBSERVATORY CODE</th>
                                <th class="sortable" data-sort-by="significance">SIGNIFICANCE</th>
                                <th class="sortable" data-sort-by="submitted">EXISTS AT MPC?</th>
                            </tr>
                        </thead>
                        <tbody>

                        </tbody>
                    </table>
                    </div>
                    <div class="page-section"></div>
                </div>

    </div>
        </div>
<%@include file="footer.jsp" %>
