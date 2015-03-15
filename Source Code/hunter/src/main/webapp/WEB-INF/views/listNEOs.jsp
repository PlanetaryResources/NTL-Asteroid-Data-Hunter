<%--
  - Author: TCSASSEMBLER
  - Version: 1.0
  - Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
  -
  - Description: The NEOs page.
--%>
<%@ page import="com.fasterxml.jackson.databind.ObjectMapper" %>
<%@taglib uri="http://www.springframework.org/tags/form" prefix="form" %>



<%@include file="header.jsp" %>

<script>
var lastImageInfo = <% out.print(new ObjectMapper().writeValueAsString(request.getAttribute("lastImageInfo"))); %>
</script>
        <div class="page-content sky-page">
            <div class="round-box sky-page-content">
                <div class="title-panel">
                    <img src='../i/search-icon.png' /><span>WHAT'S UP IN THE SKY</span>
                </div>
                <div class="seperator"></div>
                <div class="sky-search-panel">
                    <div class="sky-search-row">
                        <span>DATE:</span><input name="neoDate" id="neoDate" type="text" class="blue-input width120 dateselection width120-ie9"/>
                        <span class="spangap">TIME:</span><input name="neoTime" id="neoTime" type="text" class="blue-input width80"/>
                        <span class="spangap">RA:</span><input name="raTXT" id="raTXT" type="text" class="blue-input width120 width120-ie9"/>
                        <span class="spangap">DEC:</span><input name="decTXT" id="decTXT" type="text" class="blue-input width120 width120-ie9"/>
                        <span class="spangap">OBSERVATIONS:</span><input name="observationsTXT" id="observationsTXT" type="text" class="blue-input width200"/>
                    </div>
                    <div class="sky-search-row">
                        <span>RADIUS:</span><input name="radiusTXT" id="radiusTXT" type="text" class="blue-input width40" /><span class="blue-text">ARC MINS</span>
                        <span class="spangap"><span class="spangap"></span>LIMITING MAGNITUDE (V):</span><input name="vTXT" id="vTXT" type="text" class="blue-input width40"/>
                        <span class="spangap">OBSERVATORY CODE:</span><input name="obCodeTXT" id="obCodeTXT" type="text" class="blue-input width40"/>
                        <a href="javascript:void(0)" class="search-neo-btn"><img src="../i/search-icon.png" />SEARCH</a>
                    </div>
                    <div class="seperator"></div>
                </div>
                <div class="asteroid-table-container neo-list-table">

                    <table>
                        <colgroup>
                    		<col style="width:17%">
                    		<col style="width:15%">
                    		<col style="width:12%">
                    		<col style="width:6%">
                    		<col style="width:10%">
                    		<col style="width:11%">
                    		<col style="width:10%">
                    		<col style="width:12%">
                    		<col >
                    	</colgroup>
                        <thead>
                            <tr>
                                <th class="sortable" data-sort-by="objectDesignation">OBJECT DESIGNATION</th>
                                <th class="sortable" data-sort-by="rightAscension">RIGHT ASCENSION</th>
                                <th class="sortable" data-sort-by="declination">DECLINATION</th>
                                <th class="sortable" data-sort-by="v">V</th>
                                <th class="sortable" data-sort-by="offsetRightAscension">OFFSETS RA</th>
                                <th class="sortable" data-sort-by="offsetDeclination">OFFSETS DEC</th>
                                <th class="sortable" data-sort-by="motionRightAscension">MOTION RA</th>
                                <th class="sortable" data-sort-by="motionDeclination">MOTION DEC</th>
                                <th class="sortable" data-sort-by="orbit">ORBIT</th>
                            </tr>
                        </thead>
                        <tbody>

                        </tbody>
                    </table>

                    <div class="page-section"></div>
                </div>

    </div>
        </div>
<%@include file="footer.jsp" %>