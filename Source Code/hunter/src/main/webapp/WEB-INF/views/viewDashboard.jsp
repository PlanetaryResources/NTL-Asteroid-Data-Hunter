<%--
  - Author: TCSASSEMBLER
  - Version: 1.0
  - Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
  -
  - Description: The dash board page.
--%>
<%@ page import="com.fasterxml.jackson.databind.ObjectMapper" %>

<%@include file="header.jsp" %>

<script>
    var latestItem = <% out.print(new ObjectMapper().writeValueAsString(request.getAttribute("latestDetectionItem"))); %>
    var gSearchResult = <% out.print(new ObjectMapper().writeValueAsString(request.getAttribute("searchResult"))); %>
</script>
        <div class="page-content dashboard-page">
            <div class="upload-section ">
                <img src="../i/select-file-icon.png"/>
                <span>ASTEROID DATA: </span>
                <input type="text" class="uploadfile-input" id="fileUploadTxt" accept=".H"></input>
                <a class="browse-btn btn" href="#">
                    <img src="../i/browse-icon.png"/>BROWSE
                </a>
                <span class="obs-code-text">OBSERVATORY CODE: </span>
                <input type="text" class="observatory-code-input" value="703"></input>
                <input type="file" class="hideUploadFile" id="fileUpload" multiple/>
                <a class="start-hunt-btn btn" href="#">
                    <img src="../i/hunt-icon.png"/>START HUNTING
                </a>
            </div>
            <div>
            
            <div class="progress-section hidden">
                <div class="progress-bar" ><div class="progress-bar-bg"></div><div class="progress-bar-stripe"></div></div>
				<div class="progress-tip">10%</div>
            </div>
            </div>
            <div class="seperator">
            </div>
            <div class="round-box latest-result">
                <div class="latest-result-title"><img src="../i/latest-icon.png"/>SEARCH RESULTS
                    <div class="constrast-slider">CONTRAST&nbsp;&nbsp;<div class="slider" ></div></div>
                </div>
                
                <div class="latest-result-content">
                
					
                    <div class="result-img box">     
                    	<c:if test="${latestDetectionItem != null}">               
                        <div class="img-container">
                        	<img class="bg-img" data-frame-index="0" src="/frameimage/${latestDetectionItem.id}/frame/0/image" alt="">
                        	<img class="red-circle" src="../i/red-circle.png" />
                        </div>  
                        </c:if>  
                    </div>

                    <div class="result-info box">
                        <div class="info-item">
                            <div class="info-item-name">NEO?</div>
                            <div class="info-item-value blue-text">Asteroid</div>
                            <div class="seperator"></div>
                        </div>
                        <div class="info-item">
                            <div class="info-item-name">NEO?</div>
                            <div class="info-item-value blue-text">Asteroid</div>
                            <div class="seperator"></div>
                        </div>
                    </div>

                    <div class="info-item-template hidden">
                        <div class="info-item-name">NEO?</div>
                        <div class="info-item-value blue-text">Asteroid</div>
                        <div class="seperator"></div>
                    </div>

                    <div class="result-controls">
                        <a href="javascript:void(0)" class="play-pause-btn" title="Will play all 4 images in order"><img class="play-btn" src="../i/play-btn.png"/><img class="pause-btn hidden" src="../i/pause-btn.png"/></a>
                        <div class="frame-control">
                            <a href="javascript:void(0);" class="frame-index-btn" data-frame-index="0">1</a>
                            <a href="javascript:void(0);" class="frame-index-btn" data-frame-index="1">2</a>
                            <a href="javascript:void(0);" class="frame-index-btn" data-frame-index="2">3</a>
                            <a href="javascript:void(0);" class="frame-index-btn" data-frame-index="3">4</a>
                        </div>
                        <a href="javascript:void(0)" class="invert-btn" title="Will invert the image"><img src="../i/invert-btn.png"/></a>
                        <a href="javascript:void(0)" class="flicker-btn hidden"><img src="../i/circle-btn.png"/></a>
                    </div>
                </div>
                <div class="latest-result-no-content hidden">
                    <p>No Recent Results</p>
                </div>
            </div>
            <div class="round-box search-history">
                <div class="search-history-title">
                     <img src="../i/search-icon.png"/>
                     <span>Object Identification History</span>
                     <a class="search_more blue-text" href="../detectionItems">MORE ></a>
                </div>
                <div class="asteroid-table-container search-history-table">
                    <table>
                    	<colgroup>
                    		<col style="width:32%">
                    		<col style="width:37%">
                    		<col style="width:31%">
                    	</colgroup>
                        <thead>
                            <tr>
                                <th>IMAGE DATE/TIME</th>
                                <th>OBSERVATORY CODE</th>
                                <th>CELESTIAL POINTING</th>
                            </tr>
                        </thead>
                        <tbody>

                        </tbody>
                    </table>
                </div>
            </div>
        </div>

<%@include file="footer.jsp" %>