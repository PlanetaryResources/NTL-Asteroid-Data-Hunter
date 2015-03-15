/**
 * The javascript using in the common pages.
 * @author TCSASSEMBLER
 * @version 1.0
 */
function paddingZero(num, size) {
    var s = num+"";
    while (s.length < size) s = "0" + s;
    return s;
}

function formatDate(timestamp, formatString) {
    if (!formatString) {
        formatString = "MM/dd/yyyy hh:mma"
    }
    return $.format.date(new Date(timestamp), formatString)
}
function convertDMS(v, delims) {
    // assume v >= 0
    var d = Math.floor(v)
    var m = Math.floor((v - d) * 60)
    var s = Math.floor((v - d - m / 60) * 3600)
    return paddingZero(d, 2) + delims[0] + paddingZero(m) + delims[1] + paddingZero(s) + delims[2]
}

function formatHMS(v) {
	v = v / 15
    if (v >= 0) {
        return convertDMS(v, 'hms')
    } else {
        return '-' + convertDMS(-v, 'hms')
    }
}

function formatDMS(v) {
    if (v >= 0) {
        return '+' + convertDMS(v, '\xB0\'\"')
    } else {
        return '-' + convertDMS(-v,  '\xB0\'\"')
    }
}
function formatLatitude(v) {
    if (v >= 0) {
        return convertDMS(v, '\xB0\'\"') + "N"
    } else {
        return convertDMS(-v, '\xB0\'\"') + "S"
    }
}

function formatLogitude(v) {
    if (v >= 0) {
        return convertDMS(v, '\xB0\'\"') + "E"
    } else {
        return convertDMS(-v, '\xB0\'\"') + "W"
    }
}


function AsteroidTable(tableContainer) {
    this.enableScroll = false
    var that = this
    this.columnNum = tableContainer.find("thead tr:eq(0) th").length
    this.renderCell = function(item, columnIndex, cell) {
        // should be override
    }
    this.changePageCallBack = function(pageIndex) {
       // should be override
    }
    this.onSortColumn = function(sortBy, sortType) {

    }
    this.onScrollToAppend = function() {
         // show be override
    }
    var theTable = tableContainer.find('table')
    var pageSection = tableContainer.find('.page-section')
    this.renderTable = function(items) {
        // copy the row from template
        var tbody = theTable.find('tbody')
        tbody.html('')
        for (var i = 0; i < items.length; i++) {
            var item = items[i]
            var row = $('<tr></tr>')
            tbody.append(row)
            for (var j = 0; j < that.columnNum; j++) {
                var cell = $('<td></td>')
                cell.addClass('blue-text')
                row.append(cell)
                that.renderCell(item, j, cell)
            }

        }
        // change the bg height
       /* if (that.enableScroll) {
            var contentHeight = $('.center-content').height();
            if (contentHeight > 728) {
                $('body').css('height', contentHeight + 'px')
            } else {
                $('body').css('height', '')
            }

        }*/

    }

    this.appendTable = function(items) {
        var tbody = theTable.find('tbody')

        for (var i = 0; i < items.length; i++) {
            var item = items[i]
            var row = $('<tr></tr>')
            tbody.append(row)
            for (var j = 0; j < that.columnNum; j++) {
                var cell = $('<td></td>')
                cell.addClass('blue-text')
                row.append(cell)
                that.renderCell(item, j, cell)
            }

        }
        // change the bg height
        /*if (that.enableScroll) {
            var contentHeight = $('.center-content').height();
            if (contentHeight > 728) {
                $('body').css('height', contentHeight + 'px')
            } else {
                $('body').css('height', '')
            }

        }*/
    }

    this.renderSortColumns = function(sortBy, sortType) {
        theTable.find("thead th").each(function() {
            $(this).find("img.sort-arrow").remove()

            if ($(this).data('sort-by') == sortBy) {

                var arrow = $('<img />')
                arrow.addClass('sort-arrow')
                if (sortType == 'ASC') {
                    arrow.addClass('sort-asc')
                    arrow[0].src = '../i/up-arrow.png'
                } else {
                    arrow.addClass('sort-desc')
                    arrow[0].src = '../i/down-arrow.png'
                }
                $(this).append(arrow)
            } else {
            	 var arrow = $('<img />')
                 arrow.addClass('sort-arrow')
                 arrow[0].src = '../i/sortable-arrow.png'
                 $(this).append(arrow)
            }
        })
    }
    this.showLoading = function() {
        tableContainer.append("<div class=\"append-loading\"><img src=\"../i/loading.gif\"/></div>")
    }
    this.hideLoading = function() {
        tableContainer.find(".append-loading").remove()
    }

    this.renderPaging = function(totalPages, pageNumber) {
        if (pageSection == null) {
            return
        }
        pageSection.html('')

        if (totalPages == 0 || totalPages == 1) {
            return;
        }

        var prevPage = pageNumber - 1 >= 1 ? pageNumber - 1 : 1
        var btn = createPageBtn(1, "<")
        if (pageNumber <= 1) {
            // disable the prev button
            btn.addClass("disabled")
        }
        pageSection.append(btn)

        var btnsNum = 8
        var beginIndex = pageNumber - parseInt(btnsNum / 2)
        if (beginIndex <= 1) {
            beginIndex = 1
        }
        var endIndex = beginIndex + btnsNum
        if (endIndex > totalPages) {
            endIndex = totalPages
        }
        if (beginIndex != 1) {
            btn = createPageBtn(beginIndex - 1, '...')
            pageSection.append(btn)
        }
        for (var i = beginIndex; i <= endIndex; i++) {
            btn = createPageBtn(i, i + "")
            if (i == pageNumber) {
                btn.addClass('page-current-btn')
            }
            pageSection.append(btn)
        }
        if (endIndex != totalPages) {
            btn = createPageBtn(endIndex + 1, '...')
            pageSection.append(btn)
        }
        var nextPage = pageNumber + 1 <= totalPages ? pageNumber + 1: totalPages
        btn = createPageBtn(totalPages, '>')
        if (pageNumber >= totalPages) {
            // disable the prev button
            btn.addClass("disabled")
        }
        btn.addClass('page-next-btn')
        pageSection.append(btn)
        pageSection.addClass('hidden')
    }

    function createPageBtn(gotoPage, text) {
        btn = $('<a class="page-btn" href="javascript:void(0);"></a>').text(text)
        btn.data('goto', gotoPage)
        pageSection.append(btn)
        return btn
    }

    tableContainer.on('click', '.page-btn', function() {
        if ($(this).hasClass("disabled")) {
            return
        }
        that.changePageCallBack($(this).data('goto'))
    })

    tableContainer.on('click', 'thead th.sortable', function() {
        var sortBy = $(this).data('sort-by')
        var sortType = 'ASC'
        if ($(this).find('img.sort-arrow.sort-asc').length > 0) {
            sortType = 'DESC'
        }
        that.onSortColumn(sortBy, sortType)
    })
    function isScrolledIntoView(elem)
    {
        var $elem = $(elem);
        var $window = $(window);

        var docViewTop = $window.scrollTop();
        var docViewBottom = docViewTop + $window.height();

        var elemTop = $elem.offset().top;
        var elemBottom = elemTop + $elem.height();

        // console.log("view bottom:" + docViewBottom + " elem bottom: " + elemBottom)
        return ((elemBottom <= docViewBottom) && (elemTop >= docViewTop));
    }
    $(document).on('scroll', function() {
        if (!that.enableScroll) {
            return;
        }
        var lastRow = tableContainer.find("tbody tr:last")
        if (lastRow.length == 0) {
            // no data
            return
        }
        if (isScrolledIntoView(lastRow)) {
            // last row visible
            that.onScrollToAppend()
        }

    })
}

var ignoreSlideChange = false

function onSliderChange() {

    var contrastValue = $('.constrast-slider .slider').slider("value")
    $('.customMarker').css('width', contrastValue / 255.0 * 100 + '%');
    if (ignoreSlideChange) {
        ignoreSlideChange = false
        return
    }

    var item = $('.img-container .bg-img')
    var vjsAPI = item.data('vintageJS');
    item.data('contrastValue', contrastValue)
    vjsAPI = item.data('vintageJS');
    vjsAPI.vintage({
        contrast  : contrastValue
    })
}

$(document).ready(function() {
    // setup the headers
    $(document).on('click', '.header-tab', function() {
        var page = $(this).data('page')
        window.location.href = '/' + page
    })

    $('.slider').slider({"value": 0, "max":255, change:onSliderChange,
    	slide: function() {
    		var contrastValue = $('.constrast-slider .slider').slider("value")
    		$('.customMarker').css('width', contrastValue / 255.0 * 100 + '%');
    	}
    });
    $('.slider').prepend("<div class='customMarker'></div>");

    $('.img-container .bg-img').each(function() {
        $(this).vintage({onStop: function() {
            $(document).trigger('RESET_ELEVATE_EVENT')
        }}, null)
    })
})

function showLoading(element) {
    element.addClass('hidden')
    $("<div></div>").addClass("loading-view").appendTo(element.parent())
}

function hideLoading(element) {
    element.removeClass('hidden')
    element.parent().find(".loading-view").remove()
}

function showItemNotFound(element) {
    element.addClass('hidden')
    $("<div>No results were found.</div>").addClass("notfound-view").appendTo(element.parent())
}

function hideItemNotFound(element) {
    element.removeClass('hidden')
    element.parent().find(".notfound-view").remove()
}


