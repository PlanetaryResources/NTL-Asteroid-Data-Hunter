/**
 * The javascript using in the SKY page.
 * @author TCSASSEMBLER
 * @version 1.0
 */

$(document).ready(function() {


    // set default value
    $('#neoDate').val(formatDate(new Date(), 'MM/dd/yyyy')).change()
    $('#neoTime').val(formatDate(new Date(), 'hh:mm a')).change()
    $('#raTXT').val('14 03 08')
    $('#decTXT').val('+01 48.3')

    $('#radiusTXT').val('129')
    $('#vTXT').val('22')
    $('#obCodeTXT').val('703')
    if (lastImageInfo) {
        var dateStr = lastImageInfo.dateOBS + 'T' + lastImageInfo.timeOBS
        var obsDate = new Date(dateStr)
        // change to UTC
        obsDate = new Date(obsDate.getTime() + obsDate.getTimezoneOffset() * 60000)
        var rightAscension = (lastImageInfo.rightAscension).replace(/:/g, ' ')
        var dec = lastImageInfo.declination.replace(/:/g, ' ')
        $('#raTXT').val(rightAscension)
        $('#decTXT').val(dec)

        $('#neoDate').val(formatDate(obsDate, 'MM/dd/yyyy')).change()
        $('#neoTime').val(formatDate(obsDate, 'hh:mm a')).change()
        $('#obCodeTXT').val(lastImageInfo.observatoryCode)
    }
    var searcher = null
    var criteria = getNEOCriteria()
    if (criteria !== false) {
        searcher = new NEOSearcher(criteria)
        searcher.search()
    }




    $(document).on('click', '.search-neo-btn', function() {
        var criteria = getNEOCriteria()
        if (criteria === false) {
            return
        }
        searcher = new NEOSearcher(criteria)
        searcher.search()

    })

    var neoTable = new AsteroidTable($('.neo-list-table'))
    neoTable.enableScroll = true
    neoTable.changePageCallBack = function(pageIndex) {
        searcher.pageNumber = pageIndex
        searcher.search()
    }
    neoTable.onScrollToAppend = function() {
        searcher.scrollAppend()
    }
    neoTable.renderCell = function(item, columnIndex, cell) {
        var text = ""
        if (columnIndex == 0) {
            text = item.objectDesignation
        } else if (columnIndex == 1) {
            text = item.rightAscension
        } else if (columnIndex == 2) {
            text = item.declination
            cell.html(text)
            text = null
        } else if (columnIndex == 3) {
            text = item.v
        } else if (columnIndex == 4) {
            text = item.offsetRightAscension
        } else if (columnIndex == 5) {
            text = item.offsetDeclination
        } else if (columnIndex == 6) {
           text = item.motionRightAscension
        } else if (columnIndex == 7) {
            text = item.motionDeclination
        } else if (columnIndex == 8) {
            text = item.orbit
        }
        if (text != null) {
            cell.text(text)
        }
    }
    neoTable.onSortColumn = function(sortBy, sortType) {
        searcher.sortBy = sortBy
        searcher.sortType = sortType
        searcher.search()
    }
    function NEOSearcher(criteria) {
        this.pageSize = 50
        this.pageNumber = 1
        this.sortBy = 'objectDesignation'
        this.sortType = 'ASC'
        var that = this
        var searched = false
        var appending = false
        var nomore = false
        this.search = function() {
            searched = false
            that.pageNumber = 1
            var options = $.extend(criteria, {
                pageSize: that.pageSize,
                pageNumber: that.pageNumber,
                sortBy: that.sortBy,
                sortType: that.sortType})
            // show loading
            showLoading($('.asteroid-table-container'))
            $.ajax({
                url:'/search/neos?' + $.param(options),
                type: "GET",
                headers: {
                    'Accept': 'application/json'
                },
                dataType : 'json',
                success:function(result) {
                    searched = true
                    nomore = false
                    hideLoading($('.asteroid-table-container'))
                    hideItemNotFound($('.asteroid-table-container'))
                    if (result.values.length == 0) {
                        showItemNotFound($('.asteroid-table-container'))
                    }
                    neoTable.renderTable(result.values)
                    neoTable.renderPaging(result.totalPages, result.pageNumber)
                    neoTable.renderSortColumns(that.sortBy, that.sortType)
                },
                error: function() {
                    hideLoading($('.asteroid-table-container'))
                	// we do not alert this, because this error will promp if we manually abort the requests.
                	console.log("Error failed in request the neos.");
                }
            })
        }

        this.scrollAppend = function() {
            if (!searched || appending || nomore) {
                return
            }
            appending = true
            neoTable.showLoading()
            var options = $.extend(criteria, {
                pageSize: that.pageSize,
                pageNumber: that.pageNumber + 1,
                sortBy: that.sortBy,
                sortType: that.sortType})
            $.ajax({
                url:'/search/neos?' + $.param(options),
                type: "GET",
                headers: {
                    'Accept': 'application/json'
                },
                dataType : 'json',
                success:function(result) {
                    neoTable.hideLoading()
                    if (result.values.length == 0) {
                        nomore = true
                    }
                    that.pageNumber = that.pageNumber + 1
                    neoTable.appendTable(result.values)
                    appending = false
                },
                error: function() {
                    neoTable.hideLoading()
                    appending = false
                    alert('error failed...')
                }
            })
        }
    }


    AnyTime.picker('neoDate', {format: "%m/%d/%Y"})
    AnyTime.picker('neoTime', {format: "%h:%i %p"})

    $('.AnyTime-pkr').each(function() {
        var el = $(this).find('.AnyTime-hdr')
        /*el.contents().filter(function(){
         return (this.nodeType == 3);
         }).remove();*/
        el.find('.AnyTime-x-btn').each(function() {
            $(this).contents().filter(function(){
                return (this.nodeType == 3);
            }).remove();
        })
    })

    function getNEOCriteria() {
        var result = {}

        // validate
        var date = new Date($('#neoDate').val() + ' ' + $('#neoTime').val())

        if (isNaN(date)) {
            alert('Invalid Date or Time Format. They are required')
            return false
        }

        result['date'] = formatDate(date, 'MM/dd/yyyy hh:mm')

        var ra = $('#raTXT').val()
        if (!ra) {
            alert("Right Ascension is required.")
            return false
        }
        if (!/^[0-5][0-9] [0-5][0-9] [0-5][0-9](\.\d+)?$/.test(ra)) {
            alert("Right Ascension is not in correct format.")
            return false
        }
        result['rightAscension'] = ra

        var dec = $('#decTXT').val()
        if (!dec) {
            alert("declination is required.")
            return false
        }
        if (!/^[+\-][0-5][0-9] \d+(\.\d+)?$/.test(dec)
            && !/^[+\-][0-5][0-9] [0-5][0-9] \d+(\.\d+)?$/.test(dec)) {
            alert("DEC is not in correct format.")
            return false
        }
        result['declination'] = dec

        var radius = $('#radiusTXT').val()
        if (!radius) {
            alert("radius is required")
            return false
        }
        if (!/^\d+(\.\d+)?$/.test(radius)) {
            alert("radius is not valid")
            return false
        }
        if (radius < 5 || radius > 300) {
            alert('radius is not in the correct range.')
            return false
        }
        result['radius'] = radius

        var v = $('#vTXT').val()
        if (!v) {
            alert("v is required")
            return false
        }
        if (!/^\d+(\.\d+)?$/.test(v)) {
            alert("V is not valid")
            return false
        }
        result['v'] = v
        var observatoryCode   = $('#obCodeTXT').val()
        if (!observatoryCode ) {
            alert("observatoryCode Code is required")
            return false
        }
        result['observatoryCode'] = observatoryCode

        var observations = $('#observationsTXT').val()
        if (observations) {
            if (!/^.* .*$/.test(observations)) {
                alert("observations is not valid")
                return false
            }
            result['observations'] = observations
        }
        return result
    }




})


