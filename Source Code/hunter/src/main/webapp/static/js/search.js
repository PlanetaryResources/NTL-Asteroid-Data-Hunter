/**
 * The javascript using in the search page.
 * @author TCSASSEMBLER
 * @version 1.0
 */
$(document).ready(function() {
    // load the first page data
    var searcher = new Searcher({})

    // load the first page
    searcher.search()

    var searchTable = new AsteroidTable($('.search-list-table'))
    searchTable.enableScroll = true
    searchTable.changePageCallBack = function(pageIndex) {
        searcher.pageNumber = pageIndex
        searcher.search()
    }
    searchTable.onScrollToAppend = function() {
        searcher.scrollAppend()
    }
    searchTable.onSortColumn = function(sortBy, sortType) {
        if (true) {
            searcher.sortBy = sortBy
            searcher.sortType = sortType
            searcher.search()
            return
        }

        var rows = []
        $('.search-list-table tbody tr').each(function() {
            rows.push($(this).clone(true, true))
        })
        $('.search-list-table tbody').html('')

        var columnIndex = 0
        var index = 0
        $('.search-list-table thead').find('th').each(function() {

            if ($(this).data('sort-by') == sortBy) {
                columnIndex = index
            }
            index++
        })

        rows.sort(function(item1, item2) {
            var v1 = item1.find("td:eq(" + columnIndex + ")").text()
            var v2 = item2.find("td:eq(" + columnIndex + ")").text()
            var t = v1 < v2 ? (-1) : (v1 > v2 ? 1 : 0)

            if (sortType == 'DESC') {
                t = -1 * t
            }
            return t
        })
        for (var i = 0; i < rows.length; i++) {
            $('.search-list-table tbody').append(rows[i])
        }

        searchTable.renderSortColumns(sortBy, sortType)

    }
    searchTable.renderCell = function(item, columnIndex, cell) {
        var text = ""
        if (columnIndex == 0) {
            text = formatDate(item.timestamp, 'MM/dd/yyyy')
        } else if (columnIndex == 1) {
            text = formatDate(item.timestamp, 'hh:mma')
        } else if (columnIndex == 2) {
            text = formatHMS(item.frames[0].rightAscension)
        } else if (columnIndex == 3) {
            text = formatDMS(item.frames[0].declination)
        } else if (columnIndex == 4) {
            text = item.observatoryCode
        } else if (columnIndex == 5) {
            text = item.significance.toFixed(2)
        } else if (columnIndex == 6) {
            text = item.submitted ? "Yes" : "No"
            cell.removeClass("blue-text")
            if (item.submitted) {
                cell.addClass("green-text")
            } else {
                cell.addClass("red-text")
            }
        }
        
        cell.text(text)
        cell.closest('tr').addClass('clickable')
        cell.closest('tr').data('item-id', item.id)
    }
    $(document).on('click', '.search-list-table tbody tr', function() {
    	window.location.href = '/detectionItems/' + $(this).data('item-id')
    })
    
    function Searcher(criteria) {
        this.pageSize = 50
        this.pageNumber = 1
        this.sortBy = 'date'
        this.sortType = 'DESC'
        var searched = false
        var appending = false
        var nomore = false
        var that = this
        this.search = function() {
            that.pageNumber = 1
            searched = false
            var options = $.extend(criteria, {
                    pageSize: that.pageSize,
                    pageNumber: that.pageNumber,
                    sortBy: that.sortBy,
                    sortType: that.sortType})
            $.ajax({
                url:'/search/detectionItems?' + $.param(options),
                type: "GET",
                headers: {
                    'Accept': 'application/json'
                },
                dataType : 'json',
                success:function(result) {
                    hideItemNotFound($('.asteroid-table-container'))
                    if (result.values.length == 0) {
                        showItemNotFound($('.asteroid-table-container'))
                    }
                    searched = true
                    nomore = false
                     searchTable.renderTable(result.values)
                     searchTable.renderPaging(result.totalPages, result.pageNumber)
                     searchTable.renderSortColumns(that.sortBy, that.sortType)

                },
                error: function() {
                    alert('error failed...')
                }
            })
        }
        this.scrollAppend = function() {
            console.log("scroll append called.")
            if (!searched || appending || nomore) {
                return
            }
            appending = true
            searchTable.showLoading()
            var options = $.extend(criteria, {
                pageSize: that.pageSize,
                pageNumber: that.pageNumber + 1,
                sortBy: that.sortBy,
                sortType: that.sortType})
            $.ajax({
                url:'/search/detectionItems?' + $.param(options),
                type: "GET",
                headers: {
                    'Accept': 'application/json'
                },
                dataType : 'json',
                success:function(result) {
                    searchTable.hideLoading()
                    if (result.values.length == 0) {
                        nomore = true
                    }
                    that.pageNumber = that.pageNumber + 1
                    searchTable.appendTable(result.values)
                    appending = false
                },
                error: function() {
                    searchTable.hideLoading()
                    appending = false
                    alert('error failed...')
                }
            })
        }
    }

   $(document).on('click', '.search-filter', function() {
        var that = this
        function hideFilter() {
            $(".filter-panel").addClass("hidden")
            $(that).removeClass("filter-on")
            $(".scrollwrapper").removeClass("shorter")
        }
        function showFilter() {
            $(".filter-panel").removeClass("hidden")
            $(that).addClass("filter-on")
            $(".scrollwrapper").addClass("shorter")
        }
         if ($(this).hasClass('filter-on')) {
            hideFilter()
         } else {
            showFilter()
         }

   })

    $(document).on('click', '.do-filter-btn', function() {
        searcher = new Searcher(getFilterCriteria())
        searcher.search()
    })

    $(document).on('click', '.clear-history-btn', function() {
        var r = window.confirm("Are you sure you want to clear the search history?")
        if (r) {
            $.ajax({
                url:"/detectionItems/deleteAll",
                type:"DELETE",
                success: function() {
                    window.location.reload()
                }
            })
        } else {
            // does nothing
        }
    })
    // the filter
    AnyTime.picker('dateStart', {format: "%m/%d/%Y"})
    AnyTime.picker('dateEnd', {format: "%m/%d/%Y"})
    AnyTime.picker('timeStart', {format: "%h:%i %p"})
    AnyTime.picker('timeEnd', {format: "%h:%i %p"})

    $('.AnyTime-pkr').each(function() {
        var dashLine = $('<div class="dash-separate-line"></div>')
        var clearBtn = $('<button>Clear</button>')
        clearBtn.addClass('date-clear-btn')
        $(this).append(dashLine)
        $(this).append(clearBtn)
    })
    $('#AnyTime--dateStart').on('click', '.date-clear-btn', function() {
        $('#dateStart').val('').change()
        $('#AnyTime--dateStart .AnyTime-x-btn').trigger('click')
    })
    $('#AnyTime--dateEnd').on('click', '.date-clear-btn', function() {
        $('#dateEnd').val('').change()
        $('#AnyTime--dateEnd .AnyTime-x-btn').trigger('click')
    })
    $('#AnyTime--timeStart').on('click', '.date-clear-btn', function() {
        $('#timeStart').val('').change()
    })
    $('#AnyTime--timeEnd').on('click', '.date-clear-btn', function() {
        $('#timeEnd').val('').change()
    })
    var lastStartDate = null
    var lastEndDate = null
    $('#dateStart, #dateEnd').on('change', function() {
        var startDateStr = $('#dateStart').val()
        var endDateStr = $('#dateEnd').val()
        if (startDateStr && endDateStr) {
            var startDate = new Date(startDateStr)
            var endDate = new Date(endDateStr)
            if (startDate > endDate) {
                if (lastStartDate !== startDateStr || lastEndDate !== endDateStr) {
                    alert("Start Date must be before the End Date");
                }
                lastStartDate = startDateStr
                lastEndDate = endDateStr

                $(this).val('').change()
            }
        }
    })

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
    $('#timeStart, #timeEnd').on('change', function() {
        var datePrefix = formatDate(new Date(), 'MM/dd/yyyy ')
        var startDateStr = datePrefix + ' ' + $('#timeStart').val()
        var endDateStr = datePrefix + ' ' +  $('#timeEnd').val()
        if ($('#timeStart').val() && $('#timeEnd').val()) {
            var startDate = new Date(startDateStr)
            var endDate = new Date(endDateStr)
            if (startDate > endDate) {
                alert("Start Time must be before the End Date");
                $(this).val('').change()
            }
        }
    })

    function getFilterCriteria() {
        var criteria = {}
        if ($('#dateStart').val()) {
            criteria['dateStart'] = $('#dateStart').val()
        }
        if ($('#dateEnd').val()) {
            criteria['dateEnd'] = $('#dateEnd').val()
        }
        var datePrefix = formatDate(new Date(), 'MM/dd/yyyy ')
        if ($('#timeStart').val()) {

            var timeValue = new Date(datePrefix + $('#timeStart').val())

            criteria['hourStart'] = timeValue.getHours()
            criteria['minuteStart'] = timeValue.getMinutes()
        }
        if ($('#timeEnd').val()) {
            var timeValue = new Date(datePrefix + $('#timeEnd').val())
            criteria['hourEnd'] = timeValue.getHours()
            criteria['minuteEnd'] = timeValue.getMinutes()
        }

        var submitted = $('#submittedOptions option:selected').text()
        if (submitted == 'Yes') {
            criteria['submitted'] = true
        }  else if (submitted == 'No') {
            criteria['submitted'] = false
        }
        return criteria
    }
})