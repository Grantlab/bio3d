
/**
$(function() {
    $( "#dialog1" ).dialog({
	autoOpen: false,
    });
    $( "#opener1" ).click(function() {
	$( "#dialog1" ).dialog( "open" );
    });
});
*/

$(function() {
    $( "#dialog" ).dialog({
	autoOpen: false,
	show: {
	    effect: "bounce",
	    duration: 500
	},
	hide: {
	    effect: "blind",
	    duration: 500
	}
    });
    $( "#opener" ).click(function() {
	$( "#dialog" ).dialog( "open" );
    });
});
