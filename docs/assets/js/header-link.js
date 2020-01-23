$(function () {
    $("h1, h2, h3, h4, h5, h6").each(function () {
        var id = $(this).attr("id");
        if (id) {
            $(this).append($("<a />").addClass("header-link").attr("href", "#" + id).html('<i class="fa fa-link"></i>'));
        }
    });
});
