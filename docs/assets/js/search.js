$(function () {
    var query = getQuery(["q", "t", "a", "d"]);

    var targets;
    switch (query.key) {
        case "t":
            targets = ["tags"];
            break;
        case "a":
            targets = ["author"];
            break;
        case "d":
            targets = ["date"];
            break;
        case "q":
        default:
            targets = ["title", "tags", "author", "url", "date", "content"];
            break;
    }
    showPosts(query.words, targets);

    if (query.key == "q") {
        $("#search").val(query.query).focus();
    }
});

function getQuery(keys)
{
    var query = "";
    var key = "";
    var words = [];

    keys.forEach(function (queryKey) {
        var regex = RegExp("[?&]" + queryKey + "=([^&]+)", 'i');
        var matched;
        if (matched = window.location.search.match(regex)) {
            query = decodeURIComponent(matched[1]).replace(/(　| |\+)+/g, ' ');
            words = query.split(' ');
            key = queryKey;
            return false;  // break;
        }
        return true;  // continue;
    });

    return { query: query, key: key, words: words };
}

function showPosts(words, targets)
{
    var getJson = function () {

        var dfd = $.Deferred();
        $.ajax({
            url: baseurl + "/search.json",
            dataType: "json",
            timeout: 3000,  // 3 sec
            success: function (posts) {
                var matchedPosts = [];
                posts.forEach(function (post) {

                    // concatenate target fields as a string.
                    var searchee = "";
                    for (var i = 0; i < targets.length; i++) {
                        var target = post[targets[i]];
                        var targetString = "";
                        if (target instanceof Array) {
                            for (var j = 0; j < target.length; j++) {
                                targetString += target[j];
                            }
                        } else if (typeof target == "object") {
                            for (key in target) {
                                targetString += target[key];
                            }
                        } else {
                            targetString = target;
                        }
                        searchee += targetString;
                    }

                    // matching.
                    var matched = true;
                    words.forEach(function (word) {
                        var regex = new RegExp(word, 'i');
                        if (searchee.match(regex) == null) {
                            matched = false;
                            return false;  // break;
                        }
                        return true;  // continue;
                    });

                    if (matched) {
                        matchedPosts.push(post);
                    }
                });

                dfd.resolve(matchedPosts);
            }
        });

        return dfd.promise();
    };

    getJson().then(function (matchedPosts) {
        matchedPosts.forEach(function (post) {
            $("#search-results").find("#" + post.id).show();
        });
    });
}
