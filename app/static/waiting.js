(function(){
    var sleep = 1000;
    var paths = location.href.split('/');
    var ID = paths[paths.length - 1];

    function watch() {
        axios.get("/status/" + ID).then(function (response) {
        console.log('status ' + response.data.status)
            if (response.data.status === 1) {
                location.href = "/result/" + ID;
            } else {
                setTimeout(watch, sleep);
            }
        }).catch(function (err) {
            location.href = "/";
        });
    }

    watch();
}());
