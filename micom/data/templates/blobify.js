files = {{files | safe}};

function download_data(key) {
    var blob = new Blob([files[key]], {type: "text/csv;charset=utf-8;"}),
        filename = key + ".csv",
        url = window.URL.createObjectURL(blob);
    if (navigator.msSaveBlob) { // IE 10+
        navigator.msSaveBlob(blob, filename);
    } else {
        var link = document.createElement("a");
        if (link.download !== undefined) { // feature detection
            // Browsers that support HTML5 download attribute
            var url = URL.createObjectURL(blob);
            link.setAttribute("href", url);
            link.setAttribute("download", filename);
            link.style.visibility = 'hidden';
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
        }
    }
}
