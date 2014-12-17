    function changeImage() {
        c = document.getElementById("img").src.split("/").pop()
        if (c == "fit_over_time.png"){
            document.getElementById("img").src = "fit_over_time_logy.png";
        }
        else{
            if (c == "fit_over_time_logy.png"){ 
                document.getElementById("img").src = "fit_over_time_logx.png";
            }
            else {
                document.getElementById("img").src = "fit_over_time.png";
            }
        }
    }
