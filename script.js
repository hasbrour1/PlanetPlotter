/*******************************************
***
*** Set Up Canvas
***
*******************************************/
var a_canvas = document.getElementById("a");
var context = a_canvas.getContext("2d");

function setUpCanvas(){
//NTS Disclamer
context.fillStyle = "black";
context.font = "25px Arial";
context.fillText("*Not to Scale",840,980);
context.font = "10px Arial";

// Draw the SUN
context.fillStyle = "yellow";
context.beginPath();
context.arc(500, 500, 40, 0, 2*Math.PI);
context.closePath();
context.fill();
context.lineWidth = 2;
context.stroke();
context.fillStyle = "black";
}


/*******************************************
***
*** Drawing Functions
***
*******************************************/
function drawPlanet(x, y, color, planetStr){
    context.beginPath();
    context.fillStyle = color;
    context.arc(x, y, 5, 0, 2*Math.PI);
    context.closePath();
    context.fill();

    //Tag
    context.fillText(planetStr, x + 10, y);
}

function drawRing(radius){
    context.beginPath();
    context.arc(500, 500, radius, 0, 2*Math.PI);
    context.closePath();
    context.lineWidth = 1;
    context.stroke();
}


/*******************************************
***
***This function is called when the submit button is pressed
***It will get the date entered then fill out the table and eventually
***the graphic will change
***Recode of meeus2.m
***
*******************************************/
function processData(){
    //clear canvas
    context.clearRect(0, 0, a_canvas.width, a_canvas.height);
    setUpCanvas();
    
    //get Data (This is the Date Value)
    //Eventually Chekc what type of date it is
    var dataType = document.getElementById("dataType").value;
    var data = document.getElementById("data").value;
    var JD = "";
    var yr;
    var month;
    var day;
    var hr;
    var min;
    var sec;
    
    //Convert Date to JD
    //Converted from projectmain.m
    switch(dataType){
        case "date":
            yr = parseInt(data.substring(6, 10));
            month = parseInt(data.substring(0, 2));
            day = parseInt(data.substring(3, 5));
            hr = 12;
            min = 0; 
            sec = 0;
            JD = 367*yr-Math.floor((7*(yr+Math.floor((month+9)/12)))/(4))+Math.floor(275*month/9)+day+1721013.5+((sec/60+min)/60+hr)/24;
            break;
        case "JD":
            JD = data;
            break;
        case "doy":
            //122.99,2019
            var dayFloat = parseFloat(data.substring(0, data.indexOf(',')));         
            yr= parseInt(data.substr(data.indexOf(',') + 1, data.length - 1)); 
            day = parseInt(dayFloat);
            //find day/month
            var monthArray;
            if(yr % 4 === 0){
                monthArray = [31 , 59 , 90 , 120 , 151 , 181 , 212 , 243 , 273 , 304 , 334 , 365];
                if(yr === 2000){
                    monthArray = [31 , 60 , 91 , 121 , 152 , 182 , 213 , 244 , 274 , 305 , 335 , 366];
                }
            }else{
                monthArray = [31 , 59 , 90 , 120 , 151 , 181 , 212 , 243 , 273 , 304 , 334 , 365];
            }
        
            if(day > monthArray[0]){
                month = 1;
            }else if(day > monthArray[1]){
                month = 2;
            }else if(day > monthArray[2]){
                month = 3;
            }else if(day > monthArray[3]){
                month = 4;
            }else if(day > monthArray[4]){
                month = 5;
            }else if(day > monthArray[5]){
                month = 6;
            }else if(day > monthArray[6]){
                month = 7;
            }else if(day > monthArray[7]){
                month = 8;
            }else if(day > monthArray[8]){
                month = 9;
            }else if(day > monthArray[9]){
                month = 10;
            }else if(day > monthArray[10]){
                month = 11;
            }else if(day > monthArray[11]){
                month = 12;
            }
            
            if(month === 1){
                day = parseInt(dayFloat);   
            }else{
                day = parseInt(Math.floor(dayFloat-monthArray[month-1]));
            }
            
            hr = (dayFloat-Math.floor(dayFloat))*24; 
            min = (hr-Math.floor(hr))*60;
            sec = (min-Math.floor(min))*60;
            JD = 367*yr-Math.floor((7*(yr+Math.floor((month+9)/12)))/(4))+Math.floor(275*month/9)+day+1721013.5+((sec/60+min)/60+hr)/24;
            break;
        case "date_time":
            yr = parseInt(data.substring(6, 10));
            month = parseInt(data.substring(0, 2));
            day = parseInt(data.substring(3, 5));
            hr = parseInt(data.substring(11, 13));
            min = parseInt(data.substring(14, 16)); 
            sec = parseInt(data.substring(17, 19));
            JD = 367*yr-Math.floor((7*(yr+Math.floor((month+9)/12)))/(4))+Math.floor(275*month/9)+day+1721013.5+((sec/60+min)/60+hr)/24;
            break;
        default:
            break;
    }  
    
    document.getElementById("enteredData").innerHTML = "You Entered " + JD + " as Julian Date";

    var T_TDB = (JD-2451545)/(36525);
    var AU = 149597870.691; //conversion AU to km

    //Set Values for Each Planet

    //Mercury- Appendix D pg 1046
    var mercury = new Planet("merc");
    mercury.a = 0.387098310*AU; //(km)
    mercury.e = 0.20563175+0.000020406*T_TDB-0.0000000284*Math.pow(T_TDB,2)-0.00000000017*Math.pow(T_TDB,3); // no unit
    mercury.i = 7.004986-0.0059516*T_TDB+0.00000081*Math.pow(T_TDB,2)+0.000000041*Math.pow(T_TDB,3); // (deg)
    mercury.omega = 48.330893-0.125422*T_TDB-0.00008833*Math.pow(T_TDB,2)-0.000000196*Math.pow(T_TDB,3); //(deg)
    mercury.wt = 77.456119+0.1588643*T_TDB-0.00001343*Math.pow(T_TDB,2)+ 0.000000039*Math.pow(T_TDB,3); //(deg)
    mercury.lamda = 252.250906+149472.6746358*T_TDB-0.00000535*Math.pow(T_TDB,2)+0.000000002*Math.pow(T_TDB,3); //(deg)
    mercury.calculateResults();
    mercury.plotData();
    // Draw Mercury
    drawRing(80);
    var mercX = 80*Math.cos(mercury.w+mercury.nu)+500;
    var mercY = 80*Math.sin(mercury.w+mercury.nu)+500;
    drawPlanet(mercX, mercY, 'gray', 'Mercury');
    
    //Venus - Appendis D pg 1046
    var venus = new Planet("ven");
    venus.a = 0.723329820*AU; //(km)
    venus.e = 0.00677188-0.000047766*T_TDB+0.0000000975*Math.pow(T_TDB,2)-0.00000000044*Math.pow(T_TDB,3); // no unit
    venus.i = 3.394662-0.0008568*T_TDB-0.00003244*Math.pow(T_TDB,2)+0.000000010*Math.pow(T_TDB,3); // (deg)
    venus.omega = 76.679920-0.2780080*T_TDB-0.00014256*Math.pow(T_TDB,2)-0.000000198*Math.pow(T_TDB,3); //(deg)
    venus.wt = 131.563707+0.0048646*T_TDB-0.00138232*Math.pow(T_TDB,2)-0.000005332*Math.pow(T_TDB,3); //(deg)
    venus.lamda = 181.979801+58517.8156760*T_TDB+0.00000165*Math.pow(T_TDB,2)-0.000000002*Math.pow(T_TDB,3); //(deg)  
    venus.calculateResults();
    venus.plotData();
    // Draw Venus
    drawRing(110);
    var venX = 110*Math.cos(venus.w+venus.nu)+500;
    var venY = 110*Math.sin(venus.w+venus.nu)+500;
    drawPlanet(venX, venY, 'green', 'Venus');
             
    //Earth - Appen. D pg 1046
    var earth = new Planet("earth");
    earth.a = 1.000001018*AU; //(km)
    earth.e = 0.01670862-0.000042037*T_TDB-0.0000001236*Math.pow(T_TDB,2)+0.00000000004*Math.pow(T_TDB,3); // no unit
    earth.i = 0.0000000+0.0130546*T_TDB-0.00000931*Math.pow(T_TDB,2)-0.000000034*Math.pow(T_TDB,3); // (deg)
    earth.omega = 174.873174-0.2410908*T_TDB+0.00004067*Math.pow(T_TDB,2)-0.000001327*Math.pow(T_TDB,3); //(deg)
    earth.wt = 102.937348+0.3225557*T_TDB+0.00015026*Math.pow(T_TDB,2)+ 0.000000478*Math.pow(T_TDB,3); //(deg)
    earth.lamda = 100.466449+35999.3728519*T_TDB-0.00000568*Math.pow(T_TDB,2)+0.000000000*Math.pow(T_TDB,3); //(deg)
    earth.calculateResults();
    earth.plotData();
    // Draw Earth
    drawRing(150);
    var earX = 150*Math.cos(earth.w+earth.nu)+500;
    var earY = 150*Math.sin(earth.w+earth.nu)+500;
    drawPlanet(earX, earY, 'blue', 'Earth');
    
    //Mars - Appen. D pg 1047
    var mars = new Planet("mars");
    mars.a = 1.523679342*AU; //(km)
    mars.e = 0.09340062+0.000090483*T_TDB-0.0000000806*Math.pow(T_TDB,2)-0.00000000035*Math.pow(T_TDB,3); //no unit
    mars.i = 1.849726-0.0081479*T_TDB-0.00002255*Math.pow(T_TDB,2)-0.000000027*Math.pow(T_TDB,3); // (deg)
    mars.omega = 49.558093-0.2949846*T_TDB-0.00063993*Math.pow(T_TDB,2)-0.000002143*Math.pow(T_TDB,3); //(deg)
    mars.wt = 336.060234+0.4438898*T_TDB-0.00017321*Math.pow(T_TDB,2)+0.000000300*Math.pow(T_TDB,3); //(deg)
    mars.lamda = 355.433275+19140.2993313*T_TDB+0.00000261*Math.pow(T_TDB,2)-0.000000003*Math.pow(T_TDB,3); //(deg)
    mars.calculateResults();
    mars.plotData();
    // Draw Mars
    drawRing(200);
    var marsX = 200*Math.cos(mars.w+mars.nu)+500;
    var marsY = 200*Math.sin(mars.w+mars.nu)+500;
    drawPlanet(marsX, marsY, 'red', 'Mars');
    
    //Jupiter - Appen D pg 1047
    var jupiter = new Planet("jup");
    jupiter.a = (5.202603191+0.0000001913*T_TDB)*AU; //(km)
    jupiter.e = 0.04849485+0.000163244*T_TDB-0.0000004719*Math.pow(T_TDB,2)-0.00000000197*Math.pow(T_TDB,3); // no unit
    jupiter.i = 1.303270-0.0019872*T_TDB+0.00003318*Math.pow(T_TDB,2)+0.000000092*Math.pow(T_TDB,3); // (deg)
    jupiter.omega = 100.464441+0.1766828*T_TDB+0.00090387*Math.pow(T_TDB,2)-0.000007032*Math.pow(T_TDB,3); //(deg)
    jupiter.wt = 14.331309+0.2155525*T_TDB+0.00072252*Math.pow(T_TDB,2)-0.000004590*Math.pow(T_TDB,3); //(deg)
    jupiter.lamda = 34.351484+3034.9056746*T_TDB-0.00008501*Math.pow(T_TDB,2)+0.000000004*Math.pow(T_TDB,3); //(deg)
    jupiter.calculateResults();
    jupiter.plotData();
    // Draw Jupiter
    drawRing(250);
    var jupX = 250*Math.cos(jupiter.w+jupiter.nu)+500;
    var jupY = 250*Math.sin(jupiter.w+jupiter.nu)+500;
    drawPlanet(jupX, jupY, 'gold', 'Jupiter');
    
    //Saturn - Appen D. pg 1047
    var saturn = new Planet("sat");
    saturn.a = (9.554909596-.0000021389*T_TDB)*AU; //(km)
    saturn.e = 0.05550862-0.000346818*T_TDB-0.0000006456*Math.pow(T_TDB,2)+0.00000000338*Math.pow(T_TDB,3); //no unit
    saturn.i = 2.488878+0.0025515*T_TDB-0.00004903*Math.pow(T_TDB,2)+0.000000018*Math.pow(T_TDB,3); // (deg)
    saturn.omega = 113.665524-0.2566649*T_TDB-0.00018345*Math.pow(T_TDB,2)+0.000000357*Math.pow(T_TDB,3); //(deg)
    saturn.wt = 93.056787+0.5665496*T_TDB+0.00052809*Math.pow(T_TDB,2)+0.000004882*Math.pow(T_TDB,3); //(deg)
    saturn.lamda = 50.077471+1222.1137943*T_TDB+0.00021004*Math.pow(T_TDB,2)-0.000000019*Math.pow(T_TDB,3); //(deg)
    saturn.calculateResults();
    saturn.plotData();
    // Draw Satrun
    drawRing(300);
    var satX = 300*Math.cos(saturn.w+saturn.nu)+500;
    var satY = 300*Math.sin(saturn.w+saturn.nu)+500;
    drawPlanet(satX, satY, 'coral', 'Saturn');
       
    //Uranus - Appen D. pg 1047-1048
    var uranus = new Planet("uran");
    uranus.a = (19.218446062-0.0000000372*T_TDB+0.000000000098*Math.pow(T_TDB,2))*AU; //(km)
    uranus.e = 0.04629590-0.000027337*T_TDB+0.0000000790*Math.pow(T_TDB,2)+0.00000000025*Math.pow(T_TDB,3); //no unit
    uranus.i = 0.773196-0.0016869*T_TDB+0.00000349*Math.pow(T_TDB,2)+0.000000016*Math.pow(T_TDB,3); // (deg)
    uranus.omega = 74.005947+0.0741461*T_TDB+0.00040540*Math.pow(T_TDB,2)+0.000000104*Math.pow(T_TDB,3); //(deg)
    uranus.wt = 173.005159+0.0893206*T_TDB-0.00009470*Math.pow(T_TDB,2)+0.000000413*Math.pow(T_TDB,3); //(deg)
    uranus.lamda = 314.055005+428.4669983*T_TDB-0.00000486*Math.pow(T_TDB,2)+0.000000006*Math.pow(T_TDB,3); //(deg)
    uranus.calculateResults();
    uranus.plotData();
    // Draw Uranus
    drawRing(350);
    var uraX = 350*Math.cos(uranus.w+uranus.nu)+500;
    var uraY = 350*Math.sin(uranus.w+uranus.nu)+500;
    drawPlanet(uraX, uraY, 'purple', 'Uranus');
    
    //Neptune - Appen D pg 1048
    var neptune = new Planet("nep");
    neptune.a = (30.110386869-0.0000001663*T_TDB+0.00000000069*Math.pow(T_TDB,2))*AU; //(km)
    neptune.e = 0.00898809+0.000006408*T_TDB-0.0000000008*Math.pow(T_TDB,2); //no unit
    neptune.i = 1.769952+0.0002257*T_TDB+0.00000023*Math.pow(T_TDB,2)-0.000000000*Math.pow(T_TDB,3); // (deg)
    neptune.omega = 131.784057-0.0061651*T_TDB-0.00000219*Math.pow(T_TDB,2)-0.000000078*Math.pow(T_TDB,3); //(deg)
    neptune.wt = 48.123691+0.0291587*T_TDB+0.00007051*Math.pow(T_TDB,2)-0.000000000*Math.pow(T_TDB,3); //(deg)
    neptune.lamda = 304.348665+218.4862002*T_TDB+0.00000059*Math.pow(T_TDB,2)-0.000000002*Math.pow(T_TDB,3); //(deg   
    neptune.calculateResults();
    neptune.plotData();
    // Draw Neptune
    drawRing(400);
    var nepX = 400*Math.cos(neptune.w+neptune.nu)+500;
    var nepY = 400*Math.sin(neptune.w+neptune.nu)+500;
    drawPlanet(nepX, nepY, 'lightblue', 'Neptune');
}

/*******************************************
***
***Planet Class
***Holds the data of the planets
***
*******************************************/
function Planet(code){
    this.code = code;
    this.a = "";
    this.e = "";
    this.i = "";
    this.omega = "";
    this.lamda = "";
    this.wt = "";
    this.r = "";
    this.v = "";
    this.m = "";
    this.w = "";
    this.nu = "";
    this.mu = "";
    this.p = "";
    
    //UNITS STRINGS
    var aUnit = "km";
    var eUnit = "";
    var iUnit = "&deg";
    var omegaUnit = "&deg";
    var lamdaUnit = "&deg";
    var wtUnit = "&deg";
    var rUnit = "km";
    var vUnit = "km/s";
    
    this.plotData = function(){
        document.getElementById(code +"A").innerHTML = this.a + aUnit;
        document.getElementById(code + "E").innerHTML = this.e + eUnit;
        document.getElementById(code + "I").innerHTML = this.i + iUnit;
        document.getElementById(code + "Omega").innerHTML = this.omega + omegaUnit;
        document.getElementById(code + "Lamda").innerHTML = this.lamda + lamdaUnit;
        document.getElementById(code + "WT").innerHTML = this.wt + wtUnit;
        document.getElementById(code + "R").innerHTML = this.r + rUnit;
        document.getElementById(code + "V").innerHTML = this.v + vUnit;
    }
    
   
    //This will bound the angles and calculate the M, W, R and V
    this.calculateResults = function(){
       
        this.lamda = boundAngle(this.lamda);
 
        //Get M
        this.m = this.lamda-this.wt; //(deg)
        this.m = boundAngle(this.m);

        //Get W
        this.w = this.wt-this.omega; //(deg)
        this.w = boundAngle(this.w);
       
        //Calculate NU
        this.nu = orbitAnomaly(this.a, this.e, this.m);
       
        //calculate r and v
        this.mu = 1.32712428e11;
        this.p = this.a*(1-Math.pow(this.e, 2));  //fix this
        var rad_nu = deg2rad(this.nu);
        var rad_i = deg2rad(this.i);
        var rad_omega = deg2rad(this.omega);
        var rad_w = deg2rad(this.w);
   
        //Arrays for computing rpqw and vpqw
        var rpqw = [];
        rpqw[0] = this.p*Math.cos(rad_nu)/(1+this.e*Math.cos(rad_nu));
        rpqw[1] = this.p*Math.sin(rad_nu)/(1+this.e*Math.cos(rad_nu)) ;
        rpqw[2] = 0;
    
        var vpqw = [];
        vpqw[0] = -Math.sqrt(this.mu/this.p)*Math.sin(rad_nu);
        vpqw[1] = Math.sqrt(this.mu/this.p)*(this.e+Math.cos(rad_nu));
        vpqw[2] = 0;
  
        var ijk_pqw = ijkCalc(rad_omega, rad_w, rad_i);
        
        var rMatrix = multiplyArrays(ijk_pqw, rpqw);
        var vMatrix = multiplyArrays(ijk_pqw, vpqw);

        this.r = "(" + rMatrix[0] + " , "  + rMatrix[1] + " , " + rMatrix[2] + ")";
        this.v = "(" + vMatrix[0] + " , "  + vMatrix[1] + " , " + vMatrix[2] + ")";

    }
}

//Function for Orbit Anomaly with M
function orbitAnomaly(a, e, m){
    var mu = 1.32712428e11;
//    var n = Math.sqrt(mu/Math.pow(a, 3)); don't need anymore
    var En;    

    var flag = true;
    var capM = deg2rad(m);
    
    //calculate E
    if((-Math.PI < capM && capM<0)|| (capM > Math.PI)){
        En = capM - e;
    }else{
        En = capM + e;
    }
    
    var En1;
    var count = 0;
    while(flag){
        En1 = En + (capM-En+e*Math.sin(En))/(1-e*Math.cos(En));
        count++;
        if(Math.abs(En1 - En) < 0.0001){
            flag = false;
        }
        if(count > 10000){
            break;
        }
        En = En1;
    }
    
    
    var capE = En1; 
    
    var tempM = deg2rad(m);
    if(tempM < 0){
        tempM = tempM+2*Math.PI;
    }
    
    //calculate nu
    var nu = 2 * Math.atan2(Math.sqrt((1+e)/(1-e))*Math.tan(capE/2),1);
    if(nu < 0){
        nu = nu+2*Math.PI;
    }
    nu = rad2deg(nu);
    return nu; 
}
    

/*******************************************
***
*** Math Functions
***
***
*******************************************/
//Function for Bounding 
function boundAngle(angle){
    var aa = Math.floor(angle/360);
    var bb = aa*360;
    angle = angle-bb;
    return angle;
}

function deg2rad(angle){
    return (angle/180) *Math.PI;
}
    
function rad2deg(angle){
    return angle * 57.29577951308232; // angle / Math.PI * 180
}

//Convert to rijk and vijk
function ijkCalc(Omega, w, i){
    var matrix = Create2DArray(3);
    
    matrix[0][0] = Math.cos(Omega)*Math.cos(w)-Math.sin(Omega)*Math.sin(w)*Math.cos(i);
    matrix[0][1] = -Math.cos(Omega)*Math.sin(w)-Math.sin(Omega)*Math.cos(w)*Math.cos(i);       
    matrix[0][2] =  Math.sin(Omega)*Math.sin(i);
    
    matrix[1][0] = Math.sin(Omega)*Math.cos(w)+Math.cos(Omega)*Math.sin(w)*Math.cos(i);
    matrix[1][1] = -Math.sin(Omega)*Math.sin(w)+Math.cos(Omega)*Math.cos(w)*Math.cos(i);
    matrix[1][2] = -Math.cos(Omega)*Math.sin(i);

    matrix[2][0] = Math.sin(w)*Math.sin(i);
    matrix[2][1] = Math.cos(w)*Math.sin(i);
    matrix[2][2] = Math.cos(i);
    
    return matrix;
}

//a is matrix(3,3) b is batrix(3,1)
function multiplyArrays(a, b){
    var temp = [];
    for(var i = 0; i < 3; i++){
        temp[i] = a[i][0] * b[0] + a[i][1] * b[1] + a[i][2] * b[2];
    }
    
    return temp;
}

function Create2DArray(rows) {
  var arr = [];

  for (var i=0;i<rows;i++) {
     arr[i] = [];
  }

  return arr;
}














