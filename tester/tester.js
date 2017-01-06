/*******************************************
***
*** Set Up Test
***
*******************************************/
//Convert to rijk and vijk

var Omega = 48.30970351777715;
var i = 7.0039806469392785;
var w = 77.48136297033933-Omega;

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

var phone = Create2DArray(3);

phone[0][0] = 30;
phone[0][1] = 10;
phone[1][0] = 66;

var day = 201.99;


document.getElementById("omega").innerHTML = Omega;
document.getElementById("i").innerHTML = i;
document.getElementById("w").innerHTML = w;
document.getElementById("matrix1").innerHTML = matrix[0][0];
document.getElementById("date").innerHTML = matrix[0][0];

function Create2DArray(rows) {
  var arr = [];

  for (var i=0;i<rows;i++) {
     arr[i] = [];
  }

  return arr;
}
