////////////////////////////////////////////////////////////////////////
//            JS-CODE FOR 1D ACTIVE FUNCTION LEARNING                 //
//                AUTHORS: CHARLEY WU, ERIC SCHULZ                    //
////////////////////////////////////////////////////////////////////////

//EXPERIMENT PARAMETERS
var fullurl=document.location.href, //url of incoming MTurk worker
  trials=16, //number of trials
  trialCounter=0, //counter for current trial number
  tracker=new Array(0), //tracker array
  investigationIndex=0, //current click number
  scoretotal=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
  scorecurrent=0,
  reward = 0.00,
  gridMax = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
  envOrder = getRandomSubarray([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39], trials),
  currentEnv = envOrder[0],
  exampleNum = 0;
  scale = [randomNum(60,80),randomNum(60,80),randomNum(60,80),randomNum(60,80), randomNum(60,80),randomNum(60,80),randomNum(60,80),randomNum(60,80),randomNum(60,80),randomNum(60,80),randomNum(60,80),randomNum(60,80), randomNum(60,80),randomNum(60,80),randomNum(60,80),randomNum(60,80)], // Randomly scale max value between 60-80
  //Color parameters for heatmap
  colors = ['#fff7ec','#fee8c8','#fdd49e','#fdbb84','#fc8d59','#ef6548','#d7301f','#b30000','#7f0000'],
  heatmapColor = d3.scale.linear().domain(d3.range(0, 90, 90.0/ (colors.length - 1))).range(colors);

//data collectors for search history
var xcollect=[],
  ycollect=[],
  ycollectScaled = [],
  initcollect = [];
//Populate data collectors with empty arrays
//Todo: change max clicks to the assigned horizon for each trial
for(var i=0; i<trials; i++) {
    xcollect[i] = Array.apply(null, Array(0)).map(Number.prototype.valueOf,-99);
    ycollect[i] = Array.apply(null, Array(0)).map(Number.prototype.valueOf,-99);
    ycollectScaled[i] = Array.apply(null, Array(0)).map(Number.prototype.valueOf,-99);
    initcollect[i] = Array.apply(null, Array(0)).map(Number.prototype.valueOf,-99);
}

//Declare variables not yet assigned 
var scenarioId,
  condition,
  kernel,
  horizonOrder,
  searchHistory,
  functionList,
  trialReward,
  xout,
  yout,
  optimaGuess;

// Retrieve assignmentID and workerID from URL
var assignmentID = turkGetParam('assignmentId');
var workerID = turkGetParam('workerId');
document.getElementById('MTurkID').value = workerID; //prepopulate MTurk

//Access the MySQL database and returns scenario id with a condition numnber, adding the now() date time to the start time field and marking the specific scenario as completed
function assignScenario() {
  //DEBUG: Database functionality removed for open source release
  /*
  var ajaxRequest = new XMLHttpRequest();
  try{
      // Opera 8.0+, Firefox, Safari
      ajaxRequest = new XMLHttpRequest();
  } catch (e){
      // Internet Explorer Browsers
      try{
          ajaxRequest = new ActiveXObject("Msxml2.XMLHTTP");
      } catch (e) {
          try{
              ajaxRequest = new ActiveXObject("Microsoft.XMLHTTP");
          } catch (e){
              // Something went wrong
              alert("Your browser broke!");
              return false;
          }
      }
  }
  var queryString = "?action=" + 'assignScenario';
  ajaxRequest.open("GET", "databasecall.php"+queryString, false); //DEBUG: Database functionality removed for open source release
  ajaxRequest.send(null);
  var response = ajaxRequest.responseText;
  var jsonArray = JSON.parse(response);
  scenarioId = parseInt(jsonArray['scenarioId']);
  condition = parseInt(jsonArray['condition']);
  kernel = parseInt(jsonArray['kernel']);
  horizon = parseInt(jsonArray['horizon']);
  */

  //DEBUG: Set psuedo-randomly assigned experiment variables randomly in lieu of database assignment
  scenarioId = Math.floor(Math.random()*10000); // random number between 1 and 10000 (since it is suppose to be unique)
  condition = Math.floor(Math.random());
  kernel = Math.floor(Math.random());
  horizon = Math.floor(Math.random());
  
  //Set horizon variable
  if (horizon===0){
    horizonOrder = [5,10,5,10,5,10,5,10,5,10,5,10,5,10,5,10];
  }else if (horizon===1){
    horizonOrder = [10,5,10,5,10,5,10,5,10,5,10,5,10,5,10,5];
  }
  clicks=horizonOrder[0];//set initial number of clicks to the first entry in the horizonOrder array
  change("remaining2",  "Number of clicks remaining <b>" + clicks +"</b>");
  //Load  kernel
  var kernelFiles = ['kernel1.json', 'kernel2.json'];
  loadJSON(kernelFiles[kernel], function(response){
    functionList = JSON.parse(response);
  });
  //Switch text of experiment instructions based on condition
  if(condition===0){
    document.getElementById("instructions1").style.display = "block"; //culmulative rewards
  } else if (condition===1){
    document.getElementById("instructions2").style.display = "block"; //global optimization
  }
   //Update experiment header
  if(condition===0){
    document.getElementById("expHeader1").style.display = "block"; //culmulative rewards
  } else if (condition===1){
    document.getElementById("expHeader2").style.display = "block"; //global optimization
  }
  //update short instructions
  var instructions = document.getElementById('explainexshow').innerHTML;
  if (condition==0){
    instructions += "<p><b>VI.</b> Your reward will be based on the total points you earn, by revealing new tiles and also by reclicking previously revealed tiles.</p>";
    change('explainexshow', instructions);
  }else if (condition==1){
    instructions += "<p><b>VI.</b> Your reward will be based the largest point value that is revealed in each grid. </p>";
    change('explainexshow', instructions);
  }
}

function showExample(){
  //first click
  if (exampleNum==0){
    //provide examples of kernel functions based on the assigned kernel
    //sample 4 environments not shown to participants in the task
    var envs = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39],
       exampleEnvs = getRandomSubarray(envs.filter(function(a){return envOrder.indexOf(a)<0;}),4);
    //sample 4 new random scaling factors
    var exScale = [randomNum(60,80),randomNum(60,80),randomNum(60,80),randomNum(60,80)];
    for (i=0; i<4; i++){
      createExample(exampleEnvs[i], exScale[i], i+1);
    }
    document.getElementById('example1').style.display = "block";
    change('exampleButton', 'Show Next Example');
  }else if (exampleNum == 1 | exampleNum == 2){ //all examples have been shown
    var prev = 'example' + exampleNum,
      next = 'example' + (exampleNum+1); 
    clickStart(prev, next);
  }else if (exampleNum == 3){
    //last example
    change('exampleButton', 'Start Task');
    var prev = 'example' + exampleNum,
      next = 'example' + (exampleNum+1); 
    clickStart(prev, next);
  }else{
    //start task
    clickStart('page2', 'page3')
  }
  exampleNum = exampleNum + 1; //increment counter
}

//intialize Gaussian Process function
//'gp' as the input is the json representing a single search environment 
function intialize(gp) {
    var t=0
    for (k = 0; k <= 29; k++){
        if(gp[k].y >= 0 && t<1){
            var xout=[gp[k].x];
            var yout=[gp[k].y * 100];
            t=t+1;
        } else if(gp[k].y>=0 && t>=1){
            xout=xout.concat([gp[k].x]);
            yout=yout.concat([gp[k].y] * 100);
            t=t+1;
        }
    }
    var choose =Math.floor(Math.random() * xout.length);
 //return x-y coordinate as well as function output
 return([xout[choose], yout[choose]]);    
}

//Checkers:
//initialize GP functionList
var init=[];
function instructioncheck()
{
    //check if correct answers are provided
    //Q1 is condition specific
    if (condition===0){
      if (document.getElementById('icheck1a').checked) {var ch1=1}
    }else if (condition===1){
      if (document.getElementById('icheck1b').checked) {var ch1=1}
    }
    //none condition specific questions
    if (document.getElementById('icheck2').checked) {var  ch2 = 1}
    if (document.getElementById('icheck3').checked) {var  ch3 = 1}
    //are all of the correct
    var checksum=ch1+ch2+ch3;
    if (checksum===3){
      //DEBUG: Remove database functionality for opensource version
      /*
      //Save instructions Complete time in Database
      var ajaxRequest = new XMLHttpRequest(); //intiate AJAX request
      try{
          // Opera 8.0+, Firefox, Safari
          ajaxRequest = new XMLHttpRequest();
      } catch (e){
        // Internet Explorer Browsers
        try{
            ajaxRequest = new ActiveXObject("Msxml2.XMLHTTP");
        } catch (e) {
            try{
                ajaxRequest = new ActiveXObject("Microsoft.XMLHTTP");
            } catch (e){
                // Something went wrong
                alert("Your browser broke!");
                return false;
            }
        }}
      var queryString = "?action=" + 'completeInstructions'+ '&scenarioId=' + scenarioId;
      ajaxRequest.open("GET", "databasecall.php"+queryString, false);
      ajaxRequest.send(null);
      */
      //start next page  
      clickStart('page3','page4');
      //BUILD SEARCH SPACE
      //initialize very first revealed point
      init = intialize(functionList[envOrder[0]]);
      var noiseyValue = Math.round(init[1]+ myNorm()), //add noise to y values
        rescaledValue = Math.max(Math.round(noiseyValue/100*scale[trialCounter] + 5),0); //rescale value to randomly generated ceiling
      $('#'+init[0]).html(rescaledValue);
      $('#'+init[0]).css('text', rescaledValue); //display rescaled value
      $('#'+init[0]).css('color', 'black');
      $('#'+init[0]).css('background-color', heatmapColor(rescaledValue)); // change background color based on rescaled value
      //add to history
      $('#'+init[0]).attr('title', rescaledValue.toString());
      //store initial values
      xcollect[trialCounter][0]=parseInt(init[0],10);
      ycollect[trialCounter][0]=noiseyValue;
      ycollectScaled[trialCounter][0]=rescaledValue;
      //calculate points for both payoff conditions
      gridMax[trialCounter]= rescaledValue;
      scoretotal[trialCounter] = rescaledValue;
      //Update text for assigned condition
      if (condition===0){
        change('scoretotal', "Current Score: " + scoretotal[trialCounter]);
      }else if (condition===1){
        change('scoretotal', "Largest Reward Found: " + gridMax[trialCounter]);
      }
      //update number of clicks
      change("remaining2", "Number of clicks left: <b>"+horizonOrder[0] + "</b>");
    } else{
        //if one or more answers are wrong, raise alert box
        alert('You have answered some of the questions wrong. Please try again.');
        clickStart('page3', 'page2');
    }
}


//DO THE TASK
//get ready
$(document).ready(function(){
  //intialize grid
  var i, j, gridHTML='', WIDTH=30, HEIGHT=1, lastRevealedCellId;
  //creating grid HTML
  gridHTML += '<tr>';
  for (i = 0; i < WIDTH; i++){
    gridHTML += '<td align="center" id="' + i + '">&nbsp;</td>';
  }
  gridHTML += '</tr>';
  //append grid HTML at grid div
  $('#grid tbody').append(gridHTML);
  //on click of a particular cell, do the following   
  $('#grid * td').click(//cell has been clicked
    function() {
      var $this=$(this), id=$this.attr('id'), x, absoluteValue, noiseyValue, rescaledValue;
      // x is the coordinates of the cell in the grid, starting from left to right
      x=id;
      currentEnvNum = envOrder[trialCounter];
      //find the corresponding function output
      for (k = 0; k <= 29; k++){//search through all possible locations
        if (functionList[currentEnvNum][k].x==x){
          absoluteValue = functionList[currentEnvNum][k].y * 100; //noise free observation
          noiseyValue = Math.round(absoluteValue + myNorm()); //add normally distributed noise
          rescaledValue = Math.max(Math.round(noiseyValue/100*scale[trialCounter] + 5),0); // scale to randomly selected ceiling, add a constant to avoid negatives, and truncate distribution at 0
          }
        }
      //Check if tile has been clicked before
      if (document.getElementById(id).hasAttribute("title") == true){
        //If tile has been clicked before, add new value to history and show the new value
        var history = JSON.parse("[" + document.getElementById(id).getAttribute("title") + "]");
        history.push(rescaledValue);
        document.getElementById(id).setAttribute("title", history.toString());
        $this.html(rescaledValue);
        $this.css('text', rescaledValue);
        $this.css('color', 'black');
        $this.css('background-color', heatmapColor(Math.round(average(history)))); //average of history uses the rescaled values, so need to re-normalize to range 0 - 100
      }else{
        //save rescaledValue in history
        document.getElementById(id).setAttribute("title", rescaledValue.toString());
        // Show the rescaled value
        $this.html(rescaledValue);
        $this.css('text', rescaledValue);
        $this.css('color', 'black');
        $this.css('background-color', heatmapColor(rescaledValue));
      }
      //update trackers and counters
      investigationIndex=investigationIndex+1, tracker[investigationIndex]= id, clicks=clicks-1;
      //Update maximum reward found
      if (rescaledValue>gridMax[trialCounter]){
        gridMax[trialCounter] = rescaledValue;
      }
      //For all trials except the last one
      if (trialCounter<15){ //trials 1 - 15, no global optima guess
        if (clicks>0){//if there are still investigations available
        //keep track of chosen cell
        xcollect[trialCounter][investigationIndex]=parseInt(x,10);
        ycollect[trialCounter][investigationIndex] = noiseyValue;
        ycollectScaled[trialCounter][investigationIndex]=rescaledValue; //store absolute value
        //update number of clicks left
        change("remaining2", "Number of clicks left: <b>"+clicks + "</b>");
        //update current score
        scorecurrent=Math.round(rescaledValue);
        //update current total
        scoretotal[trialCounter]=scoretotal[trialCounter]+scorecurrent;
        if (condition==0){//Change reward total if cumulative reward
          reward = rewardCum(scoretotal);
          trialReward = singleTrialReward(scoretotal[trialCounter], scale[trialCounter], condition);
          change('scoretotal', "Current Score: " + scoretotal[trialCounter]);
        }else if(condition==1){//change reward total based on maxGrid
          reward = rewardMax(gridMax);
          trialReward = singleTrialReward(gridMax[trialCounter], scale[trialCounter], condition);
          change('scoretotal', "Largest Reward Found: " + gridMax[trialCounter]);
        }          
      } else if (clicks==0){//if out of investigations
        //keep track of chosen cell
        xcollect[trialCounter][investigationIndex]=parseInt(x,10);
        ycollect[trialCounter][investigationIndex]= noiseyValue;
        ycollectScaled[trialCounter][investigationIndex]=rescaledValue; //store absolute value
        //update number of clicks left
        change("remaining2", "Number of clicks left: <b>"+clicks + "</b>");
        //update current score
        scorecurrent=Math.round(rescaledValue);
        //update current total
        scoretotal[trialCounter]=scoretotal[trialCounter]+scorecurrent;
        if (condition==0){//Change reward total if cumulative reward
          reward = rewardCum(scoretotal);
          trialReward = singleTrialReward(scoretotal[trialCounter], scale[trialCounter], condition);
          change('scoretotal', "Current Score: " + scoretotal[trialCounter]);
        }else if(condition==1){//change reward total based on maxGrid
          reward = rewardMax(gridMax);
          trialReward = singleTrialReward(gridMax[trialCounter], scale[trialCounter], condition);
          change('scoretotal', "Largest Reward Found: " + gridMax[trialCounter]);
        }
        //move to next page
        clickStart('gridDiv', 'page5finished');
        //update trials remaining
        trials = trials - 1;
        //percentMaxReward = toFixed(trialReward / (1.5/8) * 100, 0);
        change("trials", "<h1> <font size='4'>You have finished exploring this environment and gain $" + trialReward +" as a bonus. You have currently earned a total performance bonus of $" + reward + " and have "+ trials  +" environment(s) left to complete. On the next environment, you will be allowed " + horizonOrder[trialCounter + 1] + " clicks.</font> </h1>");
      }

      }else{//Last trial, ask about global optima
        if (clicks>0){//if there are still investigations available
        //keep track of chosen cell
        xcollect[trialCounter][investigationIndex]=parseInt(x,10);
        ycollect[trialCounter][investigationIndex]=noiseyValue;
        ycollectScaled[trialCounter][investigationIndex]=rescaledValue; //store absolute value
        //update number of clicks left
        change("remaining2", "Number of clicks left: <b>"+clicks + "</b>");
        //update current score
        scorecurrent=Math.round(rescaledValue);
        //update current total
        scoretotal[trialCounter]=scoretotal[trialCounter]+scorecurrent;
        if (condition==0){//Change reward total if cumulative reward
          reward = rewardCum(scoretotal);
          trialReward = singleTrialReward(scoretotal[trialCounter], scale[trialCounter], condition);
          change('scoretotal', "Current Score: " + scoretotal[trialCounter]);
        }else if(condition==1){//change reward total based on maxGrid
          reward = rewardMax(gridMax);
          trialReward = singleTrialReward(gridMax[trialCounter], scale[trialCounter], condition);
          change('scoretotal', "Largest Reward Found: " + gridMax[trialCounter]);
        }
        }else if(clicks==0){//if out of investigations
          //keep track of chosen cell
          xcollect[trialCounter][investigationIndex]=parseInt(x,10);
          ycollect[trialCounter][investigationIndex]=noiseyValue;
          ycollectScaled[trialCounter][investigationIndex]=rescaledValue; //store absolute value
          //update number of clicks left
          change("remaining2", "Number of clicks left: <b>"+clicks + "</b>");
          //update current score
          scorecurrent=Math.round(rescaledValue);
          //update current total
          scoretotal[trialCounter]=scoretotal[trialCounter]+scorecurrent;
         if (condition==0){//Change reward total if cumulative reward
          reward = rewardCum(scoretotal);
          trialReward = singleTrialReward(scoretotal[trialCounter], scale[trialCounter], condition);
          change('scoretotal', "Current Score: " + scoretotal[trialCounter]);
        }else if(condition==1){//change reward total based on maxGrid
          reward = rewardMax(gridMax);
          trialReward = singleTrialReward(gridMax[trialCounter], scale[trialCounter], condition);
          change('scoretotal', "Largest Reward Found: " + gridMax[trialCounter]);
        }
          //Ask to identify global maxima
          //percentMaxReward = toFixed(trialReward / (1.5/8) * 100, 0);
          change("Guess", "<h1> <font size='4'>You have finished exploring all the environments, and gain $" + trialReward + " for the previous environment. You have earned a total performance bonus of $" + reward + "<br>We would like to ask you to complete one more task, by clicking the tile where you believe the largest reward is located. <br>Note: This will not affect your performance bonus. </font> </h1>");
          clickStart('gridDiv', 'page5guess');
        } else {//out of investigations
          //record global optimal guess
          optimaGuess = [parseInt(x,10), noiseyValue, rescaledValue];
          //move to next page
          clickStart('page5', 'page6');
          trials = trials-1;
        }
      }
    });
});



function nexttrial(){
	//debugData();
  //proceed only if there are more trials available
  if (trials>0){
    initcollect[trialCounter]=init;//retrieve initially revealed tile from previous trial before updating trial counter
    //update trialCounter
    trialCounter=trialCounter+1;
    //get rid of shown first output
    $('#'+init[0]).html("&nbsp");
    $('#'+init[0]).css("background-color", "");
    //get rid of revealed output                
    for (h = 0; h < tracker.length; h++){
        $('#' + tracker[h]).html("&nbsp");
        $('#' + tracker[h]).css("background-color", "");
        $('#' + tracker[h]).attr("title", "")
    }
    if (trialCounter<=envOrder.length){//initialize new GP
            init = intialize(functionList[envOrder[trialCounter]]);
      }
    var noiseyValue = Math.round(init[1] + myNorm()),
      rescaledValue = Math.max(Math.round(noiseyValue/100*scale[trialCounter] + 5), 0);
    //get first point on grid for next trial
    $('#'+init[0]).html(rescaledValue);
    $('#'+init[0]).css('text', rescaledValue);
    $('#'+init[0]).css('color', 'black');
    $('#'+init[0]).css('background-color', heatmapColor(rescaledValue));
    //add to hover text
    $('#'+init[0]).attr("title", rescaledValue.toString());
    //store initial values
    xcollect[trialCounter][0]=parseInt(init[0],10);
    ycollect[trialCounter][0]=noiseyValue;
    ycollectScaled[trialCounter][0]=rescaledValue; //store noisey value
    //update gridMax with initial tile
    gridMax[trialCounter] = rescaledValue;
    scoretotal[trialCounter] = rescaledValue;
    //Update text for assigned condition
    if (condition===0){
      change('scoretotal', "Current Score: " + scoretotal[trialCounter]);
    }else if (condition===1){
      change('scoretotal', "Largest Reward Found: " + gridMax[trialCounter]);
    }
    //go back to task
    clickStart('page5finished', 'gridDiv');
    //renew investigations
    clicks=horizonOrder[trialCounter];
    //renew investigationIndex
    investigationIndex=0;
    //update current reward, number of trials and clicks
    change("remaining1", "Number of environments left: <b>"+ trials + "</b>");
    change("remaining2", "Number of clicks left: <b>" + clicks +"</b>");
  //if out of trials go to next page
  } else{
    //Check there's no funny business
    if (reward > 1.5){
      reward = 1.50;
    }
    //move to final page
    clickStart('page5', 'page6');
}}

//create Example grids
function createExample(envNum, scale, exNum){
  var i, j, gridHTML='', WIDTH=11, HEIGHT=11, env = functionList[envNum], remainder, payoff, col;
  //beginning of gridHTML
  gridHTML += "<table class='gridEx'><tbody id='exGrid"+exNum+"'>";
 //beginning of row
  gridHTML += '<tr>';
  //loop through JSON, creating grid HTML
  for (i = 0; i < 30; i++){
    payoff = env[i]['y']; //absolute value of payoff between 0 and 100
    payoff = Math.round(payoff * scale) + 5; //rescaled payoff, rounded to an int, plus a constant of 5 to avoid any negative values
    col = heatmapColor(payoff); //color code
    gridHTML += '<td align="center" bgcolor='+col+'>'+payoff+'</td>';
   }
  gridHTML += '</tr>'; //cap end of row
  //cap end of gridHTML
  gridHTML +="</tbody> </table>";
  var divName = "example" + exNum;
  change(divName, gridHTML);//write gridHTML into place
}

function debugData(){
	console.log(xcollect)
	console.log(ycollect)
}

function senddata(){
    //DEBUG remove database saving for open source version and do nothing with the data
    //Check that all forms were filled out
    var MTurkID = document.getElementById('MTurkID').value,
        age = document.getElementById('age').value,
        gender=-1,
        processDescription=document.getElementById('processDescription').value;
    //Collect data from the last demographics form
    //gender
    if (document.getElementById('gender1').checked) {gender = 0}; //Female
    if (document.getElementById('gender2').checked) {gender = 1}; //Male
    if (MTurkID.length > 0 && age.length > 0 && gender >-1){
      //Subject has filled out all the required fields
      //Compile search history together
      searchHistory = {'xcollect':xcollect,'ycollect':ycollect, 'ycollectScaled':ycollectScaled};
      /*
      //Initiate AJAX request
      var ajaxRequest = new XMLHttpRequest();
      try{
          // Opera 8.0+, Firefox, Safari
          ajaxRequest = new XMLHttpRequest();
      } catch (e){
        // Internet Explorer Browsers
        try{
            ajaxRequest = new ActiveXObject("Msxml2.XMLHTTP");
        } catch (e) {
            try{
                ajaxRequest = new ActiveXObject("Microsoft.XMLHTTP");
            } catch (e){
                // Something went wrong
                alert("Your browser broke!");
                return false;
            }
        }}
      var queryString = "?action=" + 'completeScenario' + '&MTurkID=' + MTurkID + '&assignmentID=' + assignmentID + '&scale=' + JSON.stringify(scale) + '&envOrder=' + JSON.stringify(envOrder) + '&searchHistory=' + JSON.stringify(searchHistory) + '&reward=' + reward + '&age=' + age + '&gender=' + gender + '&processDescription=' +processDescription + '&optimaGuess=' + optimaGuess + '&scenarioId=' + scenarioId;
      ajaxRequest.open("GET", "databasecall.php"+queryString, false);
      ajaxRequest.send(null);
      */
      //display reward on final page
      var sumReward = toFixed(parseFloat(reward) + .5, 2); //Add participation cost
      var rewardText = "You earned a performance bonus of: $" + reward + ", which will be automatically assigned to your MTurk account in 1-3 days. Together with the base pay of $0.50 for this HIT, you earned $" + sumReward + " for this task. <br> <br><br>";
      change("bonus", rewardText);
      clickStart('page6','page7');
    }else{
      alert("Please enter your MTurk ID, age, and gender to complete the task");
    }
}

//*************UTILITIES***************************************


//changes from one page to another
function clickStart(hide, show)
{
        document.getElementById(hide).style.display="none";
        document.getElementById(show).style.display = "block";
        window.scrollTo(0,0);
}

//changes inner HTML of div with ID=x to y
function change (x,y){
    document.getElementById(x).innerHTML=y;
}

//Function to randomly shuffle an array:
function shuffle(o){ //v1.0
    for(var j, x, i = o.length; i; j = Math.floor(Math.random() * i), x = o[--i], o[i] = o[j], o[j] = x);
    return o;
};

//Randomly sample n values from an array
function getRandomSubarray(arr, size) {
    var shuffled = arr.slice(0), i = arr.length, temp, index;
    while (i--) {
        index = Math.floor((i + 1) * Math.random());
        temp = shuffled[index];
        shuffled[index] = shuffled[i];
        shuffled[i] = temp;
    }
    return shuffled.slice(0, size);
}

//load JSON file
function loadJSON(file, callback) {
    var rawFile = new XMLHttpRequest();
    rawFile.overrideMimeType("application/json");
    rawFile.open("GET", file, true);
    rawFile.onreadystatechange = function() {
        if (rawFile.readyState === 4 && rawFile.status == "200") {
            callback(rawFile.responseText);
        }
    }
    rawFile.send(null);
}

//Create normal noise distribution
function myNorm() {
   var x, x2, rad, c;
    do {
       x1 = 2 * Math.random() - 1;
       x2 = 2 * Math.random() - 1;
       rad = x1 * x1 + x2 * x2;
   } while(rad >= 1 || rad == 0);
    c = Math.sqrt(-2 * Math.log(rad) / rad);
    return (x1 * c);
};

//average the values in an array
function average(inputArray){
  var total = 0
  for (var i=0; i < inputArray.length; i++){
    total += inputArray[i];
  }
  var avg = total / inputArray.length;
  return avg;
};

//Convert cumulative score to reward value
function rewardCum(scoreTotal){
  var r = 0, r_i;
  for (var i=0; i<scoreTotal.length; i++){
    r_i = scoreTotal[i] / (scale[i]+5) / 120 * 1.5;
    r = r + r_i
  } 
  if (r>1.5){
    r = 1.5; //limit to max reward, in case of any funny business
  }
  return toFixed(r, 2);
}

//convert max values of each grid to reward value
function rewardMax(gridMax){
  var r=0, r_i;
  for(var i=0; i<gridMax.length; i++){
    r_i = Math.pow(gridMax[i]/(scale[i]+5),4)/16*1.5;
    if (r_i>1.5/16){
      r_i = 1.5/16; //cap reward if r_i is larger than 1/8th of the max reward
    }
    r = r + r_i;
  }
  if (r>1.5){
    r = 1.5; //limit to max reward, in case of any funny business
  }
  return toFixed(r, 2);
}

//single trial reward
function singleTrialReward(points, scale, condition){
  var r=0;
  if (condition==0){// if cumulative reward
    r = points/(scale+5)/120 * 1.5;
  }else if (condition==1){
    r = Math.pow(points/(scale+5),4)/16*1.5;
    if (r>1.5/16){
      r = 1.5/16;
    }
  }
  return toFixed(r, 2);
}

//random number generator
function randomNum(min, max){
  return Math.floor(Math.random() * (max-min+1)+min)
}

//Display a float to a fixed percision
function toFixed(value, precision) {
    var precision = precision || 0,
        power = Math.pow(10, precision),
        absValue = Math.abs(Math.round(value * power)),
        result = (value < 0 ? '-' : '') + String(Math.floor(absValue / power));

    if (precision > 0) {
        var fraction = String(absValue % power),
            padding = new Array(Math.max(precision - fraction.length, 0) + 1).join('0');
        result += '.' + padding + fraction;
    }
    return result;
}

// extract URL parameters (FROM: https://s3.amazonaws.com/mturk-public/externalHIT_v1.js)
function turkGetParam( name ) { 
    var regexS = "[\?&]"+name+"=([^&#]*)"; 
    var regex = new RegExp( regexS ); 
    var tmpURL = fullurl; 
    var results = regex.exec( tmpURL ); 
    if( results == null ) { 
      return ""; 
    } else { 
      return results[1];    
    } 
}

//END