
let mean = 0;
let sigma = 50;

let increment = sigma/100;
let spectrum = 50;
let scl = 150;

let meanslider;
let sigmaslider;

let ratio = 0.3;
let ratioslider;

let nit = 3;
let nitslider;
let cnv;


function setup() {
  cnv = createCanvas(600,400);

  // cnv.position(width/10,height/10);
  cnv.parent('sketch-holder');

  cnv.changed(redraw);

  meanslider = createSlider(0, 320, 14);
  meanslider.position(cnv.position().x + width, cnv.position().y + 33);
  meanslider.parent('sketch-holder');

  sigmaslider = createSlider(0, 100, 50);
  sigmaslider.parent('sketch-holder');
  sigmaslider.position(cnv.position().x + width, cnv.position().y + 33 + 20);

  ratioslider = createSlider(0,1,0.3, 0.01);
  ratioslider.parent('sketch-holder');
  ratioslider.position(cnv.position().x + width, cnv.position().y + 33 + 20 + 20);

  nitslider = createSlider(0,10, 3);
  nitslider.parent('sketch-holder');
  nitslider.position(cnv.position().x + width, cnv.position().y + 33 + 20 + 20 + 20);



  noLoop();
}

function draw() {

  meanslider.position(cnv.position().x + width, cnv.position().y + 33);

  sigmaslider.position(cnv.position().x + width, cnv.position().y + 33 + 20);

  ratioslider.position(cnv.position().x + width, cnv.position().y + 33 + 20 + 20);

  nitslider.position(cnv.position().x + width, cnv.position().y + 33 + 20 + 20 + 20);


  // update values

  mean = meanslider.value();
  sigma = sigmaslider.value();
  ratio = ratioslider.value();
  nit = nitslider.value();

  spectrum = constrain(sigma, 20,100)/3;

  background(255);
  push();
  translate(width/6, 3*height/4);

  stroke(0);
  strokeWeight(1);
  line(-width/2, 0, width, 0);
  line(0, -height/2, 0, height/2);

  strokeWeight(1.5);
  line(-10, -scl, +10, -scl);

  line(mean, -5, mean, 5);


  strokeWeight(3);
  stroke(0, 120, 240);

  noFill();
  beginShape();

  let fx;

  for (let i = mean - spectrum*sigma; i < mean + spectrum*sigma; i += increment) {

    if (sigma != 0) {fx = sigmoid(i)};
    if (sigma == 0) {fx = 0};


    vertex(i, -scl*fx);
  }

  endShape();


  // d/dx sigmoid

  stroke(100, 240, 150);

  beginShape();


  for (let i = mean - spectrum*sigma; i < mean + spectrum*sigma; i += increment) {

    if (sigma != 0) {fx = -dsigmoid(i)*0.1};
    if (sigma == 0) {fx = 0};


    vertex(i, scl*scl*fx);
  }

  endShape();


  let xo = mean;
  let s = 0;

  beginShape();

  for (let j = 0; j < nit; j++) {
    s = (ratio - sigmoid(xo));

    vertex(xo, -scl*sigmoid(xo));
    push();
    stroke(230, 100, 200);
    point(xo, -scl*sigmoid(xo));
    pop();
    let dx = -dsigmoid(xo);
    xo = -s/dx + xo;

    // puntopendiente(xo, -scl*sigmoid(xo), scl*dx, 100);

  }

  endShape();



  xo = mean;
  s = 0;

  for (let j = 0; j < nit; j++) {
    s = (ratio - sigmoid(xo));


    push();
    stroke(230, 100, 200);
    strokeWeight(10);
    point(xo, -scl*sigmoid(xo));
    pop();
    let dx = -dsigmoid(xo);
    xo = -s/dx + xo;

    // puntopendiente(xo, -scl*sigmoid(xo), scl*dx, 100);

  }


  // puntopendiente(14, -0.5*scl, -dsigmoid(14)*scl, 100);



  strokeWeight(2);
  stroke(255, 0, 0);
  line(-width, -scl*ratio, width, -scl*ratio);

  pop();



  text("mu = " + mean, width - 100, 50);
  text("sigma = " + sigma, width - 100, 70);
  text("ratio = " + ratio, width - 100, 90);
  text("n = " + nit, width - 100, 110);
  text("theta = " + xo, width - 180, 130);

  // stroke(0);
  // strokeWeight(1);
  // line(0,0,width,0);
  // line(0,height,width,height);
  // line(0,0,0,height);
  // line(width, 0, width, height);


}


function sigmoid(x) {
  return 1 / (1 + exp((x - mean)/sigma));
}

function dsigmoid(x) {
  // return sigmoid(x)*(1 - sigmoid(x));
  return -(exp((x - mean)/sigma))/(sigma*pow(1+exp((x - mean)/sigma),2))
}

function puntopendiente(x0, y0, m, r) {
  let x1, y1, x2, y2;
  x1 = x0 - r;
  x2 = x0 + r;
  y1 = m*(x1 - x0) + y0;
  y2 = m*(x2 - x0) + y0;
  line(x1, y1, x2, y2);
}


function mouseDragged() {
  redraw();
}

function mouseReleased() {
  redraw();
}
