

let scalex = 1;
let scaley = 1;
let tabla = [];

let s1 = function(sketch) {

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

  sketch.setup = function() {
    cnv = sketch.createCanvas(600,400);

    // cnv.position(width/10,height/10);
    cnv.parent('sketch-holder');

    cnv.changed(sketch.redraw);



    meanslider = sketch.createSlider(0, 320, 14);
    meanslider.position(cnv.position().x + sketch.width, cnv.position().y + 33);
    meanslider.parent('sketch-holder');

    sigmaslider = sketch.createSlider(0, 100, 50);
    sigmaslider.parent('sketch-holder');
    sigmaslider.position(cnv.position().x + sketch.width, cnv.position().y + 33 + 20);

    ratioslider = sketch.createSlider(0.001,0.999,0.3, 0.001);
    ratioslider.parent('sketch-holder');
    ratioslider.position(cnv.position().x + sketch.width, cnv.position().y + 33 + 20 + 20);

    nitslider = sketch.createSlider(0,10, 3);
    nitslider.parent('sketch-holder');
    nitslider.position(cnv.position().x + sketch.width, cnv.position().y + 33 + 20 + 20 + 20);



    sketch.noLoop();

    for (let i = 0; i < 5; i++) {
      tabla[i] = [];
      for (let j = 0; j < 7; j++) {
        tabla[i][j] = 0;
      }
    }

    tabla[0][1] = 14;
    tabla[3][1] = 14;
    tabla[0][2] = 9;
    tabla[3][2] = 9;

    // console.log(probabilidad(1,0,1));
  }

  sketch.draw = function() {

    for (let i = 0; i < 5; i++) {
      let s = "#p";
      if (sketch.select(s + (i+1)) != null) tabla[i][0] = sketch.float(sketch.select(s + (i+1)).value());

      s = "#mu";
      if (sketch.select(s + (i+1)) != null) tabla[i][1] = sketch.float(sketch.select(s + (i+1)).value());

      s = "#sigma";
      if (sketch.select(s + (i+1)) != null) tabla[i][2] = sketch.float(sketch.select(s + (i+1)).value());
    }


    emsrb();
    proteger();


    for (let i = 0; i < 5; i++) {
      for (let j = 0; j < 7; j++) {
        let s = "#A";

        if (sketch.select(s + (i+1) + (j+1)) != null) sketch.select(s + (i+1) + (j+1)).html(sketch.round(tabla[i][j]));
      }
    }


    meanslider.position(cnv.position().x + sketch.width, cnv.position().y + 33);

    sigmaslider.position(cnv.position().x + sketch.width, cnv.position().y + 33 + 20);

    ratioslider.position(cnv.position().x + sketch.width, cnv.position().y + 33 + 20 + 20);

    nitslider.position(cnv.position().x + sketch.width, cnv.position().y + 33 + 20 + 20 + 20);


    // update values

    mean = meanslider.value();
    sigma = sigmaslider.value();
    ratio = ratioslider.value();
    nit = nitslider.value();

    spectrum = sketch.constrain(sigma, 20,100)/3;

    sketch.background(255);
    sketch.push();
    sketch.translate(sketch.width/6, 3*sketch.height/4);

    sketch.stroke(0);
    sketch.strokeWeight(1);
    sketch.line(-sketch.width/2, 0, sketch.width, 0);
    sketch.line(0, -sketch.height/2, 0, sketch.height/2);

    sketch.strokeWeight(1.5);
    sketch.line(-10, -scl, +10, -scl);

    sketch.line(mean, -5, mean, 5);


    sketch.strokeWeight(3);
    sketch.stroke(0, 120, 240);

    sketch.noFill();
    sketch.beginShape();

    let fx;

    for (let i = mean - spectrum*sigma; i < mean + spectrum*sigma; i += increment) {

      if (sigma != 0) {fx = sigmoid(i)};
      if (sigma == 0) {fx = 0};


      sketch.vertex(i, -scl*fx);
    }

    sketch.endShape();


    // d/dx sigmoid

    sketch.stroke(100, 240, 150);

    sketch.beginShape();


    for (let i = mean - spectrum*sigma; i < mean + spectrum*sigma; i += increment) {

      if (sigma != 0) {fx = -dsigmoid(i)*0.1};
      if (sigma == 0) {fx = 0};


      sketch.vertex(i, -scl*scl*fx);
    }

    sketch.endShape();


    let xo = mean;
    let s = 0;

    sketch.beginShape();

    for (let j = 0; j < nit; j++) {
      s = (ratio - sigmoid(xo));

      sketch.vertex(xo, -scl*sigmoid(xo));
      sketch.push();
      sketch.stroke(230, 100, 200);
      sketch.point(xo, -scl*sigmoid(xo));
      sketch.pop();
      let dx = -dsigmoid(xo);
      xo = -s/dx + xo;

      // puntopendiente(xo, -scl*sigmoid(xo), scl*dx, 100);

    }

    sketch.endShape();



    xo = mean;
    s = 0;

    for (let j = 0; j < nit; j++) {
      s = (ratio - sigmoid(xo));


      sketch.push();
      sketch.stroke(230, 100, 200);
      sketch.strokeWeight(10);
      sketch.point(xo, -scl*sigmoid(xo));
      sketch.pop();
      let dx = -dsigmoid(xo);
      xo = -s/dx + xo;

      // puntopendiente(xo, -scl*sigmoid(xo), scl*dx, 100);

    }


    // puntopendiente(14, -0.5*scl, -dsigmoid(14)*scl, 100);



    sketch.strokeWeight(2);
    sketch.stroke(255, 0, 0);
    sketch.line(-sketch.width, -scl*ratio, sketch.width, -scl*ratio);
    sketch.line(xo, -sketch.height, xo, sketch.height);

    sketch.pop();



    sketch.text("mu = " + mean, sketch.width - 100, 50);
    sketch.text("sigma = " + sigma, sketch.width - 100, 70);
    sketch.text("ratio = " + ratio, sketch.width - 100, 90);
    sketch.text("n = " + nit, sketch.width - 100, 110);
    sketch.text("theta = " + xo, sketch.width - 180, 130);

    // stroke(0);
    // strokeWeight(1);
    // line(0,0,width,0);
    // line(0,height,width,height);
    // line(0,0,0,height);
    // line(width, 0, width, height);


  }


  function sigmoid(x) {
    return 1 / (1 + sketch.exp((x - mean)/sigma));
  }

  function dsigmoid(x) {
    // return sigmoid(x)*(1 - sigmoid(x));
    return -(sketch.exp((x - mean)/sigma))/(sigma*sketch.pow(1+sketch.exp((x - mean)/sigma),2))
  }

  function puntopendiente(x0, y0, m, r) {
    let x1, y1, x2, y2;
    x1 = x0 - r;
    x2 = x0 + r;
    y1 = m*(x1 - x0) + y0;
    y2 = m*(x2 - x0) + y0;
    sketch.line(x1, y1, x2, y2);
  }

  sketch.mouseDragged = function() {
    sketch.redraw();

  }

  sketch.mouseReleased = function() {
    sketch.redraw();
  }

  function emsrb() {
    for (let i = 0; i < 5; i++) {
      let s = 0;
      for (let k = 0; k < i+1; k++) {
        s += tabla[k][1];
      }
      tabla[i][4] = s;

      s = 0;
      for (let k = 0; k < i+1; k++) {
        s += tabla[k][0]*tabla[k][1];
      }
      tabla[i][3] = s/tabla[i][4];

      s = 0;
      for (let k = 0; k < i+1; k++) {
        s += tabla[k][2]*tabla[k][2];
      }
      tabla[i][5] = sketch.sqrt(s);

    }
  }

  function gauss(x, mu, sigma) {
    let a = 1/sigma/sketch.sqrt(sketch.TWO_PI);
    let b = sketch.exp(-0.5*sketch.pow(x/sigma - mu/sigma,2));

    return a*b;
  }

  function integralCDF(t,mu,sigma,n) {

    let h = (t-mu)/n;
    let k = 0;
    let x;

    for (let i = 0; i < n-1; i++) {
      x = mu + h*i;
      k += gauss(x,mu,sigma);
    }

    return 0.5*h*(gauss(t,mu,sigma) + gauss(mu,mu,sigma) + 2*k);
  }

  function probabilidad(t, mu, sigma) {
    let s = integralCDF(t, mu, sigma, 1000);

    return s + 0.5;
  }

  function proteger() {
    let tol = 0.000001;
    let s;
    for (let i = 0; i < 5-1; i++) {
      let x = tabla[i][4];
      for (let j = 0; j < 50; j++) {
        s = tabla[i+1][0]/tabla[i][3] - 1 + probabilidad(x, tabla[i][4], tabla[i][5]);


        if (sketch.abs(s) < tol) {
          // console.log("i: ", i, "  j: ", j, "  s: ", s);
          tabla[i][6] = x;
          break;
        }
        let dx = -gauss(x, tabla[i][4], tabla[i][5]);
        x = s/dx + x;
      }
      // tabla[i][6] = x;
    }
  }



};



new p5(s1);


let s2 = function(q) {

  let tasapax = 1.5;
  let tasavag = 500;
  let nplazas = 80;
  let d = [];
  let prob = [];
  let maxd = 0;
  let mind = 0;

  let nvagones;

  let tasapaxslider;
  let tasavagslider;
  let cnv;

  q.setup = function() {
    cnv = q.createCanvas(600,400);

    // cnv.position(width/10,height/10);
    cnv.parent('q-holder');

    cnv.changed(q.redraw);

    //
    //
    tasapaxslider = q.createSlider(0, 100, 1.5, 0.5);
    tasapaxslider.position(cnv.position().x + q.width, cnv.position().y + q.height/77);
    tasapaxslider.parent('q-holder');
    //
    tasavagslider = q.createSlider(0, 5000, 500);
    tasavagslider.parent('q-holder');
    tasavagslider.position(cnv.position().x + q.width, cnv.position().y + q.height/77 + 20);
    //

    //
    // nitslider = q.createSlider(0,10, 3);
    // nitslider.parent('q-holder');
    // nitslider.position(cnv.position().x + q.width, cnv.position().y + 33 + 20 + 20 + 20);



    q.noLoop();

    for (let i = 0; i < nplazas*4; i++) {
      d[i] = [];
      prob[i] = [];
      for (let j = 0; j < 2; j++) {
        d[i][j] = 0;
        prob[i][j] = 0;
      }
    }

    // delta();
  }

  q.draw = function() {


    mind = 0;
    maxd = 0;

    tasapaxslider.position(cnv.position().x + q.width, cnv.position().y + q.height/77);
    tasapax = tasapaxslider.value();

    tasavagslider.position(cnv.position().x + q.width, cnv.position().y + q.height/77 + 20);
    tasavag = tasavagslider.value();


    q.background(255);

    // q.stroke(0);
    // q.strokeWeight(5);
    // q.point(0,0);

    delta();
    q.stroke(0);
    q.rect(0,0,q.width, q.height);
    q.translate(q.width/10, 9*q.height/10);

    q.noFill();
    q.strokeWeight(0.8);
    q.stroke(100);
    q.line(0, -q.map(0, mind, maxd, 0, q.height - q.height/4), q.width - q.width/6, -q.map(0, mind, maxd, 0, q.height - q.height/4));

    q.stroke(200, 200, 100);
    for (let i = 0; i < 5; i++) {
      q.line(q.map(i*nplazas, 1, 320, 0, q.width - q.width/6, true), -q.height/10, q.map(i*nplazas, 1, 320, 0, q.width - q.width/6, true), -q.height + q.height/5);
    }

    q.stroke(255, 0, 0);
    q.strokeWeight(1);
    for (let i = 0; i < 5-1; i++) {
      q.line(q.map(tabla[i][6], 1, 320, 0, q.width - q.width/6, true), -q.height/10, q.map(tabla[i][6], 1, 320, 0, q.width - q.width/6, true), -q.height + q.height/5);
    }


    q.stroke(180, 50, 200);
    q.strokeWeight(2);



    q.beginShape();

    for (let i = 0; i < d.length; i++) {
      q.vertex(q.map(d[i][0], 1, 320, 0, q.width - q.width/6), -q.map(d[i][1], mind, maxd, 0, q.height - q.height/4));
    }

    q.endShape();

    q.stroke(70, 220, 150);

    q.beginShape();

    for (let i = 0; i < d.length; i++) {
      q.vertex(q.map(prob[i][0], 1, 320, 0, q.width - q.width/6), -q.map(prob[i][1], 0, 1, 0, q.height - q.height/4));
    }

    q.endShape();





    q.stroke(180);
    q.strokeWeight(1);
    q.line(0, -q.height + q.height/4, q.width - q.width/6, -q.height + q.height/4);
    q.text(q.round(maxd), -q.width/12, -q.height + q.height/4);

    q.text("Tasa por pasajero: " + tasapax, q.width - q.width/3, -q.height + q.height/7);
    q.text("Coste por vagón: " + tasavag, q.width - q.width/3, -q.height + q.height/7 + 20);

    q.text("N vagones: " + nvagones, q.width/2.5, 0 + q.height/15);

  }

  function delta() {
    let clase = 0;
    let s = 0;

    for (let i = 0; i < nplazas*4; i++) {
      if (i > tabla[clase][6]) {
        if (tabla[clase][6] != 0) clase++;
      }

      if (q.ceil((i+1)/nplazas) - q.ceil(i/nplazas) == 1) s -= tasavag;

      s += (tabla[clase][0] - tasapax)*(1 - probabilidad(i + 1, tabla[clase][4], tabla[clase][2]));

      // console.log(tabla[clase][6]);

      prob[i][0] = i+1;
      prob[i][1] = (1 - probabilidad(i + 1, tabla[clase][4], tabla[clase][2]));

      d[i][0] = i + 1;
      d[i][1] = s;

      if (s > maxd) {
        maxd = s;
        nvagones = q.ceil((i+1)/nplazas);
      }
      if (s < mind) mind = s;
      // console.log(d[i][0], d[i][1], );
    }
  }


    function gauss(x, mu, sigma) {
      let a = 1/sigma/q.sqrt(q.TWO_PI);
      let b = q.exp(-0.5*q.pow(x/sigma - mu/sigma,2));

      return a*b;
    }

    function integralCDF(t,mu,sigma,n) {

      let h = (t-mu)/n;
      let k = 0;
      let x;

      for (let i = 0; i < n-1; i++) {
        x = mu + h*i;
        k += gauss(x,mu,sigma);
      }

      return 0.5*h*(gauss(t,mu,sigma) + gauss(mu,mu,sigma) + 2*k);
    }

    function probabilidad(t, mu, sigma) {
      let s = integralCDF(t, mu, sigma, 1000);

      return s + 0.5;
    }

    q.mouseDragged = function() {
      q.redraw();
    }

    q.mouseReleased = function() {
      q.redraw();
    }

};

new p5(s2);
