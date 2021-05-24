

let tabla = [];
let nplazas = 80;
let tasapax = 1.5;
let tasavag = 500;

let maxd = 0;
let mind = 0;


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
    cnv = sketch.createCanvas(600,300);

    // cnv.position(width/10,height/10);
    cnv.parent('sketch-holder');

    cnv.changed(sketch.redraw);



    meanslider = sketch.select("#is1mu");

    sigmaslider = sketch.select("#is1sigma");

    ratioslider = sketch.select("#is1frac");

    nitslider = sketch.select("#is1n");



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




    // update values

    mean = meanslider.value();
    sigma = sigmaslider.value();
    ratio = ratioslider.value();
    nit = nitslider.value();

    sketch.select("#s1mu").html(mean);
    sketch.select("#s1sigma").html(sigma);
    sketch.select("#s1frac").html(ratio);
    sketch.select("#s1n").html(nit);





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

    sketch.select("#s1theta").html(sketch.round(xo));

    //
    // sketch.text("mu = " + mean, sketch.width - 100, 50);
    // sketch.text("sigma = " + sigma, sketch.width - 100, 70);
    // sketch.text("ratio = " + ratio, sketch.width - 100, 90);
    // sketch.text("n = " + nit, sketch.width - 100, 110);
    // sketch.text("theta = " + xo, sketch.width - 180, 130);

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

  let d = [];
  let prob = [];


  let nvagones;

  let tasapaxslider;
  let tasavagslider;
  let cnv;

  q.setup = function() {
    cnv = q.createCanvas(600,300);

    // cnv.position(width/10,height/10);
    cnv.parent('q-holder');

    cnv.changed(q.redraw);

    //
    //
    tasapaxslider = q.select("#is2tp");
    //
    tasavagslider = q.select("#is2tv");
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


    tasapax = tasapaxslider.value();
    tasavag = tasavagslider.value();

    q.select("#s2tp").html(tasapax);
    q.select("#s2tv").html(tasavag);

    q.background(255);

    // q.stroke(0);
    // q.strokeWeight(5);
    // q.point(0,0);

    delta();
    q.stroke(0);
    q.noFill();
    q.rect(0,0,q.width, q.height);


    let h = 50;
    let nticks = q.int(q.width/h);
    q.strokeWeight(1);
    for (let i = 1; i < nticks; i++) {
      q.line(i*h, 0, i*h, 5);
      q.line(i*h, q.height, i*h, q.height - 5);
    }

    h = 50;
    nticks = q.int(q.height/h);
    q.strokeWeight(1);
    for (let i = 1; i < nticks; i++) {
      q.line(0, i*h, 5, i*h);
      q.line(q.width, i*h, q.width - 5, i*h);
    }

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





    q.stroke(100);
    q.strokeWeight(1);
    q.line(0, -q.height + q.height/4, q.width - q.width/6, -q.height + q.height/4);
    // q.text(q.round(maxd), -q.width/12, -q.height + q.height/4);


    q.select("#s2mx").html(q.round(maxd));
    q.select("#s2nv").html(nvagones);

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


let s3 = function(q) {

  let cnv;
  let p = [];
  let trenes = [];
  let d = [];

  let alphafade = 2;

  let centros = [];


  let prevTP = tasapax;
  let prevTV = tasavag;

  q.setup = function() {
    cnv = q.createCanvas(600, 300);

    cnv.parent('monte-carlo');

    for (let i = 0; i < 50; i++) {
      p[i] = [];
      for (let j = 0; j < 5; j++) {
        p[i][j] = 0;
      }
    }

    for (let i = 0; i < 50*4; i++) {
      trenes[i] = [];
      for (let j = 0; j < 4; j++) {
        trenes[i][j] = [];
        for (let k = 0; k < 2; k++) {
          trenes[i][j][k] = 0;
        }
      }
    }

    centros[0] = q.createVector(0,0);
    centros[1] = q.createVector(0,0);
    centros[2] = q.createVector(0,0);
    centros[3] = q.createVector(0,0);

    // q.noLoop();
    // q.frameRate(10);
  }

  q.draw = function() {
    // mind = 0;
    // maxd = 0;

    // q.background(255);

    if (tasapax != prevTP) q.background(255);
    if (tasavag != prevTV) q.background(255);

    if (q.frameCount % 1800 == 0) q.background(255);

    q.stroke(0);
    q.strokeWeight(1);
    q.noFill();
    q.rect(0,0,q.width, q.height);

    let h = 50;
    let nticks = q.int(q.width/h);
    q.strokeWeight(0.05);
    for (let i = 1; i < nticks; i++) {
      q.line(i*h, 0, i*h, 5);
      q.line(i*h, q.height, i*h, q.height - 5);
    }

    h = 50;
    nticks = q.int(q.height/h);
    q.strokeWeight(0.05);
    for (let i = 1; i < nticks; i++) {
      q.line(0, i*h, 5, i*h);
      q.line(q.width, i*h, q.width - 5, i*h);
    }

    q.push();

    q.translate(q.width/10, 9*q.height/10);

    for (let g = 0; g < 10; g++) {
    montecarlo();

    q.strokeWeight(1);
    q.stroke(20, 50, 200, alphafade);
    for (let i = 0; i < 50; i++) {
      q.point(q.map(trenes[i][0][0], 0, nplazas*4, 0, q.width - q.width/6), -q.map(trenes[i][0][1], mind, maxd, 0, q.height - q.height/4));
      centros[0].add(trenes[i][0][0], trenes[i][0][1]);
    }

    q.stroke(200, 50, 20, alphafade);
    for (let i = 0; i < 50; i++) {
      q.point(q.map(trenes[i][1][0], 0, nplazas*4, 0, q.width - q.width/6), -q.map(trenes[i][1][1], mind, maxd, 0, q.height - q.height/4));
      centros[1].add(trenes[i][1][0], trenes[i][1][1]);
    }

    q.stroke(20, 180, 100, alphafade);
    for (let i = 0; i < 50; i++) {
      q.point(q.map(trenes[i][2][0], 0, nplazas*4, 0, q.width - q.width/6), -q.map(trenes[i][2][1], mind, maxd, 0, q.height - q.height/4));
      centros[2].add(trenes[i][2][0], trenes[i][2][1]);

    }

    q.stroke(220, 150, 240, alphafade);
    for (let i = 0; i < 50; i++) {
      q.point(q.map(trenes[i][3][0], 0, nplazas*4, 0, q.width - q.width/6), -q.map(trenes[i][3][1], mind, maxd, 0, q.height - q.height/4));
      centros[3].add(trenes[i][3][0], trenes[i][3][1]);
    }
    centros[0].mult(1/50);
    centros[1].mult(1/50);
    centros[2].mult(1/50);
    centros[3].mult(1/50);
    }



    q.strokeWeight(1);
    q.stroke(20, 50, 200, alphafade*10*5);
    q.point(q.map(centros[0].x, 0, nplazas*4, 0, q.width - q.width/6), -q.map(centros[0].y, mind, maxd, 0, q.height - q.height/4));

    q.stroke(200, 50, 20, alphafade*10*5);
    q.point(q.map(centros[1].x, 0, nplazas*4, 0, q.width - q.width/6), -q.map(centros[1].y, mind, maxd, 0, q.height - q.height/4));

    q.stroke(20, 180, 100, alphafade*10*5);
    q.point(q.map(centros[2].x, 0, nplazas*4, 0, q.width - q.width/6), -q.map(centros[2].y, mind, maxd, 0, q.height - q.height/4));

    q.stroke(220, 150, 240, alphafade*10*5);
    q.point(q.map(centros[3].x, 0, nplazas*4, 0, q.width - q.width/6), -q.map(centros[3].y, mind, maxd, 0, q.height - q.height/4));

    q.stroke(200, 200, 50, alphafade);
    q.beginShape();
    for (let j = 0; j < 4; j++) {
      q.vertex(q.map(centros[j].x, 0, nplazas*4, 0, q.width - q.width/6), -q.map(centros[j].y, mind, maxd, 0, q.height - q.height/4));
    }
    q.endShape();

    q.pop();


    prevTP = tasapax;
    prevTV = tasavag;

    // if (q.frameCount > 120) q.noLoop();
  }

  function montecarlo() {
    for (let index = 0 + 1; index < 4 + 1; index++) {
      for (let i = 0; i < 50; i++) {
        for (let j = 0; j < 5; j++) {
          p[i][j] = (q.randomGaussian(tabla[j][1], tabla[j][2]));
        }
      }

      let libres;
      let t;
      for (let i = 0; i < 50; i++) {
        d[i] = 0;
        libres = index*nplazas;

        for (let j = 5-1; j >= 0; j--) {
          t = 0;
          for (let k = j + 1; k < 5; k++) {
            t += p[i][k];
          }

          if (j == 0) {
            libres = index*nplazas - t;
          } else {
            libres = index*nplazas - tabla[j-1][6] - t;
          }

          if (p[i][j] > libres) {
            p[i][j] = libres;
            if (libres < 0) p[i][j] = 0;
          }

          let factor = 1;
          if (centros[index-1].x > 0) {
            // factor = 1/(nplazas*index/centros[index-1].x);
            factor = q.constrain(centros[index-1].x/(nplazas*index), 0, 1);
          }

          d[i] += p[i][j]*(tabla[j][0] - tasapax - tasavag/(nplazas*factor));
        }
      }

      for (let i = 0; i < 50; i++) {
        t = 0;
        for (let j = 0; j < 5; j++) {
          t += p[i][j];
        }

        trenes[i][index-1][0] = t;
        trenes[i][index-1][1] = d[i];

        if (q.abs(d[i])/q.abs(maxd) > 1.2) {
          maxd = d[i];
          q.background(255);
          let h = 50;
          let nticks = q.int(q.width/h);
          q.strokeWeight(0.05);
          for (let i = 1; i < nticks; i++) {
            q.line(i*h, 0, i*h, 5);
            q.line(i*h, q.height, i*h, q.height - 5);
          }

          h = 50;
          nticks = q.int(q.height/h);
          q.strokeWeight(0.05);
          for (let i = 1; i < nticks; i++) {
            q.line(0, i*h, 5, i*h);
            q.line(q.width, i*h, q.width - 5, i*h);
          }
        }
        if (q.abs(d[i])/q.abs(mind) < 0.9) {
          mind = d[i];
          q.background(255);
          let h = 50;
          let nticks = q.int(q.width/h);
          q.strokeWeight(0.05);
          for (let i = 1; i < nticks; i++) {
            q.line(i*h, 0, i*h, 5);
            q.line(i*h, q.height, i*h, q.height - 5);
          }

          h = 50;
          nticks = q.int(q.height/h);
          q.strokeWeight(0.05);
          for (let i = 1; i < nticks; i++) {
            q.line(0, i*h, 5, i*h);
            q.line(q.width, i*h, q.width - 5, i*h);
          }
        }
      }

    }



  }


}

new p5(s3);
