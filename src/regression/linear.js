/*
    Linear regressions
*/

const Decimal = require("decimal.js");

const tdist = require("../distributions/studentt.js");
const fdist = require("../distributions/f.js");
const chidist = require("../distributions/chi.js");

const M = require("../tools/matrix.js");
const Arr = require("../tools/array.js");

// Dependencies, will remove in later releases
const math = require("mathjs");
const cho = require("cholesky");

/*
    @desc Adds a column of 1's to the dataset

    @param {array} exog: object to prepend to
    @param {booleen} prepend: column should be added to the front of dataset

    @return {array} returns the exogenous array with the constant object added
*/
function addConstant(exog, prepend = true) {
  // only values, not references
  exog = [...exog];

  if (prepend) {
    exog.unshift({
      title: "intercept",
      data: new Array(exog[0]["data"].length).fill(1),
      constant: true,
    });
  } else {
    exog.push({
      title: "intercept",
      data: new Array(exog[0]["data"].length).fill(1),
      constant: true,
    });
  }

  return exog;
}

/*
    @desc Base model for all linear regressions

    @param {object} endog: endogenous variables
    @param {array} exog: exogenous variables
    @param {array} sigma: array
    @param {array} cholsigmainv: array
*/
class Model {
  constructor(endog, exog, { sigma, cholsigmainv, weights } = {}) {
    // common storage as object
    this._data = {};
    this._data.endog = endog;
    this._data.exog = exog;
    this._data.df = {};
    this.alpha = 0.05;
    this.sigma = sigma;
    this.cholsigmainv = cholsigmainv;
    this._data.weights = weights;

    // array length check
    this._data.exog.forEach((a) => {
      if (a.data.length !== this._data.endog.data.length) {
        throw Error("Input lengths are not equal");
      }
    });

    let X = new Array(this._data.endog.data.length).fill(0);
    this._data.X = X.map((v, i) => {
      let arr = [];
      this._data.exog.forEach((a) => arr.push(a.data[i]));
      return arr;
    });

    // calculate whitened values
    this._data.wendog = this.whiten(this._data.endog.data);
    this._data.wexog = this.whiten(this._data.X);

    // length of data series
    this._data.N = this._data.endog.data.length;

    // # of non constants
    this._data.K = this._data.exog.filter((v) => v.constant !== true).length;

    // # of constants
    this._data.kConstant = this._data.exog.filter(
      (v) => v.constant === true
    ).length;

    // degrees of freedom
    this._data.df.regression = this._data.K;
    this._data.df.residual = this._data.N - this._data.K - 1;
    this._data.df.total = this._data.N - 1;
  }

  /*
        @desc Residual values

        @return {array} residuals
    */
  get resid() {
    let resids = this._data.residuals;
    return resids;
  }

  /*
        @desc The beta coefficients of the fitted model

        @return {array} beta coefficients
    */
  get params() {
    return this._data.exog.map((v) => v.beta);
  }

  /*
        @desc The standard errors of the regression coefficients

        @return {array} SEs
    */
  get bse() {
    return this._data.exog.map((v) => v.sterror);
  }

  /*
        @desc The T-statistic of the regression coefficients

        @return {array} coefficient t-statistics
    */
  get btstats() {
    return this._data.exog.map((v) => v.tstat);
  }

  /*
        @desc The p-values of the regression coefficients

        @return {array} coefficient P-values
    */
  get pvalues() {
    return this._data.exog.map((v) => v.pvalue);
  }

  /*
        @desc Confidence Interval for the regression coefficients

        @return {object} CI Confidence Intervals
        @return {object} CI.lower lower 95% CI value
        @return {object} CI.higher higher 95% CI value
    */
  get confInt() {
    return this._data.exog.map((v) => {
      let obj = {
        lower: v.ciLower,
        higher: v.ciHigher,
      };

      return obj;
    });
  }

  /*
        @desc F-statistic of the model

        @return {number} F-statistic
    */
  get fvalue() {
    let mseModel = this._data.mseModel ? this._data.mseModel : this.mseModel;
    let mseRes = this._data.mseResidual
      ? this._data.mseResidual
      : this.mseResidual;

    let fvalue = Decimal(mseModel).dividedBy(mseRes).toNumber();

    this._data.fvalue = fvalue;

    return fvalue;
  }

  /*
        @desc P-value of the F-statistic

        @return {number} F-statistic p-value
    */
  get fProb() {
    return fdist.probability(
      this.fvalue,
      this._data.df.regression,
      this._data.df.residual
    );
  }

  /*
        @desc Liklihood funcion of the fitted model value

        @return {number} liklihood function value
    */
  get llf() {
    let n = this._data.N;
    let ssr = this.ssResidual;
    let sigma = this.sigma;
    let weights = this._data.weights ? this._data.weights : undefined;

    let nobs2 = n / 2;
    let llf = -nobs2 * Math.log(2 * Math.PI);
    llf -= nobs2 * Math.log(ssr / n);
    llf -= nobs2;

    // if sigma exists: do the following:
    if (sigma) {
      //calculate determinate
      let logdet = Math.log(math.det(sigma));
      llf -= 0.5 * logdet;
    }

    if (weights) {
      llf += 0.5 * Arr.sum(weights.map((i) => Math.log(i)));
    }

    this._data.llf = llf;

    return llf;
  }

  /*
        @desc Akaike's information criteria

        @return {number} AIC value
    */
  get aic() {
    let llf = this.llf ? this._data.llf : this.llf;

    return -2 * llf + 2 * (this._data.df.regression + this._data.kConstant);
  }

  /*
        @desc Bayes' information criteria

        @return {number} BIC value
    */
  get bic() {
    let llf = this.llf ? this._data.llf : this.llf;

    return (
      -2 * llf +
      Math.log(this._data.N) * (this._data.df.regression + this._data.kConstant)
    );
  }

  /*
        @desc N is the number of observations in the dataset

        @return {number} # of observations
    */
  get N() {
    return this._data.N;
  }

  /*
        @desc K is the number of parameters in the model

        @return {number} # of parameters
    */
  get K() {
    return this._data.K;
  }

  get X() {
    return this._data.X;
  }

  /*
        @desc R-squared of a model with an intercept

        @return {number} rsq
    */
  get rsq() {
    let rsq = null;
    let ssResidual = this._data.ssResidual
      ? this._data.ssResidual
      : this.ssResidual;

    if (this._data.kConstant > 0) {
      let ssTotalCentered = this._data.ssTotalCentered
        ? this._data.ssTotalCentered
        : this.ssTotalCentered;

      rsq = Decimal(1)
        .minus(Decimal(ssResidual).dividedBy(ssTotalCentered))
        .toNumber();
    } else {
      let ssTotalUncentered = this._data.ssTotalUncentered
        ? this._data.ssTotalUncentered
        : this.ssTotalUncentered;

      rsq = Decimal(1)
        .minus(Decimal(ssResidual).dividedBy(ssTotalUncentered))
        .toNumber();
    }

    this._data.rsq = rsq;

    return rsq;
  }

  /*
        @desc Adjusted R-squared

        @return {number} rsq adjusted
    */
  get rsqAdj() {
    let rsq = this._data.rsq ? this._data.rsq : this.rsq;

    return 1 - ((1 - rsq) * (this.N - 1)) / (this.N - this.K - 1);
  }

  /*
        @desc Omnibus normality test


        @return {number} Omni score
    */
  get omnibus() {
    return Arr.omnibus(this._data.residuals);
  }

  /*
        @desc Skew

        @return {number} skew
    */
  get skew() {
    return Arr.skew(this._data.residuals);
  }

  /*
        @desc Kurtosis

        @return {number} kurtosis
    */
  get kurtosis() {
    return Arr.kurtosis(this._data.residuals);
  }

  /*
        @desc Jarque-Bera test for normality

        @return {number} jb value
    */
  get jarqueBera() {
    let jb = Arr.jarqueBera(this.skew, this.kurtosis, this.N);

    let obj = {
      jb: jb,
      pval: 1 - chidist.cdf(jb, 2),
    };
    return obj;
  }

  /*
        @desc Durbin-Watson statistic

        @return {number} dw value
    */
  get durbinWatson() {
    return Arr.durbinWatson(this._data.residuals);
  }

  /*
        @desc Condition Number
            via Euclidean norm aka 2- norm

        @return {number} condition number
    */
  get conditionNumber() {
    let Asq = this._data.A.reduce((sum, col) => {
      return sum.plus(
        col.reduce((_s, i) => {
          return Decimal(_s).plus(Decimal(i).pow(2));
        }, Decimal(0))
      );
    }, Decimal(0)).sqrt();

    let Ainvsq = this._data.Ainv.reduce((sum, col) => {
      return sum.plus(
        col.reduce((_s, i) => {
          return Decimal(_s).plus(Decimal(i).pow(2));
        }, Decimal(0))
      );
    }, Decimal(0)).sqrt();

    let frac = Asq.times(Ainvsq).sqrt().toNumber();

    return frac;
  }

  /*
        @desc Sum Squared of the regression
            ---

        @return {number} tbd
    */
  get ssRegression() {
    let wMean = Decimal(this._data.wendogMean);

    let ss = Arr.sum(
      this._data.predicted.map((yh) =>
        Decimal(yh).minus(wMean).pow(2).toNumber()
      )
    );

    this._data.ssRegression = ss;

    return ss;
  }

  /*
        @desc Sum Squared of the regression
            ---

        @return {number} tbd
    */
  get ssResidual() {
    let ss = Arr.sum(
      this._data.wendog.map((y, i) =>
        Decimal(y).minus(this._data.predicted[i]).pow(2).toNumber()
      )
    );

    this._data.ssResidual = ss;

    return ss;
  }

  /*
        @desc Sum Squared total centered
            ---

        @return {number} tbd
    */
  get ssTotalCentered() {
    let ss = null;

    // is regression weighted?
    if (this._data.weights) {
      // calculate weighted average of endog and weights
      let wa = Arr.mean(this._data.endog.data, { weights: this._data.weights });

      // diff between endog and wa, squared
      let emWA = this._data.endog.data.map((y, i) =>
        Decimal(y).minus(wa).pow(2).toNumber()
      );

      // sum of weights * emWA
      ss = Arr.sum(
        emWA.map((v, i) => Decimal(v).times(this._data.weights[i]).toNumber())
      );
    } else {
      ss = Arr.sum(
        this._data.wendog.map((y, i) =>
          Decimal(y).minus(this._data.wendogMean).pow(2).toNumber()
        )
      );
    }

    this._data.ssTotalCentered = ss;

    return ss;
  }

  /*
        @desc Sum Squared total centered
            ---

        @return {number} tbd
    */
  get ssTotalUncentered() {
    let ss = Arr.sum(
      this._data.wendog.map((y) => Decimal(y).pow(2).toNumber())
    );

    this._data.ssTotalUncentered = ss;

    return ss;
  }

  /*
        @desc Mean Square error of the Regression
            ---

        @return {number} tbd
    */
  get mseRegression() {
    let ssReg = this._data.ssRegression
      ? this._data.ssRegression
      : this.ssRegression;

    let mseReg = Decimal(ssReg).dividedBy(this._data.df.regression).toNumber();

    this._data.mseRegression = mseReg;

    return mseReg;
  }

  /*
        @desc Mean Squared error of residual
            ---

        @return {number} tbd
    */
  get mseResidual() {
    let ssRes = this._data.ssResidual ? this._data.ssResidual : this.ssResidual;

    let mseRes = Decimal(ssRes).dividedBy(this._data.df.residual).toNumber();

    this._data.mseResidual = mseRes;

    return mseRes;
  }

  /*
        @desc Mean Squared error of model
            ---

        @return {number} tbd
    */
  get mseModel() {
    let ess = this._data.ess ? this._data.ess : this.ess;

    let mseModel = Decimal(ess).dividedBy(this._data.df.regression).toNumber();

    this._data.mseModel = mseModel;

    return mseModel;
  }

  /*
        @desc Explained Sum of Squares
            ---

        @return {number} ESS
    */
  get ess() {
    let ess = null;

    if (this._data.kConstant > 0) {
      let ssTotalCentered = this._data.ssTotalCentered
        ? this._data.ssTotalCentered
        : this.ssTotalCentered;
      let ssResidual = this._data.ssResidual
        ? this._data.ssResidual
        : this.ssResidual;

      ess = Decimal(ssTotalCentered).minus(ssResidual).toNumber();
    } else {
      let ssTotalUncentered = this._data.ssTotalUncentered
        ? this._data.ssTotalUncentered
        : this.ssTotalUncentered;
      let ssResidual = this._data.ssResidual
        ? this._data.ssResidual
        : this.ssResidual;

      ess = Decimal(ssTotalUncentered).minus(ssResidual).toNumber();
    }

    this._data.ess = ess;

    return ess;
  }

  /*
        Predict function
    */
  predict({ exog = null, params = null } = {}) {
    if (!params) params = this.params;
    if (!exog) exog = this._data.X;

    return Arr.dot(exog, params);
  }

  /*
        @desc Full fit of the model

        @calculation
            - coefficient: beta, standard error, t-statistic, p-value, CIs
            - Predicted values
            - Residual values
    */
  fit() {
    this._data.A = Arr.dot(M.transpose(this._data.wexog), this._data.wexog);
    let xty = Arr.dot(M.transpose(this._data.wexog), this._data.wendog);

    // copy not reference
    let qrFactor = M.householder(this._data.A.map((i) => i.slice()).slice());

    // Qt * Y
    let qty = Arr.dot(M.transpose(qrFactor.Q), xty);

    // beta coefficients
    // tranposed bc bacsolve expects matrix format
    let betas = M.backsolve(qrFactor.R, M.transpose([qty]));
    betas.forEach((beta, i) => {
      this._data.exog[i]["beta"] = beta;
    });

    // calculate predicted Y values
    this._data.predicted = this._data.wendog.map((y, i) => {
      let _sum = Decimal(0);
      for (let _b = 0; _b < betas.length; _b++) {
        _sum = _sum.plus(Decimal(this._data.wexog[i][_b]).times(betas[_b]));
      }
      return _sum.toNumber();
    });

    this._data.residuals = this._data.wendog.map((y, i) =>
      Decimal(y).minus(this._data.predicted[i]).toNumber()
    );

    // Add mean attribute to endog object
    this._data.wendogMean = Arr.mean(this._data.wendog);

    // get the inverse:
    //A^-1 = R^-1 * Q'
    let Rinv = M.invUpperTriangle(qrFactor.R);
    this._data.Ainv = Arr.dot(Rinv, M.transpose(qrFactor.Q));

    // z value
    let zval = tdist.ppf(1 - this.alpha / 2, this._data.df.residual);

    // Map Standard Errors of coefficients & T stats
    this._data.exog.map((coeff, i) => {
      coeff.sterror = Decimal(this._data.Ainv[i][i])
        .times(this.mseResidual)
        .sqrt()
        .toNumber();
      coeff.tstat = Decimal(coeff.beta).dividedBy(coeff.sterror).toNumber();
      coeff.pvalue = Decimal(1)
        .minus(
          tdist.cdf(
            Decimal(coeff.tstat).abs().toNumber(),
            this._data.df.residual
          )
        )
        .times(2)
        .toNumber();
      coeff.ciLower = Decimal(coeff.beta)
        .minus(Decimal(zval).times(coeff.sterror))
        .toNumber();
      coeff.ciHigher = Decimal(coeff.beta)
        .plus(Decimal(zval).times(coeff.sterror))
        .toNumber();
      return coeff;
    });

    // make chainable
    return this;
  }

  summaryJson() {
    const json = {};
    // title
    json["title"] = this.modelType + " Regression Results";

    // model info
    json["Dep. Variable"] = this._data.endog.title;
    json["R-squared"] = this.rsq.toFixed(3);
    json["Model"] = this.modelType;
    json["Adj. R-squared"] = this.rsqAdj.toFixed(3);
    json["Method"] = this.method;
    json["F-statistic"] = this.fvalue.toFixed(1);
    try {
      json["Prob (F-stat)"] = this.fProb.toFixed(6);
    } catch (e) {
      console.log("this.fProb failed");
      json["Prob (F-stat)"] = "-";
    }
    json["Log-Likelihood"] = this.llf.toFixed(3);
    json["No. Observations"] = this.N;
    json["AIC"] = this.aic.toFixed(2);
    json["BIC"] = this.bic.toFixed(2);
    json["Df Residuals"] = this._data.df.residual;
    json["Df Regression"] = this.K;
    json["Covariance Type"] = this.covtype;
    json["Omnibus"] = this.omnibus.omnibus.toFixed(3);
    json["Durbin-Watson"] = this.durbinWatson.toFixed(3);
    json["Prob(Omnibus)"] = this.omnibus.pval.toFixed(3);
    json["Jarque-Bera (JB)"] = this.jarqueBera.jb.toFixed(3);
    json["Skew"] = this.skew.toFixed(3);
    json["Prob(JB)"] = this.jarqueBera.pval.toFixed(3);
    json["Kurtosis"] = this.kurtosis.toFixed(3);
    json["Cond. No."] = this.conditionNumber.toFixed(1);

    json.summary = {};

    for (const c of this._data.exog) {
      json.summary[c.title] = {
        coef: c.beta.toFixed(4),
        stdErr: c.sterror.toFixed(3),
        t: c.tstat.toFixed(3),
        p: c.pvalue.toFixed(3),
        ciLower: c.ciLower.toFixed(3),
        ciHigher: c.ciHigher.toFixed(3),
      };
    }

    return json;
  }

  /*
        @desc Outputs a string respresentaiton of the regression results
            for the console

        @references
            - statsmodels - python library format style

        @return {string} terminal formatted string
    */
  summary() {
    let terminalString = "";
    let termWidth = process.stdout.columns;

    // title
    let title = this.modelType + " Regression Results";

    terminalString += " ".repeat((termWidth - title.length) / 2) + title + "\n";

    // model info
    terminalString += "=".repeat(termWidth);

    let twoColW = Math.floor((termWidth - 4) / 2);

    let A1L = "Dep. Variable: ";
    terminalString +=
      A1L +
      " ".repeat(twoColW - A1L.length - this._data.endog.title.length) +
      this._data.endog.title +
      " ".repeat(4);

    let A1R = "R-squared: ";
    terminalString +=
      A1R +
      " ".repeat(twoColW - A1R.length - this.rsq.toFixed(3).toString().length) +
      this.rsq.toFixed(3) +
      "\n";

    let A2L = "Model: ";
    terminalString +=
      A2L +
      " ".repeat(twoColW - A2L.length - this.modelType.length) +
      this.modelType +
      " ".repeat(4);

    let A2R = "Adj. R-squared: ";
    terminalString +=
      A2R +
      " ".repeat(
        twoColW - A2R.length - this.rsqAdj.toFixed(3).toString().length
      ) +
      this.rsqAdj.toFixed(3) +
      "\n";

    let A3L = "Method: ";
    terminalString +=
      A3L +
      " ".repeat(twoColW - A3L.length - this.method.length) +
      this.method +
      " ".repeat(4);

    let A3R = "F-statistic: ";
    terminalString +=
      A3R +
      " ".repeat(twoColW - A3R.length - this.fvalue.toFixed(1).length) +
      this.fvalue.toFixed(1) +
      "\n";

    let A4L = "Date:";
    let today = new Date();
    let datestamp =
      today.getFullYear() +
      "-" +
      ("0" + (today.getMonth() + 1)).slice(-2) +
      "-" +
      ("0" + today.getDate()).slice(-2);
    terminalString +=
      A4L +
      " ".repeat(twoColW - A4L.length - datestamp.length) +
      datestamp +
      " ".repeat(4);

    let A4R = "Prob (F-stat): ";
    terminalString +=
      A4R +
      " ".repeat(twoColW - A4R.length - this.fProb.toFixed(6).length) +
      this.fProb.toFixed(6) +
      "\n";

    let A5L = "Time:";
    let time = today.toLocaleTimeString();
    terminalString +=
      A5L +
      " ".repeat(twoColW - A5L.length - time.length) +
      time +
      " ".repeat(4);

    let A5R = "Log-Likelihood: ";
    terminalString +=
      A5R +
      " ".repeat(twoColW - A5R.length - this.llf.toFixed(3).length) +
      this.llf.toFixed(3) +
      "\n";

    let A6L = "No. Observations: ";
    terminalString +=
      A6L +
      " ".repeat(twoColW - A6L.length - this.N.toString().length) +
      this.N +
      " ".repeat(4);

    let A6R = "AIC: ";
    terminalString +=
      A6R +
      " ".repeat(twoColW - A6R.length - this.aic.toFixed(2).length) +
      this.aic.toFixed(2) +
      "\n";

    let A7L = "Df Residuals: ";
    terminalString +=
      A7L +
      " ".repeat(
        twoColW - A7L.length - this._data.df.residual.toString().length
      ) +
      this._data.df.residual +
      " ".repeat(4);

    let A7R = "BIC: ";
    terminalString +=
      A7R +
      " ".repeat(twoColW - A7R.length - this.bic.toFixed(2).length) +
      this.bic.toFixed(2) +
      "\n";

    let A8L = "Df Regression: ";
    terminalString +=
      A8L +
      " ".repeat(twoColW - A8L.length - this.K.toString().length) +
      this.K +
      "\n";

    let A9L = "Covariance Type: ";
    terminalString +=
      A9L +
      " ".repeat(twoColW - A9L.length - this.covtype.length) +
      this.covtype +
      "\n";

    // coefficients
    terminalString += "=".repeat(termWidth);

    // get column widths
    let b1colW = Math.max(...this._data.exog.map((a) => a.title.length));
    let b2colW = Math.max(
      ...this._data.exog.map((a) => a.beta.toFixed(4).length)
    );
    let b3colW = Math.max(
      ...this._data.exog.map((a) => a.sterror.toFixed(3).length)
    );
    let b4colW = Math.max(
      ...this._data.exog.map((a) => a.tstat.toFixed(3).length)
    );
    let b5colW = Math.max(
      ...this._data.exog.map((a) => a.pvalue.toFixed(3).length)
    );
    let b6colW = Math.max(
      ...this._data.exog.map((a) => a.ciLower.toFixed(3).length)
    );
    let b7colW = Math.max(
      ...this._data.exog.map((a) => a.ciHigher.toFixed(3).length)
    );

    // space to distribute
    let bFillSpace = Math.floor(
      Math.max(
        termWidth -
          b1colW -
          b2colW -
          b3colW -
          b4colW -
          b5colW -
          b6colW -
          b7colW -
          10,
        0
      ) / 6
    );

    b1colW += 10;
    b2colW += bFillSpace;
    b3colW += bFillSpace;
    b4colW += bFillSpace;
    b5colW += bFillSpace;
    b6colW += bFillSpace;
    b7colW += bFillSpace;

    terminalString += " ".repeat(b1colW);

    terminalString += " ".repeat(Math.floor(b2colW - " Coef".length)) + " Coef";
    terminalString +=
      " ".repeat(Math.floor(b3colW - " std err".length)) + " std err";
    terminalString += " ".repeat(Math.floor(b4colW - " t".length)) + " t";
    terminalString +=
      " ".repeat(Math.floor(b5colW - " p-value".length)) + " p-value";
    terminalString +=
      " ".repeat(Math.floor(b6colW - " lower 95%".length)) + " lower 95%";
    terminalString +=
      " ".repeat(Math.floor(b7colW - " higher 95%".length)) + " higher 95%";

    terminalString += "\n";
    terminalString += "-".repeat(termWidth);

    // get min column width
    // then add the difference to each item per column

    this._data.exog.forEach((c) => {
      terminalString +=
        c.title + " ".repeat(Math.floor(b1colW - c.title.length));
      terminalString +=
        " ".repeat(Math.floor(b2colW - c.beta.toFixed(4).toString().length)) +
        c.beta.toFixed(4);
      terminalString +=
        " ".repeat(
          Math.floor(b3colW - c.sterror.toFixed(3).toString().length)
        ) + c.sterror.toFixed(3);
      terminalString +=
        " ".repeat(Math.floor(b4colW - c.tstat.toFixed(3).toString().length)) +
        c.tstat.toFixed(3);
      terminalString +=
        " ".repeat(Math.floor(b5colW - c.pvalue.toFixed(3).toString().length)) +
        c.pvalue.toFixed(3);
      terminalString +=
        " ".repeat(
          Math.floor(b6colW - c.ciLower.toFixed(3).toString().length)
        ) + c.ciLower.toFixed(3);
      terminalString +=
        " ".repeat(
          Math.floor(b7colW - c.ciHigher.toFixed(3).toString().length)
        ) + c.ciHigher.toFixed(3);

      terminalString += "\n";
    });

    terminalString += "=".repeat(termWidth) + "\n";

    let omnibus = this.omnibus;

    let C1L = "Omnibus: ";
    terminalString +=
      C1L +
      " ".repeat(
        twoColW - C1L.length - omnibus.omnibus.toFixed(3).toString().length
      ) +
      omnibus.omnibus.toFixed(3) +
      " ".repeat(4);

    let C1R = "Durbin-Watson: ";
    terminalString +=
      C1R +
      " ".repeat(
        twoColW - C1R.length - this.durbinWatson.toFixed(3).toString().length
      ) +
      this.durbinWatson.toFixed(3) +
      "\n";

    let C2L = "Prob(Omnibus): ";
    terminalString +=
      C2L +
      " ".repeat(
        twoColW - C2L.length - omnibus.pval.toFixed(3).toString().length
      ) +
      omnibus.pval.toFixed(3) +
      " ".repeat(4);

    let jarquebera = this.jarqueBera;

    let C2R = "Jarque-Bera (JB): ";
    terminalString +=
      C2R +
      " ".repeat(
        twoColW - C2R.length - jarquebera.jb.toFixed(3).toString().length
      ) +
      jarquebera.jb.toFixed(3) +
      "\n";

    let C3L = "Skew: ";
    terminalString +=
      C3L +
      " ".repeat(
        twoColW - C3L.length - this.skew.toFixed(3).toString().length
      ) +
      this.skew.toFixed(3) +
      " ".repeat(4);

    let C3R = "Prob(JB): ";
    terminalString +=
      C3R +
      " ".repeat(
        twoColW - C3R.length - jarquebera.pval.toFixed(3).toString().length
      ) +
      jarquebera.pval.toFixed(3) +
      "\n";

    let C4L = "Kurtosis: ";
    terminalString +=
      C4L +
      " ".repeat(
        twoColW - C4L.length - this.kurtosis.toFixed(3).toString().length
      ) +
      this.kurtosis.toFixed(3) +
      " ".repeat(4);

    let C4R = "Cond. No. ";
    terminalString +=
      C4R +
      " ".repeat(
        twoColW - C4R.length - this.conditionNumber.toFixed(1).toString().length
      ) +
      this.conditionNumber.toFixed(1) +
      "\n";
    terminalString += "=".repeat(termWidth) + "\n";

    // return should be the string,
    // data accessible from post fit attributes
    console.log(terminalString);

    return this;
  }
}

/*
    @class Ols

    @params {array | object} exog
        Exog - desc here
    @params {array | object} endog
        Endog - desc here
*/
class Ols extends Model {
  constructor(endog, exog) {
    // 1st remove references and copy as values
    endog = JSON.parse(JSON.stringify(endog));
    exog = JSON.parse(JSON.stringify(exog));

    super(endog, exog);

    this.modelType = "Ols";
    this.method = "Least Squares";
    this.covtype = "nonrobust";
  }

  /*
        @desc whiten the x values

        // TODo add des here
    */
  whiten(X) {
    return X;
  }
}

/*
    @class Gls

    @params {array | object} exog
        Exog - desc here
    @params {array | object} endog
        Endog - desc here
*/
class Gls extends Model {
  constructor(endog, exog, { sigma, cholsigmainv } = {}) {
    // 1st remove references and copy as values
    endog = JSON.parse(JSON.stringify(endog));
    exog = JSON.parse(JSON.stringify(exog));

    // define sigma and & Cholesky before super call
    let rho = null;

    // sigma vs given
    if (sigma) {
      // TODO add code here
      this.sigma = sigma;
    } else {
      // calculate rho
      // aka the correlation of the residual
      let ols_resid = new Ols(endog, exog).fit().resid;
      let res_fid = new Ols({ title: "a", data: ols_resid.slice(1) }, [
        { title: "b", data: ols_resid.slice(0, -1) },
      ]).fit();
      rho = res_fid.params[0];

      // creates rho ^ toeplitz matrix
      // defines the autocorrelation structure
      let sig = new Array(ols_resid.length);

      for (let a = 0; a < ols_resid.length; a++) {
        let arr = [];
        for (let b = 0; b < ols_resid.length; b++) {
          arr.push(Math.pow(rho, Math.abs(a - b)));
        }
        sig[a] = arr;
      }

      sigma = sig;
    }

    // mathjs inverse
    let inv = math.inv(sigma);

    let n = 8;
    let tri = cho(inv);
    let triR = tri.map((row) => row.map((i) => math.round(i, 6)));
    let tri2 = tri.map((row) => row.concat(new Array(n - row.length).fill(0)));

    // get transposed Cholesky
    let tri2t = math.transpose(tri2);

    // end temp code

    super(endog, exog, { sigma: sigma, cholsigmainv: tri2t });

    this.modelType = "Gls";
    this.method = "Least Squares";
    this.covtype = "nonrobust";

    // apply to model and see if I can get whitened vals
    this.cholsigmainv = tri2t;

    this.sigma = sigma;
  }

  /*
        @desc whiten the x values

        // TODo add des here
        X is a matrix and returns a matrix
    */
  whiten(X) {
    let w = Arr.dot(this.cholsigmainv, X);

    return w;
  }
}

/*
    @class Wls

    @params {array | object} exog
        Exog - desc here
    @params {array | object} endog
        Endog - desc here
    @params {array | object} weights
        Weights - desc here
*/
class Wls extends Model {
  constructor(endog, exog, { weights } = {}) {
    // Change input to bignumber or create array of 1s if no provided weights
    weights = weights ? weights : new Array(endog.data.length).fill(1);

    super(endog, exog, { weights: weights });

    this.modelType = "Wls";
    this.method = "Least Squares";
    this.covtype = "nonrobust";
  }

  /*
        @desc whiten the x values

        // TODo add des here
        multiplies each each by the square root of Weights
        we want to modify the row level items in each column,
        weights correspond to rows
    */
  whiten(X) {
    let w = undefined;

    // if matrix
    if (Array.isArray(X[0])) {
      w = X.map((row, i) =>
        row.map((v) =>
          Decimal(v).times(Decimal(this._data.weights[i]).sqrt()).toNumber()
        )
      );
    } else {
      // array
      w = X.map((v, i) =>
        Decimal(v).times(Decimal(this._data.weights[i]).sqrt()).toNumber()
      );
    }

    return w;
  }
}

module.exports = {
  Ols: Ols,
  Gls: Gls,
  Wls: Wls,
  addConstant: addConstant,
};
