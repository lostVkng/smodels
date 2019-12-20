/*
    Discrete regressions

    Currently only Logit regression
*/

const Decimal = require('decimal.js')

const normdist = require('../distributions/norm.js')
const chidist = require('../distributions/chi.js')

const M = require('../tools/matrix.js')
const Arr = require('../tools/array.js')

// TODO delete this
const math = require('mathjs')


/*
    @desc Base model for all discrete regressions

    @param {object} endog: endogenous variables
    @param {array} exog: exogenous variables
*/
class DiscreteModel{
    constructor(endog, exog, {} = {} ){

        // common storage as object
        this._data = {}
        this._data.endog = endog
        this._data.exog = exog
        this._data.df = {}
        this.alpha = 0.05
        this._data.N = 0

        // array length check
        this._data.exog.forEach((a) => {
            if(a.data.length !== this._data.endog.data.length){
                throw Error('Input lengths are not equal')
            }
        })

        // length of data series
        this._data.N = this._data.endog.data.length

        // Create data only matrix
        // X is rows by columns
        let X = new Array(this._data.endog.data.length).fill(0)
        this._data.X = X.map((v, i) => {
            let arr = []
            this._data.exog.forEach((a) => arr.push(a.data[i]) )
            return arr
        })

        // # of non constants
        this._data.K = this._data.exog.filter((v) => v.constant !== true).length

        // # of constants
        this._data.kConstant = this._data.exog.filter((v) => v.constant === true).length

        // degrees of freedom
        this._data.df.regression = this._data.K
        this._data.df.residual = this._data.N - this._data.K - 1
        this._data.df.total = this._data.N - 1
    }



    /*
        @desc The beta coefficients of the fitted model

        @return {array} beta coefficients
    */
    get params(){
        return this._data.exog.map(v => v.beta)
    }


    /*
        @desc The standard errors of the regression coefficients

        @return {array} SEs
    */
    get bse(){
        return this._data.exog.map(v => v.sterror)
    }


    /*
        @desc The z-statistic of the regression coefficients

        @return {array} coefficient z-statistics
    */
    get bzstats(){
        return this._data.exog.map(v => v.zstat)
    }


    /*
        @desc The p-values of the regression coefficients

        @return {array} coefficient P-values
    */
    get pvalues(){
        return this._data.exog.map(v => v.pvalue)
    }


    /*
        @desc Confidence Interval for the regression coefficients

        @return {object} CI Confidence Intervals
        @return {object} CI.lower lower 95% CI value
        @return {object} CI.higher higher 95% CI value
    */
    get confInt(){
        return this._data.exog.map(v => {
            let obj = {
                'lower': v.ciLower,
                'higher': v.ciHigher
            }

            return obj
        })
    }

    /*
        @desc N is the number of observations in the dataset

        @return {number} # of observations
    */
    get N(){
        return this._data.N
    }

    /*
        @desc K is the number of parameters in the model

        @return {number} # of parameters
    */
    get K(){
        return this._data.K
    }

    get X(){
        return this._data.X
    }


    /*
        @desc Outputs a string respresentaiton of the regression results
            for the console

        @references
            - statsmodels - python library format style

        @return {string} terminal formatted string
    */
    summary(){
        let terminalString = ''
        let termWidth = process.stdout.columns

        // title
        let title = this.modelType + ' Regression Results'

        terminalString += ' '.repeat( (termWidth - title.length) / 2 ) + title + '\n'

        // model info
        terminalString += '='.repeat(termWidth)

        let twoColW = Math.floor((termWidth - 4) / 2)

        let A1L = 'Dep. Variable: '
        terminalString += A1L + ' '.repeat(twoColW-A1L.length-this._data.endog.title.length)+this._data.endog.title+' '.repeat(4)

        let A1R = 'No. Observations: '
        terminalString += A1R+ ' '.repeat(twoColW-A1R.length-this._data.N.toString().length)+this._data.N+'\n'

        let A2L = 'Model: '
        terminalString += A2L + ' '.repeat(twoColW-A2L.length-this.modelType.length)+this.modelType+' '.repeat(4)

        let A2R = 'Df Residuals: '
        let dfResid = this._data.df.residual ? this._data.df.residual : NaN
        terminalString += A2R + ' '.repeat(twoColW-A2R.length-dfResid.toString().length)+dfResid+'\n'

        let A3L = 'Method: '
        terminalString += A3L + ' '.repeat(twoColW-A3L.length-this.method.length)+this.method+' '.repeat(4)

        let A3R = 'Df Regression: '
        let dfModel = this._data.df.regression ? this._data.df.regression : NaN
        terminalString += A3R + ' '.repeat(twoColW-A3R.length-dfModel.toString().length)+dfModel+'\n'

        let A4L = 'Date:'
        let today = new Date()
        let datestamp = today.getFullYear()+'-'+('0'+(today.getMonth()+1)).slice(-2)+'-'+('0'+today.getDate()).slice(-2)
        terminalString += A4L + ' '.repeat(twoColW-A4L.length-datestamp.length)+datestamp+' '.repeat(4)

        let A4R = 'Pseudo R-squ.: '
        terminalString += A4R + ' '.repeat(twoColW-A4R.length-this.pseudoRsq.toFixed(3).length)+this.pseudoRsq.toFixed(3)+'\n'

        let A5L = 'Time:'
        let time = today.toLocaleTimeString()
        terminalString += A5L + ' '.repeat(twoColW-A5L.length-time.length)+time+' '.repeat(4)

        let A5R = 'Log-Likelihood: '
        terminalString += A5R + ' '.repeat(twoColW-A5R.length-this.llf.toFixed(3).length)+this.llf.toFixed(3)+'\n'

        let A6L = 'Converged: '
        terminalString += A6L + ' '.repeat(twoColW-A6L.length-this.converged.toString().length)+this.converged.toString()+' '.repeat(4)

        let A6R = 'LL-Null: '
        terminalString += A6R + ' '.repeat(twoColW-A6R.length-this.llnull.toFixed(3).length)+this.llnull.toFixed(3)+'\n'

        let A7L = ''
        terminalString += A7L + ' '.repeat(twoColW)+' '.repeat(4)

        let A7R = 'LR p-value: '
        terminalString += A7R + ' '.repeat(twoColW-A7R.length-this.llrPvalue.toFixed(3).length)+this.llrPvalue.toFixed(3)+'\n'

        // coefficients
        terminalString += '='.repeat(termWidth)

        // get column widths
        let b1colW = Math.max(...this._data.exog.map(a => a.title.length))
        let b2colW = Math.max(...this._data.exog.map(a => 'NaN'.length))
        let b3colW = Math.max(...this._data.exog.map(a => 'NaN'.length))
        let b4colW = Math.max(...this._data.exog.map(a => 'NaN'.length))
        let b5colW = Math.max(...this._data.exog.map(a => 'NaN'.length))
        let b6colW = Math.max(...this._data.exog.map(a => 'NaN'.length))
        let b7colW = Math.max(...this._data.exog.map(a => 'NaN'.length))

        // space to distribute
        let bFillSpace = Math.floor(Math.max( (termWidth-b1colW-b2colW-b3colW-b4colW-b5colW-b6colW-b7colW-10), 0) / 6)

        b1colW += 10
        b2colW += bFillSpace
        b3colW += bFillSpace
        b4colW += bFillSpace
        b5colW += bFillSpace
        b6colW += bFillSpace
        b7colW += bFillSpace

        terminalString += ' '.repeat(b1colW)

        terminalString += ' '.repeat( Math.floor( (b2colW - ' Coef'.length) )  )+' Coef'
        terminalString += ' '.repeat( Math.floor( (b3colW - ' std err'.length) )  )+' std err'
        terminalString += ' '.repeat( Math.floor( (b4colW - ' t'.length) )  )+' t'
        terminalString += ' '.repeat( Math.floor( (b5colW - ' p-value'.length) )  )+' p-value'
        terminalString += ' '.repeat( Math.floor( (b6colW - ' lower 95%'.length) )  )+' lower 95%'
        terminalString += ' '.repeat( Math.floor( (b7colW - ' higher 95%'.length) )  )+' higher 95%'

        terminalString += '\n'
        terminalString += '-'.repeat(termWidth)


        // get min column width
        // then add the difference to each item per column
        this._data.exog.forEach((c) => {

            terminalString += c.title+' '.repeat( Math.floor( (b1colW - c.title.length) )  )
            terminalString += ' '.repeat( Math.floor( (b2colW - c.beta.toFixed(4).toString().length) )  )+c.beta.toFixed(4)
            terminalString += ' '.repeat( Math.floor( (b3colW - c.sterror.toFixed(3).toString().length) )  )+c.sterror.toFixed(3)
            terminalString += ' '.repeat( Math.floor( (b4colW - c.zstat.toFixed(3).toString().length) )  )+c.zstat.toFixed(3)
            terminalString += ' '.repeat( Math.floor( (b5colW - c.pvalue.toFixed(3).toString().length) )  )+c.pvalue.toFixed(3)
            terminalString += ' '.repeat( Math.floor( (b6colW - c.ciLower.toFixed(3).toString().length) )  )+c.ciLower.toFixed(3)
            terminalString += ' '.repeat( Math.floor( (b6colW - c.ciHigher.toFixed(3).toString().length) )  )+c.ciHigher.toFixed(3)

            terminalString += '\n'
        })

        terminalString += '='.repeat(termWidth) + '\n'

        // return should be the string,
        // data accessible from post fit attributes
        console.log(terminalString)

        return this
    }

}




class Logit extends DiscreteModel{

    constructor(endog, exog, {startParams} = {} ){

        super(endog, exog, {})

        this.modelType = 'Logit'
        this.method = 'MLE'
        this.converged = false
    }

    get llf(){

        let lf = this.loglike(this.params)

        this._data.llf = lf

        return lf
    }

    get llnull(){
        /*
            This could ceratinly be better, but it is acceptable for a temporary fix
        */

        // fit a model with a length of 1s to equal endog
        let onesExog = new Array(this._data.endog.data.length).fill(1)

        let testOne = new Logit(this._data.endog, [{'data': onesExog}]).fit()

        this._data.llnull = testOne.llf

        return testOne.llf
    }

    get llr(){

        let llnull = this._data.llnull ? this._data.llnull : this.llnull
        let llf = this._data.llf ? this._data.llf : this.llf

        let llr = Decimal(-2).times( Decimal(llnull).minus(llf) ).toNumber()

        this._data.llr = llr

        return llr
    }

    get llrPvalue(){

        let llr = this._data.llr ? this._data.llr : this.llr

        let llrPvalue = chidist.sf(this.llr, this._data.df.regression )

        this._data.llrPvalue = llrPvalue

        return llrPvalue
    }

    get pseudoRsq(){
        /*
            Using McFadden's pseudo - R - Squared
        */

        let llf = this._data.llf ? this._data.llf : this.llf
        let llnull = this._data.llnull ? this._data.llnull : this.llnull

        let prsq = Decimal(1).minus(Decimal(llf).dividedBy(llnull)).toNumber()

        this._data.pseudoRsq = prsq

        return prsq
    }

    loglike(params){

        let q = this._data.endog.data.slice()
        q = q.map(v => Decimal(v).times(2).minus(1).toNumber())
        let X = this._data.X

        let qx = Arr.dot(X, params).map((v, i) => Decimal(v).times(q[i]).toNumber())

        let cdf = this.cdf(qx)
        let llf = cdf.reduce((sum, val, _i) => sum.plus(Decimal(val).ln()) , Decimal(0)).toNumber()

        return llf
    }

    cdf(X){
        // The logistic cumulative distribution function

        // assume x is an array
        let e = Decimal(2.718281)

        // convert X array to Decimal
        let one = Decimal(1)
        X = X.slice()
        X = X.map(v => Decimal(v))

        let cdf = X.map(x => one.dividedBy( one.plus( Decimal.pow(e, Decimal(-1).times(x)) ) ))

        return cdf.map(v => v.toNumber())
    }

    hessian(params){
        // Logit model Hessian matrix of the log-likelihood
        // this could all be cleaner, come back and revisit once it works
        let X = this._data.X.slice()

        let L = this.cdf(Arr.dot(X, params))

        // apply hessian
        L = L.map(v => Decimal(v).times(Decimal(1).minus(v)) )
        let Xt = M.transpose(X)
        let lcalcT = Xt.map(a => a.map((v, i) => Decimal(v).times(L[i]).toNumber()))

        let dot = Arr.dot(lcalcT, X)
        dot = Arr.dot(dot, -1)

        return dot
    }

    score(params){
        // Logit model (gradient) vector of log-liklihood

        let y = this._data.endog.data.slice()
        let x = this._data.X.slice()

        let L = this.cdf(Arr.dot(x, params))

        let a = y.map((v, i) => Decimal(v).minus(Decimal(L[i])).toNumber())

        return Arr.dot(a, x)
    }

    newton({maxiter=100} = {} ){
        /*
            this is the Newton-Raphson solver, it should be moved to another Something
            but is here until we need ti something/somewhere else
            move to some optimizer file

        */
        // current number of iterations
        let iteration = 0

        // tolerance for optimization
        let tol = 0.000000001

        let oldParams = new Array(this._data.X[0].length).fill(Decimal(-Infinity))
        let newParams = new Array(this._data.X[0].length).fill(Decimal(0))
        let nobs = this._data.N

        // while loop
        while(
            (iteration < maxiter) &&
            (newParams.some((v, i) => v.minus(oldParams[i]).abs() > tol))
        ){
            let H = this.hessian(newParams)
            H = H.map(a => a.map(b => Decimal(b).dividedBy(nobs).toNumber()))

            // mathjs inverse
            let Hinv = math.inv(H)

            // score
            let score = this.score(newParams)
            score = score.map(a => Decimal(a).dividedBy(nobs).toNumber())

            oldParams = newParams
            let invHS = Arr.dot(Hinv, score)
            newParams = oldParams.map((v, i) => Decimal(v).minus(invHS[i]))

            iteration += 1
        }

        return {
            'betas': newParams.map(a => a.toNumber()),
            'converged': true,
        }
    }

    /*
        Predict function
    */
    predict({exog=null, params=null, linear=false} = {}){

        if(!params) params = this.params
        if(!exog) exog = this._data.X

        if(!linear){

            return this.cdf(Arr.dot(exog, params))
        }else{
            return Arr.dot(exog, params)
        }
    }

    fit(){
        let newtonOpt = this.newton()
        let betas = newtonOpt['betas']
        betas.forEach((beta, i) => {
            this._data.exog[i]['beta'] = beta
        })

        this.converged = newtonOpt.converged

        // calculate predicted Y values
        this._data.predicted = this._data.endog.data.map((y, i) => {
            let _sum = Decimal(0)
            for(let _b=0; _b<betas.length; _b++){
                _sum = _sum.plus(Decimal(this._data.exog[_b]['data'][i]).times(betas[_b]))
            }
            return _sum.toNumber()
        })

        this._data.residuals = this._data.endog.data.map((y, i) => Decimal(y).minus(this._data.predicted[i]).toNumber())

        /*
            inverse of negative hessian retval )) / nobs
        */

        let hessOpt = this.hessian(betas).map(a => a.map(b => Decimal(b).dividedBy(this._data.N).toNumber()))
        let hInv = math.inv(Arr.dot(hessOpt, -1)).map(a => a.map(b => Decimal(b).dividedBy(this._data.N).toNumber()))

        let zval = normdist.ppf((Decimal(1).minus(this.alpha / 2).toNumber()))

        // Map Standard Errors of coefficients & T stats
        this._data.exog.map((coeff, i) => {
            coeff.sterror = Decimal(hInv[i][i]).sqrt().toNumber()
            coeff.zstat = Decimal(coeff.beta).dividedBy(coeff.sterror).toNumber()
            coeff.pvalue = normdist.sf(Math.abs(coeff.zstat)) * 2
            coeff.ciLower = Decimal(coeff.beta).minus(Decimal(zval).times(coeff.sterror)).toNumber()
            coeff.ciHigher = Decimal(coeff.beta).plus(Decimal(zval).times(coeff.sterror)).toNumber()
            return coeff
        })

        // make chainable
        return this

    }

}



module.exports = {
    Logit: Logit
}
