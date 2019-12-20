/*

    Array math, similar to numpy array

    Use decimal for math but return in number form
*/

const Decimal = require('decimal.js')

const chidist = require('../distributions/chi.js')

/*
    @desc Calculates the dimensions of an array or matrix
*/
function ndim(A){

    let rows = null
    let cols = null

    if(Array.isArray(A)){

        // 1 dim?
        if(Array.isArray(A[0])){
            rows = A.length
            cols = A[0].length
        }else{
            cols = A.length
        }
    }

    return [cols, rows]
}






/*
    @desc Dot product of two arrays

    supports 1 and 2d arrays on both sides
    Scaler to 2d multiplication
*/
function dot(A, B){

    // one of the two are a number
    if( (ndim(A)[0] == null) || (ndim(B)[0] == null) ){

        // one is an array?
        if(Array.isArray(A) || Array.isArray(B)){

            let scaler = Array.isArray(A) ? Decimal(B) : Decimal(A)
            let arr = Array.isArray(A) ? A : B
            arr = arr.slice()

            let result = []

            // multiply it out
            for(let i=0; i<arr.length; i++){

                if(Array.isArray(arr[i])){
                    result.push(arr[i].map(v => Decimal(v).times(scaler).toNumber()))
                }else{
                    result.push(Decimal(arr[i]).times(scaler).toNumber())
                }
            }

            return result

        }else {
            return Decimal(A).times(Decimal(B)).toNumber()
        }
    } else if( (ndim(A)[1] == null) && (ndim(B)[1] == null) ){

        A = A.slice()
        B = B.slice()

        A = A.map(v => Decimal(v))
        B = B.map(v => Decimal(v))

        let c = A.map((v, i) => v.times(B[i]))
        c = c.reduce((sum, v) => sum.plus(v), Decimal(0))

        return c.toNumber()
    } else{
        A = A.slice()
        B = B.slice()

        let aIsArray = !Array.isArray(A[0]) ? true : false
        let bIsArray = !Array.isArray(B[0]) ? true : false

        if(aIsArray){
            // B must be a matrix

            let xRow = A.map(v => Decimal(v))

            let newArr = []

            for(let i=0; i<B[0].length; i++){
                let col = []
                B.forEach((row) => col.push(Decimal(row[i])))

                let _xy = xRow.map((v, i) => v.times(col[i]))
                _xy = _xy.reduce((sum, val, _i) => sum.plus(val) , Decimal(0))
                newArr.push(_xy.toNumber())
            }

            return newArr

        } else if(bIsArray){
            // A must be a matrix

            let newArr = []

            let yRow = B.map(v => Decimal(v))

            for(let i=0; i<A.length; i++){

                let xRow = A[i].map(v => Decimal(v))

                let _xy = xRow.map((v, i) => v.times(yRow[i]))
                _xy = _xy.reduce((sum, val, _i) => sum.plus(val) , Decimal(0))
                newArr.push(_xy.toNumber())
            }

            return newArr
        } else{
            // both are matrices
            // create empty matrix of proper dimensions
            let newMtrx = new Array(A.length).fill(0).map(i => new Array(B[0].length))

            // Loop through each row of A
            for(let rowI=0; rowI < A.length; rowI++){

                // setup the arrays to multiply
                let x = A[rowI].map(v => Decimal(v))
                let y = []

                // loop through each column of B
                for(let colI=0; colI<B[0].length; colI++){

                    // reset and add new Y values
                    y = []
                    B.forEach((row) => y.push(Decimal(row[colI])))

                    // the tonumber after val.times is a bug, temporary but working
                    let _xyArr = x.map((_x, _i) => _x.times(y[_i]))
                    let xyVal = _xyArr.reduce((sum, val, _i) => sum.plus(val) , Decimal(0))

                    newMtrx[rowI][colI] = xyVal.toNumber()
                }
            }

            return newMtrx
        }
    }
}



/*
    @desc Sum a set of values together

    @params {array} vals
        values to be summed

    @return {number} Summed value
*/
function sum(vals){

    return vals.reduce((sum, v) => sum.plus(Decimal(v)),Decimal(0)).toNumber()
}

/*
    @desc Mean of a set of values

    @params {array} vals
        values to be used

    @return {number} Mean value
*/
function mean(vals, {weights} = {}){

    if(weights){
        let weightsSum = sum(weights)
        let percWeights = weights.map(w => Decimal(w).dividedBy(weightsSum))
        vals = vals.map((v, i) => Decimal(v).times(percWeights[i]).toNumber())
        return sum(vals)
    } else{
        return Decimal(sum(vals)).dividedBy(vals.length).toNumber()
    }
}

/*
    @desc Calculates the standard deviation

    @params {array} vals
        values to be used
    @params {string} type
        Should standard deviation be calculated as population or sample

    @return {number} Standard deviation
*/
function std(vals, type='population'){

    let avg = mean(vals)

    let num = sum(vals.map(v => Decimal(v).minus(avg).pow(2).toNumber() ))
    let variance = Decimal(num).dividedBy(vals.length)
    let std = variance.sqrt().toNumber()

    return std
}


/*
    @desc Calculates the skew of a set of values

    @params {array} vals
        Values

    @return {number} skew
*/
function skew(vals){

    let avg = mean(vals)

    let num = mean(vals.map(v => Decimal(v).minus(avg).pow(3).toNumber() ) )
    return Decimal(num).dividedBy(Decimal(std(vals)).pow(3)).toNumber()
}


/*
    @desc Skew Test
        Tests if the skew is different from the normal distribution

        Requires a sample size of at least 8

        Based around D'Agostino(1970)

    @params {array} vals
        Values

    @return {number} Z-score for the test
*/
function skewtest(vals){

    let n = vals.length

    let skw = skew(vals)

    // derived expressions
    let u1 = 0
    let u2 = (6*(n-2)) / ((n+1)*(n+3))
    let v1 = 0
    let v2 = (36*(n-7)*(Math.pow(n,2)+(2*n)-5)) / ((n-2)*(n+5)*(n+7)*(n+9))

    // Sample Skewness Transformation
    let w2 = Math.sqrt(2*v2 + 4) - 1
    let w = Math.sqrt(w2)
    let delta = 1 / Math.sqrt(Math.log(w))
    let alpha2 = 2 / (w2 - 1)
    let alpha = Math.sqrt(alpha2)

    let z = delta * Math.asinh(skw / (alpha*Math.sqrt(u2)) )

    return z
}

/*
    @desc kurtosis

    @params {array} vals
        Values

    @return {number} kurtosis
*/
function kurtosis(vals){

    let avg = mean(vals)


    let num = mean(vals.map(v => Decimal(v).minus(avg).pow(4).toNumber() ) )
    return Decimal(num).dividedBy(Decimal(std(vals)).pow(4)).toNumber()
}

/*
    @desc Kurtosis Test
        Tests if the kurtosis is different from the normal distribution

        Requires a sample size of at least 20

        Based around D'Agostino(1970)

    @params {array} vals
        Values

    @return {number} Z-score for the test
*/
function kurtosistest(vals){

    let n = vals.length

    let krt = kurtosis(vals) - 3

    // derived expressions
    let u1 = -6 / (n+1)
    let u2 = (24*n*(n-2)*(n-3)) / ((n+1)*(n+1)*(n+3)*(n+5))
    let v1 = (6*(Math.pow(n,2)-5*n+2)) / ((n+7)*(n+9))
    v1 = v1 * Math.sqrt((6*(n+3)*(n+5)) / (n*(n-2)*(n-3)))
    let v2 = 36*(15*Math.pow(n,6)-36*Math.pow(n,5)-628*Math.pow(n,4)+928*Math.pow(n,3)+5777*Math.pow(n,2)-6402*n+900)
    v2 = v2 / (n*(n-3)*(n-2)*(n+7)*(n+9)*(n+11)*(n+13))

    // Sample Skewness Transformation
    let A = 6 + (8/v1)*((2/v1)+Math.sqrt(1+(4/Math.pow(v1,2))))
    let z = Math.sqrt(9*A/2)
    z = z*(1 - (2/(9*A)) - Math.pow( (1-(2/A))/(1+((krt-u1)/Math.sqrt(u2))*Math.sqrt(2/(A-4))),(1/3)))

    return z
}


/*
    @desc omnibus
        K^2 statistic
        Able to detect deviations from normality due to skewness or Kurtosis
        based on D'agostino, Belanger & D'Agostino 1990

        https://en.wikipedia.org/wiki/D%27Agostino%27s_K-squared_test

    @params {array} vals
        Residuals

    @return {object} omnibus
    @return {object} ombnibus.omnibus K2 value
    @return {object} ombnibus.pval chi squared p-value
*/
function omnibus(resids){

    let st = skewtest(resids)
    let kt = kurtosistest(resids)

    // K2 statistic
    let k2 = Decimal(st).pow(2).plus(Decimal(kt).pow(2)).toNumber()

    let obj = {
        omnibus: k2,
        pval: 1 - chidist.cdf(k2, 2)
    }

    return obj
}


/*
    @desc Jarque Bera

    @params {array} vals
        Values

    @return {number} JB
*/
function jarqueBera(skew, kurtosis, N){

    let k = Decimal(kurtosis).minus(3).pow(2).times(0.25)
    let s = Decimal(skew).pow(2)

    return Decimal(N).dividedBy(6).times(s.plus(k)).toNumber()
}


/*
    @desc Durbin Watson

    @params {array} vals
        Values

    @return {number} Durbin Watson
*/
function durbinWatson(resids){

    let diff_resids = resids.map((v,i) => {
        if(i === 0){
            return 0
        }else{
            return Decimal(resids[i-1]).minus(Decimal(v)).pow(2).toNumber()
        }
    })

    let ssr = sum(resids.map(v => Decimal(v).pow(2).toNumber()))
    let dw = Decimal(sum(diff_resids)).dividedBy(ssr).toNumber()

    return dw
}



module.exports = {
    dot: dot,
    ndim: ndim,
    sum: sum,
    mean: mean,
    std: std,
    skew: skew,
    kurtosis: kurtosis,
    omnibus: omnibus,
    jarqueBera: jarqueBera,
    durbinWatson: durbinWatson
}
