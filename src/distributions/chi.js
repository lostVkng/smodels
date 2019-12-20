/*

    Chi-distribution
*/

const math = require('../tools/math.js')
const Decimal = require('decimal.js')


function Gcf(X,A) {

    // Good for X>A+1
    let A0 = 0
    let B0 = 1
    let A1 = 1
    let B1 = X
    let AOLD = 0
    let N = 0

    while (Math.abs((A1-AOLD)/A1)>.00001) {
        AOLD = A1
        N = N+1
        A0 = A1+(N-A)*A0
        B0 = B1+(N-A)*B0
        A1 = X*A0+N*A1
        B1 = X*B0+N*B1
        A0 = A0/B1
        B0 = B0/B1
        A1 = A1/B1
        B1 = 1
    }
    let Prob = Math.exp(A*Math.log(X)-X-math.lngamma(A))*A1
    return 1-Prob
}

function Gser(X,A) {

    // Good for X<A+1.
    let T9 = 1/A
    let G = T9
    let I = 1
    while (T9>G*.00001) {
        T9 = T9*X/(A+I)
        G = G+T9
        I = I+1
    }
    G = G*Math.exp(A*Math.log(X)-X-math.lngamma(A))

    return G
}


/*
    @desc Cumulative Density Function
        The probability that X will be less than or equal to X

    @params {array} x
        The value to get the probability at or below
    @params {array} df
        The number of degrees of freedom

    @return {number} probability value
*/
function cdf(x, df){

    if(df<=0){
        throw new Error('Degrees of freedom must be positive')
    }

    x = x/2
    df = df/2

    let GI;
    if (x <= 0) {
        GI = 0
    } else if (x <df + 1){
        GI = Gser(x,df)
    } else{
        GI = Gcf(x,df)
    }
    return GI
}

/*
    @desc Survival Function
        Expressed as 1 - CDF

    @params {array} x
        Array of values to calculate SF from
    @params {array} df
        The number of degrees of freedom

    @return {number} Survival Function
*/
function sf(x, df){

    return Decimal(1).minus(cdf(x, df)).toNumber()

}


module.exports = {
    cdf: cdf,
    sf: sf,
}
