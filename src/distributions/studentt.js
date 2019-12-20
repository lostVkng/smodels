/*

    Student t-distribution
*/

const math = require('../tools/math.js')

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

    let fac = Math.sqrt(x * x + df)

    return math.incBeta((x + fac) / (2 * fac), (df / 2), (df / 2))
}

/*
    @desc Percent Point Function
        The value at a given percentage point
        The inverse of the cumulative density function

    @params {array} p
        The percentage to get the function for
    @params {array} df
        The number of degrees of freedom

    @return {number} effectively the z value
*/
function ppf(p, df){
    let fac = math.invIncBeta(2 * Math.min(p, 1 - p), (df / 2), 0.5);
    let y = Math.sqrt(df * (1 - fac) / fac);
    return (p > 0.5) ? y : -y;
}

module.exports = {
    cdf: cdf,
    ppf: ppf
}
