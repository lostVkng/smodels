/*

    Gamma
        http://en.wikipedia.org/wiki/Lanczos_approximation

    Beta
        Adapted from:
        https://github.com/AndreasMadsen/mathfn/blob/master/functions/beta.js

*/


// p constant to be used in gamma function
const p = [
    676.5203681218851,
    -1259.1392167224028,
    771.32342877765313,
    -176.61502916214059,
    12.507343278686905,
    -0.13857109526572012,
    9.9843695780195716e-6,
    1.5056327351493116e-7
]

// p constant to be used in log gamma function
const p_ln = [
    0.99999999999999709182,
    57.156235665862923517,
    -59.597960355475491248,
    14.136097974741747174,
    -0.49191381609762019978,
    0.33994649984811888699e-4,
    0.46523628927048575665e-4,
    -0.98374475304879564677e-4,
    0.15808870322491248884e-3,
    -0.21026444172410488319e-3,
    0.21743961811521264320e-3,
    -0.16431810653676389022e-3,
    0.84418223983852743293e-4,
    -0.26190838401581408670e-4,
    0.36899182659531622704e-5
]


/*
    @desc Gamma function calculation using the lanczos approximation

    @params {number} z
        Value for calculating gamma

    @return {number} Gamma
*/
function gamma(z){

    // lanczos method invalid
    if (z < 0.5) {
        // reflection formula
        return Math.PI / (Math.sin(Math.PI * z) * gamma(1 - z))
    } else if(z > 100){
        return Math.exp(lngamma(z))
    } else {
        z -= 1
        let x = 0.99999999999980993
        for(let i = 0; i < p.length; i++){
            x += p[i] / (z + i)
        }
        let t = z + p.length - 0.5

        return Math.sqrt(2 * Math.PI) * Math.pow(t, z + 0.5) * Math.exp(-t) * x
    }
}

/*
    @desc Natural log Gamma
        Spouge approximation
        Useful for large values when lanczos isn't suitible

    @params {number} z
        Value for calculating log gamma

    @return {number} log gamma
*/
function lngamma(z) {

    if(z < 0) return Number('0/0')

    let x = p_ln[0]
    for(let i = p_ln.length - 1; i > 0; --i){
        x += p_ln[i] / (z + i)
    }

    // issue is i'm taking athe log of a negative x
    let t = z + (607/128) + 0.5;

    return 0.5 * Math.log(2*Math.PI) + (z+0.5) * Math.log(t) - t + Math.log(x) - Math.log(z)
}



function log1p(x) {
    if (x <= -1.0) {
        throw new RangeError('Argument mustbe greater than -1.0');
    }

    // x is large enough that the obvious evaluation is OK
    else if (Math.abs(x) > 1e-4) {
        return Math.log(1.0 + x);
    }

    // Use Taylor approx. log(1 + x) = x - x^2/2 with error roughly x^3/3
    // Since |x| < 10^-4, |x|^3 < 10^-12, relative error less than 10^-8
    else {
        return (-0.5*x + 1.0)*x;
    }
}

// evaluates the continued fraction for incomplete beta function by modified Lentz's method.
function betacf(x, a, b) {

    let fpmin = 1e-30
    let m = 1
    let m2, aa, c, d, del, h, qab, qam, qap;

    // These q's will be used in factors that occur in the coefficients
    qab = a + b
    qap = a + 1
    qam = a - 1
    c = 1
    d = 1 - qab * x / qap
    if (Math.abs(d) < fpmin) d = fpmin
    d = 1 / d
    h = d

    for (; m <= 100; m++) {
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        // One step (the even one) of the recurrence
        d = 1 + aa * d
        if (Math.abs(d) < fpmin) d = fpmin
        c = 1 + aa / c
        if (Math.abs(c) < fpmin) c = fpmin
        d = 1 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        // Next step of the recurrence (the odd one)
        d = 1 + aa * d
        if (Math.abs(d) < fpmin) d = fpmin
        c = 1 + aa / c
        if (Math.abs(c) < fpmin) c = fpmin
        d = 1 / d
        del = d * c
        h *= del
        if (Math.abs(del - 1.0) < 3e-7) break
    }
    return h
}



/*
    @desc the Incomplete Beta Function

    @params {number} x
        Random variable X Beta distribution
    @params {number} a
        Parameter a of beta distribution
    @params {number} b
        Parameter b of beta distribution

    @return {number} the incomplete beta value
*/
function incBeta(x, a, b) {

	if(x < 0 || x > 1) {
        throw new RangeError('First argument must be between 0 and 1.');
	}
    // Special cases
    else if (a === 1 && b === 1) return x;
    else if (x === 0) return 0;
    else if (x === 1) return 1;
    else if (a === 0) return 1;
    else if (b === 0) return 0;
    else {
        let bt = Math.exp(lngamma(a + b) - lngamma(a) - lngamma(b) +
                    a * Math.log(x) + b * log1p(-x));

        // Use continued fraction directly
        if (x < (a + 1) / (a + b + 2)){

            return bt * betacf(x, a, b) / a
        } else {
            // else use continued fraction after making
            // the symmetry transformation

            return 1 - bt * betacf(1 - x, b, a) / b
        }
    }
}


/*
    @desc the Inverse of the Incomplete Beta Function

    @params {number} x
        Random variable X Beta distribution
    @params {number} a
        Parameter a of beta distribution
    @params {number} b
        Parameter b of beta distribution

    @return {number} the incomplete beta value
*/
function invIncBeta(p, a, b) {

    if(p < 0 || p > 1) {
        throw new RangeError('First argument must be between 0 and 1.')
    }
    // Special cases
    else if (a === 1 && b === 1) return p;
    else if (p === 1) return 1;
    else if (p === 0) return 0;
    else if (a === 0) return 0;
    else if (b === 0) return 1;

    else {
        let EPS = 1e-8
        let a1 = a - 1
        let b1 = b - 1
        let j = 0
        let lna, lnb, pp, t, u, err, x, al, h, w, afac;

        if(a >= 1 && b >= 1) {
            pp = (p < 0.5) ? p : 1 - p;
            t = Math.sqrt(-2 * Math.log(pp))
            x = (2.30753 + t * 0.27061) / (1 + t* (0.99229 + t * 0.04481)) - t

            if(p < 0.5){
                x = -x
            }

            al = (x * x - 3) / 6
            h = 2 / (1 / (2 * a - 1)  + 1 / (2 * b - 1))
            w = (x * Math.sqrt(al + h) / h) - (1 / (2 * b - 1) - 1 / (2 * a - 1)) * (al + 5 / 6 - 2 / (3 * h))
            x = a / (a + b * Math.exp(2 * w))

        } else {
            lna = Math.log(a / (a + b))
            lnb = Math.log(b / (a + b))
            t = Math.exp(a * lna) / a
            u = Math.exp(b * lnb) / b
            w = t + u


            if (p < t / w){
                x = Math.pow(a * w * p, 1 / a)
            } else{
                x = 1 - Math.pow(b * w * (1 - p), 1 / b)
            }

        }

        afac = -lngamma(a) - lngamma(b) + lngamma(a + b)

        for(let j=0; j < 10; j++) {

            if(x === 0 || x === 1){
                return x
            }

            err = incBeta(x, a, b) - p

            t = Math.exp(a1 * Math.log(x) + b1 * log1p(-x) + afac)
            u = err / t
            x -= (t = u / (1 - 0.5 * Math.min(1, u * (a1 / x - b1 / (1 - x)))))

            if (x <= 0) x = 0.5 * (x + t)
            if (x >= 1) x = 0.5 * (x + t + 1)

            if (Math.abs(t) < EPS * x && j > 0) break
        }

        return x
    }
}



module.exports = {
    gamma: gamma,
    lngamma: lngamma,
    incBeta: incBeta,
    invIncBeta: invIncBeta
}
