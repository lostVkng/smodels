/*

    F-distribution
*/

function L504(a, f, b, iv){
        let q = a * f / (a * f + b)
        let sa = Math.sqrt(q)
        let sl = Math.log(sa)
        let ca = Math.sqrt(1 - q)
        let cl = Math.log(ca)
        let al = Math.atan(sa / Math.sqrt(-sa * sa + 1))
        let fp = 1 - 2 * al / Math.PI
        let r = 0.0
        if (b != 1){
            let c = Math.log(2 * sa / Math.PI)
            fp -= Math.exp(c + cl)

            if (b != 3){
                let n = Math.floor((b - 3) / 2)

                for (let i = 1; i <= n; i++){
                    let x = 2 * i + 1
                    r += Math.log((x - 1) / x)
                    let rr = r + cl * x + c
                    if (rr > -78.4){
                        fp -= Math.exp(rr)
                    }
                }
            }
        }

        if (a != 1){
            let c = r

            if (b > 1){
                c += Math.log(b - 1)
            }

            c += Math.log(2 / Math.PI) + sl + cl * b

            if (c > -78.4) {
                fp += Math.exp(c)
            }

            if (a != 3){
                let n = Math.floor((a - 3) / 2)
                r = 0
                for (let i = 1; i <= n; i++){
                    let x = i * 2 + 1
                    r += Math.log((b + x - 2) / x)
                    let rr = r + sl * (x - 1) + c
                    if (rr > -78.4){
                        fp += Math.exp(rr)
                    }
                }
            }
        }
        return fp
}

function L401(a, f, b, iv){
        let q = a * f / (a * f + b)
        let ql = Math.log(q)
        let fp = 0.0
        let c = Math.log(1 - q) * b / 2
        if (c > -78.4){
            fp = Math.exp(c)
        }

        if (a != 2){
            let n = Math.floor(a / 2 - 1)
            let r = 0.0
            for (let i = 1; i <= n; i++){
                let x = 2 * i
                r += Math.log(b + x - 2) - Math.log(x) + ql
                if (r + c > -78.4){
                    fp += Math.exp(r + c)
                }
            }
        }

        if (iv == 1){
            fp = 1 - fp
        }

        return fp
}


/*
    @desc probability
        The probability based on a F distribution

    @params {array} f
        The f value
    @params {array} dfReg
        Degrees of freedom of regression
    @params {array} dfRes
        Degrees of freedom of residuals

    @return {number} probability value
*/
function probability(f, dfReg, dfRes){
        let a = dfReg
        let b = dfRes
        let iv = 0

        if (Math.floor(a / 2) * 2 == a){
            //even numerator df
            let fp = L401(a, f, b, iv)
            return fp
        } else if (Math.floor(b / 2) * 2 != b){
            let fp = L504(a, f, b, iv)
            return fp
        }

        f = 1 / f
        a = dd
        b = dn
        iv = 1
        return L401(a, f, b, iv)
}

/*
    @desc probacumulativeProbbility
        The cumulative probability based on a F value

    @params {array} f
        The f value
    @params {array} dfReg
        Degrees of freedom of regression
    @params {array} dfRes
        Degrees of freedom of residuals

    @return {number} cumulative probability value
*/
function cumulativeProb(f, dfReg, dfRes){

    if(f > 0 & dfReg > 0 & dfRes > 0){
        let p = 1 - probability(f, dfReg, dfRes)
        return p
    } else{
        console.log('needs to be > 0, should be an error')
    }
}

module.exports = {
    cumulativeProb: cumulativeProb,
    probability: probability
}
