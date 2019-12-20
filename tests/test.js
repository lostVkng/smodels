
const stats = require('../index.js')
const assert = require('assert')



describe('Ols', () => {

    // test data
    let td = {
        endog: {'title':'Y', 'data':[150.697, 179.323, 203.212, 226.505, 249.633, 281.422, 256.2, 231.2]},
        exog: [
            {'title':'Cubed HH Size', 'data': [0, 0.04, 0.16, 0.36, 0.64, 1.00, 0.8, 0.9]},
            {'title':'HH Size', 'data': [0, 0.2, 0.4, 0.6, 0.8, 1.0, 0.8, 0.9]}
        ]
    }

    let exog = stats.addConstant(td.exog)

    let Ols = new stats.Ols(td.endog, exog)
    Ols.fit()
    //Ols.summary()

    // predict test
    //let predIS = Ols.predict(Ols.X) // In Sample
    let predOS = Ols.predict({exog:[ [1, 0.12, 0.8] ]}) // Out of Sample

    it('OLS Summary matches expected values', () => {

        assert.equal(Ols.N, 8)
        assert.equal(Ols.K, 2)
        assert.equal(Ols.rsq, 0.9261362051477364)
        assert.equal(Ols.rsqAdj, 0.896590687206831)
        assert.equal(Ols.fvalue, 31.346081222881917)
        assert.equal(Ols.fProb, 0.0014827879255655651)
        assert.equal(Ols.llf, -30.448831403056587)
        assert.equal(Ols.aic, 66.89766280611317)
        assert.equal(Ols.bic, 67.13598743115269)

        // coeff values
        let exog = [
            { beta: 150.83975116848418,
               sterror: 12.14082304666748,
               tstat: 12.424178376431243,
                pvalue: 0.0000598802456448,
                ciLower: 119.63077198563445,
                ciHigher: 182.0487303513339,
                title: 'intercept'},
             { beta:  -24.992701146227464,
               sterror: 48.401002252046304,
               tstat:  -0.5163674300808683,
               pvalue: 1.3723849626541746,
               ciLower: -149.4114383197178,
               ciHigher: 99.42603602726288 ,
               title: 'Cubed HH Size'},
             { beta: 142.32883513242822,
               sterror: 54.2164643511603,
               tstat: 2.6251958115631386,
               pvalue: 0.04680352636320961,
               ciLower: 2.96097672619588,
               ciHigher: 281.69669353866055,
               title: 'HH Size'} ]

        // coeff values
        Ols.params.every((v, i) => assert.equal(v, exog[i]['beta']))
        Ols.bse.every((v, i) => assert.equal(v, exog[i]['sterror']))
        Ols.btstats.every((v, i) => assert.equal(v, exog[i]['tstat']))
        Ols.pvalues.every((v, i) => assert.equal(v, exog[i]['pvalue']))
        Ols.confInt.every((v, i) => assert.equal(v.lower, exog[i]['ciLower']) )
        Ols.confInt.every((v, i) => assert.equal(v.higher, exog[i]['ciHigher']) )

        assert.equal(Ols.omnibus.omnibus, 7.749169238768937)
        assert.equal(Ols.omnibus.pval, 0.020762960614076342)
        assert.equal(Ols.skew, -1.1871942084390812)
        assert.equal(Ols.kurtosis, 4.047992597075074)
        assert.equal(Ols.durbinWatson, 1.5946579573660176)
        assert.equal(Ols.jarqueBera.jb, 2.245336279243115)
        assert.equal(Ols.jarqueBera.pval, 0.32541071224140594)
        assert.equal(Ols.conditionNumber, 19.352085218320827)

        assert.equal(predOS[0], 261.7036951368795)
    })
})



describe('Gls', () => {

    // test data
    let td = {
        endog: {'title':'Y', 'data':[150.697, 179.323, 203.212, 226.505, 249.633, 281.422, 256.2, 231.2]},
        exog: [
            {'title':'Cubed HH Size', 'data': [0, 0.04, 0.16, 0.36, 0.64, 1.00, 0.8, 0.9]},
            {'title':'HH Size', 'data': [0, 0.2, 0.4, 0.6, 0.8, 1.0, 0.8, 0.9]}
        ]
    }

    let exog = stats.addConstant(td.exog)

    let Gls = new stats.Gls(td.endog, exog)
    Gls.fit()
    //Gls.summary()

    // predict test
    //let predIS = Gls.predict(Gls.X) // In Sample
    let predOS = Gls.predict({exog:[ [1, 0.12, 0.8] ]}) // Out of Sample

    it('GLS Summary matches expected values', () => {

        assert.equal(Gls.N, 8)
        assert.equal(Gls.K, 2)
        assert.equal(Gls.rsq, 0.9711200457416768)
        assert.equal(Gls.rsqAdj, 0.9595680640383475)
        assert.equal(Gls.fvalue, 84.06523405952066)
        assert.equal(Gls.fProb, 0.00014173961624692008)
        assert.equal(Gls.llf, -30.179846942678548)
        assert.equal(Gls.aic, 66.3596938853571)
        assert.equal(Gls.bic, 66.5980185103966)

        // coeff values
        let exog = [
            { beta: 150.2198920910882,
               sterror: 10.071303830742226,
               tstat: 14.915635017637774,
                pvalue: 0.0000245134300308,
                ciLower: 124.33078141139134,
                ciHigher: 176.10900277078505,
                title: 'intercept'},
             { beta:  -23.014464,
               sterror: 38.723897705357196,
               tstat:  -0.5943219858741988,
               pvalue: 0.578166,
               ciLower: -122.5574117976947,
               ciHigher: 76.52848422762025,
               title: 'Cubed HH Size'},
             { beta: 143.400106,
               sterror: 44.432182971672056,
               tstat: 3.2273927581835875,
               pvalue: 0.023273255737175536,
               ciLower: 29.183543127164036,
               ciHigher: 257.61666797896095,
               title: 'HH Size'} ]

        // coeff values
        Gls.params.every((v, i) => assert.equal(v, exog[i]['beta']))
        Gls.bse.every((v, i) => assert.equal(v, exog[i]['sterror']))
        Gls.btstats.every((v, i) => assert.equal(v, exog[i]['tstat']))
        Gls.pvalues.every((v, i) => assert.equal(v, exog[i]['pvalue']))
        Gls.confInt.every((v, i) => assert.equal(v.lower, exog[i]['ciLower']) )
        Gls.confInt.every((v, i) => assert.equal(v.higher, exog[i]['ciHigher']) )

        assert.equal(Gls.omnibus.omnibus, 7.679225789102609)
        assert.equal(Gls.omnibus.pval,0.021501923245905386)
        assert.equal(Gls.skew, -1.072756197420161)
        assert.equal(Gls.kurtosis, 4.250475009514842)
        assert.equal(Gls.durbinWatson, 1.1224726436511494)
        assert.equal(Gls.jarqueBera.jb, 2.055637061944866)
        assert.equal(Gls.jarqueBera.pval, 0.3577867532515172)
        assert.equal(Gls.conditionNumber, 22.348442867688227)

        assert.equal(predOS[0], 262.17824087933667)
    })
})



describe('Wls', () => {

    // test data
    let td = {
        endog: {'title':'Y', 'data':[150.697, 179.323, 203.212, 226.505, 249.633, 281.422, 256.2, 231.2]},
        exog: [
            {'title':'Cubed HH Size', 'data': [0, 0.04, 0.16, 0.36, 0.64, 1.00, 0.8, 0.9]},
            {'title':'HH Size', 'data': [0, 0.2, 0.4, 0.6, 0.8, 1.0, 0.8, 0.9]}
        ]
    }

    let exog = stats.addConstant(td.exog)
    let weights = [1,2,3,4,5,6,7,8]

    let Wls = new stats.Wls(td.endog, exog, {weights:weights})
    Wls.fit()
    //Wls.summary()

    // predict test
    //let predIS = Wls.predict(Wls.X) // In Sample
    let predOS = Wls.predict({exog:[ [1, 0.12, 0.8] ]}) // Out of Sample

    it('WLS Summary matches expected values', () => {

        assert.equal(Wls.N, 8)
        assert.equal(Wls.K, 2)
        assert.equal(Wls.rsq, 0.7895218377773362)
        assert.equal(Wls.rsqAdj, 0.7053305728882707)
        assert.equal(Wls.fvalue, 9.377716783536256)
        assert.equal(Wls.fProb, 0.0203243941089844)
        assert.equal(Wls.llf, -33.152110412410764)
        assert.equal(Wls.aic, 72.30422082482153)
        assert.equal(Wls.bic, 72.54254544986104)

        // coeff values
        let exog = [
            { beta: 150.23392385880598,
               sterror: 27.48341453638466,
               tstat: 5.4663485739704845,
                pvalue: 0.002789159420967,
                ciLower: 79.58555769428179,
                ciHigher: 220.88229002333017,
                title: 'intercept'},
             { beta:  -31.502500,
               sterror: 75.823808,
               tstat:  -0.415470,
               pvalue: 0.695026,
               ciLower: -226.414,
               ciHigher: 163.409,
               title: 'Cubed HH Size'},
             { beta: 147.475176,
               sterror: 96.884791,
               tstat: 1.522171,
               pvalue: 0.188458,
               ciLower: -101.575,
               ciHigher: 396.525,
               title: 'HH Size'} ]

        // coeff values
        Wls.params.every((v, i) => assert.equal(v, exog[i]['beta']))
        Wls.bse.every((v, i) => assert.equal(v, exog[i]['sterror']))
        Wls.btstats.every((v, i) => assert.equal(v, exog[i]['tstat']))
        Wls.pvalues.every((v, i) => assert.equal(v, exog[i]['pvalue']))
        Wls.confInt.every((v, i) => assert.equal(v.lower, exog[i]['ciLower']) )
        Wls.confInt.every((v, i) => assert.equal(v.higher, exog[i]['ciHigher']) )

        assert.equal(Wls.omnibus.omnibus, 6.231299630154271)
        assert.equal(Wls.omnibus.pval, 0.04434967869204043)
        assert.equal(Wls.skew, -1.0438436225361825)
        assert.equal(Wls.kurtosis, 3.821541637623891)
        assert.equal(Wls.durbinWatson, 1.626482559060597)
        assert.equal(Wls.jarqueBera.jb, 1.6777895651958619)
        assert.equal(Wls.jarqueBera.pval, 0.4321881884835118)
        assert.equal(Wls.conditionNumber, 28.667621889442824)

        assert.equal(predOS[0], 264.43376510811385)
    })
})



describe('Logit', () => {

    let td = {
        endog: {'title':'Y', 'data':[1, 1, 0, 1, 0, 0, 1, 1]},
        exog: [
            {'title':'Cubed HH Size', 'data': [0, 0.04, 0.16, 0.36, 0.64, 1.00, 0.8, 0.9]},
            {'title':'HH Size', 'data': [0, 0.2, 0.4, 0.6, 0.8, 1.0, 0.8, 0.9]}
        ]
    }

    let exog = stats.addConstant(td.exog)

    let model = new stats.Logit(td.endog, exog)
    model.fit()
    //model.summary()

    // predict test
    //let predIS = model.predict(model.X) // In Sample
    let predOS = model.predict({exog:[ [1, 0.12, 0.8] ]}) // Out of Sample

    it('Logit Summary matches expected values', () => {

        assert.equal(model.N, 8)
        assert.equal(model.K, 2)
        assert.equal(model.pseudoRsq, 0.2692722699709513)
        assert.equal(model.llf, -3.867380826318793)
        assert.equal(model.llnull, -5.292505905263856)
        assert.equal(model.llrPvalue, 0.2404787187667813)

        // coeff values
        let exog = [
            { beta: 6.818155827165296,
               sterror: 6.1753483876846,
               zstat: 1.1040924979654,
                pvalue: 0.26955303293986477,
                ciLower: -5.285304373891712,
                ciHigher: 18.921616028222303,
                title: 'intercept'},
             { beta:  15.34894987454284,
               sterror: 13.858924491809349,
               zstat:  1.1075137817234013,
               pvalue: 0.2680718747835029,
               ciLower: -11.814042475911718,
               ciHigher: 42.5119422249974 ,
               title: 'Cubed HH Size'},
             { beta: -22.210055025922774,
               sterror: 19.634013406200417,
               zstat: -1.1312030080874267,
               pvalue: 0.2579696736255528,
               ciLower: -60.692013440267004,
               ciHigher: 16.27190338842146,
               title: 'HH Size'} ]

        // coeff values
        model.params.every((v, i) => assert.equal(v, exog[i]['beta']))
        model.bse.every((v, i) => assert.equal(v, exog[i]['sterror']))
        model.bzstats.every((v, i) => assert.equal(v, exog[i]['zstat']))
        model.pvalues.every((v, i) => assert.equal(v, exog[i]['pvalue']))
        model.confInt.every((v, i) => assert.equal(v.lower, exog[i]['ciLower']) )
        model.confInt.every((v, i) => assert.equal(v.higher, exog[i]['ciHigher']) )

        assert.equal(predOS[0], 0.00011076250965810184)
    })
})



// end
