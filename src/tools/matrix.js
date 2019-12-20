/*
    Matrix related functions
*/


const Decimal = require('Decimal.js')
const Arr = require('./array.js')


/*
    @desc Matrix transposition

    @params {array} matrix

    @return {array} the transposed matrix
*/
function transpose(mtrx){

    // create empty matrix of proper dimensions
    let newMtrx = new Array(mtrx[0].length).fill(0).map(i => new Array(mtrx.length))

    // populate new matrix with data
    mtrx.forEach((col, i) => {
        col.forEach((val, rowI) => {
            newMtrx[rowI][i] = val
        })
    })
    return newMtrx
}

/*
    @desc Householder transformation

    @params {array} A
        A is a matrix of big numbers

    @return {object} Q & R factorisation of A
    @return {object.Q} Q of A
    @return {object.R} R of A
*/
function householder(A){

    let R = A.slice()
    let Q = null

    // Householder Rotation Process
    for(let i=0; i<A.length; i++){

        // sub matrices represent offseting from top row
        let Asub = A.slice(i,A.length).map((col) => col.slice(i,col.length))

        // X is the column householder is based around
        let X = Asub.map(a => a[0]).slice()

        // Euclidean Norm
        let xnorm = X.reduce((sum, val) => sum.plus(Decimal(val).times(val)), Decimal(0)).sqrt().toNumber()

        // e term
        let e = new Array(X.length).fill(0)
        e[0] = 1.0

        // Householder vector
        let xe = e.map(_x => Decimal(xnorm).times(_x).toNumber())
        let v = X.map((_x, _i) => Decimal(_x).plus(Decimal(xe[_i]).times(Decimal(X[0]).s)).toNumber())

        // Householder reflection
        let c = Decimal(2).dividedBy(v.reduce((sum, val) => sum.plus(Decimal(val).pow(2)), Decimal(0))).toNumber()
        let vvt = Arr.dot(transpose([v]), [v])
        let cvvt = vvt.map(col => col.map(v => Decimal(c).times(v).toNumber()))

        let I = new Array(vvt.length).fill(0).map(_i => new Array(vvt.length).fill(0))
        I.map((col, _i) => col[_i] = 1)
        let H = cvvt.map((col, colI) => col.map((v, rowI) => Decimal(I[colI][rowI]).minus(v).toNumber()))

        let HA = Arr.dot(H, Asub)

        // zero out HA column
        HA = HA.map((row, _i) => {
            if(_i > 0) row[0] = 0
            return row
        })

        // Update Q, R
        for(let _haR=0; _haR<HA.length; _haR++){
            // ha Rows
            for(let _haC=0; _haC<HA[_haR].length; _haC++){
                // update R values
                R[_haR + i][_haC + i] = HA[_haR][_haC]
            }
        }

        if(Q){
            let ih1 = new Array(A[0].length).fill(0).map(_i => new Array(A[0].length).fill(0))
            ih1.map((col, _i) => col[_i] = 1)

            for(let _hC=0; _hC<H.length; _hC++){
                // HA Rows
                for(let _hR=0; _hR<H[_hC].length; _hR++){
                    // update R values
                    ih1[_hC + i][_hR + i] = H[_hC][_hR]
                }
            }

            Q = Arr.dot(Q, ih1)
        }else{
            Q = H
        }

    }

    return {Q:Q, R:R}
}


/*
    @desc Solve varaibles via backward substitution

    @params {array} X
        X * B values
    @params {array} Y
        Y values

    @return {array} Array of solved values
*/
function backsolve(X, Y){

    let betas = new Array(X.length).fill(0)

    for(let _row=X.length-1; _row>-1; _row--){

        let yval = Decimal(Y[_row][0])

        for(let _col=X[_row].length-1; _col>-1; _col--){

            if(_col === _row){
                yval = yval.dividedBy(X[_row][_col])
                break
            }else{
                yval = yval.minus( Decimal(X[_row][_col]).times(betas[_col]) )
            }
        }
        betas[_row] = yval.toNumber()
    }

    return betas

}


/*
    @desc Calculates the inverse upper triangle of a matrix

    @params {array} mtrx
        Matrix to get the upper triangle from

    @return {array} Inverted upper triangle
*/
function invUpperTriangle(mtrx){

    let invMtrx = mtrx.slice()

    // create identity matrix
    let I = new Array(mtrx.length).fill(0).map(_i => new Array(mtrx.length).fill(0))
    I.map((col, _i) => col[_i] = 1)

    //loop through rows and columns
    for(let row=mtrx.length-1; row>=0; row--){
        for(let col=mtrx[row].length-1; col>=0; col--){

            // matrix cell to calculate
            let Ival = I[row][col]

            if(Ival === 1){
                invMtrx[row][col] = Decimal(1).dividedBy(mtrx[row][col]).toNumber()
            }else if(mtrx[row][col] === 0){
                continue
            }else{
                let sum = Decimal(0)

                // loop backwards in row
                for(let _c=mtrx[row].length-1; _c>=0; _c--){

                    if(_c === row){
                        // Only sum the numbers after the value being calculated
                        break
                    } else{

                        sum = sum.plus( ( Decimal(mtrx[row][_c]).times(invMtrx[_c][col]) ) )
                    }
                }

                sum = Decimal(Ival).minus(sum)
                invMtrx[row][col] = sum.dividedBy(mtrx[row][row]).toNumber()
            }
        }
    }
    return invMtrx

    // loop through rows in column
    for(let col=mtrx.length-1; col>=0; col--){

        for(let row=mtrx[col].length-1; row>=0; row--){

            // matrix cell to calculate
            let Ival = I[col][row]

            if(Ival === 1){
                invMtrx[col][row] = Decimal(1).dividedBy(mtrx[col][row]).toNumber()
            }else if(mtrx[col][row].isEqualTo(0)){
                continue
            }else{
                let sum = Decimal(0)

                for(let _c=mtrx.length-1; _c>=0; _c--){
                    if(_c === row){
                        // Only sum the numbers after the value being calculated
                        break
                    } else{
                        sum = sum.plus( ( mtrx[_c][row].times(invMtrx[col][_c]) ) )
                    }
                }

                sum = Decimal(Ival).minus(sum)
                invMtrx[col][row] = sum.dividedBy(mtrx[row][row]).toNumber()
            }
        }
    }
    return invMtrx
}






module.exports = {
    transpose: transpose,
    householder: householder,
    backsolve: backsolve,
    invUpperTriangle: invUpperTriangle,
}
