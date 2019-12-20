

const linear = require('./src/regression/linear')
const discrete = require('./src/regression/discrete')

// suport funcitons
exports.addConstant = linear.addConstant

// Regressions
exports.Ols = linear.Ols
exports.Gls = linear.Gls
exports.Wls = linear.Wls
exports.Logit = discrete.Logit
