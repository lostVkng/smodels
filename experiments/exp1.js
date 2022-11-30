// require it
const stats = require("../index.js");

// sample data
// let tdExample = {
//   endog: {
//     title: "Y",
//     data: [150.697, 179.323, 203.212, 226.505, 249.633, 281.422, 256.2, 231.2],
//   },
//   exog: [
//     {
//       title: "Cubed HH Size",
//       data: [0, 0.04, 0.16, 0.36, 0.64, 1.0, 0.8, 0.9],
//     },
//     { title: "HH Size", data: [0, 0.2, 0.4, 0.6, 0.8, 1.0, 0.8, 0.9] },
//   ],
// };

let tdJamesLimitEight = {
  endog: {
    title: "Energy (kWh)",
    data: [129211, 60470, 68314, 89671, 82607, 108221, 100673, 115006],
  },
  exog: [
    {
      title: "OAT",
      data: [37, 6, 9, 17, 14, 24, 21, 27],
    },
    {
      title: "Week/Weekend",
      data: [0, 1, 1, 0, 0, 0, 0, 0],
    },
  ],
};

// let tdJamesFull = {
//   endog: {
//     title: "Energy (kWh)",
//     data: [129211, 60470, 68314, 89671, 82607, 108221, 100673, 115006, 109119],
//   },
//   exog: [
//     {
//       title: "OAT",
//       data: [37, 6, 9, 17, 14, 24, 21, 27, 29],
//     },
//     {
//       title: "Week/Weekend",
//       data: [0, 1, 1, 0, 0, 0, 0, 0, 1],
//     },
//   ],
// };
// create and fit model
// let model = new stats.Ols(tdJames.endog, tdJames.exog).fit();
// let model = new stats.Ols(td.endog, td.exog).fit();

// create and fit model
let model = new stats.Gls(
  tdJamesLimitEight.endog,
  tdJamesLimitEight.exog
).fit();

console.log(model);
