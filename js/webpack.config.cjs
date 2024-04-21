const path = require('path');

module.exports = {
  entry: './src/index.js',
  mode : 'production',
  output: {
    path: path.resolve(__dirname, 'dist'),
    filename: 'jplephem.js',
    library: {
        name: "jplephem",
        type: "umd",
    },
  },
};