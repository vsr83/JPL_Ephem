import {AssertionError, strict as assert} from 'assert';
import { checkFloat, checkFloatArray} from './common.js';
import { legendreDeriv, legendreValue, legendreAssoc, legendreAssocd } from '../src/Legendre.js';

describe('Legendre', function() {
    it('legendreValue', function() {
        checkFloat(legendreValue(0, -1.0), 1.0, 1e-10);
        checkFloat(legendreValue(0, -0.5), 1.0, 1e-10);
        checkFloat(legendreValue(0,  0.0), 1.0, 1e-10);
        checkFloat(legendreValue(0,  0.5), 1.0, 1e-10);
        checkFloat(legendreValue(0,  1.0), 1.0, 1e-10);

        checkFloat(legendreValue(1, -1.0), -1.0, 1e-10);
        checkFloat(legendreValue(1, -0.5), -0.5, 1e-10);
        checkFloat(legendreValue(1,  0.0),  0.0, 1e-10);
        checkFloat(legendreValue(1,  0.5),  0.5, 1e-10);
        checkFloat(legendreValue(1,  1.0),  1.0, 1e-10);

        checkFloat(legendreValue(2, -1.0),    1.0, 1e-10);
        checkFloat(legendreValue(2, -0.5), -0.125, 1e-10);
        checkFloat(legendreValue(2,  0.0),   -0.5, 1e-10);
        checkFloat(legendreValue(2,  0.5), -0.125, 1e-10);
        checkFloat(legendreValue(2,  1.0),    1.0, 1e-10);

        checkFloat(legendreValue(3, -1.0),    -1.0, 1e-10);
        checkFloat(legendreValue(3, -0.5),  0.4375, 1e-10);
        checkFloat(legendreValue(3,  0.0),     0.0, 1e-10);
        checkFloat(legendreValue(3,  0.5), -0.4375, 1e-10);
        checkFloat(legendreValue(3,  1.0),     1.0, 1e-10);

        checkFloat(legendreValue(4, -1.0),        1.0, 1e-10);
        checkFloat(legendreValue(4, -0.5), -0.2890625, 1e-10);
        checkFloat(legendreValue(4,  0.0),      0.375, 1e-10);
        checkFloat(legendreValue(4,  0.5), -0.2890625, 1e-10);
        checkFloat(legendreValue(4,  1.0),        1.0, 1e-10);

        checkFloat(legendreValue(5, -1.0),        -1.0, 1e-10);
        checkFloat(legendreValue(5, -0.5), -0.08984375, 1e-10);
        checkFloat(legendreValue(5,  0.0),         0.0, 1e-10);
        checkFloat(legendreValue(5,  0.5),  0.08984375, 1e-10);
        checkFloat(legendreValue(5,  1.0),         1.0, 1e-10);

        checkFloat(legendreValue(6, -1.0),          1.0, 1e-10);
        checkFloat(legendreValue(6, -0.5), 0.3232421875, 1e-10);
        checkFloat(legendreValue(6,  0.0),      -0.3125, 1e-10);
        checkFloat(legendreValue(6,  0.5), 0.3232421875, 1e-10);
        checkFloat(legendreValue(6,  1.0),          1.0, 1e-10);
    });
    it('legendreDeriv', function() {
        checkFloat(legendreDeriv(0, -1.0), 0.0, 1e-10);
        checkFloat(legendreDeriv(0, -0.5), 0.0, 1e-10);
        checkFloat(legendreDeriv(0,  0.0), 0.0, 1e-10);
        checkFloat(legendreDeriv(0,  0.5), 0.0, 1e-10);
        checkFloat(legendreDeriv(0,  1.0), 0.0, 1e-10);

        checkFloat(legendreDeriv(1, -1.0), 1.0, 1e-10);
        checkFloat(legendreDeriv(1, -0.5), 1.0, 1e-10);
        checkFloat(legendreDeriv(1,  0.0), 1.0, 1e-10);
        checkFloat(legendreDeriv(1,  0.5), 1.0, 1e-10);
        checkFloat(legendreDeriv(1,  1.0), 1.0, 1e-10);

        checkFloat(legendreDeriv(2, -1.0), -3.0, 1e-10);
        checkFloat(legendreDeriv(2, -0.5), -1.5, 1e-10);
        checkFloat(legendreDeriv(2,  0.0),  0.0, 1e-10);
        checkFloat(legendreDeriv(2,  0.5),  1.5, 1e-10);
        checkFloat(legendreDeriv(2,  1.0),  3.0, 1e-10);

        checkFloat(legendreDeriv(3, -1.0),   6.0, 1e-10);
        checkFloat(legendreDeriv(3, -0.5), 0.375, 1e-10);
        checkFloat(legendreDeriv(3,  0.0),  -1.5, 1e-10);
        checkFloat(legendreDeriv(3,  0.5), 0.375, 1e-10);
        checkFloat(legendreDeriv(3,  1.0),   6.0, 1e-10);

        checkFloat(legendreDeriv(4, -1.0),   -10.0, 1e-10);
        checkFloat(legendreDeriv(4, -0.5),  1.5625, 1e-10);
        checkFloat(legendreDeriv(4,  0.0),       0, 1e-10);
        checkFloat(legendreDeriv(4,  0.5), -1.5625, 1e-10);
        checkFloat(legendreDeriv(4,  1.0),    10.0, 1e-10);

        checkFloat(legendreDeriv(5, -1.0),       15.0, 1e-10);
        checkFloat(legendreDeriv(5, -0.5), -2.2265625, 1e-10);
        checkFloat(legendreDeriv(5,  0.0),      1.875, 1e-10);
        checkFloat(legendreDeriv(5,  0.5), -2.2265625, 1e-10);
        checkFloat(legendreDeriv(5,  1.0),       15.0, 1e-10);

        checkFloat(legendreDeriv(6, -1.0),       -21.0, 1e-10);
        checkFloat(legendreDeriv(6, -0.5),  0.57421875, 1e-10);
        checkFloat(legendreDeriv(6,  0.0),         0.0, 1e-10);
        checkFloat(legendreDeriv(6,  0.5), -0.57421875, 1e-10);
        checkFloat(legendreDeriv(6,  1.0),        21.0, 1e-10); 
    });

    it('legendreAssoc', function() {
        checkFloat(legendreAssoc(1, 0, -1.0), -1.0, 1e-10);
        checkFloat(legendreAssoc(1, 0, -0.5), -0.5, 1e-10);
        checkFloat(legendreAssoc(1, 0,  0.0),  0.0, 1e-10);
        checkFloat(legendreAssoc(1, 0,  0.5),  0.5, 1e-10);
        checkFloat(legendreAssoc(1, 0,  1.0),  1.0, 1e-10);

        checkFloat(legendreAssoc(1, 1, -1.0),  0.0, 1e-10);
        checkFloat(legendreAssoc(1, 1, -0.5), -0.866025403784439, 1e-10);
        checkFloat(legendreAssoc(1, 1,  0.0), -1.0, 1e-10);
        checkFloat(legendreAssoc(1, 1,  0.5), -0.866025403784439, 1e-10);
        checkFloat(legendreAssoc(1, 1,  1.0),  0.0, 1e-10);

        checkFloat(legendreAssoc(2, 0,  1.0),    1.0, 1e-10);
        checkFloat(legendreAssoc(2, 0, -0.5), -0.125, 1e-10);
        checkFloat(legendreAssoc(2, 0,  0.0),   -0.5, 1e-10);
        checkFloat(legendreAssoc(2, 0,  0.5), -0.125, 1e-10);
        checkFloat(legendreAssoc(2, 0,  1.0),    1.0, 1e-10);

        checkFloat(legendreAssoc(2, 1, -1.0),    0.0, 1e-10);
        checkFloat(legendreAssoc(2, 1, -0.5), 1.299038105676658, 1e-10);
        checkFloat(legendreAssoc(2, 1,  0.0),    0.0, 1e-10);
        checkFloat(legendreAssoc(2, 1,  0.5), -1.299038105676658, 1e-10);
        checkFloat(legendreAssoc(2, 1,  1.0),    0.0, 1e-10);

        checkFloat(legendreAssoc(2, 2, -1.0),    0.0, 1e-10);
        checkFloat(legendreAssoc(2, 2, -0.5),   2.25, 1e-10);
        checkFloat(legendreAssoc(2, 2,  0.0),   3.00, 1e-10);
        checkFloat(legendreAssoc(2, 2,  0.5),   2.25, 1e-10);
        checkFloat(legendreAssoc(2, 2,  1.0),    0.0, 1e-10);

        checkFloat(legendreAssoc(3, 0, -1.0),    -1.0, 1e-10);
        checkFloat(legendreAssoc(3, 0, -0.5),  0.4375, 1e-10);
        checkFloat(legendreAssoc(3, 0,  0.0),     0.0, 1e-10);
        checkFloat(legendreAssoc(3, 0,  0.5), -0.4375, 1e-10);
        checkFloat(legendreAssoc(3, 0,  1.0),     1.0, 1e-10);

        checkFloat(legendreAssoc(3, 1, -1.0),     0.0, 1e-10);
        checkFloat(legendreAssoc(3, 1, -0.5), -0.324759526419165, 1e-10);
        checkFloat(legendreAssoc(3, 1,  0.0),     1.5, 1e-10);
        checkFloat(legendreAssoc(3, 1,  0.5), -0.324759526419165, 1e-10);
        checkFloat(legendreAssoc(3, 1,  1.0),     0.0, 1e-10);

        checkFloat(legendreAssoc(3, 2, -1.0),    0.0, 1e-10);
        checkFloat(legendreAssoc(3, 2, -0.5), -5.625, 1e-10);
        checkFloat(legendreAssoc(3, 2,  0.0),    0.0, 1e-10);
        checkFloat(legendreAssoc(3, 2,  0.5),  5.625, 1e-10);
        checkFloat(legendreAssoc(3, 2,  1.0),    0.0, 1e-10);

        checkFloat(legendreAssoc(3, 3, -1.0),     0.0, 1e-10);
        checkFloat(legendreAssoc(3, 3, -0.5), -9.742785792574935, 1e-10);
        checkFloat(legendreAssoc(3, 3,  0.0),   -15.0, 1e-10);
        checkFloat(legendreAssoc(3, 3,  0.5), -9.742785792574935, 1e-10);
        checkFloat(legendreAssoc(3, 3,  1.0),     0.0, 1e-10);

        checkFloat(legendreAssoc(4, 0, -1.0),        1.0, 1e-10);
        checkFloat(legendreAssoc(4, 0, -0.5), -0.2890625, 1e-10);
        checkFloat(legendreAssoc(4, 0,  0.0),      0.375, 1e-10);
        checkFloat(legendreAssoc(4, 0,  0.5), -0.2890625, 1e-10);
        checkFloat(legendreAssoc(4, 0,  1.0),        1.0, 1e-10);

        checkFloat(legendreAssoc(4, 1, -1.0),        0.0, 1e-10);
        checkFloat(legendreAssoc(4, 1, -0.5), -1.353164693413186, 1e-10);
        checkFloat(legendreAssoc(4, 1,  0.0),        0.0, 1e-10);
        checkFloat(legendreAssoc(4, 1,  0.5), 1.353164693413186, 1e-10);
        checkFloat(legendreAssoc(4, 1,  1.0),        0.0, 1e-10);

        checkFloat(legendreAssoc(4, 2, -1.0),     0.0, 1e-10);
        checkFloat(legendreAssoc(4, 2, -0.5), 4.21875, 1e-10);
        checkFloat(legendreAssoc(4, 2,  0.0),    -7.5, 1e-10);
        checkFloat(legendreAssoc(4, 2,  0.5), 4.21875, 1e-10);
        checkFloat(legendreAssoc(4, 2,  1.0),     0.0, 1e-10);

        checkFloat(legendreAssoc(4, 3, -1.0),        0.0, 1e-10);
        checkFloat(legendreAssoc(4, 3, -0.5), 34.099750274012273, 1e-10);
        checkFloat(legendreAssoc(4, 3,  0.0),        0.0, 1e-10);
        checkFloat(legendreAssoc(4, 3,  0.5),-34.099750274012273, 1e-10);
        checkFloat(legendreAssoc(4, 3,  1.0),        0.0, 1e-10);

        checkFloat(legendreAssoc(4, 4, -1.0),      0.0, 1e-10);
        checkFloat(legendreAssoc(4, 4, -0.5),  59.0625, 1e-10);
        checkFloat(legendreAssoc(4, 4,  0.0),    105.0, 1e-10);
        checkFloat(legendreAssoc(4, 4,  0.5),  59.0625, 1e-10);
        checkFloat(legendreAssoc(4, 4,  1.0),      0.0, 1e-10);
    });

    it('legendreAssocd', function() {
        const delta = 1e-9;
        for (let degree = 1; degree < 5; degree++) {
            for (let order = 0; order <= degree; order++) {
                for (let value = -0.99; value <= 0.99; value+= 0.1) {
                    const derivExp = (legendreAssoc(degree, order, value + delta/2) 
                                   - legendreAssoc(degree, order, value - delta/2)) / delta;
                    const deriv = legendreAssocd(degree, order, value);
                    checkFloat(deriv, derivExp, 1e-4);
                }
            }
        }
    });
});
