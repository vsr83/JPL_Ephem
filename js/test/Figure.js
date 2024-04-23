import {AssertionError, strict as assert} from 'assert';
import { checkFloat, checkFloatArray} from './common.js';
import { cross, norm, vecDiff} from '../src/MathUtils.js';
import { accBody, accOblateness, librationMoon} from '../src/Figure.js';
import { constants } from '../src/Constants.js';

describe('Figure', function() {
    it('accBody', function() {
        const rPoint = [0.943722189960001, 0.382748021970811, 0.027075526655788];
        const a = 1.161781241920150e-05;
        const mu = 1;
        const accExp = [0.475697396181525e-13, 0.318936432638165e-13, 0.038407095415640e-13];
        const acc = accBody(rPoint, a, mu, constants.Jm, constants.CSnm);
        const T = cross(rPoint, acc);
        const Texp = [0.006064867916584e-13, -0.023365870665248e-13, 0.118915151222176e-13];
        checkFloatArray(acc, accExp, 1.0e-28);
        checkFloatArray(T, Texp, 1.0e-28);
    });

    it('librationMoon', function() {
        const T = [-0.015729978061302e-16, -0.388800549855321e-16, -0.054437564119517e-16];
        const libMoon = {
            phi    : 0.005128132058714,
            phi1   : 1.165507165777481e-04,
            theta  : 0.382393200523007,
            theta1 : 1.461912823858170e-05,
            psi    : 1.294168056057082,
            psi1   : 0.229836728242082
        }

        librationMoon(libMoon, T);

        const phi2Exp = 8.388986888018661e-06;
        const theta2Exp = -9.300530241211267e-06;
        const psi2Exp = -7.885683461389150e-06;
        checkFloat(libMoon.phi2, phi2Exp, 2e-20);
        checkFloat(libMoon.theta2, theta2Exp, 2e-20);
        checkFloat(libMoon.psi2, psi2Exp, 2e-20);
    });

    it('accOblateness', function() {
        const stateInitial = constants.stateInitial;
        const OSV = {
            Sun   : stateInitial.objects[0],
            Earth : stateInitial.objects[3],
            Moon  : stateInitial.objects[4]
        };
        const libMoon = stateInitial.libration;
        const JT = stateInitial.JT;

        const accOut = accOblateness(OSV, libMoon, JT, constants.Je, constants.Jm, constants.CSnm);
        const expSun = [0.006693249200862e-20, -0.048757372219467e-20, -0.21549982209912e-20];
        const expEarth = [-0.001397684554884e-12, -0.006003096556495e-12, -0.865835515062385e-12];
        const expMoon = [0.001118207727728e-10, 0.005012531494318e-10, 0.704512523742023e-10];

        checkFloatArray(accOut.Sun, expSun, 1e-34);
        checkFloatArray(accOut.Earth, expEarth, 1e-22);
        checkFloatArray(accOut.Moon, expMoon, 1e-20);
        checkFloat(libMoon.phi2, -0.213971696871218e-5, 1e-15);
        checkFloat(libMoon.theta2, 0.245634047922518e-5, 1e-15);
        checkFloat(libMoon.psi2, -0.102242383380786e-5, 1e-15);
    });
});
