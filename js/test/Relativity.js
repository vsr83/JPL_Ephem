import {AssertionError, strict as assert} from 'assert';
import { checkFloat, checkFloatArray} from './common.js';
import { objectsInitial } from '../src/Integration.js';
import { accPointMass, barycenter } from '../src/Relativity.js';
import { norm, vecDiff} from '../src/MathUtils.js';

describe('Relativity', function() {
    it('barycenter', function() {
        const baryRel = barycenter(objectsInitial, true);
        const baryClassical = barycenter(objectsInitial, false);

        const au = 149597870.691;
        const mm = 1e-6 / au;
        const mmPerYear = mm / 365.25;

        // The initial state is w.r.t. relativistic barycenter so we verify that the
        // barycenter is at the origin at 1 mm precision and is moving less than 1 mm/d.
        checkFloat(norm(vecDiff(baryRel.r, [0, 0, 0])), 0, mm);
        checkFloat(norm(vecDiff(baryRel.v, [0, 0, 0])), 0, mmPerYear);

        checkFloat(norm(vecDiff(baryClassical.r, 
            [0.246396862761073e-12, 0.055122686802596e-12, 0.016489701714104e-12]
        )), 0, mm);
        checkFloat(norm(vecDiff(baryClassical.v, 
            [-0.063962645203393e-15, 0.284219346631659e-15, 0.123357968819851e-15])), 
            0, mmPerYear);
    });

    it('accPointMass', function() {
        const aNewtExp = 
        [
            [-0.000000007030480,  -0.000000002546937,  -0.000000001039335],
            [-0.001946635811530,   0.000498825346480,   0.000468489629444],
            [-0.000466355841907,   0.000267680540236,   0.000149924625681],
            [-0.000033317995592,   0.000259427905441,   0.000112322333977],
            [ 0.000018984536242,   0.000386349145341,   0.000181690771306],
            [ 0.000010794404372,   0.000125116490059,   0.000057093067063],
            [ 0.000009828288481,   0.000001517434253,   0.000000410914142],
            [-0.000002938971744,  -0.000001711679075,  -0.000000580434488],
            [ 0.000000882739155,   0.000000056108733,   0.000000012074041],
            [ 0.000000170613041,   0.000000254479431,   0.000000099914327],
            [ 0.000000281675618,   0.000000008050719,  -0.000000082378875]            
        ];
        const aRelExp = [
           [-0.000000024253668e-9,   0.000000069390127e-9,  0.000000036034055e-9],
           [ 0.145795913913519e-9,  -0.078334511308043e-9, -0.056973076302699e-9],
           [ 0.018985452799065e-9,  -0.010958291332444e-9, -0.006130905646534e-9],
           [ 0.001013271312639e-9,  -0.007552496802565e-9, -0.003264414995136e-9],
           [-0.002075686979579e-9,  -0.015173476028597e-9, -0.007447267191877e-9],
           [-0.000544215613018e-9,  -0.002503603330035e-9, -0.001133554561092e-9],
           [-0.000054231478083e-9,  -0.000008659334278e-9, -0.000002390288444e-9],
           [ 0.000009707271894e-9,   0.000004898197807e-9,  0.000001604957487e-9],
           [-0.000001400936842e-9,  -0.000000105835574e-9, -0.000000026640241e-9],
           [-0.000000166948284e-9,  -0.000000250978377e-9, -0.000000098599461e-9],
           [-0.000000250602345e-9,   0.000000058330022e-9,  0.000000094357851e-9]
        ];

        accPointMass(objectsInitial, true);

        for (let indTarget = 0; indTarget < objectsInitial.length; indTarget++) {
            const target = objectsInitial[indTarget];

            checkFloatArray(target.accNewton, aNewtExp[indTarget], 1e-15);
            checkFloatArray(target.accRel,    aRelExp[indTarget],  1e-23);
        }
    });
});
