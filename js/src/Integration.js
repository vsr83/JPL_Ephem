import { coordJ2000Mod, coordModTod, coordTodMod } from "./Frames.js";
import {cosd, sind, norm, linComb, vecSum, vecMul, cross} from "./MathUtils.js";
import { rotateCart3 } from "./Rotations.js";
import { constants } from "./Constants.js";
import { accPointMass } from "./Relativity.js";
import { accOblateness } from "./Figure.js";

/**
 * This method performs a single integration step for the initial value
 * problem
 *    dy/dt = f(t, y) 
 *    y(t_in) = y_in * 
 * 
 * @param {*} funcIn 
 *      Function handle for f(t, y).
 * @param {*} tIn 
 *      Time before the step.
 * @param {*} yIn 
 *      DoFs before the step.
 * @param {*} h 
 *      Step size.
 * @param {*} param
 *      Parameter passed to the method funcIn. 
 * @returns Object with yOut and tOut fields for the DoFs and time after the step.
 */
export function runge4(funcIn, tIn, yIn, h, param) {
    const k1 = funcIn(tIn, yIn, param);
    const k2 = funcIn(tIn + h/2, linComb([1, h/2], [yIn, k1]), param);
    const k3 = funcIn(tIn + h/2, linComb([1, h/2], [yIn, k2]), param);
    const k4 = funcIn(tIn + h,   linComb([1, h],   [yIn, k3]), param);

    const yOut = linComb([1, h/6, h/3, h/3, h/6], [yIn, k1, k2, k3, k4]);
    const tOut = tIn + h;

    return {yOut : yOut, tOut : tOut};
}

export function adams8(funcIn, tIn, yIn, Fin, h, param) {

    const predCof = [
         h *  434241/120960, 
        -h * 1152169/120960,
         h * 2183877/120960,
        -h * 2664477/120960,
         h * 2102243/120960,
        -h * 1041723/120960,
         h *  295767/120960,
        -h *   36799/120960
    ];
    const corrCof = [
         h *  36799/120960,
         h * 139849/120960,
        -h * 121797/120960,
         h * 123133/120960,
        -h *  88547/120960,
         h *  41499/120960,
        -h *  11351/120960,
         h *   1375/120960
    ];

    const numDof = yIn.length;
    // Predictor step.
    let yNew = yIn.slice();
    let yOut = yIn.slice();

    for (let indCof = 0; indCof < predCof.length; indCof++) {
        const coeff = predCof[indCof];

        for (let indDof = 0; indDof < numDof; indDof++) {
            yNew[indDof] += coeff * Fin[indCof][indDof];
        }
    }

    let f1 = funcIn(tIn + h, yNew, param);
    // Corrector step.
    const Ftmp = [f1].concat(Fin).slice(0, 8);

    for (let indCof = 0; indCof < corrCof.length; indCof++) {
        const coeff = corrCof[indCof];

        for (let indDof = 0; indDof < numDof; indDof++) {
            yOut[indDof] += coeff * Ftmp[indCof][indDof];
        }
    }

    f1 = funcIn(tIn + h, yOut, param);
    const Fout = [f1].concat(Fin).slice(0, 8);

    const tOut = tIn + h;

    return {tOut : tOut, yOut : yOut, Fout : Fout};
}

/**
 * 
 * 
 * @param {*} state 
 */
export function stateToDof(state) {
    const {JT, objects, libration} = state;
    const numObjects = state.objects.length;
    const dof = [];

    dof.push(libration.phi);
    dof.push(libration.phi1);
    dof.push(libration.theta);
    dof.push(libration.theta1);
    dof.push(libration.psi);
    dof.push(libration.psi1);

    for (let indObject = 0; indObject < objects.length; indObject++) {
        const {r, v, mu} = objects[indObject];

        dof.push(r[0]);
        dof.push(r[1]);
        dof.push(r[2]);
        dof.push(v[0]);
        dof.push(v[1]);
        dof.push(v[2]);
    }

    return dof;
}

/**
 * Compute the vector f in the equation y' = f(y, t).
 * 
 * @param {*} tIn 
 *      Time after epoch.
 * @param {*} yNew 
 *      Degrees of freedom.
 * @param {*} JT 
 *      Julian time.
 * @returns 
 */
export function funcLibration(tIn, yNew, JT) {
    let state = dofToState(yNew, tIn, constants.stateInitial, JT);    
    accPointMass(state, true);
    accOblateness(state);

    const acc = [];
    const libration = state.libration;
    acc.push(libration.phi1);
    acc.push(libration.phi2);
    acc.push(libration.theta1);
    acc.push(libration.theta2);
    acc.push(libration.psi1);
    acc.push(libration.psi2);

    const objects = state.objects;
    for (let indObject = 0; indObject < objects.length; indObject++) {
        const indDof = 6 + indObject * 6;

        const object = objects[indObject];
        if ('accObl' in object) {
            object.acc = linComb([1, 1, 1], 
                [object.accNewton, object.accRel, object.accObl]);
        } else {
            object.acc = linComb([1, 1], 
                [object.accNewton, object.accRel]);            
        }
        acc.push(yNew[indDof + 3]);
        acc.push(yNew[indDof + 4]);
        acc.push(yNew[indDof + 5]);
        acc.push(object.acc[0]);
        acc.push(object.acc[1]);
        acc.push(object.acc[2]);
    }

    return acc;
}

/**
 * Convert degrees of freedom to a state.
 * 
 * @param {*} dof 
 *      The degrees of freedom.
 * @param {*} tIn 
 *      Time after epoch.
 * @param {*} stateOld 
 *      Old state object used to fill remaining fields.
 * @returns The new state.
 */
export function dofToState(dof, tIn, stateOld) {
    const objectsOld = stateOld.objects;
    const numObjects = objectsOld.length;

    const libration = {
        phi    : dof[0],
        phi1   : dof[1],
        theta  : dof[2],
        theta1 : dof[3],
        psi    : dof[4],
        psi1   : dof[5],
    }
    const objects = [];

    for (let indObject = 0; indObject < numObjects; indObject++) {
        const indDof = 6 + indObject * 6;
        const name = objectsOld[indObject].name;
        const mu   = objectsOld[indObject].mu;
        const r    = [dof[indDof], dof[indDof + 1], dof[indDof + 2]]; 
        const v    = [dof[indDof + 3], dof[indDof + 4], dof[indDof + 5]]; 

        objects.push({name : name, mu : mu, r : r, v : v});
    }

    return {JT : stateOld.JT, t : tIn, objects : objects, objectIndices : stateOld.objectIndices,
         libration : libration};
}