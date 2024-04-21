import { linComb, norm, vecDiff, vecMul, vecSum, dot } from "./MathUtils.js";

// Speed of light (au/d).
const c = 173.144632720536344565;
const c2 = c * c;

/**
 * Compute the classical or relativistic barycenter for the given point 
 * masses.
 * 
 * REFERENCES: 
 * [1] Newhall, Standish, Williams - DE 102: a numerically integrated
 * ephemeris of the Moon and planets spanning forty-four centuries,
 * Astronomy and Astrophysics, 125, 150-167, 1983.
 * [2] Urban, Seidelmann - Explanatory Supplement to the Astronomical
 * Almanac, 3rd edition, University Science Books, 2013.
 * 
 * @param {*} state 
 *      State of all objects with fields mu, r and v.
 * @param {*} relavistic 
 *      Flag indicating whether to compute relativistic instead of classical 
 *      barycenter.
 * @returns Object with fields r and v for the position and velocity of the 
 *          barycenter.
 */
export function barycenter(state, relavistic) {
    const numTargets = state.length;

    // Relativistic or classical standard gravitational parameter for the 
    // evaluation of the barycenter. 
    const muStar = [];
    let muStarSum = 0;

    for (let indTarget = 0; indTarget < numTargets; indTarget++) {
        const target = state[indTarget];

        if (relavistic) {
            let tmp = Math.pow(norm(target.v), 2.0);
        
            for (let indSource = 0; indSource < numTargets; indSource++) {
                if (indSource == indTarget) {
                    continue;
                }
                const source = state[indSource];
                tmp = tmp - source.mu / norm(vecDiff(target.r, source.r));
            }
            tmp = 1 - tmp * 0.5 / c2;
            muStar.push(target.mu * tmp);
        } else {
            muStar.push(target.mu);
        }
        muStarSum += muStar[muStar.length - 1];
    }

    // Compute the part of the equation for barycenter not including the
    // contribution from the Sun.
    let rBary = [0, 0, 0];
    let vBary = [0, 0, 0]; 

    for (let indTarget = 0; indTarget < numTargets; indTarget++) {
        const target = state[indTarget];
        rBary = linComb([1, muStar[indTarget]], [rBary, target.r]);
        vBary = linComb([1, muStar[indTarget]], [vBary, target.v]);
    }

    rBary = vecMul(rBary, 1 / muStarSum);
    vBary = vecMul(vBary, 1 / muStarSum);

    return {r : rBary, v : vBary};
}

/**
 * Compute the relativistic and Newtonian parts of acceleration.
 * 
 * This routine computes the relativistic and Newtonian parts of the
 * acceleration for an arbitrary number of point-masses. The total
 * accelerations are obtained as the sum of the two matrices. The
 * computation is based on the equation (1) in [1] and (8.1) in [2].
 *
 * REFERENCES: 
 *  [1] Newhall, Standish, Williams - DE 102: a numerically integrated
 *  ephemeris of the Moon and planets spanning forty-four centuries,
 *  Astronomy and Astrophysics, 125, 150-167, 1983.
 *  [2] Urban, Seidelmann - Explanatory Supplement to the Astronomical
 *  Almanac, 3rd edition, University Science Books, 2013.
 *  [3] Steve Moshier, DE118i available at 
 *  http://www.moshier.net/de118i-2.zip
 * 
 * @param {*} state 
 *      Input state with fields r, v and mu for each point mass. The
 *      results will be updated to the fields accNewton and accRel.
 * @param {*} relavistic 
 *      Flag indicating whether to compute relativistic parts of the
 *      acceleration.
 */
export function accPointMass(state, relativistic) {
    const numTargets = state.length;

    // Compute distances and third powers of distances between every pair
    // of objects.
    const Rij = []
    const Rij3 = []
    for (let indTarget = 0; indTarget < numTargets; indTarget++) {
        const target = state[indTarget];
        const RijRow = [];
        const RijRow3 = [];
        for (let indSource = 0; indSource < numTargets; indSource++) {
            const source = state[indSource];
            const distance = norm(vecDiff(target.r, source.r));

            RijRow.push(distance);
            RijRow3.push(Math.pow(distance, 3.0));
        }
        Rij.push(RijRow);
        Rij3.push(RijRow3);
    }

    // For numerical accuracy, it is very important to compute the relativistic 
    // acceleration separately. Otherwise, one has to do large amount of 
    // floating computations that involve adding small numbers to larger ones.

    // Compute the Newtonian accelerations first.
    for (let indTarget = 0; indTarget < numTargets; indTarget++) {
        const target = state[indTarget];
        let accSum = [0, 0, 0];

        for (let indSource = 0; indSource < numTargets; indSource++) {
            if (indSource == indTarget) {
                continue;
            }
            const source = state[indSource];

            const rTargetSource = vecDiff(source.r, target.r);
            const accNewton = vecMul(rTargetSource, source.mu / Rij3[indSource][indTarget]);
            accSum = vecSum(accSum, accNewton);
        }

        target.accNewton = accSum;
    }

    // Compute relativistic accelerations.
    for (let indTarget = 0; indTarget < numTargets; indTarget++) {
        const target = state[indTarget];

        target.accRel = [0, 0, 0];
        if (!relativistic) {
            continue;
        }
    
        for (let indSource = 0; indSource < numTargets; indSource++) {
            if (indSource == indTarget) {
                continue;
            }
            const source = state[indSource];

            const rTargetSource = vecDiff(source.r, target.r);
            const vTargetSource = vecDiff(source.v, target.v);
            const accNewton = vecMul(rTargetSource, source.mu / Rij3[indSource][indTarget]);

            // The first part of the acceleration formula involves
            // multiplication of the Newtonian acceleration.
            let newtonMult = 0.0;
            for (let indTarget2 = 0; indTarget2 < numTargets; indTarget2++) {
                const target2 = state[indTarget2];

                if (indTarget != indTarget2) {
                    newtonMult -= (4 / c2) * target2.mu / Rij[indTarget][indTarget2];
                }

                if (indSource != indTarget2) {
                    newtonMult -= (1 / c2) * target2.mu / Rij[indSource][indTarget2];
                }
            }

            const dist = Rij[indSource][indTarget];
            const dist3 = Rij3[indSource][indTarget];

            newtonMult += Math.pow(norm(target.v), 2.0) / c2 
                       + 2 * Math.pow(norm(source.v), 2.0) / c2 
                       - (4.0/c2) * dot(target.v, source.v) 
                       - (1.5/c2) * Math.pow(dot(rTargetSource, source.v) / dist, 2.0)
                       + (0.5/c2) * dot(rTargetSource, source.accNewton);
                       
            // Add the Newtonian part and the remaining terms.
            target.accRel = linComb([
                1.0,
                newtonMult, 
                source.mu / c2 / dist3 * dot(rTargetSource, linComb([4, -3], [target.v, source.v])),
                3.5 / c2 * source.mu / dist ],
                [target.accRel, accNewton, vTargetSource, source.accNewton]);
        }
    }
}