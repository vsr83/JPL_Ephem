import { nutationTerms } from "./Nutation.js";
import { rotateCart1, rotateCart2, rotateCart3 } from "./Rotations.js";
import { constants } from "./Constants.js";
import { cross, norm, linComb, vecMul, vecDiff, vecSum } from "./MathUtils.js";
import { legendreValue, legendreDeriv, legendreAssoc, legendreAssocd } from "./Legendre.js";
import { coordModTod, coordJ2000Mod, coordTodMod, coordModJ2000 } from "./Frames.js";

/**
 * Compute the acceleration and torque due to zonal and tesseral harmonics 
 * from an extended body.
 * 
 * This method computes the expression in equation (2) of [1] or (8.3) in
 * [2] and transforms the acceleration to body coordinates. 
 * 
 * REFERENCES: 
 *  [1] Newhall, Standish, Williams - DE 102: a numerically integrated
 *  ephemeris of the Moon and planets spanning forty-four centuries,
 *  Astronomy and Astrophysics, 125, 150-167, 1983.
 *  [2] Urban, Seidelmann - Explanatory Supplement to the Astronomical
 *  Almanac, 3rd edition, University Science Books, 2013.
 *  [3] Steve Moshier, DE118i available at 
 *  http://www.moshier.net/de118i-2.zip * 
 * 
 * @param {*} rPoint 
 *      The position of the point-mass w.r.t. the body center in body 
 *      coordinates (au).
 * @param {*} a 
 *      The equatorial radius of the extended body (au).
 * @param {*} mu 
 *      Standard gravitational parameter (au^3/d^2) or 1 if the results 
 *      are multiplied with -mu afterwards.
 * @param {*} Jn 
 *      Zonal harmonics for the extended body starting from n = 2.
 * @param {*} CSnm 
 *      Tesseral harmonics in the (n, m, C_nm, Snm) row format.
 * @returns The acceleration of the point mass in body coordinates 
 *      (au/d^2, num_targets x 3).
 */
export function accBody(rPoint, a, mu, Jn, CSnm) {
    // Distance between center of the extended body and the point mass.
    const r = norm(rPoint);

    // Latitude and longitude of the point mass w.r.t. body coordinates (rad).
    const sinLat = rPoint[2] / norm(rPoint);
    const latPoint = Math.asin(sinLat);
    const lonPoint = Math.atan2(rPoint[1], rPoint[0]);
    const cosLat = Math.cos(latPoint);

    // Number of zonal harmonics starting from n=2.
    const numberZonal = Jn.length;
    const numberTesseral = CSnm.length;
    
    let accPointZonal = [0, 0, 0];
    let accPointTesseral = [0, 0, 0];
    let accPoint = [0, 0, 0];

    // Evaluate zonal harmonics.
    for (let indZonal = 0; indZonal < numberZonal; indZonal++) {
        const n = indZonal + 2;

        // Legendre value and derivative terms.
        const Pn    = legendreValue(n, sinLat);
        const PnDot = legendreDeriv(n, sinLat);

        accPointZonal = linComb([1, Jn[indZonal] * Math.pow(a/r, n)], 
            [accPointZonal, [(n + 1) * Pn, 0, -cosLat * PnDot]]);
    }
    accPointZonal = vecMul(accPointZonal, -mu / (r * r));

    // Evaluate tesseral harmonics.
    for (let indTesseral = 0; indTesseral < numberTesseral; indTesseral++) {
        const n    = CSnm[indTesseral][0];
        const m    = CSnm[indTesseral][1];
        const Cnm  = CSnm[indTesseral][2];
        const Snm  = CSnm[indTesseral][3];
    
        const cosMlon = Math.cos(m * lonPoint);
        const sinMlon = Math.sin(m * lonPoint);

        const Pnm    = Math.pow(-1, m) * legendreAssoc(n, m, sinLat);
        const PnmDot = Math.pow(-1, m) * legendreAssocd(n, m, sinLat);

        accPointTesseral = linComb([1, Math.pow(a / r, n)],
            [accPointTesseral,
            [-(n + 1)     * Pnm    * ( Cnm * cosMlon + Snm * sinMlon), 
              (m/cosLat) * Pnm    * (-Cnm * sinMlon + Snm * cosMlon), 
              cosLat     * PnmDot * ( Cnm * cosMlon + Snm * sinMlon)]]);
    }
    accPointTesseral = vecMul(accPointTesseral, -mu / (r * r));
    accPoint = linComb([1, 1], [accPointZonal, accPointTesseral]);

    return rotateCart3(rotateCart2(accPoint, latPoint), -lonPoint);
}

/**
 * Transform position to body coordinates.
 * 
 * @param {*} r   The position vector.
 * @param {*} phi The clockwise angle along the xy-plane to the line of 
 *                nodes from the x axis (radians)
 * @param {*} theta The clockwise inclination of the body equator (radians).
 * @param {*} psi The clockwise angle from the node to the prime meridian
 *                along the body equator (radians).
 * @returns Position in body coordinates.
 */
function coordToBody(r, phi, theta, psi) {
    return rotateCart3(rotateCart1(rotateCart3(r, phi), theta), psi);
}

/**
 * Transform position from body coordinates.
 * 
 * @param {*} r   The position vector.
 * @param {*} phi The clockwise angle along the xy-plane to the line of 
 *                nodes from the x axis (radians)
 * @param {*} theta The clockwise inclination of the body equator (radians).
 * @param {*} psi The clockwise angle from the node to the prime meridian
 *                along the body equator (radians).
 * @returns Position in body coordinates.
 */
function coordFromBody(r, phi, theta, psi) {
    return rotateCart3(rotateCart1(rotateCart3(r, -psi), -theta), -phi);
}


/**
 * Compute the accelerations due to Earth and Moon figure and tides.
 * 
 * This method is heavily based on the oblate method in [2].
 * 
 * REFERENCES: 
 * [1] Newhall, Standish, Williams - DE 102: a numerically integrated
 * ephemeris of the Moon and planets spanning forty-four centuries,
 * Astronomy and Astrophysics, 125, 150-167, 1983.
 * [2] Steve Moshier, DE118i available at 
 *  http://www.moshier.net/de118i-2.zip  
 * 
 * @param {*} OSV 
 *      The fields r (au), v (au/d), mu (au^3/d^2) for "Sun", "Earth", "Moon".
 * @param {*} libMoon 
 *      Libration state with the fields phi, phi1, theta, theta1, psi, psi1 
 *      (rad or rad/day).
 * @param {*} JT 
 *      Julian time.
 * @returns Object field accelerations (au/d^2, 3) for "Sun", "Earth", "Moon".
 */
export function accOblateness(state) {
    const osvSun = state.objects[state.objectIndices["Sun"]];
    const osvEarth = state.objects[state.objectIndices["Earth"]];
    const osvMoon = state.objects[state.objectIndices["Moon"]];

    const rS = osvSun.r;
    const rE = osvEarth.r;
    const rM = osvMoon.r;
    const muS = osvSun.mu;
    const muE = osvEarth.mu;
    const muM = osvMoon.mu;
    const Je = constants.Je;
    const Jm = constants.Jm;
    const CSnm = constants.CSnm;

    // Parse libration angles.
    const libMoon = state.libration;
    const phi    = libMoon.phi;
    const theta  = libMoon.theta;
    const psi    = libMoon.psi;
    const JT = state.JT;

    const nutData = nutationTerms((JT - 2451545.0) / 36525.0);

    // The position of the Earth w.r.t. Moon body center in DE118/J2000 and 
    // body coordinates.
    const rEmJ2000 = vecDiff(rE, rM);
    const rEmBody = coordToBody(rEmJ2000, phi, theta, psi);

    // Acceleration/mu and
    const accEmBodyTmp = accBody(rEmBody, constants.aMoon, 1, Jm, CSnm);
    // Torque per unit mass.
    const Tearth = cross(rEmBody, accEmBodyTmp);

    // 1. Accelerations from the interaction between the Moon figure and Earth.
    const accEmBody     = vecMul(accEmBodyTmp, -muM);
    const accEmJ2000Fig = coordFromBody(accEmBody, phi, theta, psi);
    const accMeBody     = vecMul(accEmBodyTmp, muE);
    const accMeJ2000Fig = coordFromBody(accMeBody, phi, theta, psi);

    const rSmJ2000 = vecDiff(rS, rM);
    const rSmBody  = coordToBody(rSmJ2000, phi, theta, psi);
    const accSmBodyTmp = accBody(rSmBody, constants.aMoon, 1, Jm, CSnm);
    const Tsun = cross(rSmBody, accSmBodyTmp);

    // 2. Accelerations from the interaction between the Moon figure and Sun.
    const accSmBody     = vecMul(accSmBodyTmp, -muM);
    const accSmJ2000Fig = coordFromBody(accSmBody, phi, theta, psi);
    const accMsBody     = vecMul(accSmBodyTmp, muS);
    const accMsJ2000Fig = coordFromBody(accMsBody, phi, theta, psi);

    // 3. Libration of the Moon.

    // Compute the total torque on the Moon and the angular accelerations.
    const T = linComb([muE, muS], [Tearth, Tsun]);
    librationMoon(state.libration, T);

    // 4. Oblateness of the Earth.

    // The position of the Moon w.r.t. Earth body center in DE118/J2000.
    const rMeJ2000 = vecDiff(rM, rE);

    // The position of the Sun w.r.t. Earth body center in DE118/J2000.
    const rSeJ2000 = vecDiff(rS, rE);

    // Transform the relative position of the Moon to the True-of-Date frame.
    const osvMeMod = coordJ2000Mod({r : rMeJ2000, v : [0, 0, 0], JT : JT});
    const rMeTod = coordModTod(osvMeMod, nutData).r;

    // Transform the relative position of the Sun to the True-of-Date frame.
    const osvSeMod = coordJ2000Mod({r : rSeJ2000, v : [0, 0, 0], JT : JT});
    const rSeTod = coordModTod(osvSeMod, nutData).r;

    const accMeTodTmp = accBody(rMeTod, constants.aEarth, 1, Je, []);
    const accSeTodTmp = accBody(rSeTod, constants.aEarth, 1, Je, []);

    const accMeTod = vecMul(accMeTodTmp, -muE);
    const accSeTod = vecMul(accSeTodTmp, -muE);
    const accEmTod = vecMul(accMeTodTmp, muM);
    const accEsTod = vecMul(accSeTodTmp, muS);

    // 5. Accelerations from the interaction between Earth tides and the Moon.
    //[acc_me_tod_tides, acc_em_tod_tides] = acc_tides(r_me_tod, mu_e, mu_m);
    const {accMeTodTides , accEmTodTides} = accTides(rMeTod, muE, muM);


    // Convert accelerations from Earth oblateness and tides to J2000 frame.
    const accMeJ2000Obl = coordModJ2000(
        coordTodMod({r : accMeTod, v : [0, 0, 0], JT : JT}, nutData)).r;
    const accSeJ2000Obl = coordModJ2000(
        coordTodMod({r : accSeTod, v : [0, 0, 0], JT : JT}, nutData)).r;
    const accEmJ2000Obl = coordModJ2000(
        coordTodMod({r : accEmTod, v : [0, 0, 0], JT : JT}, nutData)).r;
    const accEsJ2000Obl = coordModJ2000(
        coordTodMod({r : accEsTod, v : [0, 0, 0], JT : JT}, nutData)).r;
    const accMeJ2000Tides = coordModJ2000(
        coordTodMod({r : accMeTodTides, v : [0, 0, 0], JT : JT}, nutData)).r;
    const accEmJ2000Tides = coordModJ2000(
        coordTodMod({r : accEmTodTides, v : [0, 0, 0], JT : JT}, nutData)).r;    

    const accSJ2000 = vecSum(accSmJ2000Fig, accSeJ2000Obl);
    const accEJ2000 = linComb([1, 1, 1, 1], 
        [accEsJ2000Obl, accEmJ2000Fig, accEmJ2000Obl, accEmJ2000Tides]);
    const accMJ2000 = linComb([1, 1, 1, 1], 
        [accMsJ2000Fig, accMeJ2000Fig, accMeJ2000Obl, accMeJ2000Tides]);

        /*
    console.log('Moon Figure <-> Earth : Earth Acceleration');
    console.log(accEmJ2000Fig);
    console.log('Moon Figure <-> Earth : Moon Acceleration');
    console.log(accMeJ2000Fig);
    console.log('Moon Figure <-> Sun : Sun Acceleration');
    console.log(accSmJ2000Fig);
    console.log('Moon Figure <-> Sun : Moon Acceleration');
    console.log(accMsJ2000Fig);
    console.log('Earth Oblateness <-> Moon : Earth Acceleration');
    console.log(accEmJ2000Obl);
    console.log('Earth Oblateness <-> Moon : Moon Acceleration');
    console.log(accMeJ2000Obl);
    console.log('Earth Oblateness <-> Sun : Earth Acceleration');
    console.log(accEsJ2000Obl);
    console.log('Earth Oblateness <-> Sun : Sun Acceleration');
    console.log(accSeJ2000Obl);
    console.log('Earth Tides <-> Moon : Earth Acceleration');
    console.log(accEmJ2000Tides);
    console.log('Earth Tides <-> Moon : Moon Acceleration');
    console.log(accMeJ2000Tides);        
    */
    osvSun.accObl   = accSJ2000;
    osvEarth.accObl = accEJ2000;
    osvMoon.accObl  = accMJ2000;
}

/**
 * Compute the second derivatives for the Moon libration angles.
 * 
 * This method computes the expression in equation (3) of [1] or (8.6)-(8.8) 
 * in [2].
 * 
 * The method fills the fields phi2, theta2 and psi2 (rad/d^2).
 * 
 * REFERENCES: 
 *  [1] Newhall, Standish, Williams - DE 102: a numerically integrated
 *  ephemeris of the Moon and planets spanning forty-four centuries,
 *  Astronomy and Astrophysics, 125, 150-167, 1983.
 *  [2] Urban, Seidelmann - Explanatory Supplement to the Astronomical
 *  Almanac, 3rd edition, University Science Books, 2013.
 *  [3] Steve Moshier, DE118i available at 
 *  http://www.moshier.net/de118i-2.zip
 *  [4] Ferrari et. al. - Geophysical Parameters of the Earth-Moon System,
 *  Journal of Geophysical Research, 1980.
 * 
 * @param {*} librationState 
 *      Libration state with the fields phi, phi1, theta, theta1, psi, psi1 
 *      (rad or rad/day).
 * @param {*} N 
 *      Torque per unit mass in body coordinates.
 */
export function librationMoon(librationState, N) {
    const phi    = librationState.phi;
    const phi1   = librationState.phi1;
    const theta  = librationState.theta;
    const theta1 = librationState.theta1;
    const psi    = librationState.psi;
    const psi1   = librationState.psi1;
    const inertiaMoon = constants.inertiaMoon;

    // Angular velocity vector.
    const omegaX = phi1 * Math.sin(theta) * Math.sin(psi) + theta1 * Math.cos(psi);
    const omegaY = phi1 * Math.sin(theta) * Math.cos(psi) - theta1 * Math.sin(psi);
    const omegaZ = phi1 * Math.cos(theta) + psi1;

    // Differential equations for the angular velocity from Euler's equations.
    const omegaX1 =  omegaY * omegaZ * (inertiaMoon.gammaL - inertiaMoon.betaL) 
                  / (1 - inertiaMoon.betaL * inertiaMoon.gammaL) + N[0] / inertiaMoon.A;
    const omegaY1 =  omegaZ * omegaX * inertiaMoon.betaL + N[1] / inertiaMoon.B;
    const omegaZ1 = -omegaX * omegaY * inertiaMoon.gammaL + N[2] / inertiaMoon.C;

    // Differential equations for the three Euler angles.
    const phi2 = (omegaX1 * Math.sin(psi) + omegaY1 * Math.cos(psi) 
               + theta1 * (psi1 - phi1 * Math.cos(theta))) / Math.sin(theta);
    const theta2 = omegaX1 * Math.cos(psi) - omegaY1 * Math.sin(psi) 
                 - phi1 * psi1 * Math.sin(theta);
    const psi2 = omegaZ1 - phi2 * Math.cos(theta) + phi1 * theta1 * Math.sin(theta);

    librationState.phi2 = phi2;
    librationState.theta2 = theta2;
    librationState.psi2 = psi2;
}

/**
 * Compute acceleration of the Moon and the Earth due to tides.
 * 
 * [1] Newhall, Standish, Williams - DE 102: a numerically integrated
 * ephemeris of the Moon and planets spanning forty-four centuries,
 * Astronomy and Astrophysics, 125, 150-167, 1983.
 * [2] Steve Moshier, DE118i available at 
 * http://www.moshier.net/de118i-2.zip
 * 
 * @param {*} rMeTod 
 *      The position of the Moon w.r.t. Earth in the ToD frame (au, 3)
 * @param {*} muE 
 *      Standard gravitational parameter (au^3/d^2) for Earth.
 * @param {*} muM 
 *      Standard gravitational parameter (au^3/d^2) for Moon.
 * @returns Objects with fields accMeTodTides and accEmTodTides for the 
 *      accelerations of the Moon and the Earth (au/d^2, 3).
 */
export function accTides(rMeTod, muE, muM) {
    // Distance between Earth and the Moon.
    const rEm = norm(rMeTod);
    const accMoon =  vecMul(
        [rMeTod[0] + constants.phase * rMeTod[1],
         rMeTod[1] - constants.phase * rMeTod[0],
         rMeTod[2]],     
        - (3 * constants.love * muM) * (1 + muM/muE) 
        * (constants.aEarth ** 5 / rEm ** 8));
    const accEarth = vecMul(accMoon, -muM / muE);

    return {accMeTodTides : accMoon, accEmTodTides : accEarth};
}