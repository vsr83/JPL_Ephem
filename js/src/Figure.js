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
        n = indZonal + 2;

        // Legendre value and derivative terms.
        const Pn    = legendreValue(n, sinLat);
        const PnDot = legendreDeriv(n, sinLat);

        accPointZonal = linComb([1, Math.pow(a/r, n)], 
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

        Pnm    = Math.pow(-1, m) * legendreAssoc(n, m, sin_lat);
        PnmDot = Math.pow(-1, m) * legendreAssocd(n, m, sin_lat);

        accPointTesseral = linComb([1, Math.pow(a / r, n)],
            [accPointTesseral,
            [-(n + 1)     * Pnm    * ( Cnm * cosMlon + Snm * sinMlon), 
              (m/cos_lat) * Pnm    * (-Cnm * sinMlon + Snm * cosMlon), 
              cos_lat     * PnmDot * ( Cnm * cosMlon + Snm * sinMlon)]]);
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

function accOblateness(OSV, mu, libMoon, JT, Je, Jm, CSnm) {
    const muS = mu[0];
    const muE = mu[1];
    const muM = mu[2];

    rS = OSV["Sun"].r;
    rE = OSV["Earth"].r;
    rM = OSV["Moon"].r;
    vS = OSV["Sun"].v;
    vE = OSV["Earth"].v;
    vM = OSV["Moon"].v;

    // Parse libration angles.
    const phi    = libMoon(1);
    const phi1   = libMoon(2);
    const theta  = libMoon(3);
    const theta1 = libMoon(4);
    const psi    = libMoon(5);
    const psi1   = libMoon(6);

    const nutData = nutationTerms((JT - 2451545.0) / 36525.0);

    // The position of the Earth w.r.t. Moon body center in DE118/J2000 and 
    // body coordinates.
    const rEmJ2000 = rE - rM;
    const rEmBody = coordToBody(rEmJ2000, phi, theta, psi);

    // Acceleration/mu and
    const accEmBodyTmp = accBody(rEmBody, aMoon, 1, Jm, CSnm);
    // Torque per unit mass.
    const Tearth = cross(rEmBody, accEmBodyTmp);

    // 1. Accelerations from the interaction between the Moon figure and Earth.
    const accEmBody     = vecMul(accEmBodyTmp, -muM);
    const accEmJ2000Fig = coordFromBody(accEmBody, phi, theta, psi);
    const accMeBody     = vecMul(accEmBodyTmp, muE);
    const accMeJ2000Fig = coordFromBody(accMeBody, phi, theta, psi);
    
    const rSmJ2000 = vecDiff(rS, rM);
    const rSmBody  = coordToBody(rSmJ2000, phi, theta, psi);
    const accSmBodyTmp = accBody(rSmBody, aMoon, 1, Jm, CSnm);
    const Tsun = cross(rSmBody, accSmBodyTmp);

    // 2. Accelerations from the interaction between the Moon figure and Sun.
    const accSmBody     = vecMul(accSmBodyTmp, -mu_m);
    const accSmJ2000Fig = coordFromBody(accSmBody, phi, theta, psi);
    const accMsBody     = vecMul(accSmBodyTmp, mu_s);
    const accMsJ2000Fig = coordFromBody(accMsBody, phi, theta, psi);

    // 3. Libration of the Moon.

    // Compute the total torque on the Moon and the angular accelerations.
    const T = linComb([mu_e, mu_s], [T_earth, T_sun]);
    const accMoon = librationMoon(phi, theta, psi, phi1, theta1, psi1, T);

    // 4. Oblateness of the Earth.

    // The position of the Moon w.r.t. Earth body center in DE118/J2000.
    const rMeJ2000 = vecDiff(rM, rE);

    // The position of the Sun w.r.t. Earth body center in DE118/J2000.
    const rSeJ2000 = vecDiff(rS, rE);

    // Transform the relative position of the Moon to the True-of-Date frame.
    const osvMeMod = coordJ2000Mod({r : rMeJ2000, v : [0, 0, 0], JT : JT});
    const rMeTod = coordModTod(osvMeMod, nutData).r;

    // Transform the relative position of the Sun to the True-of-Date frame.
    const osvSeMod = coordJ2000Mod({r : r_se_j2000, v : [0, 0, 0], JT : JT});
    const rSeTod = coordModTod(osvSeMod, nutData).r;

    const accMeTodTmp = accBody(rMeTod, aEarth, 1, Je, []);
    const accSeTodTmp = accBody(rSeTod, aEarth, 1, Je, []);

    const accMeTod = vecMul(accMeTodTmp, -muE);
    const accSeTod = vecMul(accSeTodTmp, -muE);
    const accEmTod = vecMul(accMeTodTmp, muM);
    const accEsTod = vecMul(accSeTodTmp, muS);

    // 5. Accelerations from the interaction between Earth tides and the Moon.
    //[acc_me_tod_tides, acc_em_tod_tides] = acc_tides(r_me_tod, mu_e, mu_m);

    // Convert accelerations from Earth oblateness and tides to J2000 frame.
    acc_me_j2000_obl = coordTodMod(acc_me_tod, nutData);
    acc_se_j2000_obl = coordTodMod(acc_se_tod, nutData);
    acc_em_j2000_obl = coordTodMod(acc_em_tod, nutData);
    acc_es_j2000_obl = coordTodMod(acc_es_tod, nutData);
    acc_me_j2000_tides = coordTodMod(acc_me_tod_tides, nutData);
    acc_em_j2000_tides = coordTodMod(acc_em_tod_tides, nutData);    
}