import { limitAngleDeg } from "./Angles.js";
import { nutationTerms } from "./Nutation.js";
import { sind, cosd, atand, rad2Deg } from "./MathUtils.js";
import { correlationUt1Tdb } from "./TimeCorrelation.js";

/**
 * Compute Greenwich Mean Sidereal Time (GMST).
 * 
 * References:
 *  [1] E. Suirana, J. Zoronoza, M. Hernandez-Pajares - GNSS Data Processing -
 *  Volume I: Fundamentals and Algorithms, ESA 2013. 
 * 
 * @param {*} JT 
 *      Julian time.
 */
 export function timeGmstOld(JT)
 {
    // For computation of the UT1 time.
    const JDmin = Math.floor(JT) - 0.5;
    const JDmax = Math.floor(JT) + 0.5;
    let JD0 = 0;
    if (JT > JDmin)
    {
        JD0 = JDmin;
    }
    if (JT > JDmax)
    {
        JD0 = JDmax;
    }

    // Julian time at 2000-01-01 12:00:00 UT1
    const epochJ2000 = 2451545.0;
    // UT1 time
    const H = (JT - JD0) * 24.0;
    const UT1 = H * 15.0;
    // Julian centuries of UT1 date (A.36)
    const T = (JD0 - epochJ2000) / 36525.0;

    // Greenwich Mean Sidereal Time (GMST) at 0h UT1 (A.35)
    const theta_G0 = 100.460618375 + 36000.77005360834 * T + 3.879333333333333e-04 * T*T - 2.583333333333333e-08 *T*T*T;

    // GMST(A.34)
    return limitAngleDeg(1.002737909350795 * UT1 + theta_G0, 360);
 }

 /**
  * Compute Greenwich Mean Sidereal Time (GMST).
  * 
  * References: 
  *  [1] S. Urban, K. Seidelmann - Explanatory Supplement to the Astronomical Almanac,
  *      3rd edition, 2013.
  * 
  * @param {*} JTut1 
  *      Julian time (UT1)
  * @param {*} JTtt 
  *      Julian time (TT). If undefined, time correlation data is used.
  * @returns GMST time.
  */
export function timeGmst(JTut1, JTtt)
{
    if (JTtt === undefined)
    {
        JTtt = correlationUt1Tdb(JTut1);
    }

    const DU = JTut1 - 2451545.0;
    const T = (JTtt - 2451545.0) / 36525.0;
    const T2 = T*T;
    const T3 = T2*T;
    const T4 = T3*T;
    const T5 = T4*T;

    // Equation (6.64).
    const GMST = 86400.0 * (0.7790572732640 + 0.00273781191135448 * DU + DU % 1.0)
               + 0.00096707 + 307.47710227 * T + 0.092772113 * T2 - 2.93e-8 * T3 
               - 1.99708e-5 * T4 - 2.453e-9 * T5;

    return (GMST * 360.0 / 86400.0) % 360.0;
} 
 

/**
 * Compute Greenwich Apparent Sidereal Time (GAST).
 * 
 * References:
 *  [1] E. Suirana, J. Zoronoza, M. Hernandez-Pajares - GNSS Data Processing -
 *  Volume I: Fundamentals and Algorithms, ESA 2013. 
 * 
 * @param {*} JT 
 *      Julian time.
 * @param {*} nutParams
 *      Nutation parameters. If missing, the parameters are computed from JT. 
 */
export function timeGast(JT, nutParams)
{
    if (nutParams === undefined)
    {
        const T = (JT - 2451545.0) / 36525.0;
        nutParams = nutationTerms(T);
    }

    const GMST = timeGmst(JT);

    // The equinox equation (A.37) for GAST
    //return limitAngleDeg(GMST - atand(N12 / N11));
    return limitAngleDeg(GMST + nutParams.dpsi * cosd(nutParams.eps));
}

/**
 * 
 * @param {*} year 
 *      Year as an integer.
 * @param {*} month 
 *      Month (1-12).
 * @param {*} mday 
 *      Day of the month (1-31).
 * @returns Julian date.
 */
export function dateJulianYmd(year, month, mday)
{
    if (month < 3)
    {
        year--;
        month += 12;
    }

    const A = Math.floor(year / 100.0);
    const B = Math.floor(A / 4.0);
    const C = Math.floor(2.0 - A + B);
    const E = Math.floor(365.25 * (year + 4716.0));
    const F = Math.floor(30.6001 * (month + 1.0));

    return C + mday + E + F - 1524.5;    
}

/**
 * Compute Julian time.
 * 
 * @param {*} year 
 *      Year as an integer.
 * @param {*} month 
 *      Month (1-12) integer.
 * @param {*} mday 
 *      Day of the month (1-31) integer.
 * @param {*} hour 
 *      Hour (0-23) integer.
 * @param {*} minute
 *      Minute (0-59) integer. 
 * @param {*} second 
 *      Second (0-60) floating point.
 * @returns An object with JD and JT for Julian date and time.
 */
export function timeJulianYmdhms(year, month, mday, hour, minute, second)
{
    const JD = dateJulianYmd(year, month, mday);
    const JT = JD + hour / 24.0 + minute/(24.0 * 60.0) + second/(24.0 * 60.0 * 60.0);

    return {JD : JD, JT : JT};
}

/**
 * Compute Julian time.
 * 
 * @param {Date} d 
 *      Date object.
 * @returns Object 
 */
export function timeJulianTs(d)
{
    let year = d.getUTCFullYear();
    let month = d.getUTCMonth() + 1;
    
    let mday = d.getUTCDate();
    let hour = d.getUTCHours();
    let minute = d.getUTCMinutes();
    let second = d.getUTCSeconds() + d.getUTCMilliseconds() / 1000.0;

    return timeJulianYmdhms(year, month, mday, hour, minute, second);
}

/**
 * Compute Gregorian date and time from Julian time.
 */
export function timeGregorian(JT)
{
    // Meeus - Astronomical Algorithms - Chapter 7.
    const Z = Math.floor(JT + 0.5);
    const F = JT + 0.5 - Z;
    let A = Z;
    if (Z >= 2299161) 
    {
        let alpha = Math.floor((Z - 1867216.25) / 36524.25);
        A = Z + 1 + alpha - Math.floor(alpha / 4.0);
    }
    const B = A + 1524;
    const C = Math.floor((B - 122.1) / 365.25);
    const D = Math.floor(365.25 * C);
    const E = Math.floor((B - D)/30.6001);

    const mday = Math.floor(B - D - Math.floor(30.6001 * E) + F);
    let month = E - 1;
    if (E >= 14)
    {
        month = E - 13;
    }
    let year = C - 4716;
    if (month < 3)
    {
        year = C - 4715;
    }

    let JTfrac = F;
    if (JTfrac < 0)
    {
        JTfrac += 1;
    }
    const hour = Math.floor(JTfrac * 24.0);
    JTfrac -= hour / 24.0;
    const minute = Math.floor(JTfrac * (24.0 * 60.0));
    JTfrac -= minute / (24.0 * 60.0);
    const second = JTfrac * (24.0 * 60.0 * 60.0);

    return {year : year, month : month, mday : mday, 
        hour : hour, minute : minute, second : second};
}