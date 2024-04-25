import { checkFloat, checkFloatArray} from './common.js';
import { cross, norm, vecDiff} from '../src/MathUtils.js';
import { accBody, accOblateness, librationMoon} from '../src/Figure.js';
import { constants } from '../src/Constants.js';
import { funcLibration, stateToDof, dofToState, runge4, adams8 } from '../src/Integration.js';

describe('Integration', function() {
    it('RK4 test', function() {
        const startTime = performance.now();

        const JT = constants.stateInitial.JT;
        let dof = stateToDof(constants.stateInitial);
        funcLibration(0, dof, JT);

        let t = 0;

        let F = [];
        const h = 0.05;
        for (let timestep = 0; timestep < 365.25 * 100 / h; timestep++) {

            let out = [];
            if (timestep < 8) {
                out = runge4(funcLibration, t, dof, h, JT);
                const y = funcLibration(out.tOut, out.yOut, JT);
                F.unshift(y);
            } else {
                out = adams8(funcLibration, t, dof, F, h, JT);
                F = out.Fout;
            }
            t = out.tOut;
            dof = out.yOut;
            if (timestep % 36525 == 0) console.log(t);
            //console.log(dofToState(dof, constants.stateInitial, JT).objects);
        }
        let state = dofToState(dof, t, constants.stateInitial, JT);
        console.log(state.objects);
        const endTime = performance.now();
        //console.log(state.objects);
        //console.log(state.libration);

        const objectsExp = {
            Sun : {
                r : [-0.001366042119904, 0.001086570534222, 0.000479276172850],
                v : [-0.000002485308045, -0.000004995123049, -0.000002104784454]
            },
            Mercury : {
                r : [0.063358904813830,  0.268961575670744, 0.136898565567972],
                v : [-0.033147172797133, 0.004957604762552,  0.006080247077021]
            },
            Venus : {
                r : [-0.695766559126292, 0.151549952264849, 0.112113175625333],
                v : [-0.005243103890166,  -0.018049444202374,  -0.007792810653165]
            },
            Earth : {
                r : [0.104190111711093, -0.926590528808781, -0.401545624050591],
                v : [0.016826949493706,   0.001580053036777,   0.000683919503085]
            },
            Moon : {
                r : [0.101689472301485, -0.927496497431367, -0.402025106811183],
                v : [0.017032975771069,   0.001086831563088,  0.000521534743146]
            },
            Mars : {
                r : [1.255005132172056, -0.506955213143151, -0.266315415568197],
                v : [0.006327204339345,   0.012725848217208,   0.005667439357743]
            },
            Jupiter : {
                r : [4.687209400064162, -1.522941534172532, -0.766669007188262],
                v : [0.002481977956446,   0.006871944602397,   0.002884822094047]
            },
            Saturn : {
                r : [-9.487592488267653,-1.048416554168364, -0.023772315422021],
                v : [0.000266675387653,  -0.005135806730118,  -0.002133431710838]
            },
            Uranus : {
                r : [-4.483181311793592, -17.004998451524436, -7.384184285966848],
                v : [0.003794695422726,  -0.000996179143481,  -0.000489913281920]
            },
            Neptune : {
                r : [-3.892225308217819,  27.415006485310887, 11.318010935320046],
                v : [-0.003130520679461,  -0.000387370583558,  -0.000080611010613]
            },
            Pluto : {
                r : [43.601672365271234,   5.815632782670460, -11.321978006214040],
                v : [0.000347731671945,   0.002250644480601,   0.000597609922177]
            }
        }
        console.log("Computation took " + ((endTime - startTime)/ 1000) + " seconds.");

        const names = constants.stateInitial.objectIndices;
        for (let ind = 0; ind < state.objects.length; ind++) {
            const objName = state.objects[ind].name;
            const rExp = objectsExp[objName].r;
            const vExp = objectsExp[objName].v;
            const r = state.objects[ind].r;
            const v = state.objects[ind].v;

            console.log(objName);

            const errR = norm(vecDiff(r, rExp)) * constants.au;
            const errV = norm(vecDiff(v, vExp)) * constants.au;
            console.log("  Position error " + errR + " km");
            console.log("  Velocity error " + errV + " km/d");
            
            //console.log(r);
            //console.log(rExp);
            //console.log(v);
            //console.log(vExp);
        }
    });
});