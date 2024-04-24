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
        for (let timestep = 0; timestep < 36525*20; timestep++) {
            const h = 0.05;

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
        let state = dofToState(dof, constants.stateInitial, JT);
        console.log(state.objects);
        const endTime = performance.now();
        console.log(endTime - startTime);
        //console.log(state.objects);
        //console.log(state.libration);
    });
});