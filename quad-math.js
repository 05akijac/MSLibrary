/*
 * Quadrupole Lab — independent numerical physics implementation.
 * No external libraries. Works in browsers and Node.js.
 *
 * Electrode convention: opposite rod pairs at +/-[U+V cos(Omega*t)],
 * where U and V are the DC and RF zero-to-peak voltages PER ROD PAIR.
 * Phi=(U+V cos Omega*t)(x*x-y*y)/(r0*r0).
 * a=8eU/(m r0^2 Omega^2), q=4eV/(m r0^2 Omega^2).
 * Floquet uses canonical Mathieu sign (a-2q*cos 2tau); a q sign
 * change corresponds to half-RF-cycle phase shift, leaving stability unchanged.
 */
(function (root, factory) {
  const result = factory();
  if (typeof module === 'object' && module.exports) module.exports = result;
  if (root) root.QuadMath = result;
})(typeof globalThis !== 'undefined' ? globalThis : this, function () {
  'use strict';
  const E = 1.602176634e-19;
  const AMU = 1.66053906660e-27;
  const TWO_PI = 2 * Math.PI;
  const cfg = {steps: 96};
  let apexCache = null;
  let curveCache = null;
  const validPos = x => Number.isFinite(x) && x > 0;

  function physicalFromAQ(mz, frequencyMHz, radiusMm, a, q) {
    if (![mz, frequencyMHz, radiusMm, a, q].every(Number.isFinite) ||
        !validPos(mz) || !validPos(frequencyMHz) || !validPos(radiusMm) || a < 0 || q < 0) {
      throw Error('Positive m/z, frequency and field radius are required; a and q must be nonnegative.');
    }
    const omega = TWO_PI * frequencyMHz * 1e6;
    const r = radiusMm * 1e-3;
    // mz is a numerical mass-to-charge value expressed in u per elementary charge.
    const c = mz * AMU * r * r * omega * omega / E;
    const U = a * c / 8;
    const V = q * c / 4;
    return {a, q, U, V, VppPerRod: 2 * V, differentialRFPeak: 2 * V,
      differentialRFpp: 4 * V, differentialDC: 2 * U, omega};
  }

  function aqFromPhysical(mz, frequencyMHz, radiusMm, U, V) {
    if (![mz, frequencyMHz, radiusMm, U, V].every(Number.isFinite) ||
        !validPos(mz) || !validPos(frequencyMHz) || !validPos(radiusMm)) {
      throw Error('Invalid geometry, frequency, m/z or voltage.');
    }
    const omega = TWO_PI * frequencyMHz * 1e6;
    const r = radiusMm * 1e-3;
    const f = E / (mz * AMU * r * r * omega * omega);
    return {a: 8 * U * f, q: 4 * V * f};
  }

  // RK4 integration of two fundamental solutions in one pass over one
  // coefficient period, tau = [0, pi]. Returns trace of monodromy M.
  // d²u/dtau² + [a - 2q cos(2tau)]u = 0.
  function floquetTrace(a, q, steps = cfg.steps) {
    const h = Math.PI / steps;
    let p = 1, v = 0, s = 0, w = 1;
    for (let i = 0; i < steps; i++) {
      const t = i * h;
      const c1 = 2 * q * Math.cos(2 * t) - a;
      const c2 = 2 * q * Math.cos(2 * (t + h / 2)) - a;
      const c4 = 2 * q * Math.cos(2 * (t + h)) - a;
      const k1p = v, k1v = c1 * p, k1s = w, k1w = c1 * s;
      const k2p = v + h * k1v / 2;
      const k2v = c2 * (p + h * k1p / 2);
      const k2s = w + h * k1w / 2;
      const k2w = c2 * (s + h * k1s / 2);
      const k3p = v + h * k2v / 2;
      const k3v = c2 * (p + h * k2p / 2);
      const k3s = w + h * k2w / 2;
      const k3w = c2 * (s + h * k2s / 2);
      const k4p = v + h * k3v;
      const k4v = c4 * (p + h * k3p);
      const k4s = w + h * k3w;
      const k4w = c4 * (s + h * k3s);
      p += h / 6 * (k1p + 2*k2p + 2*k3p + k4p);
      v += h / 6 * (k1v + 2*k2v + 2*k3v + k4v);
      s += h / 6 * (k1s + 2*k2s + 2*k3s + k4s);
      w += h / 6 * (k1w + 2*k2w + 2*k3w + k4w);
    }
    return p + w;
  }

  function bisect(fn, lo, hi, iterations = 33) {
    let flo = fn(lo), fhi = fn(hi);
    if (Math.abs(flo) < 1e-12) return lo;
    if (Math.abs(fhi) < 1e-12) return hi;
    if (flo * fhi > 0 || !Number.isFinite(flo * fhi)) throw Error('Root not bracketed.');
    for (let i=0; i<iterations; i++) {
      const mid=(lo+hi)/2, fm=fn(mid);
      if (flo*fm <= 0) {hi=mid;fhi=fm;} else {lo=mid;flo=fm;}
    }
    return (lo+hi)/2;
  }

  // Region I boundaries (not an arbitrary finite-time amplitude threshold):
  // x: first antiperiodic characteristic (Tr Mx = -2)
  // y: first periodic characteristic (Tr My = +2)
  function xLimit(q) {
    if (q === 0) return 1;
    return bisect(a => floquetTrace(a,q) + 2, 0, 1, 24);
  }
  function yLimit(q) {
    if (q === 0) return 0;
    return bisect(a => floquetTrace(-a,q) - 2, 0, 0.7, 24);
  }
  function apex() {
    if (apexCache) return apexCache;
    const q=bisect(q => xLimit(q)-yLimit(q),0.69,0.72,27);
    const a=(xLimit(q)+yLimit(q))/2;
    const qZero=bisect(q => floquetTrace(0,q)+2, 0.89, 0.94, 33);
    apexCache={q,a,slope:a/q,qZero};
    return apexCache;
  }

  // Strict |trace|<2 is stability of the linearized equation, excluding the
  // neutrally periodic boundary. Restrict to the first stability island.
  function isFirstStable(a,q) {
    const {qZero}=apex();
    if (!Number.isFinite(a+q) || q < 0 || a < 0 || q >= qZero || a > .45) return false;
    const tx=floquetTrace(a,q), ty=floquetTrace(-a,q);
    return Math.abs(tx) < 2 && Math.abs(ty) < 2;
  }

  // Given scan-line slope a/q, locate two crossings with the Region I edge.
  function crossings(slope) {
    const tip=apex();
    if (!(slope >= 0 && slope < tip.slope)) throw Error('Choose a/q below the first-region apex slope.');
    const qLow=slope === 0 ? 0 : bisect(q => floquetTrace(-slope*q,q)-2, 1e-7, tip.q, 34);
    const qHigh=bisect(q => floquetTrace(slope*q,q)+2, tip.q, tip.qZero+0.0005, 34);
    return {qLow,qHigh};
  }

  function centredBand(mz, slope) {
    const {qLow,qHigh}=crossings(slope);
    const qCentre=2*qLow*qHigh/(qLow+qHigh);
    const lowerMz=mz*qCentre/qHigh;
    const upperMz=qLow === 0 ? Infinity : mz*qCentre/qLow;
    return {qLow,qHigh,qCentre,lowerMz,upperMz,
      width:upperMz-lowerMz, resolution:mz/(upperMz-lowerMz)};
  }

  function solveWidth(mz, width) {
    if (!validPos(mz) || !validPos(width) || width>=2*mz) {
      throw Error('Target passband must be > 0 and less than twice the target m/z.');
    }
    const tip=apex();
    let lo=0, hi=tip.slope * (1-1e-8);
    // Width decreases monotonically with increasing operating slope.
    for (let i=0; i<33; i++) {
      const mid=(lo+hi)/2;
      const band=centredBand(mz,mid);
      if (band.width > width) lo=mid; else hi=mid;
    }
    const slope=(lo+hi)/2;
    const band=centredBand(mz,slope);
    if (Math.abs(band.width-width)/width > 0.02) {
      throw Error('Requested width is below the numerical precision of this prototype. Try a larger Δ(m/z).');
    }
    return {a:band.qCentre*slope,q:band.qCentre,slope,...band};
  }

  function stabilityCurve(n=140) {
    if (curveCache && curveCache.n===n) return curveCache.data;
    const tip=apex();
    const data=[];
    for(let i=0;i<=n;i++) {
      const q=tip.qZero*i/n;
      const a=i===n ? 0 : q<tip.q ? yLimit(q) : q>tip.q ? xLimit(q) : tip.a;
      data.push({q,a:Math.max(0,a)});
    }
    curveCache={n,data};
    return data;
  }

  function trajectory({targetMz, testMz, frequencyMHz, radiusMm, lengthMm,
      axialEVPerCharge, a, q, phaseDeg=0, x0Mm=0.12, y0Mm=0.12,
      vx0=0, vy0=0, stepsPerCycle=60}) {
    if (![targetMz,testMz,frequencyMHz,radiusMm,lengthMm,axialEVPerCharge,a,q,
      phaseDeg,x0Mm,y0Mm,vx0,vy0,stepsPerCycle].every(Number.isFinite) ||
      ![targetMz,testMz,frequencyMHz,radiusMm,lengthMm,axialEVPerCharge,stepsPerCycle].every(validPos)) {
      throw Error('Trajectory settings must be finite and positive where required.');
    }
    const omega=TWO_PI*frequencyMHz*1e6;
    // E is axial eV per elementary charge: z cancels when m/z is used.
    const vz=Math.sqrt(2*axialEVPerCharge*E/(testMz*AMU));
    const flight=(lengthMm/1000)/vz;
    const cycles=flight*frequencyMHz*1e6;
    const steps=Math.ceil(cycles*stepsPerCycle);
    if (steps>60000) throw Error('Trajectory requires too many RK4 steps; reduce length or frequency, or raise axial energy.');
    const dTau=Math.PI*cycles/steps;
    const ratio=targetMz/testMz, ap=a*ratio, qp=q*ratio;
    const r=radiusMm/1000, phase=phaseDeg*Math.PI/180;
    let x=x0Mm/1000, y=y0Mm/1000;
    let dx=vx0*2/omega, dy=vy0*2/omega;
    const samples=[];
    const stride=Math.max(1,Math.floor(steps/850));
    let transmitted=true, lostAtMm=null, maxRadius=Math.hypot(x,y);
    const acceleration=(tau,xx,yy) => {
      const factor=ap+2*qp*Math.cos(2*tau+phase);
      return [-factor*xx, factor*yy];
    };
    if(maxRadius>=r) {transmitted=false;lostAtMm=0;}
    samples.push({z:0,x:x*1000,y:y*1000});
    for(let i=0;i<steps && transmitted;i++) {
      const tau=i*dTau, h=dTau;
      const [a1x,a1y]=acceleration(tau,x,y);
      // Conventional RK4 on state [x,dx,y,dy] (not Nyström):
      const k1x=dx,k1dx=a1x,k1y=dy,k1dy=a1y;
      const k2x=dx+h*k1dx/2, k2dx=acceleration(tau+h/2,x+h*k1x/2,y+h*k1y/2)[0];
      const k2y=dy+h*k1dy/2, k2dy=acceleration(tau+h/2,x+h*k1x/2,y+h*k1y/2)[1];
      const k3x=dx+h*k2dx/2, k3y=dy+h*k2dy/2;
      const [k3dx,k3dy]=acceleration(tau+h/2,x+h*k2x/2,y+h*k2y/2);
      const k4x=dx+h*k3dx, k4y=dy+h*k3dy;
      const [k4dx,k4dy]=acceleration(tau+h,x+h*k3x,y+h*k3y);
      x+=h*(k1x+2*k2x+2*k3x+k4x)/6;
      dx+=h*(k1dx+2*k2dx+2*k3dx+k4dx)/6;
      y+=h*(k1y+2*k2y+2*k3y+k4y)/6;
      dy+=h*(k1dy+2*k2dy+2*k3dy+k4dy)/6;
      const radius=Math.hypot(x,y);
      maxRadius=Math.max(maxRadius,radius);
      if(i%stride===0 || i===steps-1) samples.push({z:lengthMm*(i+1)/steps,x:x*1000,y:y*1000});
      if(!Number.isFinite(radius) || radius>=r) {transmitted=false;lostAtMm=lengthMm*(i+1)/steps;}
    }
    return {transmitted,lostAtMm,flight,cycles,vz,ap,qp,samples,maxRadiusMm:maxRadius*1000};
  }

  return {E,AMU,physicalFromAQ,aqFromPhysical,floquetTrace,xLimit,yLimit,
    apex,isFirstStable,crossings,centredBand,solveWidth,stabilityCurve,trajectory};
});
