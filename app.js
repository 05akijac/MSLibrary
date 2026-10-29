'use strict';
(() => {
  const Q=window.QuadMath;
  const $=id=>document.getElementById(id);
  const read=id=>Number($(id).value);
  const fmt=(n,d=3)=>Number.isFinite(n) ? n.toLocaleString('en-GB',{maximumFractionDigits:d}) : '—';
  let selected=null, trace=null, scan=[];
  const colours={grid:'#26394c',text:'#9bafc1',cyan:'#72d8e8',amber:'#f4bc78',mint:'#3fc7ab',rose:'#fb8199',white:'#f9fcff',bg:'#0c1a2b'};
  const getCanvas=(id)=>{
    const canvas=$(id),rect=canvas.getBoundingClientRect();
    if(rect.width<1)return null;
    const dpi=Math.min(2,window.devicePixelRatio||1);
    canvas.width=Math.round(rect.width*dpi);canvas.height=Math.round(rect.height*dpi);
    const c=canvas.getContext('2d');c.setTransform(dpi,0,0,dpi,0,0);
    return {c,w:rect.width,h:rect.height};
  };
  function basePlot(ctx,xMax,yMax,yTicks,xTicks,xLabel,yLabel, margins={left:52,right:17,top:17,bottom:37}, tickLabels={}){
    const {c,w,h}=ctx;const {left,right,top,bottom}=margins;
    const x0=left,y0=h-bottom,pw=w-left-right,ph=h-top-bottom;
    const X=v=>x0+v/xMax*pw,Y=v=>y0-v/yMax*ph;
    c.fillStyle=colours.bg;c.fillRect(x0,top,pw,ph);
    c.strokeStyle=colours.grid;c.lineWidth=1;c.font='11px ui-sans-serif,system-ui';c.textBaseline='middle';
    for(let j=0;j<=yTicks;j++){let val=yMax*j/yTicks;let py=Y(val);c.beginPath();c.moveTo(x0,py);c.lineTo(x0+pw,py);c.stroke();c.fillStyle=colours.text;c.textAlign='right';c.fillText(tickLabels.y?tickLabels.y(val):fmt(val,yMax<1?2:0),left-10,py);}
    for(let j=0;j<=xTicks;j++){let val=xMax*j/xTicks;let px=X(val);c.beginPath();c.moveTo(px,top);c.lineTo(px,y0);c.stroke();c.fillStyle=colours.text;c.textAlign='center';c.fillText(tickLabels.x?tickLabels.x(val):fmt(val,xMax<2?2:0),px,y0+14);}
    c.textAlign='center';c.fillStyle=colours.text;c.fillText(xLabel,x0+pw/2,h-9);
    c.save();c.translate(12,top+ph/2);c.rotate(-Math.PI/2);c.fillText(yLabel,0,0);c.restore();
    c.save();c.beginPath();c.rect(x0,top,pw,ph);c.clip();
    return {c,X,Y,x0,y0,top,pw,ph,end:()=>c.restore()};
  }
  function linePath(c,data,fnX,fnY,color,weight=1.7){
    c.beginPath();let started=false;
    for(let i=0;i<data.length;i++){const point=data[i];const x=fnX(point),y=fnY(point);if(!Number.isFinite(x+y)){started=false;continue;}if(!started){c.moveTo(x,y);started=true;}else c.lineTo(x,y);}
    c.strokeStyle=color;c.lineWidth=weight;c.stroke();
  }
  function drawStability(){
    if(!selected)return;
    const ctx=getCanvas('stabilityPlot');if(!ctx)return;
    const P=basePlot(ctx,.95,.30,6,5,'Mathieu q','Mathieu a');
    const {c,X,Y}=P,cur=Q.stabilityCurve(160);
    c.beginPath();c.moveTo(X(0),Y(0));cur.forEach(p=>c.lineTo(X(p.q),Y(p.a)));c.closePath();
    c.fillStyle='rgba(63,199,171,.17)';c.fill();
    linePath(c,cur,p=>X(p.q),p=>Y(p.a),colours.mint,2);
    const s=selected.q?selected.a/selected.q:0;
    c.beginPath();c.moveTo(X(0),Y(0));c.lineTo(X(.95),Y(.95*s));c.strokeStyle=colours.amber;c.lineWidth=1.5;c.setLineDash([6,5]);c.stroke();c.setLineDash([]);
    c.beginPath();c.arc(X(selected.q),Y(selected.a),5,0,Math.PI*2);c.fillStyle=colours.white;c.fill();c.lineWidth=2;c.strokeStyle=colours.amber;c.stroke();
    P.end();
    $('apexNote').textContent=`Numerical apex: q = ${fmt(Q.apex().q,6)}, a = ${fmt(Q.apex().a,6)}`;
  }
  function drawTrajectory(){
    const ctx=getCanvas('trajectoryPlot');if(!ctx||!selected)return;
    const rad=selected.radiusMm;const P=basePlot(ctx,selected.lengthMm,2*rad,4,5,'Axial distance z (mm)','Displacement (mm)',{left:58,right:17,top:15,bottom:36},{y:val=>fmt(val-rad,1)});
    const {c,X,Y}=P;
    const ypos=n=>Y(n+rad);
    c.setLineDash([5,5]);c.strokeStyle='#6d849b';c.lineWidth=1;
    for (let n of [-rad,rad]) {c.beginPath();c.moveTo(X(0),ypos(n));c.lineTo(X(selected.lengthMm),ypos(n));c.stroke();}
    c.setLineDash([]);
    if(trace){linePath(c,trace.samples,p=>X(p.z),p=>ypos(p.x),colours.cyan,1.7);linePath(c,trace.samples,p=>X(p.z),p=>ypos(p.y),colours.rose,1.7);}
    P.end();
  }
  function drawScan(){
    const ctx=getCanvas('scanPlot');if(!ctx||!scan.length)return;
    const m0=scan[0].mz,m1=scan[scan.length-1].mz,max=Math.max(scan[scan.length-1].V,scan[scan.length-1].U)*1.08||1;
    // Use a normalized x axis internally; display exact start and stop as labels.
    const P=basePlot(ctx,1,max,4,5,'m/z (Th)','Voltage (V)',{left:65,right:22,top:16,bottom:38},{x:val=>fmt(m0+(m1-m0)*val,1)});
    const {c,X,Y}=P;
    linePath(c,scan,p=>X((p.mz-m0)/(m1-m0)),p=>Y(p.V),colours.cyan,2.2);
    linePath(c,scan,p=>X((p.mz-m0)/(m1-m0)),p=>Y(p.U),colours.amber,2.2);
    P.end();
  }
  function renderResult(){
    if(!selected)return;
    const o=selected;
    $('vRf').textContent=fmt(o.V,2);$('vDc').textContent=fmt(o.U,2);
    $('vPp').textContent=fmt(o.VppPerRod,2);$('vDiff').textContent=fmt(o.differentialRFpp,2);
    $('outQ').textContent=fmt(o.q,6);$('outA').textContent=fmt(o.a,6);
    $('outSlope').textContent=o.q?fmt(o.a/o.q,6):'—';
    const band=o.band;
    $('outBand').textContent=band&&Number.isFinite(band.width)?fmt(band.width,4):'not isolated';
    $('outBandNote').textContent=band&&Number.isFinite(band.width) ? `Th · R = ${fmt(band.resolution,1)}` : 'RF-only / outside finite ideal passband';
    const stable=Q.isFirstStable(o.a,o.q);
    $('status').textContent=o.rfOnly?'RF-only mode':stable?'Region I stable':'Outside Region I';
    $('status').classList.toggle('warn',!stable);
    $('regionNote').textContent=`${stable?'Selected point inside':'Selected point outside'} region I`;
    drawStability();drawTrajectory();drawScan();
  }
  function getBand(mz,a,q){
    if(!(q>0 && a>=0))return null;
    try{
      const slope=a/q;
      if(slope>=Q.apex().slope)return null;
      const window=Q.centredBand(mz,slope);
      // For manually chosen a,q, actual band endpoints depend on selected q,
      // not on the *centred* window operating point.
      const low=mz*q/window.qHigh,high=window.qLow===0?Infinity:mz*q/window.qLow;
      return {lowerMz:low,upperMz:high,width:high-low,resolution:mz/(high-low)};
    }catch{return null;}
  }
  function computeOperatingPoint(){
    $('error').textContent='';
    const mz=read('mz'),radiusMm=read('radius'),frequencyMHz=read('frequency');
    const rfOnly=$('rfOnly').checked;
    if(![mz,radiusMm,frequencyMHz].every(v=>Number.isFinite(v)&&v>0)) throw Error('m/z, field radius and frequency must be positive.');
    let a,q,band=null;
    if(rfOnly){a=0;q=read('guideQ');if(!(q>0&&q<0.9))throw Error('RF-only q must be between 0 and 0.9.');}
    else if($('mode').value==='manual'){
      a=read('manualA');q=read('manualQ');
      if(!Number.isFinite(a)||!Number.isFinite(q)||a<0||q<0)throw Error('Mathieu a and q must be nonnegative.');
    }else{
      const targetWidth=read('width');
      const solution=Q.solveWidth(mz,targetWidth);
      a=solution.a;q=solution.q;
    }
    const physical=Q.physicalFromAQ(mz,frequencyMHz,radiusMm,a,q);
    band=rfOnly?null:getBand(mz,a,q);
    selected={...physical,mz,radiusMm,frequencyMHz,band,rfOnly,
      lengthMm:read('length')};
    return selected;
  }
  function runTrajectory(){
    if(!selected)return;
    selected.lengthMm=read('length');
    try{
      trace=Q.trajectory({targetMz:selected.mz,testMz:read('testMz'),frequencyMHz:selected.frequencyMHz,
        radiusMm:selected.radiusMm,lengthMm:selected.lengthMm,axialEVPerCharge:read('energy'),
        a:selected.a,q:selected.q,phaseDeg:read('phase'),x0Mm:read('x0'),y0Mm:read('y0'),
        vx0:read('vx0'),vy0:read('vy0')});
      const p=$('transmissionResult');p.textContent=trace.transmitted?'TRANSMITTED · single-ion model':`LOST · at z = ${fmt(trace.lostAtMm,2)} mm`;
      p.className=trace.transmitted?'pass':'fail';
      $('cycleResult').textContent=`${fmt(trace.cycles,1)} RF cycles · vᶻ = ${fmt(trace.vz,0)} m/s`;
    }catch(err){trace=null;$('transmissionResult').textContent=err.message;$('transmissionResult').className='fail';$('cycleResult').textContent='';}
    drawTrajectory();
  }
  function updateScan(){
    if(!selected)return;
    const start=read('scanStart'),end=read('scanEnd');
    if(!(Number.isFinite(start)&&Number.isFinite(end)&&start>0&&end>start)){
      $('error').textContent='Scan end must be greater than scan start, both > 0.';return;
    }
    $('error').textContent='';scan=[];
    for(let i=0;i<=150;i++){
      const mz=start+(end-start)*i/150;
      const {U,V}=Q.physicalFromAQ(mz,selected.frequencyMHz,selected.radiusMm,selected.a,selected.q);
      scan.push({mz,U,V});
    }
    drawScan();
  }
  function downloadCsv(){
    if(!scan.length)return;
    const csv=['m_z_Th,RF_rod_peak_V,DC_rod_V,RF_differential_peak_to_peak_V,DC_differential_V,mathieu_a,mathieu_q',
      ...scan.map(p=>[p.mz,p.V,p.U,4*p.V,2*p.U,selected.a,selected.q].join(','))].join('\n');
    const blob=new Blob([csv],{type:'text/csv;charset=utf-8'});const url=URL.createObjectURL(blob);
    const a=document.createElement('a');a.href=url;a.download='quadrupole-mass-scan.csv';a.click();
    setTimeout(()=>URL.revokeObjectURL(url),1000);
  }
  function calculate(){
    try{
      computeOperatingPoint();$('testMz').value=String(selected.mz);
      renderResult();runTrajectory();updateScan();
    }catch(err){$('error').textContent=err.message;}
  }
  function setMode(){
    const rf=$('rfOnly').checked,manual=$('mode').value==='manual';
    $('guideQField').classList.toggle('hidden',!rf);
    $('widthFields').classList.toggle('hidden',rf||manual);
    $('manualFields').classList.toggle('hidden',rf||!manual);
    $('mode').disabled=rf;
  }
  $('settings').addEventListener('submit',e=>{e.preventDefault();calculate();});
  $('trajectorySettings').addEventListener('submit',e=>{e.preventDefault();runTrajectory();});
  $('updateScan').addEventListener('click',updateScan);
  $('downloadCsv').addEventListener('click',downloadCsv);
  $('mode').addEventListener('change',setMode);
  $('rfOnly').addEventListener('change',setMode);
  let resizeTask=0;
  window.addEventListener('resize',()=>{cancelAnimationFrame(resizeTask);resizeTask=requestAnimationFrame(()=>{drawStability();drawTrajectory();drawScan();});});
  setMode();calculate();
})();
