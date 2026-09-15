const fs=require('fs');
const path='/mnt/data/K4_reproducibility_raw_20260915';
function readCsv(file){
  const lines=fs.readFileSync(file,'utf8').trim().split(/\r?\n/); const h=lines[0].split(',');
  return lines.slice(1).map(line=>{const v=line.split(','); const o={}; h.forEach((k,i)=>o[k]=v[i]); return o;});
}
let all=true, details=[];
for(const m of [2,3]){
  const w=readCsv(`${path}/wave_samples_html_800_m${m}.csv`);
  const n=800, W={x0:70,x1:1010,yBase:100,yHarm:215,ySum:335,amp:48};
  for(const [key,center,amp,fn] of [
    ['base',W.yBase,W.amp,t=>Math.cos(t)],
    ['harmonic',W.yHarm,W.amp,t=>Math.cos(m*t)],
    ['composite',W.ySum,W.amp*0.72,t=>Math.cos(t)+Math.cos(m*t)]]){
      let ok=true,max=0;
      for(let i=0;i<n;i++){
        const theta=2*Math.PI*i/(n-1);
        const raw=Number(w[i][key]);
        max=Math.max(max,Math.abs(raw-fn(theta)));
        const y0=(center-amp*fn(theta)).toFixed(2);
        const y1=(center-amp*raw).toFixed(2);
        if(y0!==y1){ok=false; break;}
      }
      details.push({m,kind:`path_${key}`,coordinate_string_exact:ok,max_raw_value_diff:max}); all=all&&ok;
  }
  const s=readCsv(`${path}/complex_states_animation_exact_i.csv`).filter(r=>Number(r.m)===m && r.direction==='forward').sort((a,b)=>Number(a.step)-Number(b.step));
  for(let k=0;k<4;k++){
    const r=s[k], theta=k*Math.PI/2;
    const o={base_re:Math.round(Math.cos(theta)),base_im:Math.round(Math.sin(theta)),harmonic_re:Math.round(Math.cos(m*theta)),harmonic_im:Math.round(Math.sin(m*theta))};
    const q={base_re:Number(r.base_re),base_im:Number(r.base_im),harmonic_re:Number(r.harmonic_re),harmonic_im:Number(r.harmonic_im)};
    // coordinates in original rendering; signed zero is visually identical for coordinate arithmetic here.
    const ok=Object.keys(o).every(x=>Object.is(o[x],q[x]) || o[x]===q[x]);
    const ampOk=(Math.cos(theta)===q.base_re || Math.abs(Math.cos(theta)-q.base_re)<1e-15) && (Math.cos(m*theta)===q.harmonic_re || Math.abs(Math.cos(m*theta)-q.harmonic_re)<1e-15);
    details.push({m,k,kind:'state_vectors',exact_numeric:ok,amplitude_equivalent:ampOk}); all=all&&ok&&ampOk;
  }
}
const report={all_html_draw_coordinates_exact:all,details};
fs.writeFileSync('/mnt/data/K4_html_coordinate_validation_20260915.json',JSON.stringify(report,null,2));
console.log(JSON.stringify(report,null,2));
