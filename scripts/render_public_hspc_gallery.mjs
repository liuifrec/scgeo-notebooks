// Deterministic SVGs from compact, audited public HSPC result tables.
import { readFileSync, writeFileSync, mkdirSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { join } from "node:path";

function renderHspc(d) {
  const ink="#24384B",muted="#607386",grid="#DFE7EB",tnf="#3975A9",lps="#D18B48",neg="#B76B43",pos="#286F9E";
  const colors=[tnf,lps], short=["Quiescent HSC","Activated HSC","MPP/progenitor","Megakaryocyte/erythroid","Myeloid","Lymphoid","Ambiguous HSPC"];
  const esc=s=>String(s).replaceAll("&","&amp;").replaceAll("<","&lt;").replaceAll(">","&gt;").replaceAll('"',"&quot;");
  const r=v=>Number(v.toFixed(3));
  const txt=(x,y,s,o={})=>`<text x="${r(x)}" y="${r(y)}" font-size="${o.size||14}" fill="${o.color||ink}" text-anchor="${o.anchor||"start"}"${o.bold?' font-weight="600"':""}>${esc(s)}</text>`;
  const rect=(x,y,w,h,f,rx=0)=>`<rect x="${r(x)}" y="${r(y)}" width="${r(w)}" height="${r(h)}" fill="${f}" rx="${rx}"/>`;
  const line=(x1,y1,x2,y2,c=grid,w=1)=>`<line x1="${r(x1)}" y1="${r(y1)}" x2="${r(x2)}" y2="${r(y2)}" stroke="${c}" stroke-width="${w}"/>`;
  const dot=(x,y,rad,c)=>`<circle cx="${r(x)}" cy="${r(y)}" r="${rad}" fill="${c}" stroke="white" stroke-width="1.2"/>`;
  const begin=(w,h,title,desc,subtitle)=>`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${w} ${h}" width="${w}" height="${h}" role="img" aria-labelledby="title description"><title id="title">${esc(title)}</title><desc id="description">${esc(desc)}</desc><rect width="100%" height="100%" fill="white"/><g font-family="Arial, Helvetica, sans-serif">${txt(26,35,title,{size:24,bold:true})}${txt(26,59,subtitle,{size:14,color:muted})}`;
  const end=s=>s+"</g></svg>\n";
  const key=(x,y,c,s)=>rect(x,y-11,13,11,c,2)+txt(x+20,y,s,{size:13});
  const p=(v,lo,hi,x,w)=>x+w*(v-lo)/(hi-lo);
  const ticks=(x,y,w,h,lo,hi,vals,format=v=>String(v))=>vals.map(v=>line(p(v,lo,hi,x,w),y,p(v,lo,hi,x,w),y+h,v===0?"#91A2AD":grid,v===0?1.3:1)+txt(p(v,lo,hi,x,w),y+h+21,format(v),{anchor:"middle",size:12,color:muted})).join("");
  const sums=new Map(d.summary.map(row=>[row[0]+":"+row[1],row]));
  if(d.summary.length!==14||sums.size!==14||d.directions.length!==7||d.genes.length!==36||d.programs.length!==14||d.program_names.length!==18)throw Error("Incomplete HSPC source data");
  if(d.validation.biological_samples!==0||d.validation.bootstrap_draws!==0||d.validation.status!=="insufficient_samples")throw Error("Pooled-only contract broken");
  const out={};
  {
    let s=begin(1160,586,"Distinct responses to TNF and LPS","Robust displacement magnitudes for seven human HSPC states under TNF and LPS relative to PBS, plus approximate directional cosines between the two shifts.","GSE249479 • original verified SVD50 basis • pooled descriptive analysis only");
    s+=key(252,86,tnf,"PBS → TNF")+key(419,86,lps,"PBS → LPS");
    const x=255,w=420,cx=785,cw=300,y=133,step=52,h=7*step;
    s+=txt(x,115,"Robust displacement magnitude",{bold:true})+txt(cx,115,"Directional cosine",{bold:true});
    s+=ticks(x,y,w,h,0,4.5,[0,1,2,3,4])+ticks(cx,y,cw,h,-1,1,[-1,-.5,0,.5,1]);
    for(let i=0;i<7;i++){
      let yy=y+i*step+8;
      s+=txt(x-14,yy+21,short[i],{anchor:"end",size:14,bold:i<3});
      for(let j=0;j<2;j++){const v=sums.get(i+":"+j)[4];s+=rect(x,yy+j*17,v/4.5*w,12,colors[j],2)+txt(x+v/4.5*w+7,yy+j*17+10,v.toFixed(2),{size:11,color:muted});}
      const c=d.directions[i][0],xp=p(c,-1,1,cx,cw);
      s+=line(p(0,-1,1,cx,cw),yy+15,xp,yy+15,pos,3)+dot(xp,yy+15,5,pos)+txt(1135,yy+20,c.toFixed(3),{anchor:"end",size:12});
    }
    s+=txt(26,547,"Cosines compare the two fitted displacement directions, not the treatment labels or trajectories.",{size:12,color:muted})+txt(26,568,"No biological replicate identities: no significance tests, intervals, or population-level claims.",{size:12,color:muted});
    out["01_geometry.svg"]=end(s);
  }
  {
    let s=begin(1080,414,"Similar magnitude need not mean similar direction","Schematic pairs of vectors for quiescent HSC, activated HSC and lymphoid-labeled populations. Length and angle reflect the observed SVD shifts.","Common-origin schematics, not actual embedding coordinates or cell trajectories.");
    s+=key(30,87,tnf,"PBS → TNF")+key(230,87,lps,"PBS → LPS");
    [0,1,5].forEach((i,k)=>{
      const ox=52+k*350,oy=324,scale=45,cos=d.directions[i][0],theta=Math.acos(cos),a=sums.get(i+":0")[4],b=sums.get(i+":1")[4];
      s+=txt(ox,124,short[i],{size:17,bold:true});
      const arrow=(x2,y2,c)=>{const a=Math.atan2(y2-oy,x2-ox),b=10;return line(ox,oy,x2,y2,c,3)+`<path d="M ${r(x2-b*Math.cos(a-.45))} ${r(y2-b*Math.sin(a-.45))} L ${r(x2)} ${r(y2)} L ${r(x2-b*Math.cos(a+.45))} ${r(y2-b*Math.sin(a+.45))}" fill="none" stroke="${c}" stroke-width="3" stroke-linecap="round"/>`;};
      s+=arrow(ox+a*scale,oy,tnf)+arrow(ox+b*scale*cos,oy-b*scale*Math.sin(theta),lps)+dot(ox,oy,4,ink);
      const rad=27;
      s+=`<path d="M ${ox+rad} ${oy} A ${rad} ${rad} 0 0 0 ${r(ox+rad*cos)} ${r(oy-rad*Math.sin(theta))}" fill="none" stroke="#99A9B5" stroke-width="1.5"/>`;
      s+=txt(ox+42,oy-31,(theta*180/Math.PI).toFixed(1)+"°",{bold:true,size:17})+txt(ox,oy+23,"PBS",{size:12,color:muted})+txt(ox,oy+48,`TNF ${a.toFixed(2)} · LPS ${b.toFixed(2)}`,{size:13});
    });
    s+=txt(1050,399,"Common scale: 45 px / SVD unit",{anchor:"end",size:11,color:muted});
    out["02_direction.svg"]=end(s);
  }
  const mix=(v,lim)=>{const q=Math.min(1,Math.abs(v)/lim),a=[249,250,251],b=v<0?[183,107,67]:[40,111,158];return "#"+a.map((u,i)=>Math.round(u+(b[i]-u)*q).toString(16).padStart(2,"0")).join("");};
  const legend=(x,y,lim)=>{let s="";for(let i=0;i<80;i++)s+=rect(x+i*2.5,y,2.5,11,mix(-lim+2*lim*i/79,lim));[-lim,0,lim].forEach((v,i)=>s+=txt(x+i*100,y+29,String(v),{anchor:"middle",size:12,color:muted}));return s};
  const genes=new Map(d.genes.map(row=>[row.slice(0,3).join(":"),row]));
  {
    let s=begin(1260,445,"Molecular accounting of HSPC responses","Two heatmaps show robust expression-center changes and signed contributions for six selected genes across three focus states and two contrasts.","Retrospectively selected genes • independent scales for expression and contribution • no replicate uncertainty");
    const rowNames=d.gene_names,y=134,w=85,h=36,xs=[145,705],lims=[1,.4],titles=["Expression-center change","Signed geometric contribution"];
    for(let k=0;k<2;k++){
      const x=xs[k],lim=lims[k];s+=txt(x,87,titles[k],{size:17,bold:true});
      for(let i=0;i<3;i++){
        s+=txt(x+(2*i+1)*w,y-27,["Quiescent HSC","Activated HSC","MPP"][i],{anchor:"middle",size:13,bold:true});
        for(let j=0;j<2;j++)s+=txt(x+(2*i+j+.5)*w,y-8,["TNF","LPS"][j],{anchor:"middle",size:12,color:muted});
      }
      for(let g=0;g<6;g++){
        if(k===0)s+=txt(x-12,y+(g+.5)*h+5,rowNames[g],{anchor:"end",size:14,bold:true});
        for(let i=0;i<3;i++)for(let j=0;j<2;j++){const row=genes.get(`${i}:${j}:${g}`),v=row[k===0?3:5],xx=x+(i*2+j)*w,yy=y+g*h;s+=rect(xx+1,yy+1,w-2,h-2,mix(v,lim),2)+txt(xx+w/2,yy+h/2+4,v.toFixed(3),{anchor:"middle",size:12,color:Math.abs(v)>.65*lim?"white":ink});}
      }
      s+=legend(x+155,365,lim);
    }
    s+=txt(26,427,"Contribution = expression-center change × loading along the fitted direction. Values near zero may round to 0.000.",{size:12,color:muted});
    out["03_genes.svg"]=end(s);
  }
  {
    let s=begin(1130,388,"FOS: the same expression-change sign, different contributions","FOS decreases in all six selected state-contrast rows, yet its small TNF contributions are negative while LPS contributions are positive.","Exploratory example • two separately scaled axes • no biological replicate intervals");
    const y=115,h=190,step=60,x=256,w=300,x2=680,w2=375;
    s+=key(26,85,tnf,"TNF")+key(128,85,lps,"LPS")+txt(x,97,"Expression-center change",{bold:true})+txt(x2,97,"Signed contribution",{bold:true});
    s+=ticks(x,y,w,h,-.30,0,[-.3,-.2,-.1,0],v=>v.toFixed(1))+ticks(x2,y,w2,h,-.003,.022,[-.002,0,.005,.01,.015,.02],v=>v.toFixed(3));
    for(let i=0;i<3;i++){
      const yy=y+12+i*step;
      s+=txt(x-14,yy+22,short[i],{anchor:"end",bold:true,size:14});
      for(let j=0;j<2;j++){
        const row=genes.get(`${i}:${j}:5`),e=row[3],c=row[5],by=yy+j*22,ex=p(e,-.30,0,x,w),z=p(0,-.003,.022,x2,w2),cx=p(c,-.003,.022,x2,w2);
        s+=rect(ex,by,w-(ex-x),12,colors[j],2)+txt(x+w+8,by+10,e.toFixed(3),{size:11,color:muted});
        s+=rect(Math.min(z,cx),by,Math.max(.7,Math.abs(cx-z)),12,colors[j],2)+txt(c>=0?cx+5:cx-5,by+10,c.toFixed(4),{anchor:c>=0?"start":"end",size:11});
      }
    }
    s+=txt(26,368,"The opposing TNF terms are numerically small; a sign flip is geometric accounting, not a causal or statistically established effect.",{size:12,color:muted});
    out["04_fos.svg"]=end(s);
  }
  {
    let s=begin(1045,554,"Predefined program sums reveal different response components","Eight selected programs across quiescent HSC, activated HSC and MPP, comparing TNF and LPS.","HVG-restricted, overlapping programs • one shared signed scale • pooled descriptive estimates");
    const selected=[0,1,2,10,12,13,14,16],labels=["Study quiescence top-50","Study inflammation top-50","Study HSC-IM top-50","TNF/NF-κB response","AP-1 stress response","Interferon response","Quiescence markers","Immediate early"];
    const x=337,y=134,w=110,h=41,lim=.85, rows=new Map(d.programs.map(row=>[row[0]+":"+row[1],row]));
    for(let i=0;i<3;i++){s+=txt(x+(2*i+1)*w,y-27,["Quiescent HSC","Activated HSC","MPP"][i],{anchor:"middle",size:15,bold:true});for(let j=0;j<2;j++)s+=txt(x+(2*i+j+.5)*w,y-8,["TNF","LPS"][j],{anchor:"middle",size:12,color:muted});}
    selected.forEach((g,k)=>{
      const yy=y+k*h;
      s+=txt(x-12,yy+h/2+5,labels[k]+" ("+d.program_sizes[g]+"/"+d.program_requested[g]+")",{anchor:"end",size:14,bold:g===2});
      for(let i=0;i<3;i++)for(let j=0;j<2;j++){
        const v=rows.get(i+":"+j)[g+2],xx=x+(i*2+j)*w;
        s+=rect(xx+1,yy+1,w-2,h-2,mix(v,lim),2)+txt(xx+w/2,yy+h/2+5,v.toFixed(3),{anchor:"middle",size:13,color:Math.abs(v)>.65*lim?"white":ink});
      }
    });
    s+=legend(785,475,lim)+txt(26,521,"Parentheses: measured HVGs / requested members. Program sums overlap and do not partition the net displacement.",{size:12,color:muted})+txt(26,540,"HSC-IM denotes the original study's selected HSC-II marker subset; this is not evidence of persistent inflammatory memory.",{size:12,color:muted});
    out["05_programs.svg"]=end(s);
  }
  return out;
}

const root = fileURLToPath(new URL("../notebooks/applications/results/gse249479/", import.meta.url));
const d = JSON.parse(readFileSync(join(root,"selected_results.json"),"utf8"));
const check = process.argv.includes("--check");
const dir = join(root,"figures");
if (!check) mkdirSync(dir,{recursive:true});
for (const [name, svg] of Object.entries(renderHspc(d))) {
  const path = join(dir,name);
  if (check) { if (readFileSync(path,"utf8") !== svg) throw Error("Figure differs: " + name); }
  else writeFileSync(path,svg);
}
console.log("PASS: five public HSPC result snapshots");
