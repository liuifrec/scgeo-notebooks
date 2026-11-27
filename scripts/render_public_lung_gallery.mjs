// Deterministic SVG snapshots from compact public result tables. No external dependencies.
import { readFileSync, writeFileSync, mkdirSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { join } from "node:path";

function renderGallery(d) {
  const colors=["#3975A9","#D18B48","#458C83"], ink="#24384B",muted="#617486",grid="#DDE5EA";
  const short=["Endothelial","Epithelial","Fibroblast/stromal","Lymphoid","Myeloid","Proliferating","Col13a1-like","Col14a1-like","Myofibroblast"];
  const contrasts=["Control → early","Control → late","Early → late"];
  const out={};
  const esc=s=>String(s).replaceAll("&","&amp;").replaceAll("<","&lt;").replaceAll(">","&gt;").replaceAll('"',"&quot;");
  const n=v=>Number(v.toFixed(3));
  const txt=(x,y,s,opt={})=>`<text x="${n(x)}" y="${n(y)}" font-size="${opt.size||14}" fill="${opt.color||ink}" text-anchor="${opt.anchor||"start"}"${opt.bold?' font-weight="600"':""}>${esc(s)}</text>`;
  const line=(x1,y1,x2,y2,color=grid,w=1,dash="")=>`<line x1="${n(x1)}" y1="${n(y1)}" x2="${n(x2)}" y2="${n(y2)}" stroke="${color}" stroke-width="${w}"${dash?' stroke-dasharray="'+dash+'"':""}/>`;
  const rect=(x,y,w,h,fill,rx=0)=>`<rect x="${n(x)}" y="${n(y)}" width="${n(w)}" height="${n(h)}" fill="${fill}" rx="${rx}"/>`;
  const dot=(x,y,r,color)=>`<circle cx="${n(x)}" cy="${n(y)}" r="${r}" fill="${color}" stroke="white" stroke-width="1.4"/>`;
  const begin=(w,h,title,desc,subtitle)=>`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${w} ${h}" width="${w}" height="${h}" role="img" aria-labelledby="title description"><title id="title">${esc(title)}</title><desc id="description">${esc(desc)}</desc><rect width="100%" height="100%" fill="white"/><g font-family="Arial, Helvetica, sans-serif">${txt(28,34,title,{size:23,bold:true})}${txt(28,57,subtitle,{size:14,color:muted})}`;
  const end=s=>s+"</g></svg>\n";
  const key=(x,y,c,s)=>rect(x,y-10,13,10,c,2)+txt(x+20,y,s,{size:13});
  const ticks=(x,y,w,h,lo,hi,values,format=v=>String(v))=>values.map(v=>{const p=x+(v-lo)/(hi-lo)*w;return line(p,y,p,y+h,v===0?"#8999A6":grid,v===0?1.4:1)+txt(p,y+h+19,format(v),{anchor:"middle",size:12,color:muted});}).join("");
  const pos=(v,lo,hi,x,w)=>x+(v-lo)/(hi-lo)*w;
  const interval=(x,lo,hi,w,y,l,v,h,c)=>line(pos(l,lo,hi,x,w),y,pos(h,lo,hi,x,w),y,c,2)+line(pos(l,lo,hi,x,w),y-4,pos(l,lo,hi,x,w),y+4,c,1.5)+line(pos(h,lo,hi,x,w),y-4,pos(h,lo,hi,x,w),y+4,c,1.5)+dot(pos(v,lo,hi,x,w),y,5,c);
  const sum=new Map(d.summary.map(r=>[r[0]+":"+r[1],r]));
  const dir=new Map(d.directions.map(r=>[r[0],r]));
  const maxError=Math.max(...d.summary.map(r=>Math.abs(r[6]-r[7])));
  if(d.summary.length!==27||sum.size!==27||maxError>=2e-5||d.directions.length!==9||d.genes.length!==18||d.programs.length!==63||d.mice.length!==24) throw Error("Unexpected public result roster");
  // 01. Complete state roster, magnitudes and directional cosine.
  {
    let s=begin(1120,700,"Radiation-associated geometric responses","All nine major and refined states, three observed PCA50 displacement magnitudes per state, and cosine between control-to-early and control-to-late axes.","GSE211713 • original verified PCA50 basis • cross-sectional mice");
    contrasts.forEach((v,i)=>s+=key(238+i*180,87,colors[i],v));
    s+=txt(200,118,"Robust displacement magnitude",{bold:true})+txt(808,118,"Early–late directional cosine",{bold:true});
    const y0=142,step=55,x=200,w=505,cx=810,cw=250,h=9*step;
    s+=ticks(x,y0,w,h,0,16,[0,4,8,12,16])+ticks(cx,y0,cw,h,-1,1,[-1,-.5,0,.5,1]);
    for(let i=0;i<9;i++){
      const y=y0+i*step+10;
      if(i===6)s+=line(25,y-12,1080,y-12,"#A7B7C2",1.5);
      s+=txt(185,y+19,short[i],{anchor:"end",size:14,bold:i>=6});
      for(let j=0;j<3;j++){
        const v=sum.get(i+":"+j)[6];s+=rect(x,y+j*12,v/16*w,9,colors[j],2);
      }
      const v=dir.get(i)[1];s+=line(pos(0,-1,1,cx,cw),y+14,pos(v,-1,1,cx,cw),y+14,v<0?"#C27B48":"#3975A9",3)+dot(pos(v,-1,1,cx,cw),y+14,5,v<0?"#C27B48":"#3975A9")+txt(1090,y+19,v.toFixed(2),{anchor:"end",size:12});
    }
    s+=txt(28,679,"Top: six major states. Bottom: three refined fibroblast states. Magnitudes are not standardized effect sizes; direction estimates have no angular intervals.",{size:12,color:muted});
    out["01_geometry.svg"]=end(s);
  }
  // 02. Common-origin two-dimensional direction diagrams preserving magnitudes and angle.
  {
    let s=begin(1050,400,"Directional continuity differs among fibroblasts","Three common-origin schematics preserve the control-to-early and control-to-late vector lengths and the estimated angle between them, not their original spatial coordinates.","Not observed trajectories or UMAP coordinates. All three panels share the same length scale.");
    s+=key(30,86,colors[0],contrasts[0])+key(230,86,colors[1],contrasts[1]);
    for(let i=6;i<=8;i++){
      const ox=58+(i-6)*346,oy=320,scale=21,cos=dir.get(i)[1],theta=Math.acos(cos),e=sum.get(i+":0")[6],l=sum.get(i+":1")[6];
      s+=txt(ox,119,short[i],{size:17,bold:true});
      const arrow=(x2,y2,c)=>{
        const a=Math.atan2(y2-oy,x2-ox),b=10;
        const p1=[x2-b*Math.cos(a-.45),y2-b*Math.sin(a-.45)],p2=[x2-b*Math.cos(a+.45),y2-b*Math.sin(a+.45)];
        return line(ox,oy,x2,y2,c,3)+`<path d="M ${n(p1[0])} ${n(p1[1])} L ${n(x2)} ${n(y2)} L ${n(p2[0])} ${n(p2[1])}" fill="none" stroke="${c}" stroke-width="3" stroke-linecap="round"/>`;
      };
      s+=arrow(ox+e*scale,oy,colors[0])+arrow(ox+l*scale*cos,oy-l*scale*Math.sin(theta),colors[1])+dot(ox,oy,4,ink);
      const r=29;s+=`<path d="M ${ox+r} ${oy} A ${r} ${r} 0 0 0 ${n(ox+r*cos)} ${n(oy-r*Math.sin(theta))}" fill="none" stroke="#A2ACB7" stroke-width="1.5"/>`;
      s+=txt(ox+45,oy-36,(theta*180/Math.PI).toFixed(1)+"°",{bold:true,size:16})+txt(ox,oy+23,"Control",{size:12,color:muted})+txt(ox,oy+49,`Early ${e.toFixed(2)} · late ${l.toFixed(2)}`,{size:13});
    }
    s+=txt(1015,382,"Vector scale: 21 px / PCA unit",{anchor:"end",size:11,color:muted});
    out["02_direction.svg"]=end(s);
  }
  // 03. Opposite contribution signs despite increased expression.
  {
    let s=begin(1060,338,"Col3a1: expression change is not geometric contribution","All three fibroblast populations show increased robust normalized expression for Col3a1 in the control-to-early contrast; two have negative signed contributions and one positive contribution.","Control → early • selected exploratory example • conditional sample-bootstrap intervals");
    const x=200,w=320,x2=600,w2=345,y=105,h=152;
    s+=txt(x,90,"Expression-center change",{bold:true})+txt(x2,90,"Signed geometric contribution",{bold:true});
    s+=ticks(x,y,w,h,0,1.1,[0,.25,.5,.75,1])+ticks(x2,y,w2,h,-.06,.10,[-.05,0,.05,.10],v=>v.toFixed(2));
    d.genes.filter(r=>r[1]===0&&r[2]===2).forEach((r,i)=>{const yy=y+25+i*49,c=r[5]<0?"#C27B48":"#3975A9";s+=txt(185,yy+5,short[r[0]],{anchor:"end",bold:true})+rect(x,yy-6,r[3]/1.1*w,12,"#7896AB",2)+txt(x+r[3]/1.1*w+8,yy+4,r[3].toFixed(3),{size:12})+interval(x2,-.06,.10,w2,yy,r[6],r[5],r[7],c)+txt(1030,yy+4,r[5].toFixed(4),{anchor:"end",size:13,color:c});});
    s+=txt(28,316,"Contribution = robust expression-center change × axis loading. Interval includes zero for Col13a1-like; this is not a causal effect.",{size:12,color:muted});
    out["03_col3a1.svg"]=end(s);
  }
  // 04. Late cross-lineage markers. Keep separate expression and signed contributions.
  {
    let s=begin(1060,475,"Late response: cross-lineage marker diagnostic","Selected Spp1, Igkc and Sftpc expression-center changes and signed contributions across refined fibroblast states, early versus late.","Early → late • selected after reviewing the full results • inspect cellular origin before biological attribution");
    const x=225,w=305,x2=620,w2=340,y=108,h=305;
    s+=txt(x,91,"Expression-center change",{bold:true})+txt(x2,91,"Signed contribution (conditional interval)",{bold:true});
    s+=ticks(x,y,w,h,-2.2,2.2,[-2,-1,0,1,2])+ticks(x2,y,w2,h,-.1,1,[-.1,0,.25,.5,.75,1]);
    const rows=d.genes.filter(r=>r[1]===2);
    rows.forEach((r,i)=>{const yy=y+18+i*32+(i>=3?5:0)+(i>=6?5:0);
      if(i%3===0){if(i)s+=line(22,yy-17,1037,yy-17,"#B9C7CF");s+=txt(26,yy+4,short[r[0]],{size:13,bold:true});}
      s+=txt(206,yy+5,d.gene_names[r[2]],{anchor:"end",size:13});
      const a=pos(0,-2.2,2.2,x,w),b=pos(r[3],-2.2,2.2,x,w);s+=rect(Math.min(a,b),yy-6,Math.abs(b-a),12,r[3]<0?"#C27B48":"#7896AB",2);
      s+=interval(x2,-.1,1,w2,yy,r[6],r[5],r[7],"#3975A9")+txt(1036,yy+5,r[5].toFixed(3),{anchor:"end",size:12});
    });
    s+=txt(27,460,"Potential ambient RNA, doublets, annotation or substate composition: these selected signals are not established fibroblast-intrinsic drivers.",{size:12,color:muted});
    out["04_late_markers.svg"]=end(s);
  }
  // 05. Selected program sums, single centered scale.
  {
    let s=begin(1140,500,"Program contributions across the three contrasts","Seven selected predefined HVG-restricted programs in three fibroblast states. Every tile is the sum of signed gene contributions, not pathway activity or enrichment.","One shared diverging scale (−0.8 to +0.8 PCA units). Programs can overlap.");
    const x=302,y=126,w=88,h=39, names=["Fibroblast/stromal (7)","Col13a1-like matrix (4)","Col14a1-like matrix (4)","Myofibroblast (5)","Epithelial (10)","B cell (2)","Macrophage (6)"];
    for(let i=0;i<3;i++){s+=txt(x+(i*3+1.5)*w,y-31,short[i+6],{anchor:"middle",bold:true,size:15});for(let j=0;j<3;j++)s+=txt(x+(i*3+j+.5)*w,y-9,["C→E","C→L","E→L"][j],{anchor:"middle",size:12,color:muted});}
    const mix=(v)=>{const q=Math.min(1,Math.abs(v)/.8),a=[249,250,251],b=v<0?[191,113,64]:[40,111,158];return "#"+a.map((u,i)=>Math.round(u+(b[i]-u)*q).toString(16).padStart(2,"0")).join("");};
    for(let p=0;p<7;p++){
      s+=txt(x-15,y+(p+.5)*h+5,names[p],{anchor:"end",bold:p===4,size:14});
      for(let i=6;i<=8;i++)for(let j=0;j<3;j++){const r=d.programs.find(a=>a[0]===i&&a[1]===j&&a[2]===p),v=r[4],xx=x+((i-6)*3+j)*w,yy=y+p*h;s+=rect(xx+1,yy+1,w-2,h-2,mix(v),2)+txt(xx+w/2,yy+h/2+5,(v<0?"−":"")+Math.abs(v).toFixed(3),{anchor:"middle",color:Math.abs(v)>.53?"white":ink,size:13});}
      if(p===3)s+=line(28,y+4*h,1100,y+4*h,"#A6B5C0",1.5);
    }
    const lx=774,ly=422;for(let i=0;i<80;i++)s+=rect(lx+i*3,ly,3,10,mix(-.8+1.6*i/79));[-.8,0,.8].forEach((v,i)=>s+=txt(lx+i*120,ly+26,String(v),{anchor:"middle",size:12}));
    s+=txt(28,467,"C→E: control to early; C→L: control to late; E→L: early to late. Parentheses: measured program members.",{size:12,color:muted})+txt(28,484,"Epithelial-marker contributions warrant source and annotation checks. Color does not encode statistical significance.",{size:12,color:muted});
    out["05_programs.svg"]=end(s);
  }
  // 06. Unpaired sample projections.
  {
    let s=begin(1050,470,"Individual mouse projections: late-state heterogeneity","Four independent early and four independent late mice per refined fibroblast state, projected onto each state's fitted early-to-late direction.","Early → late • no pairing • axis fitted using these same mice");
    const x=310,w=690,y=109,h=280,lo=-2,hi=14;
    s+=key(310,85,colors[0],"Early 17 Gy")+key(480,85,colors[1],"Late 17 Gy");
    s+=ticks(x,y,w,h,lo,hi,[-2,0,2,4,6,8,10,12,14]);
    for(let i=6;i<=8;i++){
      let yy=y+24+(i-6)*91;s+=txt(28,yy+17,short[i],{bold:true,size:16});
      for(let j=0;j<2;j++){
        const cond=j===0?"17Gy_early":"17Gy_late",r=d.mice.filter(a=>a[0]===i&&a[1]===cond),v=r.map(a=>a[4]),by=yy+j*35;
        s+=txt(290,by+4,j===0?"Early":"Late",{anchor:"end",size:12});
        s+=line(pos(Math.min(...v),lo,hi,x,w),by,pos(Math.max(...v),lo,hi,x,w),by,colors[j],2);
        r.forEach((a,k)=>{s+=dot(pos(a[4],lo,hi,x,w),by+[-5,3,-2,6][k],5.5,colors[j]);});
      }
      if(i<8)s+=line(28,yy+65,1010,yy+65,"#D4DDE3");
    }
    s+=txt(655,431,"Mouse projection along state-specific fitted axis (PCA units)",{anchor:"middle",size:14})+txt(28,455,"Points are public GEO mouse samples; lines show observed range, not confidence intervals. Separation is descriptive, not predictive validation.",{size:12,color:muted});
    out["06_mouse_scores.svg"]=end(s);
  }
  return out;
}

const root = fileURLToPath(new URL("../notebooks/applications/results/gse211713/", import.meta.url));
const data = JSON.parse(readFileSync(join(root, "selected_results.json"), "utf8"));
const figures = renderGallery(data);
const check = process.argv.includes("--check");
const dir = join(root, "figures");
if (!check) mkdirSync(dir, {recursive:true});
for (const [name, svg] of Object.entries(figures)) {
  const path = join(dir, name);
  if (check) { if (readFileSync(path, "utf8") !== svg) throw Error("Figure differs: " + name); }
  else writeFileSync(path, svg);
}
console.log("PASS: " + Object.keys(figures).length + " public result snapshots");
