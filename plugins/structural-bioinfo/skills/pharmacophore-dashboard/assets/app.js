/* ===== pharmacophore-dashboard client app — target-agnostic, data-driven ===== */
const FAMILIES = ["Donor","Acceptor","PosIonizable","NegIonizable","Aromatic","LumpedHydrophobe"];
const FAMILY_COLORS = {Donor:"#2563eb",Acceptor:"#dc2626",PosIonizable:"#0d9488",
  NegIonizable:"#c026d3",Aromatic:"#d97706",LumpedHydrophobe:"#16a34a"};
const ITYPE_COLORS = {hbond:"#0891b2", salt_bridge:"#d97706"};

const $ = (s,r=document)=>r.querySelector(s);
const $$ = (s,r=document)=>Array.from(r.querySelectorAll(s));
function ce(t,c,txt){const e=document.createElement(t); if(c)e.className=c; if(txt!=null)e.textContent=txt; return e;}
function pct(a,b){return b>0?Math.round(100*a/b):0;}

const EXCLUDED = new Set(["excluded","control","negative","negative control","inactive","decoy"]);
const isProductive = e => e.n_after>0 && !EXCLUDED.has(e.category);
const RET = DATA.retention || DATA.retention_orthosteric || {};
const HAS_CONSENSUS = !!(DATA.consensus && DATA.consensus.features && DATA.consensus.features.length);
const HAS_MODES = !!(CONFIG.modes && CONFIG.modes.length);
const HAS_MAP = !!(CONFIG.interaction_map);

/* ---------------- Header + stat tiles ---------------- */
function renderHeader(){
  $("#target-name").textContent = DATA.target.name || "Pharmacophore analysis";
  const t=DATA.target, bits=[];
  if(t.gene)bits.push(`Gene <b>${t.gene}</b>`);
  if(t.uniprot)bits.push(`UniProt <a href="https://www.uniprot.org/uniprotkb/${t.uniprot}/entry" target="_blank" rel="noopener">${t.uniprot}</a>`);
  if(t.organism)bits.push(t.organism);
  $("#target-sub").innerHTML = bits.join(" · ");
  const productive = DATA.entries.filter(isProductive);
  const totalFiltered = DATA.entries.filter(e=>!EXCLUDED.has(e.category)).reduce((s,e)=>s+e.n_after,0);
  const tiles = [
    ["Structures analyzed", t.n_structures || new Set(DATA.entries.map(e=>e.pdb_id)).size],
    ["Ligand entries", DATA.entries.length],
    ["Productive poses", productive.length],
    ["Retained features", totalFiltered],
  ];
  if(HAS_CONSENSUS) tiles.push(["Consensus query", DATA.consensus.summary.total_features+" feats"]);
  else if(HAS_MODES) tiles.push(["Query modes", CONFIG.modes.length]);
  const wrap=$("#stat-tiles");
  tiles.forEach(([k,v])=>{const d=ce("div","tile");d.appendChild(ce("div","tile-v",String(v)));d.appendChild(ce("div","tile-k",k));wrap.appendChild(d);});
  $("#meta-line").textContent = CONFIG.meta || "";
}

function familyLegend(families){
  const w=ce("div","legend");
  families.forEach(f=>{const i=ce("span","leg-item");const sw=ce("span","sw");sw.style.background=FAMILY_COLORS[f];
    i.appendChild(sw);i.appendChild(ce("span",null,f));w.appendChild(i);});
  return w;
}

/* ---------------- Retention chart ---------------- */
function renderRetention(){
  const fams=FAMILIES.filter(f=>RET[f]&&RET[f].all>0);
  const host=$("#retention-chart");
  if(!fams.length){host.appendChild(ce("p","muted","No retention data available."));return;}
  const W=Math.max(560,fams.length*120),H=300,padL=44,padB=58,padT=18,padR=12;
  const maxV=Math.max(...fams.map(f=>RET[f].all));
  const ns="http://www.w3.org/2000/svg";
  const svg=document.createElementNS(ns,"svg");svg.setAttribute("viewBox",`0 0 ${W} ${H}`);svg.setAttribute("class","chart");
  const plotH=H-padB-padT,plotW=W-padL-padR,ticks=4;
  for(let i=0;i<=ticks;i++){const v=Math.round(maxV*i/ticks),y=padT+plotH-(plotH*i/ticks);
    const ln=document.createElementNS(ns,"line");ln.setAttribute("x1",padL);ln.setAttribute("x2",W-padR);ln.setAttribute("y1",y);ln.setAttribute("y2",y);ln.setAttribute("class","grid");svg.appendChild(ln);
    const tx=document.createElementNS(ns,"text");tx.setAttribute("x",padL-8);tx.setAttribute("y",y+4);tx.setAttribute("class","ax-y");tx.textContent=v;svg.appendChild(tx);}
  const band=plotW/fams.length;
  fams.forEach((f,idx)=>{const x0=padL+band*idx,bw=Math.min(34,band*0.32),gap=10;
    const allH=plotH*RET[f].all/maxV,filH=plotH*RET[f].filtered/maxV;
    const r1=document.createElementNS(ns,"rect");r1.setAttribute("x",x0+band/2-bw-gap/2);r1.setAttribute("y",padT+plotH-allH);r1.setAttribute("width",bw);r1.setAttribute("height",allH);r1.setAttribute("rx",3);r1.setAttribute("fill",FAMILY_COLORS[f]);r1.setAttribute("opacity","0.28");svg.appendChild(r1);
    const l1=document.createElementNS(ns,"text");l1.setAttribute("x",x0+band/2-bw/2-gap/2);l1.setAttribute("y",padT+plotH-allH-6);l1.setAttribute("class","ax-v");l1.textContent=RET[f].all;svg.appendChild(l1);
    const r2=document.createElementNS(ns,"rect");r2.setAttribute("x",x0+band/2+gap/2);r2.setAttribute("y",padT+plotH-filH);r2.setAttribute("width",bw);r2.setAttribute("height",filH);r2.setAttribute("rx",3);r2.setAttribute("fill",FAMILY_COLORS[f]);svg.appendChild(r2);
    const l2=document.createElementNS(ns,"text");l2.setAttribute("x",x0+band/2+bw/2+gap/2);l2.setAttribute("y",padT+plotH-filH-6);l2.setAttribute("class","ax-v");l2.textContent=RET[f].filtered;svg.appendChild(l2);
    const fl=document.createElementNS(ns,"text");fl.setAttribute("x",x0+band/2);fl.setAttribute("y",H-padB+20);fl.setAttribute("class","ax-x");fl.textContent=f;svg.appendChild(fl);
    const fp=document.createElementNS(ns,"text");fp.setAttribute("x",x0+band/2);fp.setAttribute("y",H-padB+38);fp.setAttribute("class","ax-pct");fp.textContent=pct(RET[f].filtered,RET[f].all)+"% kept";svg.appendChild(fp);});
  host.appendChild(svg);
  const lg=ce("div","legend");
  const a=ce("span","leg-item");const asw=ce("span","sw");asw.style.background="#9aa3af";asw.style.opacity="0.5";a.appendChild(asw);a.appendChild(ce("span",null,"All detected"));lg.appendChild(a);
  const b=ce("span","leg-item");const bsw=ce("span","sw");bsw.style.background="#475569";b.appendChild(bsw);b.appendChild(ce("span",null,"Retained (interaction-filtered)"));lg.appendChild(b);
  host.appendChild(lg);
  const vd=$("#retention-verdicts");
  fams.concat(FAMILIES.filter(f=>!RET[f]||RET[f].all===0)).forEach(f=>{
    const all=(RET[f]&&RET[f].all)||0,fil=(RET[f]&&RET[f].filtered)||0,keep=fil>0;
    const chip=ce("div","verdict "+(keep?"keep":"drop"));
    const dot=ce("span","vdot");dot.style.background=FAMILY_COLORS[f];chip.appendChild(dot);
    chip.appendChild(ce("b",null,f));
    chip.appendChild(ce("span","vmeta",all>0?`${fil}/${all} kept (${pct(fil,all)}%)`:"absent"));
    chip.appendChild(ce("span","vtag",keep?"INCLUDE":"EXCLUDE"));
    vd.appendChild(chip);});
}

/* ---------------- Structure overview table ---------------- */
let SELECTED=null;
function compositionBar(e){
  const ns="http://www.w3.org/2000/svg";
  const total=FAMILIES.reduce((s,f)=>s+(e.family_filtered[f]||0),0);
  const W=120,H=14;const svg=document.createElementNS(ns,"svg");svg.setAttribute("viewBox",`0 0 ${W} ${H}`);svg.setAttribute("width",W);svg.setAttribute("height",H);svg.setAttribute("class","compbar");
  if(total===0){const t=document.createElementNS(ns,"text");t.setAttribute("x",0);t.setAttribute("y",11);t.setAttribute("class","comp0");t.textContent="—";svg.appendChild(t);return svg;}
  let x=0;
  FAMILIES.forEach(f=>{const n=e.family_filtered[f]||0;if(!n)return;const w=W*n/total;
    const r=document.createElementNS(ns,"rect");r.setAttribute("x",x);r.setAttribute("y",0);r.setAttribute("width",Math.max(1,w-1));r.setAttribute("height",H);r.setAttribute("rx",2);r.setAttribute("fill",FAMILY_COLORS[f]);
    const ti=document.createElementNS(ns,"title");ti.textContent=`${f}: ${n}`;r.appendChild(ti);svg.appendChild(r);x+=w;});
  return svg;
}
function catBadge(cat){
  const map={productive:["Productive","b-ortho"],binding:["Binding","b-ortho"],orthosteric:["Orthosteric","b-ortho"],
    control:["Control","b-neg"],negative:["Negative control","b-neg"],excluded:["Excluded","b-excl"]};
  const [t,c]=map[cat]||[cat,"b-excl"];return `<span class="badge ${c}">${t}</span>`;
}
function renderStructures(){
  const tb=$("#struct-table tbody");
  DATA.entries.forEach(e=>{
    const tr=ce("tr");tr.dataset.key=e.key;if(!e.pocket_pdb)tr.classList.add("no3d");
    const ret=pct(e.n_after,e.n_before),chemo=e.chemotype_label||e.chemotype||"—";
    tr.innerHTML=
      `<td data-sort="${e.label}"><b>${e.label}</b></td>`+
      `<td data-sort="${e.chemotype||''}"><span class="chemo">${chemo}</span></td>`+
      `<td data-sort="${e.role||''}">${e.role||''}</td>`+
      `<td data-sort="${e.category}">${catBadge(e.category)}</td>`+
      `<td class="num" data-sort="${e.n_before}">${e.n_before}</td>`+
      `<td class="num" data-sort="${e.n_after}"><b>${e.n_after}</b></td>`+
      `<td class="num" data-sort="${ret}">${ret}%</td>`+
      `<td class="comp"></td><td class="act"></td>`;
    tr.querySelector(".comp").appendChild(compositionBar(e));
    if(e.pocket_pdb){const btn=ce("button","mini-btn","View 3D");btn.onclick=(ev)=>{ev.stopPropagation();selectEntry(e.key);};
      tr.querySelector(".act").appendChild(btn);tr.onclick=()=>selectEntry(e.key);tr.classList.add("clickable");}
    else{tr.querySelector(".act").innerHTML=`<span class="muted">no 3D</span>`;}
    tb.appendChild(tr);});
  makeSortable($("#struct-table"));
}

/* ---------------- Sortable tables ---------------- */
function makeSortable(table){
  const ths=$$("thead th",table);
  ths.forEach((th,ci)=>{
    if(th.classList.contains("nosort"))return;
    th.classList.add("sortable");
    th.onclick=()=>{
      const asc=!(th.dataset.asc==="true");ths.forEach(t=>{t.removeAttribute("data-asc");t.classList.remove("sorted");});
      th.dataset.asc=asc;th.classList.add("sorted");
      const rows=$$("tbody tr",table);
      rows.sort((a,b)=>{const av=a.children[ci].dataset.sort??a.children[ci].textContent;const bv=b.children[ci].dataset.sort??b.children[ci].textContent;
        const an=parseFloat(av),bn=parseFloat(bv);let r;if(!isNaN(an)&&!isNaN(bn))r=an-bn;else r=String(av).localeCompare(String(bv));return asc?r:-r;});
      const tb=$("tbody",table);rows.forEach(r=>tb.appendChild(r));
    };});
}

/* ---------------- 3D viewer ---------------- */
let viewerMain=null,viewerCons=null;
const vstate={fams:{},hbonds:true,labels:true,pocket:true,tol:false,spin:false};
function ensureViewer(id){if(typeof $3Dmol==="undefined")return null;return $3Dmol.createViewer($(id),{backgroundColor:"#0b1220"});}
function xyz(a){return {x:a[0],y:a[1],z:a[2]};}
function drawScene(viewer,pdb,features,opt){
  if(!viewer)return;
  viewer.clear();viewer.addModel(pdb,"pdb");viewer.setStyle({},{});
  if(opt.pocket){viewer.setStyle({hetflag:false},{line:{colorscheme:"whiteCarbon"}});
    const resis=[...new Set(features.filter(f=>f.partner_resnum).map(f=>parseInt(f.partner_resnum)))];
    if(resis.length)viewer.addStyle({hetflag:false,resi:resis},{stick:{radius:0.11,colorscheme:"whiteCarbon"}});}
  viewer.setStyle({hetflag:true},{stick:{radius:0.2,colorscheme:"cyanCarbon"}});
  features.forEach(f=>{if(opt.fams[f.family]===false)return;if(!f.position||f.position.length<3)return;
    viewer.addSphere({center:xyz(f.position),radius:opt.tol?(f.tolerance||1):0.45,color:FAMILY_COLORS[f.family],opacity:opt.tol?0.25:0.92});});
  if(opt.hbonds)features.forEach(f=>{if(opt.fams[f.family]===false)return;if(!f.partner_xyz)return;
    viewer.addLine({start:xyz(f.position),end:xyz(f.partner_xyz),dashed:true,color:ITYPE_COLORS[f.itype]||"#94a3b8"});});
  if(opt.labels){const seen=new Set();features.forEach(f=>{if(!f.partner_xyz||!f.partner_residue)return;if(opt.fams[f.family]===false)return;
    if(seen.has(f.partner_residue))return;seen.add(f.partner_residue);
    viewer.addLabel(f.partner_residue,{position:xyz(f.partner_xyz),fontSize:11,fontColor:"white",backgroundColor:"#1e293b",backgroundOpacity:0.85,borderThickness:0.5,borderColor:"#475569",inFront:true});});}
  viewer.zoomTo({hetflag:true});viewer.zoom(0.85);
  if(opt.spin)viewer.spin("y");else viewer.spin(false);
  viewer.render();
}
function selectEntry(key){
  const e=DATA.entries.find(x=>x.key===key);if(!e||!e.pocket_pdb)return;
  SELECTED=key;
  $$("#struct-table tbody tr").forEach(tr=>tr.classList.toggle("active",tr.dataset.key===key));
  const lig=e.ligand;
  $("#viewer-caption").innerHTML=
    `<b>${e.label}</b> — ${e.chemotype_label||e.chemotype||''} · ligand <code>${lig.resname||'?'} ${lig.chain||''}${lig.resnum||''}</code> · `+
    `${e.n_after} retained / ${e.n_before} detected · ${e.n_pocket_res} pocket residues`+
    (lig.smiles_used?`<div class="smiles">SMILES (${lig.smiles_source||'n/a'}): <code>${lig.smiles_used}</code></div>`:"");
  buildFamilyToggles(e.features);
  renderFeatureTable(e);
  if(!viewerMain)viewerMain=ensureViewer("#viewer-main");
  if(!viewerMain){$("#viewer-main").innerHTML=NO3D;return;}
  drawScene(viewerMain,e.pocket_pdb,e.features,vstate);viewerMain.resize();
}
function buildFamilyToggles(features){
  const present=[...new Set(features.map(f=>f.family))];
  const host=$("#fam-toggles");host.innerHTML="";
  present.forEach(f=>{if(vstate.fams[f]===undefined)vstate.fams[f]=true;
    const lab=ce("label","famtog");const cb=ce("input");cb.type="checkbox";cb.checked=vstate.fams[f]!==false;
    cb.onchange=()=>{vstate.fams[f]=cb.checked;redrawMain();};
    const sw=ce("span","sw");sw.style.background=FAMILY_COLORS[f];
    lab.appendChild(cb);lab.appendChild(sw);lab.appendChild(ce("span",null,f));host.appendChild(lab);});
}
function redrawMain(){if(!SELECTED)return;const e=DATA.entries.find(x=>x.key===SELECTED);drawScene(viewerMain,e.pocket_pdb,e.features,vstate);}
const NO3D=`<div class="no3d-msg"><b>3D viewer unavailable offline.</b><br>The 3Dmol.js library could not load. Reconnect to the internet and reload to view structures.</div>`;

/* ---------------- Feature detail table ---------------- */
function renderFeatureTable(e){
  const host=$("#feature-detail");
  host.innerHTML=`<h3>Retained features — ${e.label} <span class="muted">(${e.features.length})</span></h3>`;
  if(e.features.length===0){host.appendChild(ce("p","muted","No features retained after the interaction filter — this entry is a negative control."));return;}
  const tbl=ce("table","data-table");
  tbl.innerHTML=`<thead><tr><th>#</th><th>Family</th><th>Type</th><th>Partner residue</th><th>Atom</th><th>Dist (Å)</th><th>Interaction</th></tr></thead><tbody></tbody>`;
  const tb=$("tbody",tbl);
  e.features.forEach(f=>{const tr=ce("tr");
    tr.innerHTML=`<td class="num" data-sort="${f.id}">${f.id}</td>`+
      `<td data-sort="${f.family}"><span class="fdot" style="background:${FAMILY_COLORS[f.family]}"></span>${f.family}</td>`+
      `<td data-sort="${f.type}" class="mono">${f.type||''}</td>`+
      `<td data-sort="${f.partner_residue||''}"><b>${f.partner_residue||'—'}</b></td>`+
      `<td data-sort="${f.partner_atom||''}" class="mono">${f.partner_atom||'—'}</td>`+
      `<td class="num" data-sort="${f.distance==null?99:f.distance}">${f.distance!=null?f.distance.toFixed(2):'—'}</td>`+
      `<td data-sort="${f.itype||''}"><span class="itype ${f.itype||''}">${(f.itype||'').replace('_',' ')}</span></td>`;
    tb.appendChild(tr);});
  host.appendChild(tbl);makeSortable(tbl);
}

/* ---------------- Residue engagement ---------------- */
function renderResidues(){
  const all=DATA.residue_engagement||[];
  if(!all.length){$("#residues").style.display="none";return;}
  const data=all.slice(0,14),host=$("#residue-chart"),maxV=Math.max(...data.map(r=>r.total));
  const ns="http://www.w3.org/2000/svg",rowH=26,padL=104,padR=44,W=620,H=data.length*rowH+16;
  const svg=document.createElementNS(ns,"svg");svg.setAttribute("viewBox",`0 0 ${W} ${H}`);svg.setAttribute("class","chart");
  const plotW=W-padL-padR;
  data.forEach((r,i)=>{const y=i*rowH+8;
    const lbl=document.createElementNS(ns,"text");lbl.setAttribute("x",padL-8);lbl.setAttribute("y",y+rowH/2+1);lbl.setAttribute("class","ax-res");lbl.textContent=r.residue;svg.appendChild(lbl);
    let x=padL;
    FAMILIES.forEach(f=>{const n=r[f]||0;if(!n)return;const w=plotW*n/maxV;
      const rect=document.createElementNS(ns,"rect");rect.setAttribute("x",x);rect.setAttribute("y",y+3);rect.setAttribute("width",Math.max(1,w-1));rect.setAttribute("height",rowH-10);rect.setAttribute("rx",2);rect.setAttribute("fill",FAMILY_COLORS[f]);
      const ti=document.createElementNS(ns,"title");ti.textContent=`${f}: ${n}`;rect.appendChild(ti);svg.appendChild(rect);x+=w;});
    const tot=document.createElementNS(ns,"text");tot.setAttribute("x",x+6);tot.setAttribute("y",y+rowH/2+1);tot.setAttribute("class","ax-tot");tot.textContent=r.total;svg.appendChild(tot);});
  host.appendChild(svg);
  host.appendChild(familyLegend(FAMILIES.filter(f=>data.some(r=>r[f]>0))));
  const tb=$("#residue-table tbody");
  all.forEach(r=>{const tr=ce("tr");
    tr.innerHTML=`<td data-sort="${r.residue}"><b>${r.residue}</b></td>`+
      `<td class="num" data-sort="${r.total}"><b>${r.total}</b></td>`+
      `<td class="num" data-sort="${r.Donor||0}">${r.Donor||0}</td>`+
      `<td class="num" data-sort="${r.Acceptor||0}">${r.Acceptor||0}</td>`+
      `<td class="num" data-sort="${r.PosIonizable||0}">${r.PosIonizable||0}</td>`+
      `<td data-sort="${r.structures.join(',')}">${r.structures.map(s=>`<span class="pill">${s}</span>`).join(' ')}</td>`;
    tb.appendChild(tr);});
  makeSortable($("#residue-table"));
}

/* ---------------- Consensus ---------------- */
const consState={mand:true,opt:true,tol:false,spin:false};
function renderConsensus(){
  if(!HAS_CONSENSUS){return;}
  const c=DATA.consensus;$("#consensus-section").style.display="";
  const s=c.summary;
  const tiles=[["Total features",s.total_features],["Mandatory",s.mandatory],["Optional",s.optional]];
  if(c.source_entries)tiles.push(["Source PDBs",c.source_entries.length]);
  if(c.frame_pdb)tiles.push(["Reference frame",c.frame_pdb]);
  const th=$("#consensus-tiles");
  tiles.forEach(([k,v])=>{const d=ce("div","tile sm");d.appendChild(ce("div","tile-v",String(v)));d.appendChild(ce("div","tile-k",k));th.appendChild(d);});
  function ftable(host,feats){const tbl=ce("table","data-table");
    tbl.innerHTML=`<thead><tr><th>ID</th><th>Family</th><th>Position (x, y, z) Å</th><th>Tol (Å)</th><th>Sources</th></tr></thead><tbody></tbody>`;
    const tb=$("tbody",tbl);
    feats.forEach(f=>{const tr=ce("tr");
      tr.innerHTML=`<td class="num" data-sort="${f.id}">${f.id}</td>`+
        `<td data-sort="${f.family}"><span class="fdot" style="background:${FAMILY_COLORS[f.family]}"></span>${f.family}</td>`+
        `<td class="mono">[${f.position.map(x=>x.toFixed(2)).join(', ')}]</td>`+
        `<td class="num">${(f.tolerance||0).toFixed(1)}</td>`+
        `<td class="num" data-sort="${f.n_sources||''}">${f.n_sources!=null?f.n_sources:'—'}</td>`;
      tb.appendChild(tr);});
    host.appendChild(tbl);makeSortable(tbl);}
  ftable($("#cons-mand"),c.features.filter(f=>f.classification==="mandatory"));
  ftable($("#cons-opt"),c.features.filter(f=>f.classification==="optional"));
  const cl=$("#cons-legend");if(cl)cl.appendChild(familyLegend([...new Set(c.features.map(f=>f.family))]));
  if(c.pocket_pdb){if(typeof $3Dmol==="undefined")$("#viewer-cons").innerHTML=NO3D;
    else{viewerCons=$3Dmol.createViewer($("#viewer-cons"),{backgroundColor:"#0b1220"});drawConsensus();}}
  else{$("#viewer-cons").innerHTML=`<div class="no3d-msg muted">No reference pocket available for overlay.</div>`;}
  const st=$("#strategy-cards");
  (CONFIG.strategies||defaultStrategies(s)).forEach(x=>{const card=ce("div","strat");
    card.appendChild(ce("div","strat-h",x.name));card.appendChild(ce("div","strat-c",x.criteria));card.appendChild(ce("div","strat-b",x.behavior));st.appendChild(card);});
  // control labels
  const m=$("#cc-mand"),o=$("#cc-opt");if(m)m.parentNode.childNodes[1].textContent=` Mandatory (${s.mandatory})`;if(o)o.parentNode.childNodes[1].textContent=` Optional (${s.optional})`;
}
function defaultStrategies(s){return [
  {name:"Strict",criteria:`Match all ${s.mandatory} mandatory features`,behavior:"High precision, low recall."},
  {name:"Relaxed",criteria:`Match ≥${Math.max(1,Math.round(s.mandatory/2))} mandatory + ≥3 optional`,behavior:"Balanced precision / recall."},
  {name:"Scoring",criteria:"score = mandatory + 0.5 × optional; tune threshold",behavior:"Continuous ranking of hits."}];}
function drawConsensus(){
  const c=DATA.consensus,v=viewerCons;if(!v)return;
  v.clear();v.addModel(c.pocket_pdb,"pdb");v.setStyle({},{});
  v.setStyle({hetflag:false},{line:{colorscheme:"whiteCarbon"}});
  v.setStyle({hetflag:true},{stick:{radius:0.16,colorscheme:"cyanCarbon"}});
  c.features.forEach(f=>{const mand=f.classification==="mandatory";
    if(mand&&!consState.mand)return;if(!mand&&!consState.opt)return;
    v.addSphere({center:xyz(f.position),radius:consState.tol?f.tolerance:(mand?0.5:0.4),color:FAMILY_COLORS[f.family],opacity:consState.tol?0.22:(mand?0.95:0.55)});});
  v.zoomTo({hetflag:true});v.zoom(0.8);if(consState.spin)v.spin("y");else v.spin(false);v.render();
}

/* ---------------- Query modes (optional, config-supplied) ---------------- */
function renderModes(){
  if(!HAS_MODES){return;}
  $("#modes-section").style.display="";
  const host=$("#mode-cards");
  CONFIG.modes.forEach(m=>{const card=ce("div","mode-card");
    card.appendChild(ce("div","mode-h",m.name));
    if(m.confidence){const meta=ce("div","mode-meta");meta.innerHTML=`<span class="badge b-conf-${m.conf_cls||'mod'}">${m.confidence}</span> <span class="muted">${m.support||''}</span>`;card.appendChild(meta);}
    if(m.desc)card.appendChild(ce("p","mode-desc",m.desc));
    if(m.features){const tbl=ce("table","data-table sm");
      tbl.innerHTML=`<thead><tr><th>Feature</th><th>Tol</th><th>Partner</th></tr></thead><tbody>${
        m.features.map(f=>`<tr><td><span class="fdot" style="background:${FAMILY_COLORS[f.fam]}"></span>${f.fam}</td><td class="num">${f.tol}</td><td><b>${f.partner}</b></td></tr>`).join('')}</tbody>`;
      card.appendChild(tbl);}
    if(m.strategy)card.appendChild(ce("div","mode-strat",m.strategy));
    host.appendChild(card);});
}

/* ---------------- Narrative ---------------- */
function renderNarrative(){
  $("#exec-summary").innerHTML=CONFIG.summary_html||autoSummary();
  if(HAS_MAP){$("#imap").textContent=CONFIG.interaction_map;$("#imap-note").innerHTML=CONFIG.interaction_map_note||"";}
  else{$("#map").style.display="none";}
  const lim=$("#limitations-body");
  (CONFIG.limitations||autoLimitations()).forEach(l=>{const row=ce("div","lim-row");
    row.appendChild(ce("div","lim-issue",l.issue));row.appendChild(ce("div","lim-impact",l.impact));row.appendChild(ce("div","lim-mit",l.mitigation));lim.appendChild(row);});
  if(CONFIG.caveat)$("#caveat").innerHTML=CONFIG.caveat;else $("#caveat-block").style.display="none";
  $("#manifest").textContent=CONFIG.manifest||autoManifest();
  $("#footer-note").innerHTML=CONFIG.footer||defaultFooter();
}
function autoSummary(){
  const t=DATA.target,prod=DATA.entries.filter(isProductive),ctrl=DATA.entries.filter(e=>!isProductive(e));
  const topf=FAMILIES.filter(f=>RET[f]&&RET[f].filtered>0).sort((a,b)=>RET[b].filtered-RET[a].filtered);
  const topr=(DATA.residue_engagement||[]).slice(0,3).map(r=>`<b>${r.residue}</b>`);
  const drop=FAMILIES.filter(f=>RET[f]&&RET[f].all>0&&RET[f].filtered===0);
  let h=`<p>${t.n_structures||new Set(DATA.entries.map(e=>e.pdb_id)).size} structure(s) and ${DATA.entries.length} ligand pose(s) were analyzed for ${t.name||'this target'}. `;
  h+=`${prod.length} pose(s) retain at least one interaction-filtered pharmacophore feature`+(ctrl.length?`; ${ctrl.length} retain none and act as built-in negative controls.`:".")+`</p>`;
  if(topf.length)h+=`<p>The retained pharmacophore is dominated by <b>${topf.slice(0,2).join(" + ")}</b> features`+(topr.length?`, anchored mainly on ${topr.join(", ")}.`:".")+`</p>`;
  h+=`<div class="keyfind"><b>For screening:</b> build queries from `+(topf.length?topf.slice(0,3).join(", "):"the retained")+` features`+(drop.length?`; ${drop.join(", ")} are detected on the ligand(s) but never engage the receptor — exclude them.`:".")+(HAS_CONSENSUS?` A ${DATA.consensus.summary.total_features}-feature consensus query (${DATA.consensus.summary.mandatory} mandatory / ${DATA.consensus.summary.optional} optional) is provided below.`:"")+`</div>`;
  return h;
}
function autoLimitations(){
  const out=[];
  const bad=DATA.entries.filter(e=>e.ligand&&e.ligand.smiles_source&&!/ccd|rcsb|template_match|matched/i.test(e.ligand.smiles_source));
  if(bad.length)out.push({issue:"Best-effort ligand bond orders",impact:`SMILES template assignment did not fully succeed for ${bad.length} entry(ies) (e.g. ${bad.slice(0,3).map(e=>e.pdb_id).join(', ')}); some Donor/Acceptor assignments may be misclassified.`,mitigation:"Supply validated SMILES on re-run and cross-check features against visual inspection."});
  const ctrl=DATA.entries.filter(e=>e.n_after===0);
  if(ctrl.length)out.push({issue:"Non-productive poses present",impact:`${ctrl.length} pose(s) retain 0 features.`,mitigation:"Treat these as negative controls; exclude from query construction."});
  out.push({issue:"Simplified protein-side feature detection",impact:"Aromatic / hydrophobic contacts may be under-detected.",mitigation:"Add manually if specific π-stacking / hydrophobic residues are known to engage the scaffold."});
  if(DATA.provenance&&DATA.provenance.missing_pdb_for&&DATA.provenance.missing_pdb_for.length)
    out.push({issue:"Missing receptor coordinates",impact:`No clean PDB found for ${DATA.provenance.missing_pdb_for.join(', ')} — 3D disabled for those.`,mitigation:"Point --raw-dir at the directory holding the cleaned PDB receptors."});
  return out;
}
function autoManifest(){
  const lines=[(DATA.provenance&&DATA.provenance.pharmacophore_dir)||"pharmacophore/"];
  DATA.entries.forEach(e=>lines.push(`  ${e.key}   ${e.n_after}/${e.n_before} features   ${e.category}`));
  if(HAS_CONSENSUS)lines.push(`  consensus query   ${DATA.consensus.summary.total_features} features (frame ${DATA.consensus.frame_pdb})`);
  return lines.join("\n");
}
function defaultFooter(){return "Generated by the pharmacophore-dashboard skill. Feature coordinates, partner residues and distances are read directly from the per-structure filtered pharmacophore JSON files. 3D structures rendered with 3Dmol.js.";}

/* ---------------- Controls ---------------- */
function wireControls(){
  $("#c-hbonds").onchange=e=>{vstate.hbonds=e.target.checked;redrawMain();};
  $("#c-labels").onchange=e=>{vstate.labels=e.target.checked;redrawMain();};
  $("#c-pocket").onchange=e=>{vstate.pocket=e.target.checked;redrawMain();};
  $("#c-tol").onchange=e=>{vstate.tol=e.target.checked;redrawMain();};
  $("#c-spin").onchange=e=>{vstate.spin=e.target.checked;redrawMain();};
  $("#c-reset").onclick=()=>{if(viewerMain&&SELECTED){const e=DATA.entries.find(x=>x.key===SELECTED);drawScene(viewerMain,e.pocket_pdb,e.features,vstate);}};
  if(HAS_CONSENSUS){const m=$("#cc-mand"),o=$("#cc-opt"),t=$("#cc-tol"),s=$("#cc-spin");
    if(m)m.onchange=e=>{consState.mand=e.target.checked;drawConsensus();};
    if(o)o.onchange=e=>{consState.opt=e.target.checked;drawConsensus();};
    if(t)t.onchange=e=>{consState.tol=e.target.checked;drawConsensus();};
    if(s)s.onchange=e=>{consState.spin=e.target.checked;drawConsensus();};}
}

/* ---------------- Dynamic chrome: TOC + section numbering + scrollspy ---------------- */
function finalizeChrome(){
  const secs=$$("main section").filter(s=>s.style.display!=="none"&&s.dataset.toc);
  const toc=$(".toc");toc.innerHTML="";
  let n=0;
  secs.forEach(s=>{
    const numEl=s.querySelector(".sec .num");
    if(numEl){n++;numEl.textContent=n;}
    const a=ce("a",null,s.dataset.toc);a.href="#"+s.id;toc.appendChild(a);
  });
  if(typeof IntersectionObserver==="undefined")return;
  const links=$$(".toc a"),map={};links.forEach(a=>map[a.getAttribute("href").slice(1)]=a);
  const obs=new IntersectionObserver((ents)=>{ents.forEach(en=>{if(en.isIntersecting){links.forEach(l=>l.classList.remove("active"));if(map[en.target.id])map[en.target.id].classList.add("active");}});},{rootMargin:"-45% 0px -50% 0px"});
  secs.forEach(s=>obs.observe(s));
}

/* ---------------- init ---------------- */
function init(){
  renderHeader();renderNarrative();renderRetention();renderStructures();renderResidues();
  renderModes();renderConsensus();wireControls();finalizeChrome();
  const first=DATA.entries.find(e=>e.pocket_pdb&&isProductive(e))||DATA.entries.find(e=>e.pocket_pdb);
  if(first)setTimeout(()=>selectEntry(first.key),60);
}
window.__init3D=function(){
  if(typeof $3Dmol==="undefined")return;
  viewerMain=null;
  const sel=SELECTED||(DATA.entries.find(e=>e.pocket_pdb&&e.n_after>0)||DATA.entries.find(e=>e.pocket_pdb)||{}).key;
  if(sel)selectEntry(sel);
  if(HAS_CONSENSUS&&DATA.consensus.pocket_pdb){const el=$("#viewer-cons");if(el){el.innerHTML="";viewerCons=$3Dmol.createViewer(el,{backgroundColor:"#0b1220"});drawConsensus();}}
};
window.addEventListener("DOMContentLoaded",init);
