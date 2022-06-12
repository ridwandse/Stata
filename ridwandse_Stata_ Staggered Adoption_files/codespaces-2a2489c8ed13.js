"use strict";(()=>{var re=Object.defineProperty;var o=(L,E)=>re(L,"name",{value:E,configurable:!0});(globalThis.webpackChunk=globalThis.webpackChunk||[]).push([["codespaces"],{88361:(L,E,f)=>{var w=f(64463),g=f(59753),d=f(65935),p=f(57654);(0,g.on)("remote-input-error","#js-codespaces-repository-select",()=>{const e=document.querySelector("#js-codespaces-unable-load-repositories-warning");e.hidden=!1});function h(e){const n=new URL(document.location.href,window.location.origin),s=new URLSearchParams(n.search),l=["vscs_target"];for(const[u,C]of e.entries()){if(l.includes(u)||!C){s.delete(u);continue}s.set(u,C)}window.history.replaceState({},"",`?${s.toString()}`)}o(h,"updateNewCodespaceUrl"),(0,d.AC)(".js-new-codespace-form",async function(e,n){const s=e.closest("[data-replace-remote-form-target]"),l=s.querySelector(".js-new-codespace-submit-button");l instanceof HTMLInputElement&&(l.disabled=!0),e.classList.remove("is-error"),e.classList.add("is-loading");try{const u=await n.html();if(u.status!==200&&r(),s.replaceWith(u.html),s.getAttribute("data-allow-update-url")==="true"){const C=new FormData(document.querySelector("form.js-new-codespace-form"));h(C)}}catch{r()}});function r(){const e=new URL(document.location.href,window.location.origin),n=new URLSearchParams(e.search);n.set("response_error","true"),window.location.replace(`${window.location.pathname}?${n.toString()}`)}o(r,"resetNewCodespacePageWithResponseError");let c=null;function _(e){c=e,e!==null&&document.querySelector(".js-codespace-loading-steps").setAttribute("data-current-state",c)}o(_,"advanceLoadingState"),(0,w.N7)(".js-codespace-loading-steps",{constructor:HTMLElement,add:e=>{const n=e.getAttribute("data-current-state");n&&_(n)}}),(0,w.N7)(".js-codespace-advance-state",{constructor:HTMLElement,add:e=>{const n=e.getAttribute("data-state");n&&_(n)}});let t=null;function i(e){(0,d.AC)(".js-fetch-cascade-token",async function(n,s){try{t=(await s.json()).json.token}catch{}}),(0,p.Bt)(e)}o(i,"fetchCascadeToken");function m(e,n,s){if(document.querySelector(e)){const u=Date.now(),S=setInterval(o(()=>{const se=Date.now()-u;if(t||se>=s){clearInterval(S),n(t||void 0);return}},"checkToken"),50)}}o(m,"waitForCascadeTokenWithTimeout"),(0,w.N7)(".js-auto-submit-form",{constructor:HTMLFormElement,initialize:p.Bt}),(0,w.N7)(".js-workbench-form-container",{constructor:HTMLElement,add:e=>{const n=e.querySelector(".js-cascade-token");b(e,n)}});function b(e,n){if(n.value!==""){const s=e.querySelector("form");(0,p.Bt)(s)}else{const s=document.querySelector(".js-fetch-cascade-token");i(s),m(".js-workbench-form-container",a,1e4)}}o(b,"resolveCascadeToken");function a(e){const n=document.querySelector(".js-workbench-form-container form");n&&e?(v(n,e),T(n,e),(0,p.Bt)(n)):_("failed")}o(a,"insertCodespaceTokenIntoShowAuthForm");function v(e,n){const s=e.querySelector(".js-cascade-token");s&&s.setAttribute("value",n)}o(v,"insertCodespaceTokenIntoCascadeField");function T(e,n){const s=e.querySelector(".js-partner-info");if(s){let l=s.getAttribute("value");l&&(l=l.replace("%CASCADE_TOKEN_PLACEHOLDER%",n),s.setAttribute("value",l))}}o(T,"insertCodespaceTokenIntoPartnerInfo");var y=f(90420),x=f(69567),P=f(7732),M=f(5638);(0,g.on)("submit",".js-toggle-hidden-codespace-form",function(e){const n=e.currentTarget;A(n)});function A(e){const n=e.querySelectorAll(".js-toggle-hidden");for(const l of n)l.hidden=!l.hidden;const s=e.querySelectorAll(".js-toggle-disabled");for(const l of s)l.getAttribute("aria-disabled")?l.removeAttribute("aria-disabled"):l.setAttribute("aria-disabled","true")}o(A,"toggleFormSubmissionInFlight");function K(e){return e.closest("[data-replace-remote-form-target]")}o(K,"getFormTarget");async function k(){const e=document.querySelector(".js-codespaces-details-container");e&&(e.open=!1);const n=document.querySelector("new-codespace");if(n&&!n.getAttribute("data-no-submit-on-create"))try{const s=await fetch("/codespaces/new");if(s&&s.ok){const l=(0,M.r)(document,await s.text());n.replaceWith(l)}}catch{}}o(k,"createFormSubmitted"),(0,g.on)("submit",".js-create-codespaces-form-command",function(e){const n=e.currentTarget;n.classList.contains("js-open-in-vscode-form")||(k(),A(n))}),(0,g.on)("submit",".js-open-in-web-form",async function(e){const n=e.currentTarget;e.preventDefault(),n.submit()}),(0,g.on)("submit","form.js-codespaces-delete-form",async function(e){e.preventDefault();const n=e.currentTarget;try{const s=await fetch(n.action,{method:n.method,body:new FormData(n),headers:{Accept:"text/fragment+html","X-Requested-With":"XMLHttpRequest"}});if(s.ok){const l=(0,M.r)(document,await s.text());K(n).replaceWith(l)}else if(s.status===401)n.submit();else throw new Error(`Unexpected response: ${s.statusText}`)}finally{A(n)}}),(0,g.on)("submit","form.js-open-in-vscode-form",async function(e){e.preventDefault();const n=e.currentTarget;await U(n)});async function B(e,n){const s=document.querySelector(`#${e}`),l=await(0,P.W)({content:s.content.cloneNode(!0),dialogClass:"project-dialog"});return n&&n.setAttribute("aria-expanded","true"),l.addEventListener("dialog:remove",function(){n&&A(n)},{once:!0}),l}o(B,"openDialog");async function U(e){const n=await fetch(e.action,{method:e.method,body:new FormData(e),headers:{Accept:"application/json","X-Requested-With":"XMLHttpRequest"}});if(n.ok){const s=await n.json();s.codespace_url?(window.location.href=s.codespace_url,A(e),k(),V()):(e.closest("get-repo")||e.closest("new-codespace")?(e.setAttribute("data-src",s.loading_url),e.dispatchEvent(new CustomEvent("pollvscode"))):e.closest("create-button")&&(e.setAttribute("data-src",s.loading_url),e.dispatchEvent(new CustomEvent("prpollvscode"))),A(e))}else if(n.status===422){const s=await n.json();if(s.error_type==="concurrency_limit_error")await B("concurrency-error",e);else{const l=document.querySelector("template.js-flash-template"),u=s.error;l.after(new x.R(l,{className:"flash-error",message:u})),A(e)}}}o(U,"createCodespaceIntoVscode");async function V(){const e=document.querySelector(".js-codespaces-completable"),n=e&&e.getAttribute("data-src");if(!n)return;const s=await fetch(n,{method:"GET",headers:{Accept:"text/fragment+html","X-Requested-With":"XMLHttpRequest"}});if(s.ok){const l=(0,M.r)(document,await s.text());e.replaceWith(l)}else throw new Error(`Unexpected response: ${s.statusText}`)}o(V,"renderAllDone");var z=Object.defineProperty,N=Object.getOwnPropertyDescriptor,$=o((e,n,s,l)=>{for(var u=l>1?void 0:l?N(n,s):n,C=e.length-1,S;C>=0;C--)(S=e[C])&&(u=(l?S(n,s,u):S(u))||u);return l&&u&&z(n,s,u),u},"__decorateClass");let W=o(class extends HTMLElement{async connectedCallback(){B("concurrency-error")}},"ConcurrencyLimitElement");W=$([y.Ih],W);var X=Object.defineProperty,G=Object.getOwnPropertyDescriptor,R=o((e,n,s,l)=>{for(var u=l>1?void 0:l?G(n,s):n,C=e.length-1,S;C>=0;C--)(S=e[C])&&(u=(l?S(n,s,u):S(u))||u);return l&&u&&X(n,s,u),u},"new_codespace_element_decorateClass");let D=o(class extends HTMLElement{async connectedCallback(){const e=new URL(document.location.href,window.location.origin),n=new URLSearchParams(e.search);n.has("response_error")&&(n.delete("response_error"),window.history.replaceState({},"",`?${n.toString()}`))}toggleLoadingVscode(){const e=this.loadingVscode.hidden,n=this.children;for(let s=0;s<n.length;s++)n[s].hidden=e;this.loadingVscode.hidden=!e}machineTypeSelected(e){const s=e.currentTarget.getAttribute("data-sku-name");this.skuNameInput&&s&&(this.skuNameInput.value=s),this.advancedOptionsForm&&(this.advancedOptionsForm.requestSubmit?this.advancedOptionsForm.requestSubmit():this.advancedOptionsForm.submit())}pollForVscode(e){this.toggleLoadingVscode();const n=e.currentTarget.getAttribute("data-src");n&&this.vscodePoller.setAttribute("src",n)}vscsTargetUrlUpdated(e){const n=e.currentTarget;this.vscsTargetUrl.value=n.value}},"NewCodespaceElement");R([y.fA],D.prototype,"vscsTargetUrl",2),R([y.fA],D.prototype,"loadingVscode",2),R([y.fA],D.prototype,"vscodePoller",2),R([y.fA],D.prototype,"advancedOptionsForm",2),R([y.fA],D.prototype,"skuNameInput",2),D=R([y.Ih],D);var Z=Object.defineProperty,Q=Object.getOwnPropertyDescriptor,F=o((e,n,s,l)=>{for(var u=l>1?void 0:l?Q(n,s):n,C=e.length-1,S;C>=0;C--)(S=e[C])&&(u=(l?S(n,s,u):S(u))||u);return l&&u&&Z(n,s,u),u},"export_branch_element_decorateClass");let j=o(class extends HTMLElement{constructor(){super(...arguments);this.abortPoll=null}connectedCallback(){this.abortPoll=new AbortController,this.loadingIndicator.hidden||this.startPoll()}disconnectedCallback(){var e;(e=this.abortPoll)==null||e.abort()}async exportBranch(e){e.preventDefault(),this.form.hidden=!0,this.loadingIndicator.hidden=!1,(await fetch(this.form.action,{method:this.form.method,body:new FormData(this.form),headers:{Accept:"text/fragment+html","X-Requested-With":"XMLHttpRequest"}})).ok?this.startPoll():(this.form.hidden=!1,this.loadingIndicator.hidden=!0)}async startPoll(){const e=this.getAttribute("data-exported-codespace-url")||"",n=await this.poll(e);if(n)if(n.ok)this.loadingIndicator.hidden=!0,this.viewBranchLink.hidden=!1;else{const s=this.getAttribute("data-export-error-redirect-url")||"";window.location.href=s}}async poll(e,n=1e3){var s,l;if((s=this.abortPoll)==null?void 0:s.signal.aborted)return;const u=await fetch(e,{signal:(l=this.abortPoll)==null?void 0:l.signal});return u.status===202?(await new Promise(C=>setTimeout(C,n)),this.poll(e,n*1.5)):u}},"ExportBranchElement");F([y.fA],j.prototype,"form",2),F([y.fA],j.prototype,"loadingIndicator",2),F([y.fA],j.prototype,"viewBranchLink",2),j=F([y.Ih],j);var J=f(66963),Y=Object.defineProperty,ee=Object.getOwnPropertyDescriptor,O=o((e,n,s,l)=>{for(var u=l>1?void 0:l?ee(n,s):n,C=e.length-1,S;C>=0;C--)(S=e[C])&&(u=(l?S(n,s,u):S(u))||u);return l&&u&&Y(n,s,u),u},"options_popover_element_decorateClass");let I=o(class extends HTMLElement{reset(e){for(e.preventDefault(),this.dropdownDetails.hidden=!1,this.modalDetails.hidden=!0,this.exportDetails.hidden=!0,this.skuForm.hidden=!1;this.resultMessage.firstChild;)this.resultMessage.removeChild(this.resultMessage.firstChild);this.resultMessage.hidden=!0,this.errorMessage.hidden=!0}showSettingsModal(){this.dropdownDetails.hidden=!0,this.modalDetails.open=!0,this.modalDetails.hidden=!1}showExportModal(){this.dropdownDetails.hidden=!0,this.exportDetails.open=!0,this.exportDetails.hidden=!1,this.insertBranchExportComponent()}async insertBranchExportComponent(){const e=this.querySelector("[data-branch-export-url]");if(!e)return;const n=e.getAttribute("data-branch-export-url");if(!n)return;const s=await(0,J.a)(document,n);!s||(e.innerHTML="",e.appendChild(s))}},"OptionsPopoverElement");O([y.fA],I.prototype,"dropdownDetails",2),O([y.fA],I.prototype,"modalDetails",2),O([y.fA],I.prototype,"settingsModal",2),O([y.fA],I.prototype,"skuForm",2),O([y.fA],I.prototype,"resultMessage",2),O([y.fA],I.prototype,"errorMessage",2),O([y.fA],I.prototype,"exportDetails",2),I=O([y.Ih],I);var te=Object.defineProperty,ne=Object.getOwnPropertyDescriptor,q=o((e,n,s,l)=>{for(var u=l>1?void 0:l?ne(n,s):n,C=e.length-1,S;C>=0;C--)(S=e[C])&&(u=(l?S(n,s,u):S(u))||u);return l&&u&&te(n,s,u),u},"vscode_forwarder_element_decorateClass");let H=o(class extends HTMLElement{async connectedCallback(){this.closeDetailsPopover();const e=this.getAttribute("data-codespace-url");e&&(window.location.href=e)}closeDetailsPopover(){const e=document.querySelector(".js-codespaces-details-container");e&&(e.open=!1)}},"VscodeForwarderElement");q([y.fA],H.prototype,"vscodeLink",2),H=q([y.Ih],H);var ie=f(10965),oe=f(37505)},63621:(L,E,f)=>{f.d(E,{H:()=>d,v:()=>g});var w=f(59753);function g(){const p=document.getElementById("ajax-error-message");p&&(p.hidden=!1)}o(g,"showGlobalError");function d(){const p=document.getElementById("ajax-error-message");p&&(p.hidden=!0)}o(d,"hideGlobalError"),(0,w.on)("deprecatedAjaxError","[data-remote]",function(p){const h=p.detail,{error:r,text:c}=h;p.currentTarget===p.target&&(r==="abort"||r==="canceled"||(/<html/.test(c)?(g(),p.stopImmediatePropagation()):setTimeout(function(){p.defaultPrevented||g()},0)))}),(0,w.on)("deprecatedAjaxSend","[data-remote]",function(){d()}),(0,w.on)("click",".js-ajax-error-dismiss",function(){d()})},7732:(L,E,f)=>{f.d(E,{W:()=>g});var w=f(59753);async function g(d){const h=document.querySelector("#site-details-dialog").content.cloneNode(!0),r=h.querySelector("details"),c=r.querySelector("details-dialog"),_=r.querySelector(".js-details-dialog-spinner");d.detailsClass&&r.classList.add(...d.detailsClass.split(" ")),d.dialogClass&&c.classList.add(...d.dialogClass.split(" ")),d.label?c.setAttribute("aria-label",d.label):d.labelledBy&&c.setAttribute("aria-labelledby",d.labelledBy),document.body.append(h);const t=await d.content;return _.remove(),c.prepend(t),r.addEventListener("toggle",()=>{r.hasAttribute("open")||((0,w.f)(c,"dialog:remove"),r.remove())}),c}o(g,"dialog")},75488:(L,E,f)=>{f.d(E,{C:()=>g,x:()=>w});const w=function(){return document.readyState==="interactive"||document.readyState==="complete"?Promise.resolve():new Promise(d=>{document.addEventListener("DOMContentLoaded",()=>{d()})})}(),g=function(){return document.readyState==="complete"?Promise.resolve():new Promise(d=>{window.addEventListener("load",d)})}()},66963:(L,E,f)=>{f.d(E,{D:()=>p,a:()=>d});var w=f(99997),g=f(5638);async function d(h,r,c){const _=new Request(r,c);_.headers.append("X-Requested-With","XMLHttpRequest");const t=await self.fetch(_);if(t.status<200||t.status>=300)throw new Error(`HTTP ${t.status}${t.statusText||""}`);return(0,w.t)((0,w.P)(h),t),(0,g.r)(h,await t.text())}o(d,"fetchSafeDocumentFragment");function p(h,r,c=1e3){return o(async function _(t){const i=new Request(h,r);i.headers.append("X-Requested-With","XMLHttpRequest");const m=await self.fetch(i);if(m.status<200||m.status>=300)throw new Error(`HTTP ${m.status}${m.statusText||""}`);if(m.status===200)return m;if(m.status===202)return await new Promise(b=>setTimeout(b,t)),_(t*1.5);throw new Error(`Unexpected ${m.status} response status from poll endpoint`)},"poll")(c)}o(p,"fetchPoll")},57654:(L,E,f)=>{f.d(E,{Bt:()=>h,DN:()=>_,KL:()=>m,Se:()=>c,qC:()=>b,sw:()=>t});var w=f(59753),g=f(2077),d=f(63621);(0,w.on)("click",".js-remote-submit-button",async function(a){const T=a.currentTarget.form;a.preventDefault();let y;try{y=await fetch(T.action,{method:T.method,body:new FormData(T),headers:{Accept:"application/json","X-Requested-With":"XMLHttpRequest"}})}catch{}y&&!y.ok&&(0,d.v)()});function p(a,v,T){return a.dispatchEvent(new CustomEvent(v,{bubbles:!0,cancelable:T}))}o(p,"fire");function h(a,v){v&&(r(a,v),(0,g.j)(v)),p(a,"submit",!0)&&a.submit()}o(h,"requestSubmit");function r(a,v){if(!(a instanceof HTMLFormElement))throw new TypeError("The specified element is not of type HTMLFormElement.");if(!(v instanceof HTMLElement))throw new TypeError("The specified element is not of type HTMLElement.");if(v.type!=="submit")throw new TypeError("The specified element is not a submit button.");if(!a||a!==v.form)throw new Error("The specified element is not owned by the form element.")}o(r,"checkButtonValidity");function c(a,v){if(typeof v=="boolean")if(a instanceof HTMLInputElement)a.checked=v;else throw new TypeError("only checkboxes can be set to boolean value");else{if(a.type==="checkbox")throw new TypeError("checkbox can't be set to string value");a.value=v}p(a,"change",!1)}o(c,"changeValue");function _(a,v){for(const T in v){const y=v[T],x=a.elements.namedItem(T);(x instanceof HTMLInputElement||x instanceof HTMLTextAreaElement)&&(x.value=y)}}o(_,"fillFormValues");function t(a){if(!(a instanceof HTMLElement))return!1;const v=a.nodeName.toLowerCase(),T=(a.getAttribute("type")||"").toLowerCase();return v==="select"||v==="textarea"||v==="input"&&T!=="submit"&&T!=="reset"||a.isContentEditable}o(t,"isFormField");function i(a){return new URLSearchParams(a)}o(i,"searchParamsFromFormData");function m(a,v){const T=new URLSearchParams(a.search),y=i(v);for(const[x,P]of y)T.append(x,P);return T.toString()}o(m,"combineGetFormSearchParams");function b(a){return i(new FormData(a)).toString()}o(b,"serialize")},99997:(L,E,f)=>{f.d(E,{P:()=>w,t:()=>d});function w(p){const h=[...p.querySelectorAll("meta[name=html-safe-nonce]")].map(r=>r.content);if(h.length<1)throw new Error("could not find html-safe-nonce on document");return h}o(w,"getDocumentHtmlSafeNonces");class g extends Error{constructor(h,r){super(`${h} for HTTP ${r.status}`);this.response=r}}o(g,"ResponseError");function d(p,h,r=!1){const c=h.headers.get("content-type")||"";if(!r&&!c.startsWith("text/html"))throw new g(`expected response with text/html, but was ${c}`,h);if(r&&!(c.startsWith("text/html")||c.startsWith("application/json")))throw new g(`expected response with text/html or application/json, but was ${c}`,h);const _=h.headers.get("x-html-safe");if(_){if(!p.includes(_))throw new g("response X-HTML-Safe nonce did not match",h)}else throw new g("missing X-HTML-Safe nonce",h)}o(d,"verifyResponseHtmlSafeNonce")},10965:(L,E,f)=>{var w=f(90420),g=Object.defineProperty,d=Object.getOwnPropertyDescriptor,p=o((r,c,_,t)=>{for(var i=t>1?void 0:t?d(c,_):c,m=r.length-1,b;m>=0;m--)(b=r[m])&&(i=(t?b(c,_,i):b(i))||i);return t&&i&&g(c,_,i),i},"__decorateClass");let h=o(class extends HTMLElement{connectedCallback(){this.control&&(this.storedInput=Array(this.control.children.length).fill("")),this.addEventListener("input",this.relayInput.bind(this)),this.addEventListener("keydown",this.relayKeydown.bind(this));const r=this.closest("details");r&&r.addEventListener("toggle",()=>{r.open&&this.source.focus()})}relayKeydown(r){if((this.isControlTab(r.target)||r.target===this.source)&&(r.key==="ArrowDown"||r.key==="Tab"))r.preventDefault(),r.stopPropagation(),this.routeCustomEvent(new CustomEvent("focus-list"));else if(r.key==="Escape"){const c=this.closest("details");c&&c.removeAttribute("open")}}isControlTab(r){return!r||!this.control?!1:Array.from(this.control.children).includes(r)}relayInput(r){if(!r.target)return;const _=r.target.value;this.routeCustomEvent(new CustomEvent("input-entered",{detail:_}))}routeCustomEvent(r){this.sinks[this.selectedIndex].dispatchEvent(r)}get selectedIndex(){if(!this.control)return 0;const r=this.control.querySelector('[aria-selected="true"]');return r?Array.from(this.control.children).indexOf(r):0}storeInput(){this.storedInput[this.selectedIndex]=this.source.value}updateInput(r){this.source.value=this.storedInput[this.selectedIndex];const c=r.detail.relatedTarget.getAttribute("data-filter-placeholder");this.source.placeholder=c,this.source.setAttribute("aria-label",c),this.notifySelected()}notifySelected(){const r=this.sinks[this.selectedIndex],c=new CustomEvent("tab-selected");r.dispatchEvent(c)}},"InputDemuxElement");p([w.fA],h.prototype,"source",2),p([w.GO],h.prototype,"sinks",2),p([w.fA],h.prototype,"control",2),h=p([w.Ih],h)},5638:(L,E,f)=>{f.d(E,{r:()=>w});function w(g,d){const p=g.createElement("template");return p.innerHTML=d,g.importNode(p.content,!0)}o(w,"parseHTML")},70290:(L,E,f)=>{f.d(E,{Z:()=>w});function w(g){var d,p;const h=(p=(d=g.head)==null?void 0:d.querySelector('meta[name="expected-hostname"]'))==null?void 0:p.content;if(!h)return!1;const r=h.replace(/\.$/,"").split(".").slice(-2).join("."),c=g.location.hostname.replace(/\.$/,"").split(".").slice(-2).join(".");return r!==c}o(w,"detectProxySite")},37505:(L,E,f)=>{var w=f(41891),g=f(69567),d=f(90420),p=f(17945),h=Object.defineProperty,r=Object.getOwnPropertyDescriptor,c=o((t,i,m,b)=>{for(var a=b>1?void 0:b?r(i,m):i,v=t.length-1,T;v>=0;v--)(T=t[v])&&(a=(b?T(i,m,a):T(a))||a);return b&&a&&h(i,m,a),a},"__decorateClass");let _=o(class extends HTMLElement{constructor(){super(...arguments);this.isCurrentVisible=!1,this.currentSelectionIndex=null,this.handleWindowResize=()=>{if(!this.virtualizedList)return;const t=this.isMobileViewport,i=this.windowHeight;this.updateViewportSize();const m=t!==this.isMobileViewport,b=i!==this.windowHeight;if(m){this.virtualizedList.destroy(),this.setupVirtualizedList();return}!this.isMobileViewport||!b||(this.listContainer.style.maxHeight=`${this.listHeight}px`,this.virtualizedList.resize(this.listHeight))},this.windowResized=()=>{this.resizeAnimationRequest&&cancelAnimationFrame(this.resizeAnimationRequest),this.resizeAnimationRequest=requestAnimationFrame(this.handleWindowResize)}}connectedCallback(){window.addEventListener("resize",this.windowResized),this.refType=this.getRequiredAttr("type")==="branch"?w.r.Branch:w.r.Tag;const t=this.getAttribute("current-committish");this.currentCommittish=t?atob(t):null,this.input=this.hasAttribute("initial-filter")&&this.currentCommittish||"",this.defaultBranch=atob(this.getRequiredAttr("default-branch")),this.nameWithOwner=atob(this.getRequiredAttr("name-with-owner")),this.canCreate=this.hasAttribute("can-create"),this.prefetchOnMouseover=this.hasAttribute("prefetch-on-mouseover");const i=this.getRequiredAttr("query-endpoint"),m=this.getRequiredAttr("cache-key");this.index=new w.W(this.refType,this,i,m,this.nameWithOwner),this.updateViewportSize(),this.setupFetchListeners()}disconnectedCallback(){this.resizeAnimationRequest&&cancelAnimationFrame(this.resizeAnimationRequest),window.removeEventListener("resize",this.windowResized)}updateViewportSize(){this.isMobileViewport=window.innerWidth<544,this.windowHeight=window.innerHeight}inputEntered(t){this.input=t.detail,this.render()}tabSelected(){this.index.fetchData()}renderTemplate(t,i){return new g.R(t,i,g.XK)}renderRow(t){const i=this.index.currentSearchResult[t];if(!i&&t>=this.listLength)return document.createElement("span");if(this.index.fetchFailed)return this.renderTemplate(this.fetchFailedTemplate,{index:t,refName:this.input});if(!i){const T=this.input===this.currentCommittish;return this.isCurrentVisible||(this.isCurrentVisible=T),this.renderTemplate(this.noMatchTemplate,{index:t,isCurrent:T,refName:this.input})}const m=this.input.length>0,b=m?"is-filtering":"",a=i===this.currentCommittish;this.isCurrentVisible||(this.isCurrentVisible=a);const v=this.renderTemplate(this.itemTemplate,{refName:i,index:t,isFilteringClass:b,urlEncodedRefName:this.urlEncodeRef(i),isCurrent:a,isNotDefault:i!==this.defaultBranch});if(m){const T=v.querySelector("span");T.textContent="";const y=i.split(this.input),x=y.length-1;for(let P=0;P<y.length;P++){const M=y[P];if(T.appendChild(document.createTextNode(M)),P<x){const A=document.createElement("b");A.textContent=this.input,T.appendChild(A)}}}return v}urlEncodeRef(t){return encodeURIComponent(t).replaceAll("%2F","/").replaceAll("%3A",":").replaceAll("%2B","+")}render(){if(this.currentSelectionIndex=null,!this.index.isLoading){if(!this.virtualizedList){this.index.search(this.input),this.setupVirtualizedList();return}this.listContainer.scrollTop=0,this.index.search(this.input),this.virtualizedList.setRowCount(this.listLength)}}get listHeight(){return this.isMobileViewport?this.windowHeight*.3:330}get listLength(){const t=this.index.currentSearchResult.length;return this.showCreateRow?t+1:t||1}get showCreateRow(){return!this.index.fetchFailed&&!this.index.exactMatchFound&&this.input!==""&&this.canCreate}getRequiredAttr(t,i=this){const m=i.getAttribute(t);if(!m)throw new Error(`Missing attribute for ${i}: ${t}`);return m}setupFetchListeners(){const t=this.closest("details");let i=!1;const m=o(()=>{i||(this.index.fetchData(),i=!0)},"fetch");if(!t||t.open){m();return}this.prefetchOnMouseover&&t.addEventListener("mouseover",m,{once:!0}),this.addEventListener("keydown",this.keydown),this.addEventListener("change",this.updateCurrent);const b=t.querySelector("input[data-ref-filter]");b&&(b.addEventListener("input",()=>{this.input=b.value,this.render()}),b.addEventListener("keydown",a=>{if(a.key==="ArrowDown"||a.key==="Tab")a.preventDefault(),a.stopPropagation(),this.focusFirstListMember();else if(a.key==="Enter"){let v=this.index.currentSearchResult.indexOf(this.input);if(v===-1)if(this.showCreateRow)v=this.listLength-1;else return;t.querySelector(`[data-index="${v}"]`).click(),a.preventDefault()}}))}focusFirstListMember(){!this.virtualizedList||(this.currentSelectionIndex=0,this.focusItemAtIndex(this.currentSelectionIndex))}updateCurrent(t){t.target instanceof HTMLInputElement&&t.target.checked&&t.target.value&&(this.currentCommittish=t.target.value)}keydown(t){if(this.currentSelectionIndex!==null){if(t.key==="Enter"){const i=document.activeElement;if(!i)return;i.click(),t.preventDefault();return}if(!(t.key==="Tab"&&t.shiftKey)&&t.key!=="Escape")switch(t.preventDefault(),t.stopPropagation(),t.key){case"ArrowUp":{this.currentSelectionIndex--,this.currentSelectionIndex<0&&(this.currentSelectionIndex=this.listLength-1),this.focusItemAtIndex(this.currentSelectionIndex);break}case"Home":{this.currentSelectionIndex=0,this.focusItemAtIndex(this.currentSelectionIndex);break}case"End":{this.currentSelectionIndex=this.listLength-1,this.focusItemAtIndex(this.currentSelectionIndex);break}case"Tab":case"ArrowDown":{this.currentSelectionIndex++,this.currentSelectionIndex>this.listLength-1&&(this.currentSelectionIndex=0),this.focusItemAtIndex(this.currentSelectionIndex);break}}}}focusItemAtIndex(t){this.virtualizedList.scrollToIndex(t,"center"),setTimeout(()=>{const i=this.listContainer.querySelector(`[data-index="${t}"]`);i&&i.focus()},20)}setupVirtualizedList(){this.listContainer.innerHTML="",this.listContainer.style.maxHeight=`${this.listHeight}px`,this.virtualizedList=new p.Z(this.listContainer,{height:this.listHeight,rowCount:this.listLength,renderRow:this.renderRow.bind(this),rowHeight:t=>{const i=this.isMobileViewport?54:33;return this.showCreateRow&&t===this.listLength-1?51:i},onRowsRendered:()=>{var t;this.hiddenCurrentElement&&(this.listContainer.removeChild(this.hiddenCurrentElement),delete this.hiddenCurrentElement),this.isCurrentVisible?this.isCurrentVisible=!1:this.hiddenCurrentItemTemplate&&(this.hiddenCurrentElement=document.createElement("div"),(t=this.hiddenCurrentElement)==null||t.appendChild(this.renderTemplate(this.hiddenCurrentItemTemplate,{refName:this.currentCommittish})),this.listContainer.appendChild(this.hiddenCurrentElement))},initialIndex:0,overscanCount:6}),this.virtualizedList.resize.bind(this.virtualizedList)}},"RefSelectorElement");c([d.fA],_.prototype,"listContainer",2),c([d.fA],_.prototype,"itemTemplate",2),c([d.fA],_.prototype,"noMatchTemplate",2),c([d.fA],_.prototype,"fetchFailedTemplate",2),c([d.fA],_.prototype,"hiddenCurrentItemTemplate",2),_=c([d.Ih],_)},41891:(L,E,f)=>{f.d(E,{W:()=>_,r:()=>r});var w=f(31579),g=f(77552);const{getItem:d,setItem:p,removeItem:h}=(0,w.Z)("localStorage",{throwQuotaErrorsOnSet:!0});var r=(t=>(t.Branch="branch",t.Tag="tag",t))(r||{});const c=o(class{constructor(t,i,m,b,a){this.knownItems=[],this.currentSearchResult=[],this.exactMatchFound=!1,this.searchTerm="",this.isLoading=!0,this.fetchInProgress=!1,this.fetchFailed=!1,this.refType=t,this.selector=i,this.refEndpoint=m,this.cacheKey=b,this.nameWithOwner=a}render(){this.selector.render()}async fetchData(){try{if(!this.isLoading||this.fetchInProgress)return;if(!this.bootstrapFromLocalStorage()){this.fetchInProgress=!0,this.fetchFailed=!1;const t=await fetch(`${this.refEndpoint}?type=${this.refType}`,{headers:{Accept:"application/json"}});await this.processResponse(t)}this.isLoading=!1,this.fetchInProgress=!1,this.render()}catch{this.fetchInProgress=!1,this.fetchFailed=!0}}async processResponse(t){if(this.emitStats(t),!t.ok){this.fetchFailed=!0;return}const i=t.clone(),m=await t.json();this.knownItems=m.refs,this.cacheKey=m.cacheKey,this.flushToLocalStorage(await i.text())}emitStats(t){if(!t.ok){(0,g.b)({incrementKey:"REF_SELECTOR_BOOT_FAILED"},!0);return}switch(t.status){case 200:{(0,g.b)({incrementKey:"REF_SELECTOR_BOOTED_FROM_UNCACHED_HTTP"});break}case 304:{(0,g.b)({incrementKey:"REF_SELECTOR_BOOTED_FROM_HTTP_CACHE"});break}default:(0,g.b)({incrementKey:"REF_SELECTOR_UNEXPECTED_RESPONSE"})}}search(t){if(this.searchTerm=t,t===""){this.currentSearchResult=this.knownItems;return}const i=[],m=[];this.exactMatchFound=!1;let b;for(const a of this.knownItems)if(b=a.indexOf(t),!(b<0)){if(b===0){t===a?(m.unshift(a),this.exactMatchFound=!0):m.push(a);continue}i.push(a)}this.currentSearchResult=[...m,...i]}bootstrapFromLocalStorage(){const t=d(this.localStorageKey);if(!t)return!1;const i=JSON.parse(t);return i.cacheKey!==this.cacheKey||!("refs"in i)?(h(this.localStorageKey),!1):(this.knownItems=i.refs,this.isLoading=!1,(0,g.b)({incrementKey:"REF_SELECTOR_BOOTED_FROM_LOCALSTORAGE"}),!0)}async flushToLocalStorage(t){try{p(this.localStorageKey,t)}catch(i){if(i.message.toLowerCase().includes("quota")){this.clearSiblingLocalStorage(),(0,g.b)({incrementKey:"REF_SELECTOR_LOCALSTORAGE_OVERFLOWED"});try{p(this.localStorageKey,t)}catch(m){m.message.toLowerCase().includes("quota")&&(0,g.b)({incrementKey:"REF_SELECTOR_LOCALSTORAGE_GAVE_UP"})}}else throw i}}clearSiblingLocalStorage(){for(const t of Object.keys(localStorage))t.startsWith(c.LocalStoragePrefix)&&h(t)}get localStorageKey(){return`${c.LocalStoragePrefix}:${this.nameWithOwner}:${this.refType}`}},"_SearchIndex");let _=c;_.LocalStoragePrefix="ref-selector"},2077:(L,E,f)=>{f.d(E,{j:()=>w,u:()=>g});function w(d){const p=d.closest("form");if(!(p instanceof HTMLFormElement))return;let h=g(p);if(d.name){const r=d.matches("input[type=submit]")?"Submit":"",c=d.value||r;h||(h=document.createElement("input"),h.type="hidden",h.classList.add("is-submit-button-value"),p.prepend(h)),h.name=d.name,h.value=c}else h&&h.remove()}o(w,"persistSubmitButtonValue");function g(d){const p=d.querySelector("input.is-submit-button-value");return p instanceof HTMLInputElement?p:null}o(g,"findPersistedSubmitButtonValue")},31579:(L,E,f)=>{f.d(E,{Z:()=>g});class w{getItem(){return null}setItem(){}removeItem(){}clear(){}key(){return null}get length(){return 0}}o(w,"NoOpStorage");function g(d,p={throwQuotaErrorsOnSet:!1},h=window){let r;try{r=h[d]}catch{r=new w}const{throwQuotaErrorsOnSet:c}=p;function _(m){try{return r.getItem(m)}catch{return null}}o(_,"getItem");function t(m,b){try{r.setItem(m,b)}catch(a){if(c&&a.message.toLowerCase().includes("quota"))throw a}}o(t,"setItem");function i(m){try{r.removeItem(m)}catch{}}return o(i,"removeItem"),{getItem:_,setItem:t,removeItem:i}}o(g,"safeStorage")},77552:(L,E,f)=>{f.d(E,{b:()=>p});var w=f(70290),g=f(75488);let d=[];function p(t,i=!1){t.timestamp===void 0&&(t.timestamp=new Date().getTime()),t.loggedIn=_(),d.push(t),i?c():r()}o(p,"sendStats");let h=null;async function r(){await g.C,h==null&&(h=window.requestIdleCallback(c))}o(r,"scheduleSendStats");function c(){var t,i;if(h=null,!d.length||(0,w.Z)(document))return;const m=(i=(t=document.head)==null?void 0:t.querySelector('meta[name="browser-stats-url"]'))==null?void 0:i.content;if(!m)return;const b=JSON.stringify({stats:d});try{navigator.sendBeacon&&navigator.sendBeacon(m,b)}catch{}d=[]}o(c,"flushStats");function _(){var t,i;return!!((i=(t=document.head)==null?void 0:t.querySelector('meta[name="user-login"]'))==null?void 0:i.content)}o(_,"isLoggedIn"),document.addEventListener("pagehide",c),document.addEventListener("visibilitychange",c)}},L=>{var E=o(w=>L(L.s=w),"__webpack_exec__");L.O(0,["vendors-node_modules_selector-observer_dist_index_esm_js","vendors-node_modules_virtualized-list_es_index_js-node_modules_github_template-parts_lib_index_js","vendors-node_modules_github_remote-form_dist_index_js-node_modules_delegated-events_dist_inde-1424531"],()=>E(88361));var f=L.O()}]);})();

//# sourceMappingURL=codespaces-03a8b81bb938.js.map