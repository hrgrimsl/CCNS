if("object"==typeof CE2&&CE2.uid)throw"CE: multiple userscripts installed";"undefined"==typeof CE2&&(CE2={}),CE2.uid=115927,"undefined"==typeof CE2&&(CE2={}),CE2.deviceType=function($){var t,e,n=$.toLowerCase(),r=0;if(-1==(t=$.indexOf("(")))return 1;if(t++,-1!=(e=$.indexOf("Android",t))){if(e+=8,$.length>e&&(r=$.charAt(e)))switch(r){case"2":if(-1!=$.indexOf("BNTV",e))return 3;if(-1!=n.indexOf("nook",e))return 3;if(-1!=$.indexOf("Kindle",e))return 3;if(-1!=$.indexOf("Touchpad",e))return 3;break;case"3":return 3;case"4":if(-1!=$.indexOf("Silk",e))return 3}return-1!=n.indexOf("tablet",e)?3:-1!=$.indexOf("Mobi",e)?2:3}if(-1!=(e=$.indexOf("iP",t)))switch(r=$.charAt(e+2)){case"a":return 3;case"h":case"o":return 2}return-1!=(e=$.indexOf("BlackBerry",t))?-1!=$.indexOf("Tablet",e+10)?3:2:-1!=$.indexOf("Windows Phone",t)?2:-1!=$.indexOf("BB10",t)?2:"M"!=$.charAt(0)&&-1!=$.indexOf("Opera Mini",t)?2:1},"undefined"==typeof CE2&&(CE2={}),CE2.ignoredElements=[],CE2.clickCaptors=[],CE2.d=document,CE2.w=window,CE2.n=navigator,CE2.p={},function(){var $=CE2.n.userAgent;/\bMSIE\b/.test($)&&(CE2.ie=1,CE2.ieVersion=parseInt(/MSIE (\d+)\.\d+/.exec($)[1],10),CE2.ieQuirksMode="BackCompat"==CE2.d.compatMode)}(),CE2.ignore=function($){$&&(CE2.ignoredElements.push($),CE2.tracker&&CE2.tracker.ignoredElements.push($))},CE2.capture=function($){CE2.clickCaptors.push($),CE2.tracker&&CE2.tracker.clickCaptors.push($)},CE2.findMatchingSnapshot=function($,t,e,n){var r,i,o,s;if($&&$.length){for(i=0;o=$[i++];)r=Math.floor((new Date).getTime()/1e3),o.e&&o.e<=r||e&&!/n/.test(o.o||"")||CE2.isMatchingSnapshot(o,t,e,n)&&(o.s&&o.s>r?CE2.p[o.id]=o:s||(s=o));return s}},CE2.isMatchingSnapshot=function($,t,e,n){return n?n==$.vid:!$.vid&&CE2.matchURL($.u,e||t,$.o,$.d,CE2.n.userAgent)},CE2.findMatchingLiveSessions=function($,t){var e,n,r=[];if($&&$.length){for(e=0;n=$[e++];)CE2.matchURL(n.u,t,n.o,n.d,CE2.n.userAgent)&&r.push(n.id);return r.length?(r.sort(),r):void 0}},CE2.sameSessions=function($,t){var e,n;if(!$||!t)return!1;if($.length!=t.length)return!1;for(e=0,n=$.length;e<n;e++)if($[e]!=t[e])return!1;return!0},CE2.startTracking=function($,t){var e,n,r;if($){if(CE2.sampleVisit($))CE2.testID=$.id,CE2.testVersion=$.v||1;else if(!t)return;8<=$.v&&(CE2.tracker&&CE2.tracker.cleanup(),CE2.tracker=new CE2.SupernovaTracker(CE2.testID,CE2.testVersion))}t&&(CE2.sessionIDs=t),!t&&$&&8<=$.v||(e=CE2.d.createElement("script"),n="https:"==CE2.w.location.protocol?CE2.TRACKING_SCRIPT_SECURE:CE2.TRACKING_SCRIPT,r="https:"==CE2.w.location.protocol?CE2.TRACKING_DEST_NEW_SECURE:CE2.TRACKING_DEST_NEW,$&&6<=$.v?(img=CE2.d.createElement("img"),img.src=r+"t.js?s="+$.id+"&t="+(new Date).getTime(),e.src=n):e.src=n+($?"?s="+$.id+"&":"?")+"t="+(new Date).getTime(),e.type="text/javascript",e.async=!0,CE2.d.body.appendChild(e))},CE2.unescape=function(t){try{return decodeURIComponent(t)}catch($){return unescape(t)}},CE2.qs2obj=function($,t){if(null==$||/^\s*$/.test($))return null;var e,n,r={},i=null,o=$.replace(/\+/g," ").split(t||"&");for(e=0,n=o.length;e<n;e++)(i=o[e].split("="))[0]&&(r[CE2.unescape(i[0])]=null==i[1]?null:CE2.unescape(i[1]));return r},CE2.each=function($,t,e){var n,r;if($)if("number"==typeof $.length&&"function"==typeof $.concat)for(var i=0,o=$.length;i<o&&(n=$[i],!1!==t.call(e,n,i));i++);else for(r in $)if((n=$[r])!==Object.prototype[r]&&!1===t.call(e,n,r))break},CE2.indexOf=function($,t,e){var n,r;for(n=e||0,r=$.length;n<r;n++)if($[n]===t)return n;return-1},CE2.listen=CE2.addListener=function($,t,e){$.addEventListener?$.addEventListener(t,e,!0):$.attachEvent("on"+t,e)},CE2.removeListener=function($,t,e){$.removeEventListener?$.removeEventListener(t,e,!0):$.detachEvent("on"+t,e)},CE2.userData={},CE2.set=function($,t){1<=($=parseInt($,10))&&$<=5&&(CE2.userData[$]=t+"")},CE2.click=function(){if(CE2.tracker)return CE2.tracker.click.apply(CE2.tracker,arguments)},CE2.getBox=function(){},CE2.sampleVisit=function($){return null==$.r||(!1!==$.r&&!0!==$.r&&(Math.random()>=1/$.r?$.r=!1:$.r=!0),$.r)},CE2.dontTrack=function($,t,e,n){if(n&&void 0!==$.external)try{if(!0===$.external.InPrivateFilteringEnabled())return!0}catch($){}var r=t.doNotTrack||e.doNotTrack||e.msDoNotTrack||$.doNotTrack;return"1"==r||"yes"==r},CE2.cookies=function(){try{return CE2.qs2obj(document.cookie,/;\s*/g)||{}}catch($){return{}}}(),CE2.parseJSON=function(src){return void 0!==JSON&&"function"==typeof JSON.parse?JSON.parse(src):eval("("+src+")")},CE2.isBot=function($){return/alexa|bot|crawl(er|ing)|facebookexternalhit|feedburner|google web preview|nagios|postrank|pingdom|slurp|spider|yahoo!|yandex/i.test($)},CE2.convertToFormData=function($){for(var t=new FormData,e=Object.keys($),n=0;n<e.length;n++)t.append(e[n],$[e[n]]);return t},"undefined"==typeof CE2&&(CE2={}),CE2.READY_STATE_PATTERN=CE2.ie?/complete/:/complete|interactive/,CE2.autoStart="undefined"==typeof CE_MANUAL_START||!CE_MANUAL_START,CE2.domReady=document.readyState&&CE2.READY_STATE_PATTERN.test(document.readyState),CE2.domReadyListeners=[],CE2.onDOMReady=function($){if(CE2.domReady)return setTimeout($,1);CE2.domReadyListeners.push($)},CE2.domReadySetup=function(){var $=function($){for(var t=CE2.domReadyListeners;0<t.length;)t.pop().call();CE2.domReady=!0};if(CE2.domReady&&$(),CE2.listen(window,"load",$),document.addEventListener&&CE2.listen(document,"DOMContentLoaded",$),document.readyState){var t=CE2.READY_STATE_PATTERN;!function(){t.test(document.readyState)?$():setTimeout(arguments.callee,10)}()}},CE2.autoStart&&CE2.domReadySetup(),"undefined"==typeof CE2&&(CE2={}),CE2.matchURL=function($,t,e,n,r){var i,o,s,a,E,C,c,u,f,d,h,p,l,g,_,m,S=/(default|index)($|\..*)/i,v=!1;if(!$||!t)return!1;if(n&&CE2.indexOf(n,CE2.deviceType(r))<0)return!1;if(/n/.test(e=e||""))return $.trim()===t.trim();if(/[re]/.test(e))try{return RegExp($,"i").test(t)}catch($){return!1}if($=new CE2.URI($),i=new CE2.URI(t.toLowerCase()),/h/.test(e)&&$.protocol!=i.protocol)return!1;if(o=(s=i.host).replace(/^www\./,""),u=$.host,f=$.ihost,/w/.test(e)&&s!=u&&s!=f)return!1;if(o!=u.replace(/^www\./,"")&&o!=(f&&f.replace(/^www\./,"")))return!1;if((d=$.path?$.path:"/")!=(a=i.path)){if(/\//.test(e))return!1;for(h=d.split("/"),E=a.split("/"),_=0,m=Math.max(h.length,E.length);_<m;_++)if(h[_]||(h[_]=""),E[_]||(E[_]=""),_==m-1&&(h[_]=h[_].replace(S,""),E[_]=E[_].replace(S,"")),h[_]!=E[_])return!1}return C=i.qs,g=/\?/.test(e),p=$.qs||"",!(g&&C&&!p||!C&&p)&&(CE2.each(p,function($,t){if(C[t]!==$)return!(v=!0)}),!v&&((!g||(CE2.each(C,function($,t){if($!=p[t])return v=!0}),!v))&&(l=$.hash||"",c=i.hash||"",!(g=/#/.test(e))&&!l||l==c)))},"undefined"==typeof CE2&&(CE2={}),void 0===CE2.URI&&(CE2.URI=function($){if(this.src=$,this.protocol=this.host=this.port=this.path=this.qs=this.hash=this.query=null,$){var t=typeof $;"string"==t?this.initWithString($):"object"==t&&this.initWithURI($)}},CE2.URI.pattern=/^\s*([\S]+?:\/\/)?([^\s\/]+?@)?([^:\/\?\#]+)?(\:\d+)?(\/?[^#\?\s]*)?([\?][^#\s]*)?([#]\S+)?/i,CE2.URI.prototype={initWithString:function($){var t,e,n,r,i,o=CE2.URI.pattern.exec($);o[1]||"/"==$.charAt(0)?((t=o[1])&&(this.protocol=t.substr(0,t.indexOf(":"))),this.host=o[3]||null,(e=o[4])&&(this.port=+e.substr(1)),(n=o[5])?this.path=CE2.unescape(n):this.host&&(this.path="/")):this.path=CE2.unescape((o[3]||"")+(o[5]||"")),(r=o[6])&&(this.qs=CE2.qs2obj(r.substr(1)),this.query=r.substr(1)),(i=o[7])&&(this.hash=CE2.unescape(i.substr(1)))},initWithURI:function($){CE2.each($,function($,t){this[t]=$},this)},isAbsolute:function(){return this.isURL()||this.path&&"/"==this.path.charAt(0)},isURL:function(){return this.protocol&&this.host},getDomain:function(){return this.host&&this.host.replace(/^www\./,"")}}),CE2.userMain=function(){try{var r,i=CE2.snapshots,o=CE2.sessions,$=CE2.sites,s=CE2.liveBootstrap||function(){};if(CE2.isBot(CE2.n.userAgent))return;if(CE2.dontTrack(CE2.w,CE2.d,CE2.n,CE2.ie))return;CE2.testID=CE2.testVersion=CE2.sessionIDs=null,CE2.initPageEdits&&CE2.initPageEdits(),CE2.initFlowTracking&&(r=CE2.initFlowTracking()),CE2.showWebsite();var t=function(){var $,t="!$%&()*+,-.0123456789;<=>?@[]^_`{|}~",e={};for($=0;$<36;$++)e[t.charAt($)]=$.toString(36);return e}(),a=function($){return parseInt($.replace(/./g,function($){return t[$]}),36)},e=function($){for(var t,e="",n=/(![^:\/a-z])|([^:\/a-z]{2})|(:[^:\/a-z]{3})|(\/[^:\/a-z]{4})/gi,r=String.fromCharCode;null!=(t=n.exec($));)t[1]||t[2]?e+=r(a(t[0])):t[3]?e+=r(a(t[3].substr(1))):t[4]&&(e+=r(a(t[4].substr(1))));return e};"string"==typeof i&&(i=CE2.parseJSON(e(i))),"string"==typeof o&&(o=CE2.parseJSON(e(o))),"string"==typeof $&&($=CE2.parseJSON(e($))),CE2.recording&&CE2.recording.main&&CE2.recording.main($);var n=function(){try{var $=CE2.w.location.href,t=r&&r.flow&&r.flow.trackByVariant&&r.variant.variantId,e=CE2.findMatchingSnapshot(i,$,"string"==typeof CE_SNAPSHOT_NAME&&CE_SNAPSHOT_NAME.trim(),t),n=CE2.findMatchingLiveSessions(o,$);if(!e&&t&&(e=CE2.findMatchingSnapshot(i,$,"string"==typeof CE_SNAPSHOT_NAME&&CE_SNAPSHOT_NAME.trim())),s())return;if(!e&&!n)return CE2.testID=CE2.testVersion=CE2.sessionIDs=null,void(CE2.tracker&&(CE2.tracker.cleanup(),CE2.tracker=null));(e&&e.id!=CE2.testID||n&&!CE2.sameSessions(n,CE2.sessionIDs))&&(CE2.startTracking(e,n),CE2.badge&&CE2.badge())}catch($){}};n(),CE2.autoStart&&(CE2.monitorInterval=setInterval(n,1e3))}finally{CE2.showWebsite()}},CE2.showWebsite=function(){CE2.bh&&(CE2.bh.parentElement.removeChild(CE2.bh),CE2.bh=null)},CE2.autoStart&&CE2.onDOMReady(CE2.userMain),"function"==typeof CE_READY?CE2.onDOMReady(CE_READY):"object"==typeof CE_READY&&CE2.onDOMReady(function(){CE2.each(CE_READY,function($){"function"==typeof $&&$()})}),CE2.TRACKING_SCRIPT="http://trk.cetrk.com/e/t.js",CE2.TRACKING_SCRIPT_SECURE="https://s3.amazonaws.com/trk.cetrk.com/e/t.js",CE2.TRACKING_DEST="http://trk.cetrk.com/",CE2.TRACKING_DEST_SECURE="https://s3.amazonaws.com/trk.cetrk.com/",CE2.TRACKING_DEST_NEW="https://user-event-tracker.crazyegg.com/",CE2.TRACKING_DEST_NEW_SECURE="https://user-event-tracker.crazyegg.com/",CE2.__CE_HOST__="http://app.crazyegg.com",CE2.sites="&%&-&!&!",CE2.snapshots="%8&4!}%|%]!}$<$4$1$8$2$9$9$1$,!}&.!}$<$7$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}&(%{&2&+%|%[&+&,&&%]%?&2$.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%?%[&,%|&&&%$0&+%{&&&0&(&*%^%_%^&*%^&%%[%^&+!}$,!}&)&+!}$<&4!}&$%^&%&-&,%?%@!}$<!}%?%[%[&&&-&%&,%|&%%_&&!}&6&6&6$,&4!}%|%]!}$<$4$2$1$6$9$4$1$,!}&.!}$<$7$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}&(%{&2&+%|%[&+&,&&%]%?&2$.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%]&&$0$2$1$.$2$1$7$4$0&(&,$.$7$.$4$.$3$1$2$;$1$3$2$6%?$0%_&-&!&!$0!}&6&6$,&4!}%|%]!}$<$4$2$8$7$7$7$5$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}&(%{&2&+%|%[&+&,&&%]%?&2$.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%]&&$0$2$1$.$2$1$7$4$0&(&,$.$7$.$2$.$3$1$2$;$1$5$1$9%?$0%_&-&!&!$0!}&6&6$,&4!}%|%]!}$<$4$2$8$5$1$9$1$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}&(&-%@&!%|&+%{%|&%%`$.%?%|&($.&&&*%`!}$,!}&(%?&,%{!}$<!}$0&(&-%@&!%|%[%?&,%|&&&%&+$0%_%|&%%]$-&,%{%^$-&*%|%`%{&,$-%}&&&-&*&%%?&!$0!}&6&6$,&4!}%|%]!}$<$4$2$8$7$4$6$5$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}&(&-%@&!%|&+%{%|&%%`$.%?%|&($.&&&*%`!}$,!}&(%?&,%{!}$<!}$0&*%^&+&&&-&*%[%^&+$0&!%|%@&*%?&*%|%?&%&+$0&,&&&&&!&+$0&(&*&&&$&&&,%^$0!}&6&6$,&4!}%|%]!}$<$4$1$8$1$1$4$;$,!}&.!}$<$7$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}&(%{&2&+%|%[&+&,&&%]%?&2$.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%]&&%|$0$2$1$.$2$1$7$4$0&(&,$.$4$.$5$2$4$3!}&6&6$,&4!}%|%]!}$<$4$2$8$7$3$6$3$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}&(%{&2&+%|%[&+&,&&%]%?&2$.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%]&&%|$0$2$1$.$2$1$7$4$0&(&,$.$4$.$5$2$9$4!}&6&6$,&4!}%|%]!}$<$4$2$8$7$3$8$1$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}&(%{&2&+%|%[&+&,&&%]%?&2$.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%]&&%|$0$2$1$.$2$1$7$4$0&(&,$.$4$.$5$2$9$1!}&6&6$,&4!}%|%]!}$<$4$2$8$7$3$8$4$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}&(%{&2&+%|%[&+&,&&%]%?&2$.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%]&&%|$0$2$1$.$2$1$7$4$0&(&,$.$4$.$5$2$9$2!}&6&6$,&4!}%|%]!}$<$4$2$8$7$4$1$1$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}&(%{&2&+%|%[&+&,&&%]%?&2$.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%]&&%|$0$2$1$.$2$1$7$4$0&(&,$.$4$.$5$2$9$3!}&6&6$,&4!}%|%]!}$<$4$2$8$7$4$2$4$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}&(%{&2&+%|%[&+&,&&%]%?&2$.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%]&&%|$0$2$1$.$2$1$7$4$0&(&,$.$4$.$5$2$8$6!}&6&6$,&4!}%|%]!}$<$4$2$8$7$7$6$1$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?&.&+$.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0&,&&%[$0%}&.%@$0%[&&&!&!%^%[&,%|&&&%$0%[&-&*&*%^&%&,!}&6&6$,&4!}%|%]!}$<$4$2$9$2$8$6$1$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?%|&($.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0&,&&%[$0%?&(&!$0%[&&&!&!%^%[&,%|&&&%$0%[&-&*&*%^&%&,!}&6&6$,&4!}%|%]!}$<$4$2$9$2$8$6$7$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?%|&($.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0&,&&%[$0%}%[&($0%[&&&!&!%^%[&,%|&&&%$0%[&-&*&*%^&%&,!}&6&6$,&4!}%|%]!}$<$4$2$9$2$8$7$3$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?%|&($.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0&,&&%[$0%}%?&($0%[&&&!&!%^%[&,%|&&&%$0%[&-&*&*%^&%&,!}&6&6$,&4!}%|%]!}$<$4$1$7$;$;$3$6$,!}&.!}$<$7$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}&(%{&2&+%|%[&+&,&&%]%?&2$.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0&,&&%[$0&(&,&&$0%[&-&*&*%^&%&,!}&6&6$,&4!}%|%]!}$<$4$2$7$9$8$5$6$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?%|&($.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%}%?&($0%?&-&,%{&&&*&+$0&$%?&%&-&+%[&*%|&(&,!}&6&6$,&4!}%|%]!}$<$4$2$7$;$6$;$9$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?%|&($.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%?&(&!$0%?&-&,%{&&&*&+$0&$%?&%&-&+%[&*%|&(&,!}&6&6$,&4!}%|%]!}$<$4$2$8$5$2$3$3$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}&(&-%@&!%|&+%{%|&%%`$.%?%|&($.&&&*%`!}$,!}&(%?&,%{!}$<!}$0&*%^&+&&&-&*%[%^&+$0&!%|%@&*%?&*%|%?&%&+$0!}&6&6$,&4!}%|%]!}$<$4$2$8$5$2$4$;$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?%|&($.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%]&&%|$0$2$1$.$2$1$7$4$0$2$.$6$1$;$9$5$7$6!}&6&6$,&4!}%|%]!}$<$4$2$8$5$2$;$8$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?%|&($.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%]&&%|$0$2$1$.$2$1$7$4$0$2$.$6$1$;$9$5$8$4!}&6&6$,&4!}%|%]!}$<$4$2$8$4$5$;$8$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}&(%{&2&+%|%[&+&,&&%]%?&2$.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%}&&&-&*&%%?&!$0&(&,&&!}&6&6$,&4!}%|%]!}$<$4$2$8$3$;$8$5$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%@%|&&%^&%%`%|&%%^%^&*%|&%%`&,&&%]%?&2$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%^&$%^&*%`%|&%%`$-&,%^%[%{!}&6&6$,&4!}%|%]!}$<$4$2$8$4$;$2$2$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?%|&($.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%?%[&,%|&&&%$0&*%^%`%|&+&,&*%?&,%|&&&%!}&6&6$,&4!}%|%]!}$<$4$2$7$4$;$4$4$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&*!}$<$3$6$,!}&-!}$<!}%<%{&,&,&(&+$<$0$0%?&.&+%9%9$.&+%[%|&,%?&,%|&&&%%9%9$.&&&*%`$0%]&&%|$0%_&-&!&!$0$2$1%9%9$.$2$2$2$7$0$.$*$!!}$,!}&&!}$<!}&*!}&6$,&4!}%|%]!}$<$4$2$7$9$8$3$2$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&*!}$<$3$6$,!}&-!}$<!}%<%{&,&,&(&+$<$0$0%?&+%?%9%9$.&+%[%|&,%?&,%|&&&%%9%9$.&&&*%`$0%]&&%|$0%_&-&!&!$0$2$1%9%9$.$2$2$3$2$0$.$*$!!}$,!}&&!}$<!}&*!}&6$,&4!}%|%]!}$<$4$2$8$5$1$2$;$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&*!}$<$3$6$,!}&-!}$<!}%<%{&,&,&(&+$<$0$0%?%|&(%9%9$.&+%[%|&,%?&,%|&&&%%9%9$.&&&*%`$0%]&&%|$0%_&-&!&!$0$2$1%9%9$.$2$1$7$4$0$.$*$!!}$,!}&&!}$<!}&*!}&6$,&4!}%|%]!}$<$4$2$9$5$9$6$9$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<!}%<%{&,&,&(&+$<$0$0&!%|%?%9%9$.&+%[%|&,%?&,%|&&&%%9%9$.&&&*%`$0%]&&%|$0%_&-&!&!$0$2$1%9%9$.$3$4$6$2$0$.$*$!!}$,!}&&!}$<!}&*!}&6$,&4!}%|%]!}$<$4$2$9$5$9$9$5$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<!}%<%{&,&,&(&+$<$0$0&+&&&*%9%9$.&+%[%|&,%?&,%|&&&%%9%9$.&&&*%`$0%]&&%|$0%_&-&!&!$0$2$1%9%9$.$2$2$3$3$0$.$*$!!}$,!}&&!}$<!}&*!}&6$,&4!}%|%]!}$<$4$2$7$9$5$6$5$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&*!}$<$2$1$,!}&-!}$<!}%<%{&,&,&(&+$<$0$0%?&.&+%9%9$.&+%[%|&,%?&,%|&&&%%9%9$.&&&*%`$0%]&&%|$0%?%@&+$0$2$1%9%9$.$2$2$2$7$0$.$*$!!}$,!}&&!}$<!}&*!}&6$,&4!}%|%]!}$<$4$2$7$;$6$;$3$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?%|&($.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%}%?&($0%|&%%_&&$0&(&&&!%|%[%|%^&+!}&6&6$,&4!}%|%]!}$<$4$2$7$;$7$1$4$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?%|&($.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%?&(&!$0%|&%%_&&$0&(&&&!%|%[%|%^&+!}&6&6$,&4!}%|%]!}$<$4$2$8$5$1$1$6$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&*!}$<$3$6$,!}&-!}$<!}%<%{&,&,&(&+$<$0$0%?%|&(%9%9$.&+%[%|&,%?&,%|&&&%%9%9$.&&&*%`$0%]&&%|$0%?%@&+$0$2$1%9%9$.$2$1$7$4$0$.$*$!!}$,!}&&!}$<!}&*!}&6$,&4!}%|%]!}$<$4$2$8$7$6$9$4$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&*!}$<$3$1$,!}&-!}$<!}%<%{&,&,&(&+$<$0$0%?&+%?%9%9$.&+%[%|&,%?&,%|&&&%%9%9$.&&&*%`$0%]&&%|$0%?%@&+$0$2$1%9%9$.$2$2$3$2$0$.$*$!!}$,!}&&!}$<!}&*!}&6$,&4!}%|%]!}$<$4$2$8$7$7$8$4$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%^&+&(%^%[&,&*%?$.%?%|&($.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%?&(&($0!}$,!}%{%?&+%{!}$<!}$0%?%[%[%^&+&+&&&(&,%|&&&%!}&6&6$,&4!}%|%]!}$<$4$2$9$5$9$;$1$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&*!}$<$2$1$,!}&-!}$<!}%<%{&,&,&(&+$<$0$0&+&&&*%9%9$.&+%[%|&,%?&,%|&&&%%9%9$.&&&*%`$0%]&&%|$0%?%@&+$0$2$1%9%9$.$2$2$3$3$0$.$*$!!}$,!}&&!}$<!}&*!}&6$,&4!}%|%]!}$<$4$2$7$9$5$8$7$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?%|&($.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%?&(&!$0%|&%%_&&$0%[&&&%&,%?%[&,!}&6&6$,&4!}%|%]!}$<$4$2$7$9$6$4$1$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&*!}$<$6$1$,!}&-!}$<!}%<%{&,&,&(&+$<$0$0%?%[%?%9%9$.&+%[%|&,%?&,%|&&&%%9%9$.&&&*%`$0%]&&%|$0%?%@&+$0$2$1%9%9$.$2$1$7$4$.$*$!!}$,!}&&!}$<!}&*!}&6$,&4!}%|%]!}$<$4$2$7$9$8$4$8$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?%|&($.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%}%?&($0%|&%%_&&$0%[&&&%&,%?%[&,!}&6&6$,&4!}%|%]!}$<$4$2$8$4$2$2$1$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?&+%?$.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0&,&&%[$0%}%?&+$0%[&-&*&*%^&%&,!}&6&6$,&4!}%|%]!}$<$4$2$8$5$2$3$7$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}&(&-%@&!%|&+%{%|&%%`$.%?%|&($.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%?%@&&&-&,$0%[%?&*%^%^&*&+$0!}&6&6$,&4!}%|%]!}$<$4$2$8$5$2$4$5$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}&(&-%@&!%|&+%{%|&%%`$.%?%|&($.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%?%@&&&-&,$0%[&&&%&,%?%[&,$0!}&6&6$,&4!}%|%]!}$<$4$2$8$7$6$2$2$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%^&+&(%^%[&,&*%?$.%?%|&($.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%?&(&($0!}$,!}%{%?&+%{!}$<!}$0&-&+%^&*&$%?&%&-%?&!!}&6&6$,&4!}%|%]!}$<$4$2$9$5$9$7$2$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}&+&&&*$.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0&,&&%[$0%}&&&*$0%[&-&*&*%^&%&,!}&6&6$,&4!}%|%]!}$<$4$2$7$7$3$2$6$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?%|&($.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%}%[&($0%|&%%_&&$0%_&&%[&-&+!}&6$,!}%]!}$<%8$2%;&6$,&4!}%|%]!}$<$4$2$7$9$5$8$2$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?%|&($.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%?&(&!$0%|&%%_&&$0%_&&%[&-&+!}&6&6$,&4!}%|%]!}$<$4$2$7$;$6$;$7$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?%|&($.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%}%?&($0%|&%%_&&$0%_&&%[&-&+!}&6&6$,&4!}%|%]!}$<$4$2$8$3$;$3$8$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%@%|&&%^&%%`%|&%%^%^&*%|&%%`&,&&%]%?&2$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%]%|&+%^%?&+%^&+!}&6&6$,&4!}%|%]!}$<$4$2$7$9$9$1$1$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?%|&($.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%}%?&($0%|&%%_&&$0&%%^&0&+!}&6&6$,&4!}%|%]!}$<$4$2$8$3$6$;$9$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?%|&($.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%?&(&!$0%|&%%_&&$0&%%^&0&+!}&6&6$,&4!}%|%]!}$<$4$2$8$3$8$;$8$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%@%|&&%^&%%`%|&%%^%^&*%|&%%`&,&&%]%?&2$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%[&&&%&,%?%[&,!}&6&6$,&4!}%|%]!}$<$4$2$8$3$;$2$5$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%@%|&&%^&%%`%|&%%^%^&*%|&%%`&,&&%]%?&2$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0&,%{%^&*%?&(&2!}&6&6$,&4!}%|%]!}$<$4$2$8$3$;$9$2$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%@%|&&%^&%%`%|&%%^%^&*%|&%%`&,&&%]%?&2$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%|&$%?%`%|&%%`!}&6&6$,&4!}%|%]!}$<$4$2$8$7$5$2$3$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%^&+&(%^%[&,&*%?$.%?%|&($.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%?&(&($0!}$,!}%{%?&+%{!}$<!}$0&*%^%`%|&+&,%^&*!}&6&6$,&4!}%|%]!}$<$4$2$7$7$3$4$6$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?%?&(&,$.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%}&&&-&*&%%?&!$0&(&,%^!}&6$,!}%]!}$<%8$2%;&6$,&4!}%|%]!}$<$4$2$7$7$3$4$7$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?%?&(&,$.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%}&&&-&*&%%?&!$0%?%}&(!}&6$,!}%]!}$<%8$2%;&6$,&4!}%|%]!}$<$4$2$8$3$8$5$6$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%@%|&&%^&%%`%|&%%^%^&*%|&%%`&,&&%]%?&2$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0&+%?%_%^&,&2!}&6&6$,&4!}%|%]!}$<$4$2$8$4$;$9$3$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}&0&0&0$.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0&(&-%@&!%|%[%?&,%|&&&%&+!}&6&6$,&4!}%|%]!}$<$4$2$9$2$7$7$9$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%^&+&(%^%[&,&*%?$.%?%|&($.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%?&(&($0!}$,!}%{%?&+%{!}$<!}$0&!%?&%%]%|&%%`!}&6&6$,&4!}%|%]!}$<$4$2$7$4$;$1$7$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?%[%?$.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%}&&&-&*&%%?&!$0&+%]&2!}&6&6$,&4!}%|%]!}$<$4$2$7$7$2$1$3$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&*!}$<$6$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?%|&($.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%}&&&-&*&%%?&!$0%}%[&(!}&6$,!}%]!}$<%8$2%;&6$,&4!}%|%]!}$<$4$2$7$7$2$8$9$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&*!}$<$3$6$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?%|&($.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%}&&&-&*&%%?&!$0%}&$&(!}&6$,!}%]!}$<%8$2%;&6$,&4!}%|%]!}$<$4$2$7$7$2$9$4$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&*!}$<$3$6$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?%|&($.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%}&&&-&*&%%?&!$0&(%{&(!}&6$,!}%]!}$<%8$2%;&6$,&4!}%|%]!}$<$4$2$7$7$2$;$8$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&*!}$<$3$6$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?%|&($.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%}&&&-&*&%%?&!$0&(%{%_!}&6$,!}%]!}$<%8$2%;&6$,&4!}%|%]!}$<$4$2$7$7$3$1$3$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?%|&($.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%}&&&-&*&%%?&!$0&*&+%^!}&6$,!}%]!}$<%8$2%;&6$,&4!}%|%]!}$<$4$2$7$7$3$2$2$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?%|&($.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%}&&&-&*&%%?&!$0%@&$%_!}&6$,!}%]!}$<%8$2%;&6$,&4!}%|%]!}$<$4$2$7$7$3$2$4$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?%|&($.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%}&&&-&*&%%?&!$0%[%{%?!}&6$,!}%]!}$<%8$2%;&6$,&4!}%|%]!}$<$4$2$8$3$6$;$;$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?&.&+$.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%}&&&-&*&%%?&!$0&+&+&+!}&6&6$,&4!}%|%]!}$<$4$2$8$3$7$2$6$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?&.&+$.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%}&&&-&*&%%?&!$0%}&.%?!}&6&6$,&4!}%|%]!}$<$4$2$8$3$7$3$8$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?&.&+$.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%}&&&-&*&%%?&!$0%}&.%@!}&6&6$,&4!}%|%]!}$<$4$2$8$3$7$7$3$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?&.&+$.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%}&&&-&*&%%?&!$0%@%|&(!}&6&6$,&4!}%|%]!}$<$4$2$8$3$;$1$2$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%@%|&&%^&%%`%|&%%^%^&*%|&%%`&,&&%]%?&2$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%]&*&-%`&+!}&6&6$,&4!}%|%]!}$<$4$2$8$3$;$7$9$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%@%|&&%^&%%`%|&%%^%^&*%|&%%`&,&&%]%?&2$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%?%@&&&-&,!}&6&6$,&4!}%|%]!}$<$4$2$8$4$1$1$3$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?&+%?$.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%}&&&-&*&%%?&!$0&+&&&-!}&6&6$,&4!}%|%]!}$<$4$2$8$4$1$6$5$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?&+%?$.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%}&&&-&*&%%?&!$0&(&$%?!}&6&6$,&4!}%|%]!}$<$4$2$8$4$1$7$5$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?&+%?$.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%}&&&-&*&%%?&!$0%}%^&!!}&6&6$,&4!}%|%]!}$<$4$2$8$4$1$;$4$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?%|&($.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%}&&&-&*&%%?&!$0&$&$&$!}&6&6$,&4!}%|%]!}$<$4$2$8$4$2$1$6$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?&+%?$.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%}&&&-&*&%%?&!$0&%&&%[!}&6&6$,&4!}%|%]!}$<$4$2$8$4$2$4$3$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?&+%?$.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%}&&&-&*&%%?&!$0%?&*&!!}&6&6$,&4!}%|%]!}$<$4$2$8$4$6$1$1$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?%|&($.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%}&&&-&*&%%?&!$0%?&(&$!}&6&6$,&4!}%|%]!}$<$4$2$8$4$6$1$3$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%@%|&&%^&%%`%|&%%^%^&*%|&%%`&,&&%]%?&2$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%@&*%?%|&%!}&6&6$,&4!}%|%]!}$<$4$2$8$4$6$6$6$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?%|&($.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%}&&&-&*&%%?&!$0%?&(&(!}&6&6$,&4!}%|%]!}$<$4$2$8$4$6$7$5$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?%|&($.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%}&&&-&*&%%?&!$0%?&(%@!}&6&6$,&4!}%|%]!}$<$4$2$8$4$6$8$7$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?%|&($.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%}&&&-&*&%%?&!$0%?%]&.!}&6&6$,&4!}%|%]!}$<$4$2$8$4$6$9$1$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?%|&($.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%}&&&-&*&%%?&!$0%?&*%^!}&6&6$,&4!}%|%]!}$<$4$2$8$4$7$2$6$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?%|&($.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%}&&&-&*&%%?&!$0%}%?&(!}&6&6$,&4!}%|%]!}$<$4$2$8$4$7$7$5$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?%|&($.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%}&&&-&*&%%?&!$0%}&(&*!}&6&6$,&4!}%|%]!}$<$4$2$8$4$7$7$6$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?%|&($.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%}&&&-&*&%%?&!$0&*&+%|!}&6&6$,&4!}%|%]!}$<$4$2$8$4$7$7$8$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?%|&($.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%}&&&-&*&%%?&!$0&!&,&(!}&6&6$,&4!}%|%]!}$<$4$2$8$4$9$9$3$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?&+%?$.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%}&&&-&*&%%?&!$0%}%?&+!}&6&6$,&4!}%|%]!}$<$4$2$8$4$9$9$7$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?%|&($.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%}&&&-&*&%%?&!$0%?&(&!!}&6&6$,&4!}%|%]!}$<$4$2$8$7$6$6$2$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%^&+&(%^%[&,&*%?$.%?%|&($.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%?&(&($0!}$,!}%{%?&+%{!}$<!}$0&+%^%?&*%[%{!}&6&6$,&4!}%|%]!}$<$4$2$7$9$6$2$6$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&*!}$<$6$1$,!}&-!}$<!}%<%{&,&,&(&+$<$0$0%?%[%?%9%9$.&+%[%|&,%?&,%|&&&%%9%9$.&&&*%`$0%]&&%|$0%_&-&!&!$0$.$*$!!}$,!}&&!}$<!}&*!}&6$,&4!}%|%]!}$<$4$2$8$4$6$5$7$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%@%|&&%^&%%`%|&%%^%^&*%|&%%`&,&&%]%?&2$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%@&&%]&2!}&6&6$,&4!}%|%]!}$<$4$2$7$4$;$2$7$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&*!}$<$3$6$,!}&-!}$<!}%<%{&,&,&(&+$<$0$0%[&(&+%9%9$.&+%[%|&,%?&,%|&&&%%9%9$.&&&*%`$0%]&&%|$0%?%@&+$0$.$*$!!}$,!}&&!}$<!}&*!}&6$,&4!}%|%]!}$<$4$2$8$7$4$9$5$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%^&+&(%^%[&,&*%?$.%?%|&($.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%?&(&($0!}$,!}%{%?&+%{!}$<!}$0%{&&&$%^!}&6&6$,&4!}%|%]!}$<$4$2$8$4$;$8$;$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}&0&0&0$.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0&(&*%|&.%?%[&2!}&6&6$,&4!}%|%]!}$<$4$2$8$5$3$2$8$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}&(&-%@&!%|&+%{%|&%%`$.%?%|&($.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%[%{%|&%%?$0!}&6&6$,&4!}%|%]!}$<$4$2$8$4$9$;$2$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%@%|&&%^&%%`%|&%%^%^&*%|&%%`&,&&%]%?&2$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0!}&6&6$,&4!}%|%]!}$<$4$2$7$4$;$7$5$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&*!}$<$2$1$,!}&-!}$<!}%<%{&,&,&(&+$<$0$0%?&.&+%9%9$.&+%[%|&,%?&,%|&&&%%9%9$.&&&*%`$0&,&&%[$0$.$*$!!}$,!}&&!}$<!}&*!}&6$,&4!}%|%]!}$<$4$2$7$9$5$2$5$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<!}%<%{&,&,&(&+$<$0$0%[&(&+%9%9$.&+%[%|&,%?&,%|&&&%%9%9$.&&&*%`$0&,&&%[$0$.$*$!!}$,!}&&!}$<!}&*!}&6$,&4!}%|%]!}$<$4$2$7$9$6$4$5$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&*!}$<$6$1$,!}&-!}$<!}%<%{&,&,&(&+$<$0$0%?%[%?%9%9$.&+%[%|&,%?&,%|&&&%%9%9$.&&&*%`$0&,&&%[$0$.$*$!!}$,!}&&!}$<!}&*!}&6$,&4!}%|%]!}$<$4$2$8$4$;$2$7$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}&0&0&0$.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0&,%^&*&$&+!}&6&6$,&4!}%|%]!}$<$4$2$8$5$1$6$4$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}&0&0&0$.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%?%@&&&-&,!}&6&6$,&4!}%|%]!}$<$4$2$7$4$;$1$1$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&*!}$<$3$6$,!}&-!}$<!}%<%{&,&,&(&+$<$0$0%?%|&(%9%9$.&+%[%|&,%?&,%|&&&%%9%9$.&&&*%`$0&,&&%[$.$*$!!}$,!}&&!}$<!}&*!}&6$,&4!}%|%]!}$<$4$2$8$4$;$7$9$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}&0&0&0$.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0%{%^&!&(!}&6&6$,&4!}%|%]!}$<$3$;$5$4$2$6$4$,!}&.!}$<$7$,!}%^!}$<$2$6$8$6$1$1$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}&,%^&+&,$.&+&,%|&1%_&&&%&,&+$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0!}&6&6$,&4!}%|%]!}$<$4$2$8$5$2$4$8$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}&(&-%@&!%|&+%{%|&%%`$.%?%|&($.&&&*%`!}$,!}&(%?&,%{!}$<!}$0!}&6&6$,&4!}%|%]!}$<$4$2$8$3$7$6$;$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?&.&+$.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0!}&6&6$,&4!}%|%]!}$<$4$2$8$4$1$;$1$,!}&.!}$<$8$,!}%^!}$<$2$6$7$2$;$6$4$7$1$1$,!}&-!}$<&4!}&(&*&&&,&&%[&&&!!}$<!}%{&,&,&(&+!}$,!}%{&&&+&,!}$<!}%?&+%?$.&+%[%|&,%?&,%|&&&%$.&&&*%`!}$,!}&(%?&,%{!}$<!}$0!}&6&6%;",CE2.sessions=null,CE2.PAGE_VIEWS_LIMIT_REACHED=!1;