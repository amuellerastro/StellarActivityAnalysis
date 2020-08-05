function compute_FAP_GLS, jd, f0, f1, FAPlevel

  M = (max(jd)-min(jd)) * abs(f1-f0)

  powlevel = dblarr(n_elements(FAPlevel))

;   for i=0,n_elements(FAPlevel)-1 do begin
; 
;     if (FAPlevel[i] gt 0.01d0) then powlevel[i] = 1.d0-(1.d0-FAPlevel[i])^(1.d0/M) $   ;powlevel[i] = 1.d0 - (1.d0 - (1.d0 - FAPlevel[i])^(1.d0/M))^(2.d0/(n_elements(jd)-3.d0)) $
;       else powlevel[i] = 1.d0-(FAPlevel[i]/M)^(2.d0/(n_elements(jd)-3.d0))
; 
;   endfor
    powlevel[i] = 1.d0-(FAPlevel[i]/M)^(2.d0/(n_elements(jd)-3.d0))

  return, powlevel

end

pro RVSPY_compute_periodograms, resdir, ofn, jd, data, edata, f0, f1, ofac, weight, GLS=GLS, CLEAN=CLEAN, BGLS1=BGLS1

readcol, 'setup.ini', tmp, format='a', comment='#', /silent
pathGLS = tmp[2]
epstopdf = tmp[3]+'epstopdf.pl '

; epstopdf = '/home/amueller/work/IDLlibs/epstopdf.pl '

;=========================================================================

; GLS

if keyword_set(GLS) then begin

  pathGLS = '/home/amueller/work/IDLlibs/SpecRoutines/GLS_v2.3.02/'
;   ofac = 10
;   weight = 2
  mass = 0.8	;dummy

  fn = ofn+'.txt'
  openw, lun, pathGLS+fn, width=1400, /get_lun
    for i=0,n_elements(jd)-1 do printf, lun, jd[i], data[i], edata[i], format='(f17.9,f23.15,f23.15)'
  close, lun
  free_lun, lun

  openw, lun, 'GLS.par', /get_lun
    printf, lun, fix(ofac)
    printf, lun, fix(weight)
    printf, lun, '0 0.0 0.0'
    printf, lun, '0 359'
  close, lun
  free_lun, lun
  
	spawn, pathGLS+'./GLS '+pathGLS+fn+' '+strcompress(f0, /rem)+' '+strcompress(f1, /rem)+' '+strcompress(mass, /rem)

  readcol, 'GLS.plt', freq, pow, win, pls, format='d,d,d,d', /silent

  spawn, 'mv GLS.plt '+resdir+'GLS_'+ofn+'.txt'
  spawn, 'rm GLS.log GLS.par RVSinRes.plt GLSFit.plt'

  FAPlevel = [1.d-3]	;10%, 1%, 0.1%
  FAPlevel_st = ['10!U-3!N']
;   powlevel = compute_FAP_GLS(jd, f0, f1, FAPlevel)

;   if (round(1./f0) ge 10000.) then posfap = 7000
;   if (round(1./f0) eq 1000.) then posfap = 700
;   if (round(1./f0) eq 100.) then posfap = 70
  if (f1 eq 10.) then posfapfreq = 9.5
  if (f1 eq 1.) then posfapfreq = 0.95

  yplotrange = [0,1.]
  period = 1./freq

  set_plot, 'ps'
    fn = resdir+'GLS_'+ofn+'.ps'
    device, filename=fn,/color,XSIZE=20, YSIZE=14
    !P.Font=0
    !p.thick=4
    !x.thick=3
    !y.thick=3
    !p.multi=[0,1,2]

    
  ;---Periodogram LOG scaling----------------------
    plot, /nodata, /XLOG, period, pow, ytit='GLS Power', yst=1, xr=[1./f1,1./f0], xst=1, charsize=1.3, yminor=2, xticklen=0.07, xtickformat="(A1)", pos=[0.11,0.465,0.95,0.97], yr=[0,1];, yr=[0,round(max(pow)*100.)/100.]
    oplot, period, pow
    ;overplot FAP level
    off = 0.01
;     for i=0,n_elements(FAPlevel)-1 do begin
; 
; ; 			if (powlevel[i] lt max(pow)) then begin
;     
; 				plots, [1./f1,1./f0], [powlevel[i], powlevel[i]], linestyle=2, color=cgcolor('dark gray')
; ; 				xyouts, posfap, powlevel[i]+off, FAPlevel_st[i], align=1, /data, charsize=1.3
; 
; ; 			endif
; 			
;     endfor
    legend, [ofn], box=0, margin=0, /top, /right, charsize=1.5

	;---window function----------------------
    plot, /nodata, /XLOG, period, win, ytit='Window Power',xtitle='Period [days]', yst=1, xr=[1./f1,1./f0], xst=1, yr=[0,round(max(win)*100.)/100.], $	;, xtickname=ticknames, yr=[0,1], yr=[0,1]
		background=cgcolor('white'),Color=cgcolor('black'), charsize=1.3, yminor=1, xticklen=0.09, pos=[0.11,0.15,0.95,0.455];, yr=[0,0.5];, yticks=4
    oplot, period, win
   
    
  device,/close
  set_plot,'x'
  !p.thick=1
  !x.thick=1
  !y.thick=1
  !p.multi=[0,1,0]

  ;spawn, 'epstopdf '+'GLS/EW/GLS_'+line[xx]+'.ps'
  spawn, epstopdf+resdir+'GLS_'+ofn+'.ps'
  spawn, 'rm '+resdir+'GLS_'+ofn+'.ps'
  
endif

;=========================================================================


;=========================================================================


end
