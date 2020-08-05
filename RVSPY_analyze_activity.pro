@mrdfits
@fxpar
@valid_num
@mrd_skip
@ts_diff
@linspace
@rebinw
@findex
@RVSPY_compute_periodograms.pro
@showsym
@readcol
@remchar
@gettok
@strsplit
@strnumber
@interpol
@closest2
@tsum
@robust_poly_fit
@poly_fit
@rob_checkfit
@robust_sigma
@poly
@proceeding_text
@ploterror
@setdefaultvalue
@cgplot
@cgsetcolorstate
@cgbitget
@convert_to_type
@cgcheckforsymbols
@cgdefaultcolor
@colorsareidentical
@cgdefcharsize
@ps_background
@cgsymcat
@robust_linefit
@linear_corrcoeff_probability
@correlate
@qsimp
@zparcheck
@trapzd
@gamma
@sigfig
@legend
@cgsymbol
@cgcolor24
@resistant_mean
@curvefit
@headfits
@fxposit
@mrd_hread
@get_eso_keyword
@cgcolor
@cggetcolorstate
@cgsnapshot
@setdecomposedstate
@decomposedcolor
@cgcolorbar
@cgcolorfill
@cgaxis
@combigen.pro
@loadct
; @filepath
; @path_sep
@uniq
; @crsprod_weightedCCF.pro
; @crstrim.pro
@bgls.pro
@reverse.pro
@mpfitpeak.pro
@mpfit.pro
@mpfitfun.pro
@fxmove.pro
@match.pro
@mrd_struct.pro
@mean.pro

;Here, A0 is the height of the Gaussian, A1 is the center of the Gaussian, A2 is the width (the standard deviation) of the Gaussian, A3 is the constant term, A4 is the linear term, and A5 is the quadratic term. 
function fit_ccf_gauss_3terms, x, p
  fit = p[0]*exp(-0.5d0*((x-p[1])/p[2])^2.)
  return, fit
end
function fit_ccf_gauss_4terms, x, p
  fit = p[0]*exp(-0.5d0*((x-p[1])/p[2])^2.)+p[3]
  return, fit
end
function fit_ccf_gauss_5terms, x, p
  fit = p[0]*exp(-0.5d0*((x-p[1])/p[2])^2.)+p[3]+(p[4]*x)
  return, fit
end
function fit_ccf_gauss_6terms, x, p
  fit = p[0]*exp(-0.5d0*((x-p[1])/p[2])^2.)+p[3]+(p[4]*x)+(p[5]*x^2.)
  return, fit
end

function poly_2deg, x, p

  return, p[0]+p[1]*x+p[2]*x*x

end

pro poly_2deg_v2, x, A, F, pder

  F = A[0]+A[1]*x+A[2]*x^2.
  pder = [[replicate(1.0, n_elements(x))],[x],[x^2.]]

end

function n_Edlen, l

  sigma = 1d4 / l
  sigma2 = sigma*sigma
	n = 1 + 1d-8 * (8342.13 + 2406030.d0 / (130.d0-sigma2) + 15997.d0/(38.9-sigma2))
  return, n

end

function ToAir, l

	return, (l / n_Edlen(l))

end

function toVacuum, l

  cond = 1
  l_prev = l
;   while(cond eq 1) do begin
  l_new = n_Edlen(l_prev) * l
;   if (max(abs(l_new - l_prev)) lt 1e-10) then cond = 0
;   endwhile
  l_prev = l_new
  return, l_prev

end

pro RVSPY_analyze_activity;, prepare=prepare, twhya=twhya

;ToDO:
;correlate individual spectrum with mean spectrum to get better overlap?
;get polynomial coefficients for LDR to convert to Teff and average them, eg Biazzo2007, Catalano 2002
;line index for Ha, Ca IRT, e.g. Fuhrmeister 2019 (Carmenes), Robertson2016,832,112
;photospheric band indices can only be used with ceres39 but no robust wl calibration

;for FWHM interpretation see Santos2010

readcol, 'setup.ini', tmp, format='a', comment='#', /silent
base = tmp[0]

readcol, tmp[1], id, format='a', /silent

pathGLS = tmp[2]
epstopdf = tmp[3]+'epstopdf.pl '

f0 = double(tmp[6])
f1 = double(tmp[7])
ofac = double(tmp[5])
weight = double(tmp[8])
mass = double(tmp[9])

qprepare = ''
read, 'Prepare data? y/n: ', qprepare
if (qprepare eq 'n') then prepare = 0 else prepare = 1

; phi = '!Mf!X'
alpha = '!Ma!X'
; beta = '!Mb!X'
; gamma = '!Mg!X'
; delta = '!Md!X'
ldelta = '!MD!X'
; angstrom = '!3' + STRING(197B) + '!X'

st1 = mrdfits(tmp[4], 1, /silent)
st2 = mrdfits(tmp[4], 2, /silent)
st3 = mrdfits(tmp[4], 3, /silent)

ewwin = st1.window_half_width
winvcorr = st3.WINDOW_HALF_WIDTH

; ewwin = [7.,7.,5.,1.3,10.];,4.,5.]	;0.5*lambda size for EW
; winvcorr = [380., 535., 490., 435., 260.,100., 500.];, 150., 180.]	;velocity correlation, line profile variance,, 0.5*size


for iid=0,n_elements(id)-1 do begin

    resdir = base+id[iid]+'/Results/'
    spawn, 'mkdir -p '+resdir

;     epstopdf = '/home/amueller/work/IDLlibs/epstopdf.pl '

    ;for periodograms
;     pathGLS = '/home/amueller/work/IDLlibs/SpecRoutines/GLS_v2.3.02/'
;     f0 = 2.d-5	;2.d-5  ;1.d-3
;     f1 = 1
;     ofac = 20
;     weight = 2	;weightening exponent (variance: wexp=0, chi2: wexp=2)
;     mass = 1.0  ;dummy

    resdirrv = resdir+'RV/'
    spawn, 'mkdir -p '+resdirrv
    resdirldr = resdir+'LDR/'
    spawn, 'mkdir -p '+resdirldr
    resdirew = resdir+'EW/'
    spawn, 'mkdir -p '+resdirew
    resdirrf = resdir+'ResidualFlux/'
    spawn, 'mkdir -p '+resdirrf
;     resdirpbi = resdir+'PBI/'
;     spawn, 'mkdir -p '+resdirpbi
    resdirvcorr = resdir+'VelCorr/'
    spawn, 'mkdir -p '+resdirvcorr
    resdirlpv = resdir+'LineProfileVariance/'
    spawn, 'mkdir -p '+resdirlpv
		resdirindex = resdir+'LineIndex/'
    spawn, 'mkdir -p '+resdirindex

    ;for EW
    ewline = strcompress(st1.line,/rem)
    ewcl = st1.wave_center
;     ewline = ['Hgamma','Hbeta','He','OI','Halpha'];,'Ca8498', 'Ca8662']
;     ewcl = [4340.47, 4861.34, 5875.618, 6300.3, 6562.81];, 8498., 8662.]
    ewjitter = 5.	;move +/-n pixels to et EW 
    ewnlines = n_elements(ewline)
    
    ;for velocity correlation and line profile variance
    vcorrline = strcompress(st3.line,/rem)
    vcorrcl = st3.wave_center
    
;     vcorrline = ['CaH', 'CaK', 'Hgamma','Hbeta','He', 'OI', 'Halpha'];,'Ca8498', 'Ca8662']
;     vcorrcl = [3968.47, 3933.66, 4340.47, 4861.34, 5875.618, 6300.3, 6562.81];, 8498., 8662.]
    vcorrnlines = n_elements(vcorrline)
    

    ;for LDR
    ;Catalano2002: '6199_6200','6211_6215','6216_6215','6243_6246','6252_6253','6266_6265','6275_6270'
    ;Gray&Brown, 2001PASP..113..723G: ,'6224_6225','6233_6233','6256_6257'	;giants
    ;Strassmeier&Schordan2000,AN,321: 6400-6460	;giants
    ; ldrline = ['6199_6200','6211_6215','6216_6215','6243_6246','6243_6247','6246_6247','6252_6253','6266_6265', '6269_6270', '6275_6270', $	;Catalano2002
    ; 				'6225_6224', '6233_6232', '6242_6244', '6252_6253', '6257_6255', $	;Gray&Brown2001
    ; 				'6411_6413','6413_6432','6413_6456','6416_6413','6419_6413','6430_6432','6435_6432','6435_6456','6439_6421','6452_6416','6452_6456','6455_6456']	;Strassmeier&Schordan2000
    ; ldrlinename = ['VI_FeI', 'ScI_FeI+TiI', 'VI_FeI+TiI', 'VI_FeI', 'VI_FeII', 'FeI_FeII', 'VI_FeI', 'VI_FeI', 'VI_FeI', 'VI_FeI', $
    ; 				'VI_NiI', 'VI_FeI', 'VI_SiI', 'VI_FeI', 'VI_FeI', $
    ; 				'FeI_ScI', 'ScI_FeII', 'ScI_FeII', 'FeI_ScI', 'FeI_ScI', 'FeI_FeII', 'YI_FeII', 'YI_FeII', 'CaI_FeI', 'SiI_FeI', 'SiI_FeII', 'CaI_FeII']
    ; ldrl1 = [6199.19,6210.67,6216.36,6243.11,6243.11,6246.33,6251.83,6266.33,6268.87,6274.66, $
    ; 				6224.51,6233.20,6242.84,6251.83,6256.89, $
    ; 				6411.647,6413.324,6413.324,6416.919,6419.942,6430.844,6435.004,6435.004,6439.075,6452.296,6452.296,6455.598]
    ; ldrl2 = [6200.32,6215.15,6215.22,6246.33,6247.56,6247.56,6252.57,6265.14,6270.23,6270.23, $
    ; 				6223.99,6232.65,6243.82,6252.57,6255.95, $
    ; 				6413.324,6432.680,6456.383,6413.324,6413.324,6432.680,6432.680,6456.383,6421.349,6416.919,6456.383,6456.383]

    ;Biazzo 2007
    
    ldrline = strcompress(st2.line_index,/rem)
    ldrl1 = st2.lambda_1
    ldrl2 = st2.lambda_2
    ldrlinename = strcompress(st2.line_pair,/rem)
    
;     ldrline = ['6199_6200','6210_6215','6213_6213','6216_6215','6224_6223','6233_6232','6242_6244','6243_6246','6243_6247','6251_6252','6257_6255','6257_6256','6266_6265','6268_6270','6275_6270']
;     ldrl1 = [6199.19,6210.67,6213.83,6216.36,6224.51,6233.20,6242.84,6243.11,6243.11,6251.83,6256.89,6256.89,6266.33,6268.87,6274.66]
;     ldrl2 = [6200.32,6215.15,6213.44,6215.15,6223.99,6232.65,6243.82,6246.33,6247.56,6252.57,6255.95,6256.35,6265.14,6270.23,6270.23]
;     ldrlinename = ['VI_FeI','ScI_FeI+TiI','VI_FeI','VI_FeI+TiI','VI_NiI','VI_FeI','VI_SiI','VI_FeI','VI_FeII','VI_FeI','VI_FeI','VI_NiI+FeI','VI_FeI','VI_FeI','VI_FeI']

    ldrnlines = n_elements(ldrline)

;     ;Photospheric Band Indices - PBI
;     ;Schoefer 2019, A&A, The Carmenes search...
;     pbiband = ['TiO7050', 'VO7436', 'VO7942', 'TiO8430']
;     numer = [[7056.0,7060.0],[7435.9,7436.9],[7941.7,7943.7],[8436.0,8440.0]]
;     denom = [[7046.0,7050.0],[7434.3,7435.3],[7936.0,7938.0],[8430.5,8434.5]]
;     npbi = n_elements(pbiband)

    
    c = 299792.458d0
    
    ;=====================================================================

    ;begin prepare
    
    ;prepare spectra for analysis
    ;merge orders and shift to rest wavelength

    if keyword_set(prepare) then begin

        file = file_search(base+id[iid]+'/ceres/*fits', count=nfiles)
        
        jd = dblarr(nfiles)
        snr = dblarr(nfiles)
        rv = dblarr(nfiles)
        rve = dblarr(nfiles)
        bs = dblarr(nfiles)
        bse = dblarr(nfiles)

        
;         resistant_mean, rv, 3, t1, t2, nbad, /double, goodvec=idxgood, badvec=idxbad
; 
;         print, 'Removed '+strcompress(nbad,/rem)+' bad values.'
;         window, 0
;         plot, jd, rv, psym=sym(1), /yn
;         oplot, jd[idxgood], rv[idxgood], psym=4, color=cgcolor('green')
;         if(idxbad[0] ne -1) then oplot, jd[idxbad], rv[idxbad], psym=1, color=cgcolor('red')
;         
;         ngood = n_elements(idxgood)
;         jd = dblarr(ngood)
;         snr = dblarr(ngood)
;         rv = dblarr(ngood)
;         rve = dblarr(ngood)
;         bs = dblarr(ngood)
;         bse = dblarr(ngood)
;         
;         file = file[idxgood]
;         filervo = filervo[idxgood]
;         file = file[idxgood]
;         nfiles = n_elements(idxgood)
;         nfiles = n_elements(idxgood)
        
        for xx=0,nfiles-1 do begin

            ;extract file name
            pos1 = strpos(file[xx], '/', /reverse_search)
            pos2 = strpos(file[xx], '_', /reverse_search)
            fn = strmid(file[xx], pos1+1, pos2-pos1-1)

            ;read data
            d = mrdfits(file[xx], 0, hdr, /silent)

            ;use RV from regular ceres reduction
            hdr = headfits(file[xx], exten=0, /silent)
            jd[xx] = double(get_eso_keyword(hdr, 'HIERARCH MBJD'))+2400000.5d0
            snr[xx] = get_eso_keyword(hdr, 'SNR')
            rv[xx] = get_eso_keyword(hdr, 'RV')
            rve[xx] = get_eso_keyword(hdr, 'RV_E')
            bs[xx] = get_eso_keyword(hdr, 'BS')
            bse[xx] = get_eso_keyword(hdr, 'BS_E')
            
            l = d[*,*,0]
            ;f = d[*,*,5]
            fe = d[*,*,4]
            n = d[*,*,7]
            f3 = d[*,*,3]

            ;correct orders for indivifual RV
            for i=0,n_elements(l[0,*])-1 do begin
            
                l[*,i] = toAir(l[*,i])
                logl = alog(l[*,i])
                logl = (c*logl - rv[xx]) / c
                l[*,i] = exp(logl)

            endfor
            
            ;f = f*n

            dx = abs(median(ts_diff(l[*],1), /even))
            minl = min(l)
            maxl = max(l)

            binl = linspace(minl, maxl, (maxl-minl)/dx)
            newl = dblarr(n_elements(binl)-1)
            for i=0L,n_elements(newl)-1 do newl[i] = binl[i+1]-(binl[i+1]-binl[i])/2.

            ;rebinning
            ;norder = n_elements(f[0,*])
            norder = n_elements(f3[0,*])
            ;stmp = dblarr(norder, n_elements(binl)-1)
            s3tmp = dblarr(norder, n_elements(binl)-1)
            setmp = dblarr(norder, n_elements(binl)-1)
            ctmp = dblarr(norder, n_elements(binl)-1)

            for i=0,norder-1 do begin

                s3tmp[i,*] = rebinw(reform(f3[*,i]), reform(l[*,i]), binl)
                ;stmp[i,*] = rebinw(reform(f[*,i]), reform(l[*,i]), binl)
                setmp[i,*] = rebinw(reform(fe[*,i]), reform(l[*,i]), binl)
                ctmp[i,*] = rebinw(reform(n[*,i]), reform(l[*,i]), binl)

            endfor

            nf3 = total(s3tmp,1)
            ;nf = total(stmp,1)
            nfe = total(setmp,1)
            nc = total(ctmp,1)

            ;merge = nf/nc
            merge = nf3/nc

            nl = newl
            nf = merge
            nfe = nfe

            idx = where(finite(nf) eq 1)
            nl = nl[idx]
            nf = nf[idx]
            nfe = nfe[idx]

            ;shift spectrum to rest wavelength

    ; 		nl = toAir(nl)
    ; 
    ; 		lognl = alog(nl)
    ; 		lognl = (c*lognl - rv[xx]) / c
    ; 		nl = exp(lognl)
            
            openw, lun, base+id[iid]+'/ceres/Spectrum_'+fn+'.txt', width=1400, /get_lun
                
                for i=0L,n_elements(nl)-1 do printf, lun, nl[i], nf[i], nfe[i], format='(f15.10,f30.17,f22.17)'
                        
            close, lun
            free_lun, lun

            proceeding_text, loop=nfiles, i=xx, prompt='> Spectrum        '+string(xx+1,form='(I4)')
                
        endfor
        
        openw, lun, base+id[iid]+'/ceres/'+'ceres_results.txt', width=1400, /get_lun
        
            for i=0,nfiles-1 do printf, lun, jd[i], rv[i], rve[i], bs[i], bse[i], snr[i], format='(f17.9, 4f10.4, f8.1)'
        
        close, lun
        free_lun, lun
        ;create average spectrum
        file = file_search(base+id[iid]+'/ceres/Spectrum_'+id[iid]+'*txt', count=nfiles)
        ;read first spectrum to get lambda referencce
        openr,lun,file[0],/get_lun
            rows = file_lines(file[0])
            data = dblarr(3, rows)  
            readf, lun, data
            refl = reform(data[0,*])
        close, lun
        free_lun, lun
        
        allf = dblarr(nfiles,n_elements(refl))
        allfe = dblarr(nfiles,n_elements(refl))
        for xx=0,nfiles-1 do begin
        
						openr,lun,file[xx],/get_lun
                rows = file_lines(file[xx])
                data = dblarr(3, rows)  
                readf, lun, data
                wl = reform(data[0,*])
                f = reform(data[1,*])
                fe = reform(data[2,*])
            close, lun
            free_lun, lun
            
            fint = interpol(f,wl,refl)
            finte = interpol(fe,wl,refl)

            allf[xx,*] = fint
            allfe[xx,*] = finte
            
        endfor
        
    ; 	meanf = mean(allf, dim=1)
    ; 	meanfe = mean(allfe, dim=1)
        meanf = median(allf, dim=1)	;less sensitive to outliers when RV wrong
        meanfe = median(allfe, dim=1)
        
        openw, lun, base+id[iid]+'/ceres/'+'mean_spectrum.txt', width=1400, /get_lun

            for i=0L,n_elements(refl)-1 do printf, lun, refl[i], meanf[i], meanfe[i], format='(f15.10,f30.17,f22.17)'
        
        close, lun
        free_lun, lun
        
        print, 'END PREPARE'
        
    endif

    ;end prepare

    ;=====================================================================

    readcol, base+id[iid]+'/ceres/'+'ceres_results.txt', jd, rv, rve, bs, bse, snr, format='d,d,d,d,d,d', /silent

    ;mean spectrum of all epochs
    openr,lun,base+id[iid]+'/ceres/'+'mean_spectrum.txt',/get_lun
        mrows = file_lines(base+id[iid]+'/ceres/'+'mean_spectrum.txt')
        data = dblarr(3, mrows)  
        readf, lun, data
        ml = reform(data[0,*])
        mf = reform(data[1,*])
        mfe = reform(data[2,*])
    close, lun
    free_lun, lun

    file = file_search(base+id[iid]+'/ceres/Spectrum_'+id[iid]+'*txt', count=nfiles)
    fileccf = file_search(base+id[iid]+'/ceres/'+id[iid]+'*CCF.txt', count=nfilesccf)

    if (nfiles ne nfilesccf) then begin
    
			print, 'Number of spectra different from CCF files. Stop.'
			stop    
    
    endif

    s_wilson = dblarr(nfiles)
    s_feros = dblarr(nfiles)
    ew = dblarr(nfiles,ewnlines)
    ewe = dblarr(nfiles,ewnlines)
    resflux = dblarr(nfiles,ewnlines)
    resfluxe = dblarr(nfiles,ewnlines)
    ldr  = dblarr(nfiles,ldrnlines)
    ldre = dblarr(nfiles,ldrnlines)
;     pbiindex = dblarr(nfiles,npbi)
;     pbiindexe = dblarr(nfiles,npbi)
		halpha_index = dblarr(nfiles)
		halpha_index_err = dblarr(nfiles)
		halpha10width = dblarr(nfiles)
		he1d3_index = dblarr(nfiles)
		he1d3_index_err = dblarr(nfiles)
		na1d1_index = dblarr(nfiles)
		na1d1_index_err = dblarr(nfiles)
		na1d2_index = dblarr(nfiles)
		na1d2_index_err = dblarr(nfiles)

    resf = dblarr(nfiles,mrows)
    intf = dblarr(nfiles,mrows)
    intfe = dblarr(nfiles,mrows)
    
    rows = file_lines(fileccf[0])
    bis_am = dblarr(nfiles,2,rows+50)	;to give some margin
    bvs_am = dblarr(nfiles)
    bvc_am = dblarr(nfiles)
    bvd_am = dblarr(nfiles)
    ccf_fwhm = dblarr(nfiles)
    ccf_fwhme = dblarr(nfiles)

    for xx=0,nfiles-1 do begin
    
    ;=====================================================================
    
    ;CCF analysis, bisector
    
			bot_i = 0.15
			bot_f = 0.4
			top_i = 0.6
			top_f = 0.9
			dt = 0.01
    
			openr,lun,fileccf[xx],/get_lun
				rows = file_lines(fileccf[xx])
				data = dblarr(2, rows)  
				readf, lun, data
				velccf = reform(data[0,*])
				fccf = reform(data[1,*])
			close, lun
			free_lun, lun
			
			;----------------------------
			
			idxl = where(velccf le rv[xx])
			idxr = where(velccf ge rv[xx])
			
			vleft = velccf[idxl]
			ccfleft = fccf[idxl]
			vright = velccf[idxr]
			ccfright = fccf[idxr]
			
			ivright = interpol(vright, ccfright, reverse(ccfleft), /spline)
			iccfright = reverse(ccfleft)
			
			bis = (ivright+reverse(vleft))/2.-rv[xx]
			bis_am[xx,0,0:n_elements(bis)-1] = bis
			bis_am[xx,1,0:n_elements(bis)-1] = iccfright
			
			;---------------------------
			
			pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},4)
			pi[3].fixed = 1
			estimates = [min(fccf-1.), rv[xx], 5., 1.]
            yfit = mpfitpeak(velccf, fccf, coef, nterms=4, estimates=estimates, /gaussian, parinfo=pi, perror=perror, chisq=chisq, bestnorm=bestnorm, dof=dof)
			;yfit = gaussfit(velccf, fccf, A, estimates=estimates,sigma=sigma,nterms=4,chisq=chisq)
			;coef = A
			
			ccf_fwhm[xx] = 2.*sqrt(2.*alog(2.))*coef[2]
			ccf_fwhme[xx] = perror[2]*sqrt(bestnorm/dof);sigma[2];2.*sqrt(2.*alog(2.))*perror[2]

			idx1 = where(velccf gt rv[xx]-3.*coef[2] and velccf lt rv[xx])
			idx2 = where(velccf lt rv[xx]+3.*coef[2] and velccf gt rv[xx])
			idx3 = where(velccf lt rv[xx]-4.*coef[2])
			idx4 = where(velccf gt rv[xx]+4.*coef[2])
			idx34 = [idx3,idx4]
			base34 = median(fccf[idx34])
			fccf = base34-fccf
			mx = max(fccf)
			fccf = fccf/max(fccf)
			
			v1 = velccf[idx1]
			x1 = fccf[idx1]
			v2 = velccf[idx2]
			x2 = fccf[idx2]
			
			iv2 = interpol(v2, x2, reverse(x1))
			ix2 = reverse(x1)

			bis = (iv2+reverse(v1))/2.-rv[xx]
			bis_am[xx,0,0:n_elements(bis)-1] = bis
			bis_am[xx,1,0:n_elements(bis)-1] = 1.-ix2
			idx1 = where(ix2 le top_f and ix2 ge top_i)
			idx2 = where(ix2 le top_i and ix2 ge bot_f)
			idx3 = where(ix2 ge bot_i and ix2 le bot_f)
			bvs_am[xx] = median(bis[idx3])-median(bis[idx1])	;BVS
			bvc_am[xx] = (median(bis[idx3])-median(bis[idx2]))-(median(bis[idx2])-median(bis[idx1]))	;BC
			bvd_am[xx] = (median(bis[idx3])+median(bis[idx2])+median(bis[idx1]))/3.	;BVD

; 			;ceres pipeline
; 			v2 = reverse(v2)
; 			x2 = reverse(x2)
; 			
; 			dp = top_f
; 			vect = dblarr(1)
; 			while dp ge top_i do begin
; 			
; 				lb = (where(x1 gt dp))[0]
; 				m = (v1[lb]-v1[lb-1])/(x1[lb]-x1[lb-1])
; 				n = v1[lb]-m*x1[lb]
; 				bs1 = m*dp+n
; 				
; 				lb = (where(x2 gt dp))[0]
; 				m = (v2[lb]-v2[lb+1])/(x2[lb]-x2[lb+1])
; 				n = v2[lb]-m*x2[lb]
; 				bs2 = m*dp+n
; 				vect = [vect,0.5*(bs2+bs1)]
; 				dp = dp-dt
; 				
; 			endwhile  
; 			vect = vect[1:*]
; 			
; 			dp = bot_f
; 			vecb = dblarr(1)
; 			while dp ge bot_i do begin
; 			
; 				lb = (where(x1 gt dp))[0]
; 				m = (v1[lb]-v1[lb-1])/(x1[lb]-x1[lb-1])
; 				n = v1[lb]-m*x1[lb]
; 				bs1 = m*dp+n
; 				
; 				lb = (where(x2 gt dp))[0]
; 				m = (v2[lb]-v2[lb+1])/(x2[lb]-x2[lb+1])
; 				n = v2[lb] - m*x2[lb]
; 				bs2 = m*dp+n
; 				vecb = [vecb,0.5*(bs2+bs1)]
; 				dp = dp-dt
; 				
; 			endwhile
; 			vecb = vecb[1:*]
; 			print, median(vecb)-median(vect)
; 			
; 			window, 0, xs=600, ys=800
; 			plot, velccf, fccf, /yn, psym=2;, xr=[25,35]
; 			oplot, velccf, yfit
; 			oplot, (rv[xx]-3.*coef[2])*[1,1], !y.crange
; 			oplot, (rv[xx]+3.*coef[2])*[1,1], !y.crange
; 			oplot, (rv[xx]-4.*coef[2])*[1,1], !y.crange
; 			oplot, (rv[xx]+4.*coef[2])*[1,1], !y.crange
; 			oplot, vleft, ccfleft, psym=1, color=cgcolor('blue')
; 			oplot, vright, ccfright, psym=1, color=cgcolor('green')
; 			oplot, ivright, iccfright, psym=sym(2), color=cgcolor('red')
; 			oplot, bis*80.+rv[xx], iccfright, psym=2, color=cgcolor('green')

    ;=====================================================================

        ; readcol, file[xx], wl, f, fe, format='d,d,d', /silent
        openr,lun,file[xx],/get_lun
            rows = file_lines(file[xx])
            data = dblarr(3, rows)  
            readf, lun, data
            wl = reform(data[0,*])
            f = reform(data[1,*])
            fe = reform(data[2,*])
        close, lun
        free_lun, lun

    ;=====================================================================

    ;Ca II H+K
    
        CaK = 3933.66
        CaH = 3968.47
        
        ;Duncan 1991, ApJS, 76, 383
        Rpass_win = [3991.07,4011.07]
        Vpass_win = [3891.07,3911.07]
        fwhm = 1.09
        
        ;triangular bandpass
        CaK_trif = dblarr(n_elements(wl))
        idx = where(wl ge (CaK-fwhm) and wl le CaK)
        tmpl = wl[idx]
        tmpf = tmpl*(1./1.09)
        tmpf = tmpf-max(tmpf)+1
        CaK_trif[idx] = tmpf
        idx = where(wl le (CaK+fwhm) and wl ge CaK)
        tmpl = wl[idx]
        tmpf = tmpl*(-1./1.09)
        tmpf = tmpf-max(tmpf)+1
        CaK_trif[idx] = tmpf
        CaH_trif = dblarr(n_elements(wl))
        idx = where(wl ge (CaH-fwhm) and wl le CaH)
        tmpl = wl[idx]
        tmpf = tmpl*(1./1.09)
        tmpf = tmpf-max(tmpf)+1
        CaH_trif[idx] = tmpf
        idx = where(wl le (CaH+fwhm) and wl ge CaH)
        tmpl = wl[idx]
        tmpf = tmpl*(-1./1.09)
        tmpf = tmpf-max(tmpf)+1
        CaH_trif[idx] = tmpf
        
        Rpass = dblarr(n_elements(wl))
        Rpass[where(wl ge Rpass_win[0] and wl le Rpass_win[1])] = 1.
        Vpass = dblarr(n_elements(wl))
        Vpass[where(wl ge Vpass_win[0] and wl le Vpass_win[1])] = 1.
        
        s_wilson[xx] = (total(f*CaH_trif) + total(f*CaK_trif))/(total(f*Rpass) + total(f*Vpass))
        
        ;Mohler2013, PhD thesis, p33
        fc = [3933.,3935.]
        fb = [3930.,3933.]
        fr = [3935.,3938.]

        s_feros[xx] = total(f[where(wl ge fc[0] and wl le fc[1])])/(total(f[where(wl ge fb[0] and wl le fb[1])]) +	total(f[where(wl ge fr[0] and wl le fr[1])]))
        
    ;=======================================================================

    ;EW and line residuals
        
        ;residual flux
        intf[xx,*] = interpol(f, wl, ml)
        intfe[xx,*] = interpol(fe, wl, ml)
        resf[xx,*] = intf[xx,*]-mf
        resfe = interpol(fe, wl, ml)
    
        for i=0,ewnlines-1 do begin	
        
            w0 = reform(ewcl[i])
            idxl = closest2(ewcl[i]-ewwin[i],wl)
            idxm = closest2(ewcl[i],wl)
            idxr = closest2(ewcl[i]+ewwin[i],wl)
            
    ; 		window, 0
    ; 		plot, wl[idxl-100:idxr+100], f[idxl-100:idxr+100], xst=1, title=ewline[i]
    ; 		oplot, ewcl[i]*[1,1], !y.crange, linestyle=2
    ; 		oplot, ewl1[i]*[1,1], !y.crange, color=cgcolor('green')
    ; 		oplot, ewl2[i]*[1,1], !y.crange, color=cgcolor('green')
            
        ;feature_MOD_2, wl, f, [wl[idxl],f[idxl]], [wl[idxm],f[idxm]], [wl[idxr],f[idxr]], w0
        ;ew[xx,i] = w0[12]
        ;;print, w0[12]
        ;tmpf = f+fe
        ;feature_MOD_2, wl, tmpf, [wl[idxl],tmpf[idxl]], [wl[idxm],tmpf[idxm]], [wl[idxr],tmpf[idxr]], w0
        ;ewe[xx,i] = abs(ew[xx,i]-w0[12])
        
    ;     ;idea is to randomly select points for EW measurement
    ;     tmpew = dblarr(100)
    ;     for j=0,n_elements(tmpew)-1 do begin
    ;     
    ; 			idx1 = round(-ewjitter+2.*ewjitter*randomu(seed,1))	;randomly select a value between -ewjitter and ewjitter
    ; 			idx2 = round(-ewjitter+2.*ewjitter*randomu(seed,1))
    ; 			feature_MOD_2, wl, f, [wl[idxl+idx1],f[idxl+idx1]], [wl[idxm],f[idxm]], [wl[idxr+idx2],f[idxr+idx2]], w0
    ; 			tmpew[j] = w0[12]
    ;     
    ;     endfor
    ;     ew[xx,i] = median(tmpew)
    ;     ewe[xx,i] = stddev(tmpew)

    ;     include random and then error  and try cont fixed
            tmpew = dblarr(1000)
            tmpewe = dblarr(1000)
            for j=0,n_elements(tmpew)-1 do begin
            
                idx1 = round(-ewjitter+2.*ewjitter*randomu(seed,1))	;randomly select a value between -ewjitter and ewjitter
                idx2 = round(-ewjitter+2.*ewjitter*randomu(seed,1))
            
                fcont = (wl-(wl[idxl+idx1])[0])*((f[idxr+idx2])[0]-(f[idxl+idx1])[0])/((wl[idxr+idx2])[0]-(wl[idxl+idx1])[0])+(f[idxl+idx1])[0]
                tmpew[j] = tsum(wl, (1.-f/fcont), idxl+idx1, idxr+idx2)
                ;https://arxiv.org/pdf/astro-ph/0606341.pdf
                tmpewe[j] = sqrt(1.+(mean(fcont[idxl:idxr])/mean(f[idxl:idxr])))*((wl[idxr]-wl[idxl])-tmpew[j])/SNR[xx]

            endfor
            
            ;weighted mean
            ew[xx,i] = (total((1./tmpewe^2.)*tmpew))/total(1./tmpewe^2.)	;weighted mean
            ewe[xx,i] = sqrt(n_elements(tmpew)/total(1./tmpewe^2.))
            
    ; 		print, 'Median'
    ; 		print, median(tmpew, /even)
    ; 		print, 'error fixed cont.'
    ; 		print, sqrt(1.+(1./mean(f[idxl:idxr])))*((wl[idxr]-wl[idxl])-ew[xx,i])/SNR[xx]	;fixed cont.
    ; 		
    ; 		print, 'bootstrap'
    ;       ;bootstrap mean
    ;       nboot = 10000.d0
    ;       boot_mean, tmpew, nboot, tmp
    ;       print, mean(tmp)
    ;       print, stddev(tmp)

    ; 		;classic single measurement
    ; 		fcont = (wl-(wl[idxl])[0])*((f[idxr])[0]-(f[idxl])[0])/((wl[idxr])[0]-(wl[idxl])[0])+(f[idxl])[0]
    ; 		ew[xx,i] = tsum(wl, (1.-f/fcont), idxl, idxr)
    ; 		;https://arxiv.org/pdf/astro-ph/0606341.pdf
    ; 		ewe[xx,i] = sqrt(1.+(mean(fcont[idxl:idxr])/mean(f[idxl:idxr])))*((wl[idxr]-wl[idxl])-ew[xx,i])/SNR[xx]

            resflux[xx,i] = tsum(ml, resf[xx,*], idxl, idxr)
            tmp1 = tsum(ml, resf+resfe, idxl, idxr)
            resfluxe[xx,i] = abs(resflux[xx,i]-tmp1)
        
    ;      print, resflux[xx,i], resfluxe[xx,i]
    ; 	   window, 0
    ; 	   plot, ml[idxl-100:idxr+100], resf[idxl-100:idxr+100], xst=1, title=ewline[i]
    ; 	   oplot, ewcl[i]*[1,1], !y.crange, linestyle=2
    ; 	   oplot, ewl1[i]*[1,1], !y.crange, color=cgcolor('green')
    ; 	   oplot, ewl2[i]*[1,1], !y.crange, color=cgcolor('green')


        endfor

    ;=======================================================================
    
    ;Halpha index
    ;Boisse+2009,A&A,495,959
    ;computation followin Zechmeister+2017
    
    l0 = 6562.808
    wl0 = 0.678
    win1 = [6545.495,6556.245]
    win2 = [6575.934,6584.684]
    
;     idxl = closest2(l0-wl0/2., wl)
;     idxr = closest2(l0+wl0/2., wl)
;     tmp1 = tsum(wl[idxl:idxr], f[idxl:idxr])
;     idxl1 = closest2(win1[0],wl)
;     idxr1 = closest2(win1[1],wl)
;     tmp2 = tsum(wl[idxl1:idxr1], f[idxl1:idxr1])
;     idxl2 = closest2(win2[0],wl)
;     idxr2 = closest2(win2[1],wl)
;     tmp3 = tsum(wl[idxl2:idxr2], f[idxl2:idxr2])
;     ;etmp1 = tsum(wl[idxl:idxr], fe[idxl:idxr])
;     ;etmp2 = tsum(wl[idxl1:idxr1], fe[idxl1:idxr1])
;     ;etmp3 = tsum(wl[idxl2:idxr2], fe[idxl2:idxr2])
;     
;     halpha_index[xx] = tmp1/(tmp2+tmp3)
;     ;halpha_index_err = etmp1/(etmp2+etmp3)
;     
;     tmp1e = abs(tmp1-tsum(wl[idxl:idxr], f[idxl:idxr]+fe[idxl:idxr]))
;     tmp2e = abs(tmp2-tsum(wl[idxl1:idxr1], f[idxl1:idxr1]+fe[idxl1:idxr1]))
; 		tmp3e = abs(tmp3-tsum(wl[idxl2:idxr2], f[idxl2:idxr2]+fe[idxl2:idxr2]))
; 	
; 		halpha_index_err[xx] = sqrt((1./(tmp2+tmp3))^2.*tmp1e^2. + (tmp1/(tmp2+tmp3)^2.)^2.*tmp2e^2. + (tmp1/(tmp2+tmp3)^2.)^2.*tmp3e^2.)

		;Zechmeister 2017, Serval
    idxl = closest2(l0-wl0/2., wl)
    idxr = closest2(l0+wl0/2., wl)
    tmp1 = mean(f[idxl:idxr])
    idxl1 = closest2(win1[0],wl)
    idxr1 = closest2(win1[1],wl)
    tmp2 = mean(f[idxl1:idxr1])
    idxl2 = closest2(win2[0],wl)
    idxr2 = closest2(win2[1],wl)
    tmp3 = mean(f[idxl2:idxr2])
    
    halpha_index[xx] = tmp1/(0.5*(tmp2+tmp3))
    
    tmp1e = sqrt(total(fe[idxl:idxr]^2.))/double(n_elements(idxr-idxl)+1.)
    tmp2e = sqrt(total(fe[idxl1:idxr1]^2.))/double(n_elements(idxr1-idxl1)+1.)
    tmp3e = sqrt(total(fe[idxl2:idxr2]^2.))/double(n_elements(idxr2-idxl2)+1.)
    
    halpha_index_err[xx] = halpha_index[xx]*sqrt(tmp1e^2./tmp1^2. + (tmp2e^2.+tmp3e^2.)/(tmp2^2.+tmp3^2.))
		
		;---------------------------	
    
    ;HeID3 index
    ;Boisse+2009,A&A,495,959
    ;using values of Gomes da Silva 2011, 534, A30
    l0 = 5875.62
    wl0 = 0.4
    win1 = [5866.5,5871.5]
    win2 = [5878.5,5883.5]
    
;     idxl = closest2(l0-wl0/2., wl)
;     idxr = closest2(l0+wl0/2., wl)
;     tmp1 = tsum(wl[idxl:idxr], f[idxl:idxr])
;     idxl1 = closest2(win1[0],wl)
;     idxr1 = closest2(win1[1],wl)
;     tmp2 = tsum(wl[idxl1:idxr1], f[idxl1:idxr1])
;     idxl2 = closest2(win2[0],wl)
;     idxr2 = closest2(win2[1],wl)
;     tmp3 = tsum(wl[idxl2:idxr2], f[idxl2:idxr2])
;     
;     he1d3_index[xx] = tmp1/(tmp2+tmp3)
    
;     tmp1e = abs(tmp1-tsum(wl[idxl:idxr], f[idxl:idxr]+fe[idxl:idxr]))
;     tmp2e = abs(tmp2-tsum(wl[idxl1:idxr1], f[idxl1:idxr1]+fe[idxl1:idxr1]))
; 		tmp3e = abs(tmp3-tsum(wl[idxl2:idxr2], f[idxl2:idxr2]+fe[idxl2:idxr2]))
	
; 		he1d3_index_err[xx] = sqrt((1./(tmp2+tmp3))^2.*tmp1e^2. + (tmp1/(tmp2+tmp3)^2.)^2.*tmp2e^2. + (tmp1/(tmp2+tmp3)^2.)^2.*tmp3e^2.)

    idxl = closest2(l0-wl0/2., wl)
    idxr = closest2(l0+wl0/2., wl)
    tmp1 = mean(f[idxl:idxr])
    idxl1 = closest2(win1[0],wl)
    idxr1 = closest2(win1[1],wl)
    tmp2 = mean(f[idxl1:idxr1])
    idxl2 = closest2(win2[0],wl)
    idxr2 = closest2(win2[1],wl)
    tmp3 = mean(f[idxl2:idxr2])
    
    he1d3_index[xx] = tmp1/(0.5*(tmp2+tmp3))
    
    tmp1e = sqrt(total(fe[idxl:idxr]^2.))/double(n_elements(idxr-idxl)+1.)
    tmp2e = sqrt(total(fe[idxl1:idxr1]^2.))/double(n_elements(idxr1-idxl1)+1.)
    tmp3e = sqrt(total(fe[idxl2:idxr2]^2.))/double(n_elements(idxr2-idxl2)+1.)
    
    he1d3_index_err[xx] = he1d3_index[xx]*sqrt(tmp1e^2./tmp1^2. + (tmp2e^2.+tmp3e^2.)/(tmp2^2.+tmp3^2.))
    
		;---------------------------	
    
    ;NaID1 index
    
    ;using values of Gomes da Silva 2011, 534, A30
    l0 = 5895.92
    wl0 = 0.5
    win1 = [5800.,5810.]
    win2 = [6080.,6100.]
    
    idxl = closest2(l0-wl0/2., wl)
    idxr = closest2(l0+wl0/2., wl)
    tmp1 = mean(f[idxl:idxr])
    idxl1 = closest2(win1[0],wl)
    idxr1 = closest2(win1[1],wl)
    tmp2 = mean(f[idxl1:idxr1])
    idxl2 = closest2(win2[0],wl)
    idxr2 = closest2(win2[1],wl)
    tmp3 = mean(f[idxl2:idxr2])
    
    na1d1_index[xx] = tmp1/(0.5*(tmp2+tmp3))
    
    tmp1e = sqrt(total(fe[idxl:idxr]^2.))/double(n_elements(idxr-idxl)+1.)
    tmp2e = sqrt(total(fe[idxl1:idxr1]^2.))/double(n_elements(idxr1-idxl1)+1.)
    tmp3e = sqrt(total(fe[idxl2:idxr2]^2.))/double(n_elements(idxr2-idxl2)+1.)
    
    na1d1_index_err[xx] = na1d1_index[xx]*sqrt(tmp1e^2./tmp1^2. + (tmp2e^2.+tmp3e^2.)/(tmp2^2.+tmp3^2.))
    
		;---------------------------	

    ;NaID2 index
    
    ;using values of Gomes da Silva 2011, 534, A30
    l0 = 5889.95
    wl0 = 0.5
    win1 = [5800.,5810.]
    win2 = [6080.,6100.]
    
    idxl = closest2(l0-wl0/2., wl)
    idxr = closest2(l0+wl0/2., wl)
    tmp1 = mean(f[idxl:idxr])
    idxl1 = closest2(win1[0],wl)
    idxr1 = closest2(win1[1],wl)
    tmp2 = mean(f[idxl1:idxr1])
    idxl2 = closest2(win2[0],wl)
    idxr2 = closest2(win2[1],wl)
    tmp3 = mean(f[idxl2:idxr2])
    
    na1d2_index[xx] = tmp1/(0.5*(tmp2+tmp3))
    
    tmp1e = sqrt(total(fe[idxl:idxr]^2.))/double(n_elements(idxr-idxl)+1.)
    tmp2e = sqrt(total(fe[idxl1:idxr1]^2.))/double(n_elements(idxr1-idxl1)+1.)
    tmp3e = sqrt(total(fe[idxl2:idxr2]^2.))/double(n_elements(idxr2-idxl2)+1.)
    
    na1d2_index_err[xx] = na1d2_index[xx]*sqrt(tmp1e^2./tmp1^2. + (tmp2e^2.+tmp3e^2.)/(tmp2^2.+tmp3^2.)) 

    ;=======================================================================
    
    if keyword_set(twhya) then begin
    
			;10% Halpha width
			idx = where(wl ge 6550. and wl le 6574.)
			tmpl = wl[idx]
			tmpf = f[idx]
			tmpf = tmpf-1.
			tmp = max(f,idx)
			;maxf = median(f[idx-5:idx+5],/even)
			maxf = mean(f[idx-10:idx+10])

			vel_halpha10width = c*(tmpl-6562.81)/tmpl

			idxl = where(vel_halpha10width le 0.)
			idx = closest2(0.1*maxf, tmpf[idxl])
			v1 = (vel_halpha10width[idxl[idx]])[0]
			idxr = where(vel_halpha10width ge 0.)
			idx = closest2(0.1*maxf, tmpf[idxr])
			v2 = (vel_halpha10width[idxr[idx]])[0]
			halpha10width[xx] = v2-v1
			
		endif else begin
		
			halpha10width = 0.
			
		endelse

    ;=======================================================================
    
    ;LDR

        for i=0,ldrnlines-1 do begin

    ; 		window, 0
    ; 		if (ldrl1[i] lt ldrl2[i]) then xr=[ldrl1[i]-3,ldrl2[i]+3] else xr=[ldrl2[i]-3,ldrl1[i]+3]
    ; 		plot, wl[*], f[*], xst=1, xr=xr, yr=[0.4,1.05], yst=1
    ; 		oplot, ldrl1[i]*[1,1], !y.crange, color=rgb(255,0,0), linestyle=1
    ; 		oplot, ldrl2[i]*[1,1], !y.crange, color=rgb(255,0,0), linestyle=1

            ;measure LDR
            idx1 = closest2(ldrl1[i], wl)
            tl1 = wl[idx1-2:idx1+2]
            tf1 = f[idx1-2:idx1+2]
            tfe1 = fe[idx1-2:idx1+2]
            d1 = 1.-min(tf1,idxmin)	;only to get d1e
            d1e = tfe1[idxmin]

    ; 		startval = [3d8,-100000.,8.]
    ; 		res = mpfitfun('poly_2deg', tl,tf,tfe, startval,bestnorm=bestnorm,dof=dof,perror=perror,/quiet)
            res1 = robust_poly_fit(tl1,tf1,2,yfit1,sig1)
            A1 = res1
    ; 		dumweights = dblarr(n_elements(tl1))+1.
    ; 		yfit = curvefit(tl1, tf1, dumweights, A1, sigma1, function_name='poly_2deg_v2', yerror=ye1);/noderivative
            ye1 = sqrt(total((yfit1-tf1)^2 )/2.)	;yerror from curvefit
            d1 = 1.-(A1[0]-A1[1]^2./(4.*A1[2]))

            idx2 = closest2(ldrl2[i], wl)
            tl2 = wl[idx1-2:idx1+2]
            tf2 = f[idx2-2:idx2+2]
            tfe2 = fe[idx2-2:idx2+2]
            d2 = 1.-min(tf2,idxmin)	;only to get d2e
            d2e = tfe2[idxmin]
            res2 = robust_poly_fit(tl2,tf2,2,yfit2,sig2)
            A2 = res2
    ; 		dumweights = dblarr(n_elements(tl2))+1.
    ; 		yfit = curvefit(tl2, tf2, dumweights, A2, sigma2, function_name='poly_2deg_v2', yerror=ye2);/noderivative
            ye2 = sqrt(total((yfit2-tf2)^2 )/2.)
            ;startval = res2
            ;resm2 = mpfitfun('poly_2deg', tl2,tf2,tfe2, startval,bestnorm=bestnorm,dof=dof,perror=perror);,/quiet)
            d2 = 1.-(A2[0]-A2[1]^2./(4.*A2[2]))
    ; 		window, 0, xs=800, ys=800
    ; 		!p.multi=[0,1,2]
    ; 		ploterror, tl1, tf1, tfe1, /yn
    ; 		plotwl1 = linspace(min(tl1),max(tl1),100)
    ; 		oplot, plotwl1, poly_2deg(plotwl1,res1), color=rgb(255,0,0)
    ; 		ploterror, tl2, tf2, tfe2, /yn
    ; 		plotwl2 = linspace(min(tl2),max(tl2),100)
    ; 		oplot, plotwl2, poly_2deg(plotwl2,res2), color=rgb(255,0,0)
    ; 		oplot, plotwl2, poly_2deg(plotwl2,A2), color=rgb(0,255,0)
    ; 		;oplot, plotwl2, poly_2deg(plotwl2,resm2), color=rgb(0,0,255)
    ; 		!p.multi=[0,1,0]
    ; 		stop
      
            ldr[xx,i] = d1/d2
            ldre[xx,i] = sqrt(ye1^2.+ye2^2.)	;ldr[xx,i]*sqrt((d1e/d1)^2. + (d2e/d2)^2.)
            ;print ,ldre[xx,i], sqrt(ye1^2.+ye2^2.)
            ;I don't know for the moment how to get the error right. From the feeling using sdev from the fits gives a good order of magnitude. Using the formula in Catalano or Biazzo give sunreasonable large errors. Using the uncertainty from the fit gives unreasonable large errors.

            ; sp1 = (1.-d1)
            ; e1 = sp1*sqrt((1./sp1)+1.)
            ; sp2 = (1.-d2)
            ; e2 = sp2*sqrt((1./sp2)+1.)
            ; err = ldr[xx,i]*sqrt((e1/d1)^2.+(e2/d2)^2.)
            ; err = sqrt(sigma2[0]^2. + (-2.*A2[1]/(4.*A2[2]))^2.*sigma2[1]^2. + (A2[1]^2./(4.*A2[2]^2.))*sigma2[2]^2.)
            ; 
            ; ran = randomu(seed,100000)
            ; t0 = res2[0]+(-sigma2[0]+2.*sigma2[0]*ran)
            ; t1 = res2[1];+(-sigma2[1]+2.*sigma2[1]*ran)
            ; t2 = res2[2]+(-sigma2[2]+2.*sigma2[2]*ran)
            ; x = 1.-(t0-t1^2./(4.*t2))
            ; 		stop

        endfor

    ;=======================================================================

;     ;Photospheric Band Indices - PBI
;     
;         for i=0,npbi-1 do begin
;     
;             idx1 = where(wl ge numer[0,i] and wl le numer[1,i])
;             idx2 = where(wl ge denom[0,i] and wl le denom[1,i])
;     
; 						if (idx1[0] ne -1 and idx2[0] ne -1) then begin
;     
; 							pbiindex[xx,i] = mean(f[idx1])/mean(f(idx2))
; 							
; 							numere = (1./n_elements(idx1))*sqrt(total(fe[idx1]^2.))
; 							denome = (1./n_elements(idx2))*sqrt(total(fe[idx2]^2.))
; 							pbiindexe[xx,i] = pbiindex[xx,i]*sqrt((numere/mean(f[idx1]))^2.+(denome/mean(f[idx2]))^2.)
;     
; 						endif
;     
;         endfor
        

    ;=======================================================================
    
        proceeding_text, loop=nfiles, i=xx, prompt='> Spectrum        '+string(xx+1,form='(I4)')
        
    endfor

    ;=======================================================================

    ;velocity correlation
    for xx=0,vcorrnlines-1 do begin
    
        twl = ml - vcorrcl[xx]
        vel = c*twl/wl
        
        idx = where(vel ge -winvcorr[xx] and vel le winvcorr[xx])
        vel = vel[idx]
        tf = intf[*,idx]
        
        corr = dblarr(n_elements(idx), n_elements(idx))
        for i=0L,n_elements(idx)-1 do begin
        
            for j=n_elements(idx)-1,0,-1 do begin

                corr[i,j] = (correlate(tf[*,i], tf[*,j], /double))

            endfor
            
        endfor
        
        save, corr, vel, filename=resdirvcorr+'Values_VelocityCorrelation_'+vcorrline[xx]+'.sav'

        set_plot, 'ps'
        device, filename=resdirvcorr+'VelocityCorrelation_'+vcorrline[xx]+'.ps',/color, xsize=18, ysize=15

        !P.Font=1
        !p.thick=4
        !x.thick=3
        !y.thick=3
        !p.multi=[0,1,0]

        nlevels = 20
        loadct, 33, NColors=nlevels, Bottom=1, /silent
        maxcorr = 1.d0
        mincorr = -1.d0
        step = (maxcorr-mincorr) / nlevels
        levels = IndGen(nlevels) * step+mincorr
        SetDecomposedState, 0, CurrentState=currentState

        contour, corr, vel, vel, /fill, levels=levels, color=cgcolor('black'), C_Colors=IndGen(nlevels)+1, $
        title=vcorrline[xx]+'!C!C', xtitle=ldelta+'V [km/s]', ytitle=ldelta+'V [km/s]', $
    ; 	Position=[0.145, 0.14, 0.95, 0.80], $
        Position=[0.17, 0.165, 0.95, 0.75], $
        xr=[min(vel),max(vel)], yr=[min(vel),max(vel)], xst=9, yst=9, $
        charsize=2., xticklen=-0.04, yticklen=-0.03;xminor=2, yminor=2

        SetDecomposedState, currentState
        tn = ['-1.0', ' ', '-0.8', ' ', '-0.6', ' ', '-0.4', ' ', '-0.2', ' ', '0.0', ' ', '0.2', ' ', '0.4', ' ', '0.6', ' ', '0.8', ' ', '1.0']
        cgColorbar, Range=[mincorr, maxcorr], Divisions=20, XTicklen=1, XMinor=0, $
        AnnotateColor='black', ncolors=20, bottom=1, Charsize=1, $
      ;Position=[0.145, 0.915, 0.95, 0.95], $
        Position=[0.17, 0.85, 0.95, 0.87], $
        title='Correlation', ticknames=tn;, xtickformat="(A1)";, , format="(f4.1)"


        !p.multi=[0,1,0]
        !P.Font=0
        !p.thick=1
        !x.thick=1
        !y.thick=1
        
        device,/close
        set_plot,'x'
        spawn, epstopdf+resdirvcorr+'VelocityCorrelation_'+vcorrline[xx]+'.ps'
        spawn, 'rm '+resdirvcorr+'VelocityCorrelation_'+vcorrline[xx]+'.ps'

        proceeding_text, loop=vcorrnlines, i=xx, prompt='> Vcorr Line        '+string(xx+1,form='(I4)')

    endfor
    
    ;=======================================================================

    ;line profile variance
    
    ;Johns, Christopher M., Basri, Gibor, 1995, AJ, 109, 2800 Hamilton Echelle Spectra of Young Stars. II. Time Series Analysis of H(alpha) Variations

    for xx=0,vcorrnlines-1 do begin
    
        twl = ml - vcorrcl[xx]
        vel = c*twl/wl
        
        idx = where(vel ge -winvcorr[xx] and vel le winvcorr[xx])
        vel = vel[idx]
        tf = intf[*,idx]
    
        varprofile = dblarr(n_elements(idx))
        for i=0,n_elements(idx)-1 do varprofile[i] = sqrt(total((tf[*,i] - mf[idx[i]])^2.,/nan) / (double(n_elements(idx)-1.d0)))
        
        varprofile = varprofile/mf[idx]
        
        save, varprofile, vel, filename=resdirlpv+'Values_LineProfileVariance_'+vcorrline[xx]+'.sav'
        
        set_plot, 'ps'
        device, filename=resdirlpv+'LineProfileVariance_'+vcorrline[xx]+'.ps',/color, xsize=18, ysize=15

        !P.Font=1
        !p.thick=4
        !x.thick=3
        !y.thick=3
        !p.multi=[0,1,0]

        maxval = max(mf[idx])
        if (maxval le 1.05) then yr = [0,1.05] else yr=[0,1.3*maxval]
        cgplot, /nodata, vel, mf[idx], color=cgcolor('black'), background=cgcolor('white'), xtitle=ldelta+'V [km/s]', ytitle='Flux Average Line Profile', title=vcorrline[xx], xminor=2, yminor=2, xst=1, yst=9, pos=[0.125,0.145,0.85,0.9], yr=yr, charsize=2;, xr=xr[*,xx]
        cgplot, vel, mf[idx], color=cgcolor('black'), /overplot
        plots, [0,0], [!y.crange], linestyle=1
        plots, [!x.crange], [1,1], linestyle=1

        maxval = max(varprofile, min=minval)
        cgplot, /noerase, /nodata, vel, varprofile , xst=5, yst=5, yr=[0,1.5*maxval], pos=[0.125,0.145,0.85,0.9], charsize=2;, xr=xr[*,xx]
        cgplot, vel, varprofile, color=cgcolor('blue'), /overplot, psym=sym(1), symsize=0.5
    ;polyfill, refvf, varprofile, color=cgcolor('dark gray')
        cgaxis, xaxis=0, yaxis=1, yminor=2, ytitle='Normalized Variance Profile', color=cgcolor('blue')

        !p.multi=[0,1,0]
        !P.Font=0
        !p.thick=1
        !x.thick=1
        !y.thick=1
        
        device,/close
        set_plot,'x'
        spawn, epstopdf+resdirlpv+'LineProfileVariance_'+vcorrline[xx]+'.ps'
        spawn, 'rm '+resdirlpv+'LineProfileVariance_'+vcorrline[xx]+'.ps'

        proceeding_text, loop=vcorrnlines, i=xx, prompt='> Line Profile Variance        '+string(xx+1,form='(I4)')

    endfor
    
    ;standard deviation of fractional changes
    
    for xx=0,vcorrnlines-1 do begin
    
        twl = ml - vcorrcl[xx]
        vel = c*twl/wl
        
        idx = where(vel ge -winvcorr[xx] and vel le winvcorr[xx])
        vel = vel[idx]
        tf = intf[*,idx]
        
        ;every possible combination of epochs
        combi = combigen(nfiles,2)
        ncombi = n_elements(combi[*,0])
        
        dt = dblarr(ncombi)
        d_lambda = dblarr(ncombi,n_elements(idx))
        
        for i=0,ncombi-1 do begin
        
            d_lambda[i,*] = reform((tf[combi[i,1],*]-tf[combi[i,0],*])/(0.5*(tf[combi[i,1],*]+tf[combi[i,0],*])))
            dt[i] = round(jd[combi[i,1]]-jd[combi[i,0]])

        endfor

        idxs = sort(dt)
        dt = dt[idxs]
        d_lambda = d_lambda[idxs,*]
    
        idxu = uniq(dt) ;uniq time bins
        lag = round(dt[idxu])
        nu = n_elements(idxu)
        mean_sdev_d_lambda = dblarr(nu)
        for i=0,nu-1 do begin
        
            idx = where(dt eq dt[idxu[i]])
            d_lambda_u = d_lambda[idx,*]
            mean_sdev_d_lambda[i] = mean(stddev(d_lambda_u, dim=1, /double))

        endfor
        
        ;if there is only 1 epoch in the bin then sdev = nan
        idxnan = where(finite(mean_sdev_d_lambda) eq 1)
        lag = lag[idxnan]
        mean_sdev_d_lambda = mean_sdev_d_lambda[idxnan]
        
        save, lag, mean_sdev_d_lambda, filename=resdirlpv+'Values_SdevFractionalChange_'+vcorrline[xx]+'.sav'
        
        set_plot, 'ps'
        device, filename=resdirlpv+'SdevFractionalChange_'+vcorrline[xx]+'.ps',/color, xsize=18, ysize=15

        !P.Font=1
        !p.thick=4
        !x.thick=3
        !y.thick=3
        !p.multi=[0,1,0]
        
        plot, lag, mean_sdev_d_lambda, xr=[min(lag)-1,max(lag)+1], thick=3, xtitle='Lag [days]', ytitle='Standard Deviation of Fractional Changes', title=vcorrline[xx], charsize=2, xst=1, pos=[0.16,0.145,0.98,0.9]
        oplot, lag, mean_sdev_d_lambda, psym=sym(1)

        !p.multi=[0,1,0]
        !P.Font=0
        !p.thick=1
        !x.thick=1
        !y.thick=1
        
        device,/close
        set_plot,'x'
        spawn, epstopdf+resdirlpv+'SdevFractionalChange_'+vcorrline[xx]+'.ps'
        spawn, 'rm '+resdirlpv+'SdevFractionalChange_'+vcorrline[xx]+'.ps'

        proceeding_text, loop=vcorrnlines, i=xx, prompt='> Sdev Fractional Changes        '+string(xx+1,form='(I4)')
        
    endfor
    
    ;2D periodograms of lines
    
    for xx=0,vcorrnlines-1 do begin
    
        twl = ml - vcorrcl[xx]
        vel = c*twl/wl
        
        idx = where(vel ge -winvcorr[xx] and vel le winvcorr[xx])
        vel = vel[idx]
        tf = intf[*,idx]
        tfe = intfe[*,idx]
        
        nf = n_elements(idx)
        
        ;get the number of GLS frequencies
        tbase = max(jd)-min(jd)
        step = 1./tbase/ofac
        fsteps = long((f1-f0)/step+1)
        power = dblarr(nf,fsteps)


        for j=0,nf-1 do begin
        
            fn = 'tmp.txt'
            openw, lun, pathGLS+fn, width=1400, /get_lun
                for i=0,n_elements(jd)-1 do printf, lun, jd[i], tf[i,j], tfe[i,j], format='(f17.9,f23.15,f23.15)'
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
;             spawn, pathGLS+'./GLS '+pathGLS+fn+' '+strcompress(f0, /rem)+' '+strcompress(f1, /rem)+' '+strcompress(mass, /rem)

            readcol, 'GLS.plt', freq, pow, win, pls, format='d,d,d,d', /silent
            power[j,*] = pow

            spawn, 'rm GLS.log GLS.par RVSinRes.plt GLSFit.plt GLS.plt'        

        endfor
        power = transpose(power)
        
        save, power, freq, vel, filename=resdirlpv+'Values_2DGLS_'+vcorrline[xx]+'.sav'
        
        set_plot, 'ps'
        device, filename=resdirlpv+'2DGLS_'+vcorrline[xx]+'.ps',/color, xsize=18, ysize=15

        !P.Font=1
        !p.thick=4
        !x.thick=3
        !y.thick=3
        !p.multi=[0,1,0]

        nlevels = 20
        ;nlevels = 24
        loadct, 33, NColors=nlevels, Bottom=1, /silent
        maxcorr = 1.d0
        mincorr = 0.d0
        step = (maxcorr-mincorr) / nlevels
        levels = IndGen(nlevels) * step+mincorr
        SetDecomposedState, 0, CurrentState=currentState

        contour, power, 1./freq, vel, /fill, levels=levels, color=cgcolor('black'), C_Colors=IndGen(nlevels)+1, $
        title=vcorrline[xx]+'!C!C', xtitle='Period [days]', ytitle=ldelta+'V [km/s]', $
        Position=[0.17, 0.165, 0.95, 0.75], $
        xr=[1./f1,1./f0], yr=[min(vel),max(vel)], xst=9, yst=9, $
        charsize=2., /xlog, xticklen=-0.04, yticklen=-0.03;, xticklen=0.02, yticklen=0.02, xminor=2, yminor=2, 

        SetDecomposedState, currentState
        tn = ['0.0', ' ', '0.1', ' ', '0.2', ' ', '0.3', ' ', '0.4', ' ', '0.5', ' ', '0.6', ' ', '0.7', ' ', '0.8', ' ', '0.9', ' ', '1.0']
        cgColorbar, Range=[mincorr, maxcorr], Divisions=20, XTicklen=1, XMinor=0, $
        AnnotateColor='black', ncolors=20, bottom=1, Charsize=1, $
        Position=[0.17, 0.85, 0.95, 0.87], $
        title='GLS Power', ticknames=tn;, xtickformat="(A1)";, , format="(f4.1)"
        
        !p.multi=[0,1,0]
        !P.Font=0
        !p.thick=1
        !x.thick=1
        !y.thick=1
        
        device,/close
        set_plot,'x'
        spawn, epstopdf+resdirlpv+'2DGLS_'+vcorrline[xx]+'.ps'
        spawn, 'rm '+resdirlpv+'2DGLS_'+vcorrline[xx]+'.ps'        

        proceeding_text, loop=vcorrnlines, i=xx, prompt='> 2D Periodograms Lines        '+string(xx+1,form='(I4)')

    endfor
    
		;=======================================================================
		
		;shorten BIS array
		nval = dblarr(nfiles)
		for i=0,nfiles-1 do begin
		
			idx = where(bis_am[i,0,*] ne 0)
			nval[i] = n_elements(idx)
		
		endfor
		mx = max(nval)
		bis_am = bis_am[*,*,0:mx-1]
		
    !p.font=1
    !p.thick=4
    !x.thick=3
    !y.thick=3

    ;BIS
    
    xr = dblarr(nfiles,2)
    for i=0,nfiles-1 do begin
    
			idx = where(bis_am[i,0,*] ne 0. and bis_am[i,1,*] le top_f and bis_am[i,1,*] ge bot_i)
			xr[i,*] = [min(bis_am[i,0,idx]), max(bis_am[i,0,idx])]
			
		endfor
    
    
    print, 'Plot BIS'
    set_plot, 'ps'
    device, filename=resdirrv+'BIS.ps',/color, xsize=20, ysize=14
; 				idx = where(bis_am[0,0,*] ne 0. and bis_am[0,1,*] le top_f and bis_am[0,1,*] ge bot_i)
        plot, bis_am[0,0,*], bis_am[0,1,*], xtitle='Velocity [km/s]', ytitle='Normalized Intensity', charsize=2, position=[0.15,0.16,0.95,0.9], xr=[min(xr[*,0])+0.2*min(xr[*,0]),max(xr[*,1])+0.2*max(xr[*,1])], yr=[0,1], /nodata, xst=1, yst=1
        for i=0,nfiles-1 do begin
					
					idx = where(bis_am[i,0,*] ne 0. and bis_am[i,1,*] le top_f and bis_am[i,1,*] ge bot_i)
					oplot, bis_am[i,0,idx], bis_am[i,1,idx]
					
				endfor
				
    device,/close
    set_plot,'x'
    spawn, epstopdf+resdirrv+'BIS.ps'
    spawn, 'rm '+resdirrv+'BIS.ps'
    
    print, 'Plot BVSCD'
    set_plot, 'ps'
    device, filename=resdirrv+'BVSCD.ps',/color, xsize=20, ysize=30
    !p.multi=[0,1,3]

			plot, rv, bvs_am, xtitle='RV [km/s]', ytitle='Bisector Span [km/s]', charsize=3.5, /yn, psym=sym(1)
			plot, rv, bvc_am, xtitle='RV [km/s]', ytitle='Bisector Curvature [km/s]', charsize=3.5, /yn, psym=sym(1)
			plot, rv, bvd_am, xtitle='RV [km/s]', ytitle='Bisector Displacement [km/s]', charsize=3.5, /yn, psym=sym(1)
				
    device,/close
    set_plot,'x'
    spawn, epstopdf+resdirrv+'BVSCD.ps'
    spawn, 'rm '+resdirrv+'BVSCD.ps'
    !p.multi=[0,1,0]
    
    
    print, 'Plot FWHM vs RV'
    coeff = robust_linefit(rv, ccf_fwhm, YFIT, SIG, COEF_SIG)
    set_plot, 'ps'
    device, filename=resdirrv+'RV_FWHM.ps',/color, xsize=20, ysize=14

        ploterror, rv, ccf_fwhm, rve, ccf_fwhme, xtitle='Radial Velocity [km/s]', ytitle='FWHM [km/s]', charsize=2, position=[0.15,0.16,0.95,0.9], psym=sym(1), /yn, /nohat
;         ploterror, rv, ccf_fwhm, rve, ccfleft, /nohat, xtitle='Velocity [km/s]', ytitle='FWHM [km/s]', charsize=2, position=[0.15,0.16,0.95,0.9], psym=sym(1), /yn
        oplot, rv, coeff[1]*rv+coeff[0], color=cgcolor('red')

    device,/close
    set_plot,'x'
    spawn, epstopdf+resdirrv+'RV_FWHM.ps'
    spawn, 'rm '+resdirrv+'RV_FWHM.ps'


    ;=======================================================================

    ;save data
    save, ml, mf, intf, intfe, filename=resdir+'Spectra.sav'
    save, jd, rv, rve, bs, bse, bis_am, bvs_am, bvd_am, bvc_am, ccf_fwhm, ccf_fwhme, filename=resdirrv+'Values_RV.sav'
    save, jd, ewline, ewcl, ewwin, resflux, resfluxe, filename=resdirrf+'Values_RF.sav'
    save, jd, s_wilson, s_feros, filename=resdirew+'Values_Sindex.sav'
;     save, jd, pbiband, numer, denom, pbiindex, pbiindexe, filename=resdirpbi+'Values_PBI.sav'
    save, jd, ewline, ewcl, ewwin, ew, ewe, filename=resdirew+'Values_EW.sav'
    save, jd, ldrline, ldrlinename, ldrl1, ldrl2, ldr, ldre, filename=resdirldr+'Values_LDR.sav'
    save, jd, halpha_index, halpha_index_err, he1d3_index, he1d3_index_err, na1d1_index, na1d1_index_err, na1d2_index, na1d2_index_err, halpha10width, filename=resdirindex+'Values_Index.sav'
    

    ;=======================================================================

;     openw, lun, resdirew+'S_Wilson.txt', width=1400, /get_lun
;         for i=0,nfiles-1 do printf, lun, jd[i], s_wilson[i], format='(f17.9,f20.15)'
;     close, lun
;     free_lun, lun
; 
;     openw, lun, resdirew+'S_Feros.txt', width=1400, /get_lun
;         for i=0,nfiles-1 do printf, lun, jd[i], s_feros[i], format='(f17.9,f20.15)'
;     close, lun
;     free_lun, lun

    ;=======================================================================

    ;periodograms
    
	RVSPY_compute_periodograms, resdirindex, 'HalphaIndex', jd, halpha_index, halpha_index_err, f0, f1, ofac, weight, /GLS
	RVSPY_compute_periodograms, resdirindex, 'HeID3Index', jd, he1d3_index, he1d3_index_err, f0, f1, ofac, weight, /GLS
	RVSPY_compute_periodograms, resdirindex, 'NaID1Index', jd, na1d1_index, na1d1_index_err, f0, f1, ofac, weight, /GLS
	RVSPY_compute_periodograms, resdirindex, 'NaID1Index', jd, na1d2_index, na1d2_index_err, f0, f1, ofac, weight, /GLS
	if keyword_set(twhya) then begin
	
		dummyerr = halpha10width
		dummyerr[*] = 1.
		RVSPY_compute_periodograms, resdirindex, 'Halpha10width', jd, halpha10width, dummyerr, f0, f1, ofac, weight, /GLS
		
	endif


    dummyerr = dblarr(nfiles)
    dummyerr[*] = 1.

    RVSPY_compute_periodograms, resdirrv, 'RV', jd, rv, rve, f0, f1, ofac, weight, /GLS;, /CLEAN, /BGLS1

    idxbs = where(bs ne -999.)
    RVSPY_compute_periodograms, resdirrv, 'BS', jd[idxbs], bs[idxbs], bse[idxbs], f0, f1, ofac, weight, /GLS;, /CLEAN, /BGLS1

    RVSPY_compute_periodograms, resdirew, 'S_Wilson', jd, s_wilson, dummyerr, f0, f1, ofac, weight, /GLS;, /CLEAN, /BGLS1
    RVSPY_compute_periodograms, resdirew, 'S_Feros', jd, s_feros, dummyerr, f0, f1, ofac, weight, /GLS;, /CLEAN, /BGLS1

    for i=0,ldrnlines-1 do begin

        RVSPY_compute_periodograms, resdirldr, 'LDR_'+ldrlinename[i]+'_'+ldrline[i], jd, ldr[*,i], ldre[*,i], f0, f1, ofac, weight, /GLS;, /CLEAN, /BGLS1

    endfor

    for i=0,ewnlines-1 do begin

        RVSPY_compute_periodograms, resdirew, 'EW_'+ewline[i], jd, ew[*,i], ewe[*,i], f0, f1, ofac, weight, /GLS;, /CLEAN, /BGLS1

    endfor

    for i=0,ewnlines-1 do begin

        RVSPY_compute_periodograms, resdirrf, 'RF_'+ewline[i], jd, resflux[*,i], resfluxe[*,i], f0, f1, ofac, weight, /GLS;, /CLEAN, /BGLS1

    endfor

;     for i=0,npbi-1 do begin
; 
;         RVSPY_compute_periodograms, resdirpbi, 'PBI_'+pbiband[i], jd, pbiindex[*,i], pbiindexe[*,i], f0, f1, ofac, weight, /GLS;, /CLEAN, /BGLS1
; 
;     endfor


    ;=========================================================================

    ;plots

    print, ''
    print, 'Creating Plots'
    
		;-------------------------------------------------------------------------
		
		if keyword_set(twhya) then begin
		
			print, 'GLS periodogram evolution'
			
			;GLS periodogram evolution
			
			weight = 2
			f0tmp = 0.01
			f1tmp = 1.

			start = 20

			;get the number of GLS frequencies
			tbase = max(jd)-min(jd)
			step = 1./tbase/ofac
			fsteps = long((f1tmp-f0tmp)/step+1)
			power = dblarr(n_elements(jd)-start+1,fsteps)
			freq = dblarr(n_elements(jd)-start+1,fsteps)
			ndata = dindgen(n_elements(jd)-start+1)+start

			for i=start-1,n_elements(jd)-1 do begin

			; 	RVSPY_compute_periodograms, '/home/amueller/Downloads/', 'test', jd[0:i], rv[0:i], rve[0:i], 0.01, 1, 20, 2, /GLS
				
				fn = 'tmp.txt'
				tmpjd = jd[0:i]
				tmprv = rv[0:i]
				tmprve = rve[0:i]
				openw, lun, pathGLS+fn, width=1400, /get_lun
					for j=0,n_elements(tmpjd)-1 do printf, lun, tmpjd[j], tmprv[j], tmprve[j], format='(f17.9,f23.15,f23.15)'
				close, lun
				free_lun, lun

				openw, lun, 'GLS.par', /get_lun
					printf, lun, fix(ofac)
					printf, lun, fix(weight)
					printf, lun, '0 0.0 0.0'
					printf, lun, '0 359'
				close, lun
				free_lun, lun
				
				spawn, pathGLS+'./GLS '+pathGLS+fn+' '+strcompress(f0tmp, /rem)+' '+strcompress(f1tmp, /rem)+' '+strcompress(mass, /rem)

				readcol, 'GLS.plt', tmpfreq, pow, win, pls, format='d,d,d,d', /silent
				power[i-start+1,0:n_elements(pow)-1] = pow
				freq[i-start+1,0:n_elements(pow)-1] = tmpfreq

			;   spawn, 'mv GLS.plt '+resdir+'GLS_'+ofn+'.txt'
				spawn, 'rm GLS.plt GLS.log GLS.par RVSinRes.plt GLSFit.plt'


			endfor

			power = transpose(power)
			
			;interpolate frequency
			reffreq = freq[n_elements(freq[*,0])-1,*]
			ipow = power
			for i=0,n_elements(ndata)-1 do begin

				idx = where(freq[i,*] ne 0.)
				tmpfreq = freq[i,idx]
				tmppow = power[idx,i]
				ipow[*,i] = interpol(tmppow, tmpfreq, reffreq)

			endfor
			
			set_plot, 'ps'
			device, filename=resdirrv+'GLS_RV_evolution.ps',/color, xsize=18, ysize=15

			!P.Font=1
			!p.thick=4
			!x.thick=3
			!y.thick=3
			!p.multi=[0,1,0]

			nlevels = 20
			loadct, 33, NColors=nlevels, Bottom=1, /silent
			maxcorr = 1.d0
			mincorr = 0.d0
			step = (maxcorr-mincorr) / nlevels
			levels = IndGen(nlevels) * step+mincorr
			SetDecomposedState, 0, CurrentState=currentState

			contour, ipow, 1./reffreq, ndata, /fill, levels=levels, color=cgcolor('black'), C_Colors=IndGen(nlevels)+1, xtitle='Period [days]', ytitle='N Observations', $
			Position=[0.17, 0.170, 0.95, 0.85], $
			xr=[3.3,3.8], yr=[min(ndata),max(ndata)], xst=9, yst=9, $
			charsize=2., xticklen=-0.04, yticklen=-0.03;, xticklen=0.02, yticklen=0.02, xminor=2, yminor=2, 

			SetDecomposedState, currentState
			tn = ['0.0', ' ', '0.1', ' ', '0.2', ' ', '0.3', ' ', '0.4', ' ', '0.5', ' ', '0.6', ' ', '0.7', ' ', '0.8', ' ', '0.9', ' ', '1.0']
			cgColorbar, Range=[mincorr, maxcorr], Divisions=20, XTicklen=1, XMinor=0, $
			AnnotateColor='black', ncolors=20, bottom=1, Charsize=1, $
			Position=[0.17, 0.95, 0.95, 0.97], $
			title='GLS Power', ticknames=tn;, xtickformat="(A1)";, , format="(f4.1)"

			!p.multi=[0,1,0]
			!P.Font=0
			!p.thick=1
			!x.thick=1
			!y.thick=1

			device,/close
			set_plot,'x'
			spawn, epstopdf+resdirrv+'GLS_RV_evolution.ps'
			spawn, 'rm '+resdirrv+'GLS_RV_evolution.ps' 
			
			power = ipow
			freq = reffreq
			
			save, power, freq, ndata, filename=resdirrv+'Values_GLS_RV_evolution.sav'
			
		endif

;-------------------------------------------------------------------------

    !p.font=1
    !p.thick=4
    !x.thick=3
    !y.thick=3

    ;RV
    print, 'Plot RV'
    set_plot, 'ps'
    device, filename=resdirrv+'RV.ps',/color, xsize=20, ysize=14
        ploterror, jd-2450000., rv, rve, psym=sym(1), errthick=3, /nohat, /yn, xtitle='Time [BJD-2450000]', ytitle='RV [km/s]', charsize=2, position=[0.15,0.16,0.95,0.9]
    device,/close
    set_plot,'x'
    spawn, epstopdf+resdirrv+'RV.ps'
    spawn, 'rm '+resdirrv+'RV.ps'

    ;BS
    print, 'Plot BS'
    set_plot, 'ps'
    device, filename=resdirrv+'BS.ps',/color, xsize=20, ysize=14
        ploterror, jd[idxbs]-2450000., bs[idxbs], bse[idxbs], psym=sym(1), errthick=3, /nohat, /yn, xtitle='Time [BJD-2450000]', ytitle='Bisector Span [km/s]', charsize=2, position=[0.15,0.16,0.95,0.9]
    device,/close
    set_plot,'x'
    spawn, epstopdf+resdirrv+'BS.ps'
    spawn, 'rm '+resdirrv+'BS.ps'

    ;RV vs BS
    print, 'Plot RV vs BS'
    ; FITEXY, rv, bs, A, B, X_SIG=rve, Y_SIG=bse, sigma_A_B
    ; sixlin, rv, bs, A, Ae, B, Be
    coeff = robust_linefit(rv, bs, YFIT, SIG, COEF_SIG)
    ; linear_corrcoeff_probability, rv, bs, corr, prob_corr
    set_plot, 'ps'
    device, filename=resdirrv+'RV_BS.ps',/color, xsize=20, ysize=14
        ploterror, rv[idxbs], bs[idxbs], rve[idxbs], bse[idxbs], psym=sym(1), errthick=3, /nohat, /yn, xtitle='RV [km/s]', ytitle='Bisector Span [km/s]', charsize=2, position=[0.15,0.16,0.95,0.9]
        oplot, rv, coeff[1]*rv+coeff[0], color=cgcolor('red')	
    ; 	legend, ['r='+sigfig(corr,2), 'p='+strcompress(string(prob_corr, format='(e7.1)'),/rem)], margin=0, /right, /top, box=0, charsize=1.5
    device,/close
    set_plot,'x'
    spawn, epstopdf+resdirrv+'RV_BS.ps'
    spawn, 'rm '+resdirrv+'RV_BS.ps'

    ;EW
    print, 'Plot EW'
    for i=0,ewnlines-1 do begin
        set_plot, 'ps'
        device, filename=resdirew+'EW_'+ewline[i]+'.ps',/color, xsize=20, ysize=14
        ploterror, jd-2450000., ew[*,i], ewe[*,i], psym=sym(1), errthick=3, /nohat, /yn, xtitle='Time [BJD-2450000]', ytitle='EW ['+cgSymbol('Angstrom')+']', title=ewline[i], charsize=2, position=[0.15,0.16,0.95,0.9]
        device,/close
        set_plot,'x'
        spawn, epstopdf+resdirew+'EW_'+ewline[i]+'.ps'
        spawn, 'rm '+resdirew+'EW_'+ewline[i]+'.ps'
    endfor
    
    ;10% Halpha
    if keyword_set(twhya) then begin
    
			print, 'Plot RV vs 10% Halpha Line Width'
			; FITEXY, rv, bs, A, B, X_SIG=rve, Y_SIG=bse, sigma_A_B
			; sixlin, rv, bs, A, Ae, B, Be
			coeff = robust_linefit(rv, halpha10width, YFIT, SIG, COEF_SIG)
			; linear_corrcoeff_probability, rv, bs, corr, prob_corr
			set_plot, 'ps'
			device, filename=resdirindex+'RV_Halpha10width.ps',/color, xsize=20, ysize=14
					plot, rv, halpha10width, psym=sym(1), /yn, xtitle='RV [km/s]', ytitle='H'+alpha+' 10% width [km/s]', charsize=2, position=[0.15,0.16,0.95,0.9]
					oplot, rv, coeff[1]*rv+coeff[0], color=cgcolor('red')	
			; 	legend, ['r='+sigfig(corr,2), 'p='+strcompress(string(prob_corr, format='(e7.1)'),/rem)], margin=0, /right, /top, box=0, charsize=1.5
			device,/close
			set_plot,'x'
			spawn, epstopdf+resdirindex+'RV_Halpha10width.ps'
			spawn, 'rm '+resdirindex+'RV_Halpha10width.ps'
			
		endif

		;line index
    print, 'Plot RV vs Halpha Line Index'
    ; FITEXY, rv, bs, A, B, X_SIG=rve, Y_SIG=bse, sigma_A_B
    ; sixlin, rv, bs, A, Ae, B, Be
    coeff = robust_linefit(rv, halpha_index, YFIT, SIG, COEF_SIG)
    ; linear_corrcoeff_probability, rv, bs, corr, prob_corr
    set_plot, 'ps'
    device, filename=resdirindex+'RV_HalphaIndex.ps',/color, xsize=20, ysize=14
        ploterror, rv, halpha_index, rve, halpha_index_err, psym=sym(1), errthick=3, /nohat, /yn, xtitle='RV [km/s]', ytitle='H'+alpha+' Index', charsize=2, position=[0.15,0.16,0.95,0.9]
        oplot, rv, coeff[1]*rv+coeff[0], color=cgcolor('red')	
    ; 	legend, ['r='+sigfig(corr,2), 'p='+strcompress(string(prob_corr, format='(e7.1)'),/rem)], margin=0, /right, /top, box=0, charsize=1.5
    device,/close
    set_plot,'x'
    spawn, epstopdf+resdirindex+'RV_HalphaIndex.ps'
    spawn, 'rm '+resdirindex+'RV_HalphaIndex.ps'    
    
    print, 'Plot Halpha Line Index'
    set_plot, 'ps'
    device, filename=resdirindex+'HalphaIndex.ps',/color, xsize=20, ysize=14
        ploterror, jd-2450000., halpha_index, halpha_index_err, psym=sym(1), errthick=3, /nohat, /yn, xtitle='Time [BJD-2450000]', ytitle='H'+alpha+' Index', charsize=2, position=[0.15,0.16,0.95,0.9]
    device,/close
    set_plot,'x'
    spawn, epstopdf+resdirindex+'HalphaIndex.ps'
    spawn, 'rm '+resdirindex+'HalphaIndex.ps'
    
    
		print, 'Plot RV vs HeID3 Line Index'
    ; FITEXY, rv, bs, A, B, X_SIG=rve, Y_SIG=bse, sigma_A_B
    ; sixlin, rv, bs, A, Ae, B, Be
    coeff = robust_linefit(rv, he1d3_index, YFIT, SIG, COEF_SIG)
    ; linear_corrcoeff_probability, rv, bs, corr, prob_corr
    set_plot, 'ps'
    device, filename=resdirindex+'RV_HeID3Index.ps',/color, xsize=20, ysize=14
        ploterror, rv, he1d3_index, rve, he1d3_index_err, psym=sym(1), errthick=3, /nohat, /yn, xtitle='RV [km/s]', ytitle='HeI D3 Index', charsize=2, position=[0.15,0.16,0.95,0.9]
        oplot, rv, coeff[1]*rv+coeff[0], color=cgcolor('red')	
    ; 	legend, ['r='+sigfig(corr,2), 'p='+strcompress(string(prob_corr, format='(e7.1)'),/rem)], margin=0, /right, /top, box=0, charsize=1.5
    device,/close
    set_plot,'x'
    spawn, epstopdf+resdirindex+'RV_HeID3Index.ps'
    spawn, 'rm '+resdirindex+'RV_HeID3Index.ps'
    
    print, 'Plot HeI D3 Line Index'
    set_plot, 'ps'
    device, filename=resdirindex+'HeID3Index.ps',/color, xsize=20, ysize=14
        ploterror, jd-2450000., he1d3_index, he1d3_index_err, psym=sym(1), errthick=3, /nohat, /yn, xtitle='Time [BJD-2450000]', ytitle='HeI D3 Index', charsize=2, position=[0.15,0.16,0.95,0.9]
    device,/close
    set_plot,'x'
    spawn, epstopdf+resdirindex+'HeID3Index.ps'
    spawn, 'rm '+resdirindex+'HeID3Index.ps'   

		print, 'Plot RV vs NaID1 Line Index'
    ; FITEXY, rv, bs, A, B, X_SIG=rve, Y_SIG=bse, sigma_A_B
    ; sixlin, rv, bs, A, Ae, B, Be
    coeff = robust_linefit(rv, na1d1_index, YFIT, SIG, COEF_SIG)
    ; linear_corrcoeff_probability, rv, bs, corr, prob_corr
    set_plot, 'ps'
    device, filename=resdirindex+'RV_NaID1Index.ps',/color, xsize=20, ysize=14
        ploterror, rv, na1d1_index, rve, na1d1_index_err, psym=sym(1), errthick=3, /nohat, /yn, xtitle='RV [km/s]', ytitle='NaI D1 Index', charsize=2, position=[0.15,0.16,0.95,0.9]
        oplot, rv, coeff[1]*rv+coeff[0], color=cgcolor('red')	
    ; 	legend, ['r='+sigfig(corr,2), 'p='+strcompress(string(prob_corr, format='(e7.1)'),/rem)], margin=0, /right, /top, box=0, charsize=1.5
    device,/close
    set_plot,'x'
    spawn, epstopdf+resdirindex+'RV_NaID1Index.ps'
    spawn, 'rm '+resdirindex+'RV_NaID1Index.ps'
    
    print, 'Plot NaID1 Line Index'
    set_plot, 'ps'
    device, filename=resdirindex+'NaID1Index.ps',/color, xsize=20, ysize=14
        ploterror, jd-2450000., na1d1_index, na1d1_index_err, psym=sym(1), errthick=3, /nohat, /yn, xtitle='Time [BJD-2450000]', ytitle='NaI D1 Index', charsize=2, position=[0.15,0.16,0.95,0.9]
    device,/close
    set_plot,'x'
    spawn, epstopdf+resdirindex+'NaID1Index.ps'
    spawn, 'rm '+resdirindex+'NaID1Index.ps' 
    
		print, 'Plot RV vs NaI D2 Line Index'
    ; FITEXY, rv, bs, A, B, X_SIG=rve, Y_SIG=bse, sigma_A_B
    ; sixlin, rv, bs, A, Ae, B, Be
    coeff = robust_linefit(rv, na1d2_index, YFIT, SIG, COEF_SIG)
    ; linear_corrcoeff_probability, rv, bs, corr, prob_corr
    set_plot, 'ps'
    device, filename=resdirindex+'RV_NaID2Index.ps',/color, xsize=20, ysize=14
        ploterror, rv, na1d2_index, rve, na1d2_index_err, psym=sym(1), errthick=3, /nohat, /yn, xtitle='RV [km/s]', ytitle='NaI D2 Index', charsize=2, position=[0.15,0.16,0.95,0.9]
        oplot, rv, coeff[1]*rv+coeff[0], color=cgcolor('red')	
    ; 	legend, ['r='+sigfig(corr,2), 'p='+strcompress(string(prob_corr, format='(e7.1)'),/rem)], margin=0, /right, /top, box=0, charsize=1.5
    device,/close
    set_plot,'x'
    spawn, epstopdf+resdirindex+'RV_NaID2Index.ps'
    spawn, 'rm '+resdirindex+'RV_NaID2Index.ps'
    
    print, 'Plot NaI D2 Line Index'
    set_plot, 'ps'
    device, filename=resdirindex+'NaID2Index.ps',/color, xsize=20, ysize=14
        ploterror, jd-2450000., na1d2_index, na1d2_index_err, psym=sym(1), errthick=3, /nohat, /yn, xtitle='Time [BJD-2450000]', ytitle='NaI D2 Index', charsize=2, position=[0.15,0.16,0.95,0.9]
    device,/close
    set_plot,'x'
    spawn, epstopdf+resdirindex+'NaID2Index.ps'
    spawn, 'rm '+resdirindex+'NaID2Index.ps' 
    
    print, 'Plot Halpha Spectrum'
    set_plot, 'ps'
    device, filename=resdirindex+'Spectrum_Halpha.ps',/color, xsize=20, ysize=14
    
    idx = where(ml ge 6532. and ml le 6592.)
    yr = [0, max(intf[*,idx])+0.1*max(intf[*,idx])]
    
    plot, ml, intf[0,*], yr=yr, xr=[6532.,6592.], yst=1, xtitle='Wavelength ['+cgSymbol('Angstrom')+']', ytitle='Normalized Flux', /nodata, title='H'+alpha, charsize=2, position=[0.15,0.16,0.95,0.9]
    for i=0,nfiles-1 do oplot, ml, intf[i,*]
    oplot, (6562.808-0.678/2.)*[1,1], !y.crange, color=cgcolor('dark gray')
    oplot, (6562.808+0.678/2.)*[1,1], !y.crange, color=cgcolor('dark gray')
    oplot, 6545.495*[1,1], !y.crange, linestyle=1, color=cgcolor('dark gray')
    oplot, 6556.245*[1,1], !y.crange, linestyle=1, color=cgcolor('dark gray')
		oplot, 6575.934*[1,1], !y.crange, linestyle=1, color=cgcolor('dark gray')
    oplot, 6584.684*[1,1], !y.crange, linestyle=1, color=cgcolor('dark gray')
    oplot, ml, mf, color=cgcolor('gray')
    device,/close
    set_plot,'x'
    spawn, epstopdf+resdirindex+'Spectrum_Halpha.ps'
    spawn, 'rm '+resdirindex+'Spectrum_Halpha.ps'

    
    print, 'Plot HeI D3 Spectrum'
    set_plot, 'ps'
    device, filename=resdirindex+'Spectrum_He1D3.ps',/color, xsize=20, ysize=14
    
    idx = where(ml ge 5860. and ml le 5890.)
    yr = [0, max(intf[*,idx])+0.1*max(intf[*,idx])]
    
    plot, ml, intf[0,*], yr=yr, xr=[5860.,5890.], yst=1, xtitle='Wavelength ['+cgSymbol('Angstrom')+']', ytitle='Normalized Flux', /nodata, title='HeI D3', charsize=2, position=[0.15,0.16,0.95,0.9]
    for i=0,nfiles-1 do oplot, ml, intf[i,*]
    oplot, (5875.62-0.4/2.)*[1,1], !y.crange, color=cgcolor('dark gray')
    oplot, (5875.62+0.4/2.)*[1,1], !y.crange, color=cgcolor('dark gray')
    oplot, 5866.5*[1,1], !y.crange, linestyle=1, color=cgcolor('dark gray')
    oplot, 5871.5*[1,1], !y.crange, linestyle=1, color=cgcolor('dark gray')
		oplot, 5878.5*[1,1], !y.crange, linestyle=1, color=cgcolor('dark gray')
    oplot, 5883.5*[1,1], !y.crange, linestyle=1, color=cgcolor('dark gray')
    oplot, ml, mf, color=cgcolor('gray')
    device,/close
    set_plot,'x'
    spawn, epstopdf+resdirindex+'Spectrum_He1D3.ps'
    spawn, 'rm '+resdirindex+'Spectrum_He1D3.ps'   

    
    print, 'Plot NaI Spectrum'
    set_plot, 'ps'
    device, filename=resdirindex+'Spectrum_NaI.ps',/color, xsize=20, ysize=14
    
    idx = where(ml ge 5790. and ml le 6110.)
    yr = [0, max(intf[*,idx])+0.1*max(intf[*,idx])]
    
    plot, ml, intf[0,*], yr=yr, xr=[5790.,6110.], yst=1, xtitle='Wavelength ['+cgSymbol('Angstrom')+']', ytitle='Normalized Flux', /nodata, title='NaI D1 D2', charsize=2, position=[0.15,0.16,0.95,0.9]
    for i=0,nfiles-1 do oplot, ml, intf[i,*]
    oplot, (5895.92-0.5/2.)*[1,1], !y.crange, color=cgcolor('dark gray')
    oplot, (5895.92+0.5/2.)*[1,1], !y.crange, color=cgcolor('dark gray')
    oplot, (5889.95-0.5/2.)*[1,1], !y.crange, color=cgcolor('dark gray')
    oplot, (5889.95+0.5/2.)*[1,1], !y.crange, color=cgcolor('dark gray')
    oplot, 5800.*[1,1], !y.crange, linestyle=1, color=cgcolor('dark gray')
    oplot, 5810.*[1,1], !y.crange, linestyle=1, color=cgcolor('dark gray')
		oplot, 6080.*[1,1], !y.crange, linestyle=1, color=cgcolor('dark gray')
    oplot, 6100.*[1,1], !y.crange, linestyle=1, color=cgcolor('dark gray')
    oplot, ml, mf, color=cgcolor('gray')
    device,/close
    set_plot,'x'
    spawn, epstopdf+resdirindex+'Spectrum_NaI.ps'
    spawn, 'rm '+resdirindex+'Spectrum_NaI.ps'  
    

    ;LDR
    print, 'Plot LDR'
    for i=0,ldrnlines-1 do begin
        set_plot, 'ps'
        device, filename=resdirldr+'LDR_'+ldrlinename[i]+'_'+ldrline[i]+'.ps',/color, xsize=20, ysize=14
        ploterror, jd-2450000., ldr[*,i], ldre[*,i], psym=sym(1), errthick=3, /nohat, /yn, xtitle='Time [BJD-2450000]', ytitle='LDR', title=ldrlinename[i]+'_'+ldrline[i], charsize=2, position=[0.15,0.16,0.95,0.9]
        device,/close
        set_plot,'x'
        spawn, epstopdf+resdirldr+'LDR_'+ldrlinename[i]+'_'+ldrline[i]+'.ps'
        spawn, 'rm '+resdirldr+'LDR_'+ldrlinename[i]+'_'+ldrline[i]+'.ps'
    endfor

    ;Residual Flux
    print, 'Plot Residual Flux'
    for i=0,ewnlines-1 do begin
        set_plot, 'ps'
        device, filename=resdirrf+'RF_'+ewline[i]+'.ps',/color, xsize=20, ysize=14
        ploterror, jd-2450000., resflux[*,i], resfluxe[*,i], psym=sym(1), errthick=3, /nohat, /yn, xtitle='Time [BJD-2450000]', ytitle='Residual Flux', title=ewline[i], charsize=2, position=[0.15,0.16,0.95,0.9]
        device,/close
        set_plot,'x'
        spawn, epstopdf+resdirrf+'RF_'+ewline[i]+'.ps'
        spawn, 'rm '+resdirrf+'RF_'+ewline[i]+'.ps'
    endfor

;     ;PBI
;     print, 'Plot PBI'
;     for i=0,npbi-1 do begin
;         set_plot, 'ps'
;         device, filename=resdirpbi+'PBI_'+pbiband[i]+'.ps',/color, xsize=20, ysize=14
;         ploterror, jd-2450000., pbiindex[*,i], pbiindexe[*,i], psym=sym(1), errthick=3, /nohat, /yn, xtitle='Time [BJD-2450000]', ytitle='PBI', title=pbiband[i], charsize=2, position=[0.15,0.16,0.95,0.9]
;         device,/close
;         set_plot,'x'
;         spawn, epstopdf+resdirpbi+'PBI_'+pbiband[i]+'.ps'
;         spawn, 'rm '+resdirpbi+'PBI_'+pbiband[i]+'.ps'
;     endfor

    ;Ca H&K
    print, 'Plot Ca H&K Spectrum'
    set_plot, 'ps'
    device, filename=resdirew+'Spectrum_CaHK.ps',/color, xsize=20, ysize=14
    
    idx = where(ml ge 3925. and ml le 3980.)
    yr = [0, max(intf[*,idx])+0.1*max(intf[*,idx])]
    
    plot, ml, intf[0,*], yr=yr, xr=[3920,3980], yst=1, xtitle='Wavelength ['+cgSymbol('Angstrom')+']', ytitle='Normalized Flux', /nodata, title='Ca H&K', charsize=2, position=[0.15,0.16,0.95,0.9]
    for i=0,nfiles-1 do oplot, ml, intf[i,*]
    oplot, CaK*[1,1], !y.crange, linestyle=1, color=cgcolor('dark gray')
    oplot, CaH*[1,1], !y.crange, linestyle=1, color=cgcolor('dark gray')
    oplot, ml, mf, color=cgcolor('gray')
    device,/close
    set_plot,'x'
    spawn, epstopdf+resdirew+'Spectrum_CaHK.ps'
    spawn, 'rm '+resdirew+'Spectrum_CaHK.ps'

    ;S_Wilson
    print, 'Plot S_Wilson'
    set_plot, 'ps'
    device, filename=resdirew+'S_Wilson.ps',/color, xsize=20, ysize=14
        plot, jd-2450000., s_wilson, psym=sym(1), /yn, xtitle='Time [BJD-2450000]', ytitle='S!DWilson!N', charsize=2, position=[0.15,0.16,0.95,0.9]
    device,/close
    set_plot,'x'
    spawn, epstopdf+resdirew+'S_Wilson.ps'
    spawn, 'rm '+resdirew+'S_Wilson.ps'

    ;S_Feros
    print, 'Plot S_Feros'
    set_plot, 'ps'
    device, filename=resdirew+'S_Feros.ps',/color, xsize=20, ysize=14
        plot, jd-2450000., s_feros, psym=sym(1), /yn, xtitle='Time [BJD-2450000]', ytitle='S!DFeros!N', charsize=2, position=[0.15,0.16,0.95,0.9]
    device,/close
    set_plot,'x'
    spawn, epstopdf+resdirew+'S_Feros.ps'
    spawn, 'rm '+resdirew+'S_Feros.ps'

    ;EW spectrum
    print, 'Plot EW Spectrum'
    for xx=0,ewnlines-1 do begin

        set_plot, 'ps'
        device, filename=resdirew+'Spectrum_EW_'+ewline[xx]+'.ps',/color, xsize=20, ysize=14
        
        xr = [ewcl[xx]-ewwin[xx]-10.,ewcl[xx]+ewwin[xx]+10.]
        idx = where(ml ge xr[0] and ml le xr[1])
        yr = [0, max(intf[*,idx])+0.1*max(intf[*,idx])]
        
        plot, ml, intf[0,*], yr=yr, xr=xr, yst=1, xtitle='Wavelength ['+cgSymbol('Angstrom')+']', ytitle='Normalized Flux', /nodata, title=ewline[xx], charsize=2, position=[0.15,0.16,0.95,0.9]
        for i=0,nfiles-1 do oplot, ml, intf[i,*]
        oplot, (ewcl[xx]-ewwin[xx])*[1,1], !y.crange, linestyle=1, color=cgcolor('dark gray')
        oplot, (ewcl[xx]+ewwin[xx])*[1,1], !y.crange, linestyle=1, color=cgcolor('dark gray')
        oplot, ewcl[xx]*[1,1], !y.crange, linestyle=1
        oplot, ml, mf, color=cgcolor('gray')
        
        device,/close
        set_plot,'x'
        spawn, epstopdf+resdirew+'Spectrum_EW_'+ewline[xx]+'.ps'
        spawn, 'rm '+resdirew+'Spectrum_EW_'+ewline[xx]+'.ps'

    endfor

    ;LDR spectrum
    print, 'Plot LDR Spectrum'
    for xx=0,ldrnlines-1 do begin

        set_plot, 'ps'
        device, filename=resdirldr+'Spectrum_LDR_'+ldrlinename[xx]+'_'+ldrline[xx]+'.ps',/color, xsize=20, ysize=14
        
        if (ldrl1[xx] lt ldrl2[xx]) then xr = [ldrl1[xx]-5.,ldrl2[xx]+5.] else xr = [ldrl2[xx]-5.,ldrl1[xx]+5.]
        plot, ml, intf[0,*], yr=[0,1.1], xr=xr, yst=1, xtitle='Wavelength ['+cgSymbol('Angstrom')+']', ytitle='Normalized Flux', /nodata, title=ldrlinename[xx]+'_'+ldrline[xx], charsize=2, position=[0.15,0.16,0.95,0.9]
        for i=0,nfiles-1 do oplot, ml, intf[i,*]
        oplot, ldrl1[xx]*[1,1], !y.crange, linestyle=1, color=cgcolor('dark gray')
        oplot, ldrl2[xx]*[1,1], !y.crange, linestyle=1, color=cgcolor('dark gray')
        oplot, ml, mf, color=cgcolor('gray')
        
        device,/close
        set_plot,'x'
        spawn, epstopdf+resdirldr+'Spectrum_LDR_'+ldrlinename[xx]+'_'+ldrline[xx]+'.ps'
        spawn, 'rm '+resdirldr+'Spectrum_LDR_'+ldrlinename[xx]+'_'+ldrline[xx]+'.ps'
            
    endfor

    ;residual flux spectrum
    print, 'Plot Residual Flux Spectrum'
    for xx=0,ewnlines-1 do begin

        set_plot, 'ps'
        device, filename=resdirrf+'Spectrum_RF_'+ewline[xx]+'.ps',/color, xsize=20, ysize=14

        xr = [ewcl[xx]-10.,ewcl[xx]+10.]
        
        ;get yrange
        idx = where(ml ge xr[0] and ml le xr[1])
        tmpf = dblarr(nfiles,n_elements(idx))
        for i=0,nfiles-1 do tmpf[i,*] = intf[i,idx]-mf[idx]
        yr = [min(tmpf)-0.1*abs(min(tmpf)),max(tmpf)+0.1*max(tmpf)]
        
        plot, ml, intf[0,*]-mf, xr=xr, yr=yr, yst=1, xst=1, xtitle='Wavelength ['+cgSymbol('Angstrom')+']', ytitle='Residual Flux', /nodata, title=ewline[xx], charsize=2, position=[0.15,0.16,0.95,0.9]
        for i=0,nfiles-1 do oplot, ml, intf[i,*]-mf
        oplot, ewcl[xx]*[1,1], !y.crange, linestyle=1, color=cgcolor('dark gray')
				oplot, (ewcl[xx]-ewwin[xx])*[1,1], !y.crange, linestyle=1, color=cgcolor('dark gray')
        oplot, (ewcl[xx]+ewwin[xx])*[1,1], !y.crange, linestyle=1, color=cgcolor('dark gray')

        device,/close
        set_plot,'x'
        spawn, epstopdf+resdirrf+'Spectrum_RF_'+ewline[xx]+'.ps'
        spawn, 'rm '+resdirrf+'Spectrum_RF_'+ewline[xx]+'.ps'

    endfor

    ;photospheric band indices
;     print, 'Plot PBI Spectrum'
;     for xx=0,npbi-1 do begin
; 
;         set_plot, 'ps'
;         device, filename=resdirpbi+'Spectrum_PBI_'+pbiband[xx]+'.ps',/color, xsize=20, ysize=14
; 
;         tmp = [numer[*,xx], denom[*,xx]]
;         xr = [min(tmp)-2,max(tmp)+2]
;         
;         ;get yrange
;         idx = where(ml ge xr[0] and ml le xr[1])
;         tmpf = dblarr(nfiles,n_elements(idx))
;         for i=0,nfiles-1 do tmpf[i,*] = intf[i,idx]-mf[idx]
;         
;         plot, ml, intf[0,*], xr=xr, yr=[0,1.2], yst=1, xtitle='Wavelength ['+cgSymbol('Angstrom')+']', ytitle='Normalized Flux', /nodata, title=pbiband[xx], charsize=2, position=[0.15,0.16,0.95,0.9]
;         for i=0,nfiles-1 do oplot, ml, intf[i,*]
;         oplot, ml, mf, color=cgcolor('gray')
;         plots, [numer[0,xx],numer[1,xx]], [1.1,1.1]
;         plots, [numer[0,xx],numer[0,xx]], [1.08,1.12]
;         plots, [numer[1,xx],numer[1,xx]], [1.08,1.12]
;         plots, [denom[0,xx],denom[1,xx]], [1.1,1.1], linestyle=1
;         plots, [denom[0,xx],denom[0,xx]], [1.08,1.12], linestyle=1
;         plots, [denom[1,xx],denom[1,xx]], [1.08,1.12], linestyle=1
;     ; 	oplot, ewcl[xx]*[1,1], !y.crange, linestyle=1, color=cgcolor('dark gray')
; 
;         device,/close
;         set_plot,'x'
;         spawn, epstopdf+resdirpbi+'Spectrum_PBI_'+pbiband[xx]+'.ps'
;         spawn, 'rm '+resdirpbi+'Spectrum_PBI_'+pbiband[xx]+'.ps'
; 
;     endfor

endfor

stop
end
