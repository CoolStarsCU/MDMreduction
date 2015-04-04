;+
; NAME: 
;           FAN
;
; PURPOSE:
;           Take the outer product of the input ARRAY and a
;           UNIT_VECTOR to "fan out" a 1D vector into an array 
;           comprised of the vector repeated row-wise NFAN times.
;           Useful for array-wise mathematics (Look Ma, no FOR loops!)
;
; CALLING SEQUENCE:
;           result = fan(array [,nfan, /transpose])
;
; INPUTS:
;           ARRAY - 1D array, input vector
;           NFAN  - number of times to repeat the input vector,
;                   default is N_ELEMENTS(ARRAY)
;
; KEYWORD PARAMETERS:
;                     
;           TRANSPOSE - Repeat the input vector column-wise
;
; OUTPUTS:
;           A 2D array with N_ELEMENTS(ARRAY) columns and NFAN
;           rows.
;
; EXAMPLE:
;           Fan a FINDGEN of 3 elements, twice.
;
;           IDL> a = findgen(3)
;           IDL> print,fan(a,2)
;                 0.00000      1.00000      2.00000
;                 0.00000      1.00000      2.00000
;
; MODIFICATION HISTORY:
;           Created sometime in ought-2 by JohnJohn
; 06 Dec 2002 JohnJohn- Added some error handling at the beginning
;-
function fan,array,nfan,transpose=transpose
   on_error,2  ;if broke then return to sender
   if n_params() lt 1 then begin 
       message,'Syntax: f = fan(array [,nfan, /transpose])',/info
       return,-1
   endif

   if n_elements(nfan) eq 0 then nfan = n_elements(array)
   unit_vector = replicate(1d,nfan)   ;dblarr(nfan)+1.
   if keyword_set(transpose) then new = array##unit_vector $
     else new = unit_vector##array
   return,new
end

;------ 
;+
; NAME: 
;       FILLARR 
;
;
; PURPOSE:
;       This function generates an array from MIN to MAX with
;       step size DEL. If an integer number of steps cannot be
;       fit between MIN and MAX, then MAX will be adjusted to
;       be as close as the specified maximum as possible.
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE:
;       f = fillarr(n, min, max [,fan=, transfan=, /double])
;
;
; INPUTS:
;       DEL:  The desired step size
;       MIN:  The value of the first array element in F
;       MAX:  The value of the last array element in F if
;             (MAX-MIN)/DEL is an integer. Adjusted otherwise.
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
;       FANNNED:    Number of times the array is to be repeated.
;                   The final dimensions of F  will be 
;                   fix((MAX-MIN)/DEL) + 1 columns by FANNED ows.
;
;       /TRANSPOSE  Final dimensions of F wil be FAN columns by 
;                   fix((MAX-MIN)/DEL) + 1 rows if FAN is specified. 
;
; OUTPUTS:
;
;       F:    Final array. If input parameters are double precision,
;             then F will be double as well. F is float otherwise.
;
; RESTRICTIONS:
;
;       You'll need FAN.PRO to use the fan= keyword. 
;       http://astron.berkeley.edu/~johnjohn/idl.html#FAN
;
; EXAMPLE:
;
;         For an array that runs from 2 to 5 in steps of .7
;
;         IDL> f = fillarr(.7,2,5)
;         IDL> print, f
;            2.00000      2.70000      3.40000     4.10000    4.80000
;         
; MODIFICATION HISTORY:
; Written by John "JohnJohn" Johnson 21-Feb-2002
; 22-Feb-2002 JohnJohn- Fixed precision bug
; 23-Feb-2002 JohnJohn- Calculations performed in double precision. 
;                       Output in double precision if input is 
;                       double.
; 01-Mar-2002 JohnJohn- Much props to Tim Robishaw (Tdogg) for helping
;                       me understand machine precision and finally fixing
;                       the precision bug.
; 23 Apr 2002 JohnJohn- Modified the /FAN operation to match my new
;                       FAN procedure. Default direction of the
;                       fanning process is in the column direction,
;                       i.e. a 5-element array with FAN=2 will yeild a
;                       5x2 array rather than the other way around.
; 06 Dec 2002 JohnJohn- Modified the /FAN operation again to run using
;                       the actuall FAN proceedure which is faster
;                       than doing two separate operations for fanning
;                       and taking the transpose. duh.
; 14 Apr 2005 JohnJohn- Fixed bug where if n_params() eq 2, then MIN
;                       was being modified. Input variable now
;                       protected by renaming MIN as MININ.
;-
function fillarr,del,minin,max,fanned=fanned,transpose=transpose
;DEAL WITH HUMANS
on_error,2		;Return to caller if an error occurs
if n_params() lt 2 then message,'INCORRECT NUMBER OF INPUTS. Syntax: f = fillarr(del,min,max)',/ioerror

if n_params() eq 2 then begin
    max = minin[1]
    min = minin[0]
endif else min = minin

if max lt min then message,'MIN must be less than MAX',/ioerror
if del eq 0 then message,'DEL cannot equal 0',/ioerror

;if all of the input parameters are double, the return the answer in
;double precision.
doub = (size(del,/type) eq 5) and (size(min,/type) eq 5) and (size(max,/type) eq 5) or keyword_set(double)
del = double(del)
min = double(min)
max = double(max)
;ARG will go into A later. These are the only real calculations performed.
arg = (max-min)/del
;test for and correct rounding errors
rnd = round(arg)
eps = (machar(/double)).eps
if abs(rnd-arg) lt rnd*eps then arg = rnd else arg = fix(arg,type=3)

a = dindgen(arg+1)*del+min      ;can you believe there's all this code just to do this?

if n_elements(fanned) ne 0 then begin
    nfan = fanned
    if keyword_set(transpose) then a = fan(a,nfan,/transpose) $
       else a = fan(a,nfan)
endif

if not doub then a = float(a)

return,a 
end
;----- 

;+
; NAME: 
;       POLYFIT
;
; PURPOSE:
;	Fit a polynomial to a function using linear least-squares
; NOTE:
;       Twice as fast as POLY_FIT.PRO as tested by Robishaw's 
;       BENCHMARK.PRO due to the miracle of loopless IDL code.
;
;       Uses JohnJohn's FAN.PRO 
;       http://astron.berkeley.edu/~johnjohn/idlprocs/fan.pro
;
; CALLING SEQUENCE:
;       coeffs = polyfit(t, y_of_t, degree [, yfit, dum, covariance=])
; 
; KEYWORDS
; MODIFICATION HISTORY:
;       Written by JohnJohn long ago in ought 1 after attending Carl
;       Heiles' lecture on line fitting.
;-
function polyfit,t,y,deg,yfit,yfit1,covariance=cov,weight=w
;on_error,2
n = n_elements(t)
pow = indgen(deg+1)
powarr = fan(pow,n,/trans)
x =  fan(double(t),deg+1)
xarr = x^powarr
xarrt = transpose(xarr)
if keyword_set(w) then xarr = fan(double(w),deg+1)*x^powarr 

alpha = xarr##xarrt
beta = xarr##(double(y))
cov = invert(alpha)
a = cov##beta
if n_params() eq 4 then yfit = poly(t,a)
if n_params() eq 5 then yfit1 = poly(t,a) ;;Legacy code, kept for compatibility
return,a
end


;---- 

;+
; NAME:
;           CCPEAK
;
; PURPOSE:
;       Locates the precise location of the peak in the
;       cross-correlation function between two vectors.
;       (Locates LAG_max)
;
; CALLING SEQUENCE:
;
;       LAG_max = CCPEAK(VEC1, VEC2 [, RADIUS ])
;
; INPUTS:
;
;       VEC1, VEC2 - Functions to be cross-correlated

; OPTIONAL INPUTS:
;
;       RADIUS -  How many array elements around the 
;                 nominal peak where polynomial fit 
;                 should be performed.
;
; OUTPUTS:
;
;       LAG_max - Lag at which cross-correlation 
;                 function is maximized
;
; RESTRICTIONS:
;
;       Uses my POLYFIT procedure, not IDL's
;       POLY_FIT
;
;       Uses C_CORRELATE to perform cross-correlation
;
; MODIFICATION HISTORY:
; Written sometime in December 2002 by JohnJohn
; 07 June 2003 JJ - Fixed bug where CF was being indexed out of
; range. Also limited the minimum and maximum lag returned to be
; min(lag) and max(lag), respectively.
; 26 June 2003 JJ - Default radius is now 50 rather than 10
;-

function ccpeak,arr1, arr2, radius, ccf=cf, lag=lag
on_error,2		;Return to caller if an error occurs
n = n_elements(arr1)
if n_elements(radius) eq 0 then radius = 50
lag = fillarr(1,-radius,radius)
cf = c_correlate(arr1,arr2,lag)
dum = max(cf, ind)

srad = 3
sublag = lag[(ind-srad) > 0:(ind+srad) < (2*radius)]
subcf  = cf[(ind-srad) > 0:(ind+srad) < (2*radius)]
a = polyfit(sublag,subcf,2)
maxlag = - a[1]/(2.*a[2])
nlag = n_elements(lag)
if maxlag lt lag[0] then maxlag = lag[0]
if maxlag gt lag[nlag-1] then maxlag = lag[nlag-1]
return,maxlag
end

;-------

PRO OIshift_corr, wavecal_table

!P.MULTI=[0,1,1]

ref_sky = ''
READ, ref_sky, PROMPT = 'Name of the reference sky (leave off wavecal/sky.): '

READCOL, 'wavecal/sky.'+ref_sky, ref_wavelength, ref_flux, /SILENT

;now read in all the rest of the spectra (including the ref.spectrum) 
READCOL, wavecal_table, spectra, obj_code, name, f='a,a,a', /SILENT
N_skies = N_ELEMENTS(spectra)

!P.MULTI = [0,1,1]

CLOSE, 1
OPENW, 1, 'OI_shifts.tbl'

;now loop through all the skies and cross correlate w/ the reference
;sky from above.
FOR i=0, N_skies-1 DO BEGIN

   ;make sure the spectrum isn't a lamp before trying to correlate
   IF STRTRIM(obj_code[i]) NE 'lamp' THEN BEGIN

      ;read in the sky spectrum
      READCOL, 'wavecal/sky.'+spectra[i], wavelength, flux, /SILENT

      ;cross correlate w/ reference spectrum 
      pixel_shift = CCPEAK(flux, ref_flux, 10, CCF = ccf, LAG = lag)

      ;PRINT, pixel_shift  

      ;PLOT, lag, ccf
      max_ccf = MAX(ccf)
      median_ccf = MEDIAN(ccf)
      ;stddev_ccf = STDDEV(ccf)
      quality_factor = (max_ccf - median_ccf) / median_ccf
      ;oplot, [pixel_shift, pixel_shift], [-100,100], linestyle = 4

      PLOT, wavelength, flux, XRANGE = [5552., 5592.], PSYM = 10
      ;OPLOT, wavelength, flux, PSYM = 10, LINESTYLE = 4
      OPLOT, [5577.399,5577.399], [0,100000000], LINESTYLE = 5
      ;LOADCT, 13, /SILENT
      ;OPLOT, wavelength + dispersion*pixel_shift, flux, PSYM = 10, COLOR = 255
      ;PRINT, spectra[i], pixel_shift, max_ccf, quality_factor

      relevant = WHERE(wavelength GT 5552. AND wavelength LT 5592., n_relevant)

      IF n_relevant GT 2 AND quality_factor GT 0.075 THEN BEGIN

         gauss = GAUSSFIT(wavelength[relevant], flux[relevant], fit_values, ESTIMATES = [ MAX(flux[relevant],max_loc), wavelength[relevant[max_loc]] ,3.0 ,MEDIAN(flux[relevant]) ], NTERMS = 4, SIGMA = shift_sigmas)
         ;PLOT, wavelength[relevant], ref_flux[relevant]
         ;OPLOT, [5577.399,5577.399], [0,100000000], LINESTYLE = 5
         ;OPLOT, [0,10000], [MEDIAN(flux[relevant]),MEDIAN(flux[relevant])], LINESTYLE = 5
         OPLOT, [wavelength[relevant[max_loc]],wavelength[relevant[max_loc]]],[0,1000000], LINESTYLE = 6, COLOR = 85
         
         OH_center = fit_values[1]
         shift = 5577.34 - OH_center
         OH_center_err = shift_sigmas[1]
         OPLOT, wavelength+shift, flux, LINESTYLE = 6, COLOR = 120, PSYM = 10
         ;dispersion = ref_wavelength[relevant[1]] - ref_wavelength[relevant[0]]

;         PRINT, ref_wavelength[relevant[1]], ref_wavelength[relevant[0]], dispersion
      ENDIF ELSE BEGIN
         shift = 0
         OH_center_err = 10.
    ENDELSE

   ENDIF ELSE BEGIN
         shift = 0
         OH_center_err = 10.
   ENDELSE
   PRINT, spectra[i], quality_factor
      PRINT, 'Shift & error:', shift, OH_center_err
      Wait = GET_KBRD(1)


      printf,1, FORMAT = '(A15,x,A5,x,A20,x,F10.4,x,F8.5,x,F8.5)', spectra[i], obj_code[i], name[i], shift, OH_center_err, quality_factor

ENDFOR
CLOSE,1

END
