function ghfunc,x,par,ipro=ipro,offset=offset,plot_key=plot_key

common psfstuff,param

; Gaus-Hermite-based function to model "arbitrarily wingy" instrumental
; profiles. 
;
; x (input vector) pixel offsets at which to evaluate function.
; par (input vector) profile function parameters as described below
;
; Returns the normalized profile function evaluated at each point in x.
;
; Modified 10/10/2002 by JohnJohn: Now uses make_herm to initialize
; GH coeffs. Other speed-based modifications made as well.
; up to 15 Gauss Hermite components used plus par[19] = width = 1/beta.

if n_params() lt 2 then begin
  print,'syntax: ghfunc(x,par)'
  retall
endif

;param.x = x ;This is here only because to preserve the interface that
            ;GPFUNC had. With GHFUNC, x (param.x) is global
amp = fan([1d,par[1:10],par[15:18]],121)

;No need for fancy mathematics if the width parameter doesn't change.
;SET_PARAM makes changes to param.p. Amplitudes can be adjusted
;without SET_PARAM
;if param.plotpsf then print,str(par[0])+' '+str(param.oldwid)
  if par[0] ne param.oldwid then begin
      param.wid = par[0]
      set_param
      param.oldwid = par[0]
  end 

  ipro = total(amp*param.p,1)

;PLOTTER is a movie routine if needed for testing or demo purposes
;if param.plotpsf then plotter,ipro,par

  return,ipro / int_tabulated(param.x,ipro)	;return normalized profile
end





