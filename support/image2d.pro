pro MAIN, SOURCENAME=sourcename

image2d, SOURCENAME=sourcename
skybrightness, SOURCENAME=sourcename

end

;;=============================================================================
;;
;; Plots a 2D image of the source from output from wuparam/wucut
;;
;;=============================================================================
PRO image2d, SOURCENAME=sourcename

!p.font = -1
!p.color=0

;; SET UP OUTPUT

path = GETENV('PWD')+"/"
subdir = "Totals/"
OUTPUTFILE = path + subdir + 'image2d-total.eps'
set_plot,'ps'
DEVICE, /ENCAPSULATED,/COLOR, FILENAME = OUTPUTFILE,/INCHES,XSIZE = 7.0, YSIZE = 7.0 

contlevels = INDGEN(16)
contlabels = [1, 1 , 1,  1,  1,  1,  1,  1 ]
clr    = [255,255,255,255,255,255,255,255,255,255,255]
nrows = 39
ncols = 39

loadImage2D, subdir+"ra-mesh.im2d", xcoord,ycoord,ara
loadImage2D, subdir+"dec-mesh.im2d", xcoord,ycoord,adec
loadImage2D, subdir+"total-signif.im2d", xcoord, ycoord, signif
loadImage2D, subdir+"total-excess.im2d", xcoord, ycoord, rex

;; PLOT THE 2D IMAGE AND SIGNIFICANCE CONTOURS:

LoadCT, 5, /SILENT
smooth=4
rex = transpose(rex)
rex = CONGRID( rex, nrows*smooth,ncols*smooth, CUBIC=-0.5, /MINUS_ONE )

plotImage2d, signif, rex,xcoord,ycoord, TITLE=sourcename

print, "MAX SIGNIF=", max(signif)
print, "MIN SIGNIF=", min(signif)
print, "MAX EXCESS=",max(rex)
print, "MIN EXCESS=",min(rex)

;; PLOT THE RA/DEC CONTOURS

ara = ara * 12/!PI
adec = adec * 180/!PI

contour, ara, xcoord, ycoord, $   
  C_LABELS = contlabels, $
  XSTYLE = 1, YSTYLE = 1, $
  LEVELS = FINDGEN(20) * 0.05 + floor(min(ara)*10)/10.0, $
  COLOR = 200,$
  /NOERASE, /OVERPLOT

contour, adec, xcoord, ycoord,$   
  C_LABELS = contlabels,$
  XSTYLE = 1, YSTYLE = 1, $
  LEVELS = FINDGEN(20) * 0.5  + min(adec), $
  COLOR = 200, $
  /NOERASE,/OVERPLOT

;; Overplot the stars

plotStars

;; FINISH UP

DEVICE, /CLOSE

END

;;=============================================================================
;;
;; Plot the sky brightness map!
;;
;;=============================================================================
PRO skyBrightness, SOURCENAME=sourcename 

!p.font = -1
path = GETENV('PWD') + "/"
subdir = "Totals/"
OUTPUTFILE = path + subdir + 'skybrightness-total.eps'
set_plot,'ps'
DEVICE, /ENCAPSULATED,/COLOR, FILENAME = OUTPUTFILE,$
  /INCHES,XSIZE = 7.0, YSIZE = 7.0 
loadCT,1, /SILENT

PX = !X.WINDOW * !D.X_VSIZE
PY = !Y.WINDOW * !D.Y_VSIZE

; Plot the sky brightness map 

loadImage2D, subdir+"skybrightness.im2d", xcoord, ycoord, skybright

xmin = xcoord[0]
ymin = ycoord[0]
xmax = xcoord[N_ELEMENTS(xcoord)-1]
ymax = ycoord[N_ELEMENTS(ycoord)-1]
dx = xcoord[1] - xcoord[0]
dy = ycoord[1] - ycoord[0]

if keyword_set(sourcename) then title="!6"+sourcename+" - Sky Brightness" $
  else title = "!6Sky Brightness"
contour, skybright,xcoord,ycoord, NLEVELS=100,/FILL,/DATA,TITLE=title,/NOERASE

; overplot the tubeoffness contours 
loadImage2D, subdir+"tubeoffness.im2d", xcoord, ycoord, tubeoff
print, "MAX= ", max(tubeoff), " MIN=",min(tubeoff)
levs = reverse(-(findgen(5))*(max(tubeoff)-min(tubeoff))/5)
print,levs
contour, -tubeoff, xcoord, ycoord, LEVELS=levs, /DATA, /NOERASE, /DOWNHILL

plotStars

DEVICE, /CLOSE


; Now plot the tubeoffness map alone

OUTPUTFILE = path + subdir + 'tubeoffness-total.eps'
DEVICE, /ENCAPSULATED,/COLOR, FILENAME = OUTPUTFILE,$
  /INCHES,XSIZE = 7.0, YSIZE = 7.0 

loadImage2D, subdir+"tubeoffness.im2d", xcoord, ycoord, tubeoff

if keyword_set(sourcename) then title="!6"+sourcename+" - Tubeoffness" $
  else title = "!6Tubeoffness"
contour, tubeoff, xcoord, ycoord, LEVELS=findgen(20)*max(tubeoff)/20, $
  /FILL,/DATA,TITLE=title
plotStars
DEVICE, /CLOSE



END

;;=============================================================================
;;
;; Overplot stars onto the current plot
;;
;;=============================================================================
PRO plotStars

!p.font = -1
path = GETENV('PWD') + "/"
subdir = "Totals/"
maxmag = 8


loadCT,1, /SILENT
loadNearbyStars, subdir+"nearbystars-on.dat", sx_on,sy_on,smag_on,$
  sname_on,starcatid_on
loadNearbyStars, subdir+"nearbystars-off.dat", sx_off,sy_off,smag_off, $
  sname_off,starcatid_off

;;---------------------------------------------------------
;; Just for debugging, theta should be 0! 
theta = 0;28*!PI/180
rsx_on = sx_on*cos(theta) + sy_on*sin(theta)
rsy_on = -sx_on*sin(theta) + sy_on*cos(theta)

rsx_off = sx_off*cos(theta) + sy_off*sin(theta)
rsy_off = -sx_off*sin(theta) + sy_off*cos(theta)
;;---------------------------------------------------------

a = FINDGEN(17)*(!PI*2/16)
USERSYM, cos(a), sin(a), /FILL

for i=0,N_ELEMENTS(sx_on)-1 do begin
    oplot, rsx_on[i:i],  rsy_on[i:i],  PSYM=8, SYMSIZE=(MAXMAG-smag_on[i])/4,  COLOR=255
    if starcatid_on[i] eq 0 then $
      xyouts, rsx_on[i], rsy_on[i]-0.07, "!6"+sname_on[i], ALIGNMENT=0.5, CHARSIZE=0.4,COLOR=255 

endfor

for i=0,N_ELEMENTS(sx_off)-1 do begin
    oplot, rsx_off[i:i], rsy_off[i:i], PSYM=7, SYMSIZE=(MAXMAG-smag_off[i])/4, COLOR=255
    if starcatid_off[i] eq 0 then $
      xyouts, rsx_off[i], rsy_off[i]-0.07, "!6"+sname_off[i], ALIGNMENT=0.5, CHARSIZE=0.4,COLOR=255

endfor


END

;;=============================================================================
;;
;; Read an .im2d file and return the image as a 2d array
;;
;;=============================================================================
PRO loadImage2D, filename, xcoord, ycoord, matrix

path = GETENV('PWD')
line = { x:0., y:0., val:0. }

fullname = path + "/" + filename
print, "Loading from '", fullname,"'"

close,2
openr, 2, fullname

data = replicate( line, 39*39 )
readf, 2, data

xcoord = reform( data.x, 39,39 )
ycoord = reform( data.y, 39,39 ) 

matrix = reform( data.val, 39,39 )
close, 2

end

;;=============================================================================
;;
;; Read list of nearby stars for plotting
;;
;;=============================================================================
PRO loadNearbyStars, filename, starx,stary,starmag,starname,starcatid

path=GETENV('PWD')
line = { x:0., y:0., ra:0., dec:0., mag:0., id:0., name:'unknown'}

fullname = path+"/"+filename
;print, "Loading nearby stars from '",filename,"'"

close,3
openr,3, fullname
readf,3, numstars

starx = fltarr(numstars)
stary = fltarr(numstars)
starmag = fltarr(numstars)
starname = strarr(numstars)
starcatid=fltarr(numstars)

i=0
while (not EOF(3)) do begin
    readf, 3, line
    starx[i] = line.x
    stary[i] = line.y
    starmag[i] = line.mag
    starcatid[i] = line.id
    starname[i] = line.name
    i = i+1
endwhile
close,3  

;print, "Read ",i," stars."
  
end


pro plotImage2d, sig,exc,xcoord,ycoord, WINDOW_SCALE = window_scale, ASPECT = aspect, INTERP = interp, TITLE = title

on_error,2                      ;Return to caller if an error occurs
sz = size(exc)			;Size of image
if sz[0] lt 2 then message, 'Parameter not 2D'

	;set window used by contour
contour,[[0,0],[1,1]],/nodata, xstyle=4, ystyle = 4

px = !x.window * !d.x_vsize	;Get size of window in device units
py = !y.window * !d.y_vsize
swx = px[1]-px[0]		;Size in x in device units
swy = py[1]-py[0]		;Size in Y
six = float(sz[1])		;Image sizes
siy = float(sz[2])
aspi = six / siy		;Image aspect ratio
aspw = swx / swy		;Window aspect ratio
f = aspi / aspw			;Ratio of aspect ratios

if (!d.flags and 1) ne 0 then begin	;Scalable pixels?
  if keyword_set(aspect) then begin	;Retain aspect ratio?
				;Adjust window size
	if f ge 1.0 then swy = swy / f else swx = swx * f
	endif

  tvscl,exc,px[0],py[0],xsize = swx, ysize = swy, /device

endif else begin	;Not scalable pixels	
   if keyword_set(window_scale) then begin ;Scale window to image?
	tvscl,exc,px[0],py[0]	;Output image
	swx = six		;Set window size from image
	swy = siy
    endif else begin		;Scale window
	if keyword_set(aspect) then begin
		if f ge 1.0 then swy = swy / f else swx = swx * f
		endif		;aspect
	tv,poly_2d(bytscl(exc),$	;Have to resample image
		[[0,0],[six/swx,0]], [[0,siy/swy],[0,0]],$
		keyword_set(interp),swx,swy), $
		px[0],py[0]
	endelse			;window_scale
  endelse			;scalable pixels

mx = !d.n_colors-1		;Brightest color
smax = floor(max(sig))
colors = fltarr(255)
if smax le 255 and smax gt 3 then colors[smax-3:smax] = mx
contlevels = FINDGEN(255)+1
contlabels = fltarr(255)+1

; set up title

if keyword_set(title) then titlestring = "!6"+title+" - Excess and Significance" $
  else titlestring = "!6Excess and Significance"

; plot

if !d.name eq 'PS' then colors = mx - colors ;invert line colors for pstscrp
contour,sig,xcoord,ycoord,/noerase,/xst,/yst,$ ;Do the contour
  pos = [px[0],py[0], px[0]+swx,py[0]+swy],/dev,$
  c_color =  colors,  c_labels=contlabels, levels=contlevels, TITLE=titlestring
return
end


