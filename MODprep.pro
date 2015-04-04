PRO MODprep

;filename = ''
;READ, filename, PROMPT = 'Name of the input file: '

SPAWN, 'ls *.fit > MODprep.lis'

newfile = ''
READ, newfile, PROMPT = 'Name of the output file: '
reference_lamp = ''
READ, reference_lamp, PROMPT = 'Specify the reference lamp: '

CLOSE, 1
OPENW, 1, newfile

READCOL, 'MODprep.lis', longname, f='a'

rawtype = ''
final_type = ''

FOR i=0, N_ELEMENTS(longname)-1 DO BEGIN

   a = READFITS(longname[i], h, /SILENT)

   rawtype = STRUPCASE(STRTRIM(SXPAR(h,'IMAGETYP'),2))
   dataregion = STRTRIM(SXPAR(h,'DATASEC'),2)
   biasregion = STRTRIM(SXPAR(h,'BIASSEC'),2)

   final_type = 'no_type'
   IF rawtype EQ 'FOCUS' THEN final_type = 'lamp'
   IF rawtype EQ 'LAMP' THEN final_type = 'lamp'
   IF rawtype EQ 'LAMPS' THEN final_type = 'lamp'
   IF rawtype EQ 'BIAS' THEN final_type = 'bias'
   IF rawtype EQ 'OBJECT' THEN final_type = 'obj'
   IF rawtype EQ 'FLATFIELD' THEN final_type = 'flat'
   IF rawtype EQ 'FLAT' THEN final_type = 'flat'
   IF rawtype EQ 'STANDARD' THEN final_type = 'std'

   raw_obj = SXPAR(h, 'OBJECT')
   ;PRINT, 'xx'+raw_obj+'xx'
   IF raw_obj EQ '        ' THEN raw_obj = 'no_name'

   ;remove the .fits extension from the filename
   name_length = STRLEN(longname[i])
   ;PRINT, longname[i], name_length
   where_fits = STRMID(longname[i],0,name_length-4) 

   printf,1, FORMAT = '(A20,x,A6,x,A15,x,A16,x,A16,x,A20)', where_fits, final_type, STRCOMPRESS(raw_obj,/REMOVE_ALL), '[1:300,1:2048]', '[305:364,1:2048]',reference_lamp

ENDFOR
CLOSE,1
END
