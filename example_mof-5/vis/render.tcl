#set TachyonPath "/home/nimish/vmd-1.9.4a12/lib/tachyon/tachyon_LINUXAMD64"
set TachyonPath "/home/nimish/Applications/vmd/lib/tachyon_LINUXAMD64"

set x 1200
set y 800
set scale 4
set xres [expr $x*$scale]
set yres [expr $y*$scale]

render Tachyon ${file_prefix} "${TachyonPath} -aasamples 12 %s -format TARGA -res $xres $yres -o %s.tga"

#exec /usr/bin/convert -trim $file_prefix.tga $file_prefix.jpeg
