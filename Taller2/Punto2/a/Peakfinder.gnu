### a simple gnuplot peak finder
reset session


set terminal jpeg enhanced
set output 'perfil0Maxima.jpg'
set xlabel "ix" font "TimesNewRoman,14"
set ylabel "iy" font "TimesNewRoman,14"

FILE = "Perfil0.dat"
colX = 1
colY = 2

# extract all peaks
set table $Temp
    plot y2=y1=x2=x1=NaN FILE u (x0=x1,x1=x2,x2=column(colX), \
                                 y0=y1,y1=y2,y2=column(colY), \
                                 y0<y1 && y1>=y2 ? sprintf("%g %g",x1,y1) : '') w table
set table $Peaks        # remove empty lines and store first values in x0,y0 for later use
    plot $Temp u (column(-2)==0 && $0==0?x0=$1:$1):(column(-2)==0 && $0==0?y0=$2:$2) w table
unset table

isNaN(v) = v!=v         # check if value is NaN
min(a,b) = isNaN(a) && isNaN(b) ? NaN : isNaN(a) ? b : isNaN(b) ? a  : a<b ? a : b   # get minimum incl. NaN

# create prominence table
set print $Prominence
    do for [n0=0:|$Peaks|-1]  {
        xn = yn = xp = yp = NaN
        stats $Peaks u ($0<n0 && $2>y0 ? (xn=$1,yn=$2) : 0, \
                        $0>n0 && $2>y0 && yp!=yp ? (xp=$1,yp=$2) : 0, \
                        $0==n0+1 ? (x1=$1, y1=$2) : 0 ) nooutput
        print sprintf("%g %g %g",x0,y0, min(x0-xn,xp-x0))
        x0=x1
        y0=y1
    }
set print

set key noautotitle
set offsets 0,0,1,0
set grid y

# filter Peaks
stats $Prominence u 3 name "P" nooutput      # get min, max
Filter(col,t) = (!valid(col) ? 0 : (1-(column(col)-P_min)/(P_max-P_min))*99+1)<=t ? $2 : NaN

Threshold = 100
set label 1 at graph 0.02,0.95 sprintf("Threshold: %g",Threshold)

plot FILE u colX:colY w l lc rgb "blue" ti "Spectrum", \
     $Prominence u 1:(Filter(3,Threshold)) w impulses lc rgb "red", \
     $Prominence u 1:(Filter(3,Threshold)):1 w labels offset 0,1
### end of script
