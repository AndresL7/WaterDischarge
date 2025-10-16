lim = 4000
retardo = 0.001
f = 0.01

set size square
set xrange [-5.6e-5:1e-3+5.6e-5+1e-3]
set yrange [-5.6e-5:1e-3+5.6e-5+1e-3]

set palette defined (0.0 "blue", 0.01 "green", 0.02 "yellow", 0.03 "red")  # Define a color palette
set cbrange [0:0.03] 

do for [i=0:lim]{
  titulo = sprintf("paso = %.4d - tiempo = %.5f",i,i*5e-5)  
  set title titulo
  file = sprintf("state_%.4d",i)
    plot file every ::0::1599 u 2:3:(sqrt($4**2 + $5**2)) w p ps 1 pt 7 palette not,\
	 "" every ::1600::1680 u 2:3 w p ps 1 pt 5 lc rgb "black" not,\
	 "" every ::1681::1759 u 2:3 w p ps 1 pt 5 lc rgb "black" not,\
	 "" every ::1760::1840 u 2:3 w p ps 1 pt 5 lc rgb "black" not,\
	 "" every ::1841::1929 u 2:3 w p ps 1 pt 5 lc rgb "black" not
	 #"" every ::1920::1960 u 2:3 w lp ps 1 pt 7 lc rgb "black" not

  pause retardo
}

do for [i=0:lim]{
  titulo = sprintf("paso = %.4d - tiempo = %.5f",i,i*5e-5)  
  set title titulo
  file = sprintf("state_%.4d",i)
    plot file every ::0::1599 u 2:3:($4*f):($5*f):(sqrt($4**2 + $5**2)) w vectors palette not,\
	 "" every ::1600::1680 u 2:3 w p ps 1 pt 5 lc rgb "black" not,\
	 "" every ::1681::1759 u 2:3 w p ps 1 pt 5 lc rgb "black" not,\
	 "" every ::1760::1840 u 2:3 w p ps 1 pt 5 lc rgb "black" not,\
	 "" every ::1841::1929 u 2:3 w p ps 1 pt 5 lc rgb "black" not
	 #"" every ::1920::1960 u 2:3 w lp ps 1 pt 7 lc rgb "black" not

  pause retardo
}
