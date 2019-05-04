set output 'plot.png'
set xlabel "x"
set ylabel "y"
set grid
m1 = "./spline.txt"
m2 = "./pont.txt" 
set title 'The Spline' 
plot m1 using 1:2 w l, m2 using 1:2 w p 
