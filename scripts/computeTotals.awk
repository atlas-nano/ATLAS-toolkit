{
  x1 += $2
  x2 += $2*$2
  if ($2 > xmax) xmax = $2
  if ($2 < xmin) xmin = $2
}
END {
  x1 = x1/NR
  x2 = x2/NR
  sigma = 0
  if((x2 - x1*x1) > 0) {
    sigma = sqrt(x2 - x1*x1)
  }
  if (NR > 1) std_err = sigma/sqrt(NR - 1)
  print "Number of points = " NR
  print "Min = " xmin
  print "Max = " xmax
  print "Mean = " x1 " Standard Deviation = " sigma
  print "Standard Deviation = " sigma
  print "Standard Error = " std_err
  print "Total = " (x1 * NR)
}

