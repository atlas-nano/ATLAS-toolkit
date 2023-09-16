{
  x1 += $3
  x2 += $3*$3
}
END {
  x1 = x1/NR
  x2 = x2/NR
  sigma = sqrt(x2 - x1*x1)
  if (NR > 1) std_err = sigma/sqrt(NR - 1)
  print "Number of points = " NR
  print "Mean = " x1
  print "Standard Deviation = " sigma
  print "Standard Error = " std_err
}

