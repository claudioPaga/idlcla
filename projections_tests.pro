pro projections_tests, lon0, lat0, x, y

  lon0 = lon0
  lat0 = lat0

  lon0=lon0*!PI/180.0
  lat0=lat0*!PI/180.0

  x=x*!PI/180.0
  y=y*!PI/180.0

  d = atan(sqrt(x^2+y^2))
  b = atan(-x, y)
  xx = sin(lat0) * sin(d) * cos(b) + cos(lat0) * cos(d)
  yy = sin(d) * sin(b)

  lon1 = lon0 + atan(yy, xx)
  lat1 = asin(sin(lat0) * cos(d) - cos(lat0) * sin(d) * cos(b))

  print, lon1, lat1
  print, lon1*180./!PI,lat1*180./!PI
  stop

  
  c = sqrt(x^2+y^2)
  
  lat1 = 1./sin(cos(c)*sin(lat0)+ y * sin(c) * cos(lat0) * 1./c)

  lon1 = lon0 + 1./tan(x*sin(c)/(c*cos(lat0)*cos(c)-y*sin(lat0)*sin(c)))
  

  stop
end
