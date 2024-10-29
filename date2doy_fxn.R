date2doy <- function(yr, mo, dy){
  yr <- as.character(yr)
  mo <- as.character(mo)
  dy <- as.character(dy)
  
  dt <- paste(c(yr,'-',mo,'-',dy), collapse='')
  doy <- strftime(dt, format='%j')
  return(doy)
}

