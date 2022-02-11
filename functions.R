chroma <- function(R, G, B) {
  require(colorscience)
  
  Xi <- 2.7689*R + 1.7517*G + 1.1302*B
  Yi <- 1.0000*R + 4.5907*G + 0.0601*B
  Zi <- 0.0565*G + 5.5943*B
  
  x <-  Xi / (Xi + Yi +  Zi)
  y <-  Yi / (Xi + Yi +  Zi)
  z <-  Zi / (Xi + Yi +  Zi)
  
  alpha <- atan2( (x - 0.33), (y - 0.33)) * 180/pi
  
  cie <- colorscience::cccie31 %>%
    mutate(a = atan2( (x - 0.33), (y - 0.33)) * 180/pi) %>%
    dplyr::filter(wlnm <= 700) %>%
    dplyr::filter(wlnm >= 380) %>%
    dplyr::select(a, wlnm) %>%
    arrange(a)
  
  wl <- cie[as.vector(sapply(alpha,function(x) which.min(abs(x - cie$a)))), 'wlnm']
  
  return(wl)
}
chroma2 <- function(R, G, B) {
  # same as chroma but using interpolation instead of nearest neighbor to assign dw from the lookup table
  require(colorscience)
  
  Xi <- 2.7689*R + 1.7517*G + 1.1302*B
  Yi <- 1.0000*R + 4.5907*G + 0.0601*B
  Zi <- 0.0565*G + 5.5943*B
  
  x <-  Xi / (Xi + Yi +  Zi)
  y <-  Yi / (Xi + Yi +  Zi)
  z <-  Zi / (Xi + Yi +  Zi)
  
  alpha <- atan2((x - 0.33), (y - 0.33)) * 180/pi
  
  cie <- colorscience::cccie31 %>%
    mutate(a = atan2( (x - 0.33), (y - 0.33)) * 180/pi) %>%
    dplyr::filter(wlnm <= 700) %>%
    dplyr::filter(wlnm >= 380) %>%
    dplyr::select(a, wlnm) %>%
    arrange(a)
  
  wl = approx(x = cie$a, y = cie$wlnm, xout = alpha, method = "linear")$y
  # wl <- cie[as.vector(sapply(alpha,function(x) which.min(abs(x - cie$a)))), 'wlnm']
  
  return(wl)
}
chromaLehmann = function(uB, B, G, R) {
  X = 11.053 * uB + 6.950 * B + 51.135 * G + 34.457 * R
  Y = 1.320 * uB + 21.053 * B + 66.023 * G + 18.034 * R
  Z = 58.038 * uB + 34.931 * B + 2.606 * G + 0.016 * R
  
  XYZ = X + Y + Z
  x = X / XYZ
  y = Y / XYZ
  z = Z / XYZ
  
  #  %% (2 * pi))
  alpha0 = (atan2(x - 1 / 3, y - 1/3)) * 180 / pi
  
  ## angle shift (hue angle in *** and *** have different starting location and range)
  alpha = 270 - (alpha0  + 180)
  a = alpha / 100
  # correction only for 37º to 230º
  deltaAlpha = -52.16 * a^5 + 373.81 * a^4 - 981.83 * a^3 + 1134.19 * a^2 - 533.61 * a + 76.72
  deltaS = -0.0099 * a^5 + 0.1199 * a^4 - 0.4594 * a^3 + 0.7515 * a^2 - 0.5095 * a + 0.1222
  alphaCorrected = alpha + deltaAlpha
  
  alpha = 270 - alphaCorrected - 180 ## shift back to -180 to 180 range (cw from neg y axis)
  
  oneThird = 1 / 3
  
  s = sqrt((x - oneThird)^2 + (y - oneThird)^2)
  s = s + deltaS
  
  cie <- colorscience::cccie31 %>%
    mutate(a = atan2((x - oneThird), (y - oneThird)) * 180/pi,
           d = sqrt((x - oneThird)^2 + (y - oneThird)^2)) %>%
    dplyr::filter(wlnm <= 700) %>%
    dplyr::filter(wlnm >= 380) %>%
    select(a, wlnm) %>%
    arrange(a)
  
  wl = approx(x = cie$a, y = cie$wlnm, xout = alpha, method = "linear")$y
  d = approx(x = cie$a, y = cie$d, xout = alpha, method = "linear")$y
  
  tibble(
    alphaXPosUncorrected = a * 100, 
    alphaXPos = alphaCorrected, 
    deltaAlpha, 
    alphaYNeg = alpha, 
    s, wl, x, y, sp = s / d)
}

### functions to assign FUI level and color given dw value
get_FUI = function() {
  thresholds = c(470, 475, 480, 485, 489, 495, 509, 530, 549, 559, 564, 567, 568, 569, 570, 571, 573, 575, 577, 579, 581, 600)
  levels = paste("FUI", 1:21)
  fui_colors = c(
    "#2158bc", "#316dc5", "#327cbb", "#4b80a0", "#568f96", "#6d9298", "#698c86", 
    "#759e72", "#7ba654", "#7dae38", "#94b660","#94b660", "#a5bc76", "#aab86d", 
    "#adb55f", "#a8a965", "#ae9f5c", "#b3a053", "#af8a44", "#a46905", "#9f4d04")
  
  return(list(thresholds = thresholds, color_table = tibble(levels, fui_colors)))
}
dw2FUI = function(dw) {
  FUI_table = get_FUI()
  color_rgb = tibble(dw = dw) %>% 
    mutate(fui = cut(dw, breaks = FUI_table[[1]], labels = FUI_table[[2]]$fui_colors)) %>% 
    pull(fui) %>% 
    as.character()
  
  return(color_rgb)
}