# library(ggplot2)
# library(sf)
# library(hexSticker)
# sfc = st_sfc(st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,1), c(0,0)))))
# hexs <- st_make_grid(sfc, cellsize = .1, square = FALSE)
# p <- ggplot(data = hexs) +
#   geom_sf(size=.1, aes(col = sf.colors(12, categorical = TRUE))) +
#   theme_bw()
# st <- sticker(p, white_around_sticker = TRUE,
#               asp = 200,
#               s_width = 3,
#               p_x = 1,
#               p_y = 1,
#               p_family = "arial",
#               p_fontface = "bold",
#               p_color = "white",
#               s_height=3.5,
#               p_size= 11,
#               package="spqdep",
#               filename="hex.png",
#               h_fill="orange",
#               h_size = 2,
#               h_color = "orange")
# plot(st)

