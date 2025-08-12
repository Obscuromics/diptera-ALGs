library(rphylopic)

diptera <- get_uuid(name = "diptera", filter = 'by', n = 100)

# imgs <- sapply(diptera, get_phylopic)
fruitfly_uuid <- "be776840-35e4-45a8-87be-b0a39d22bd56"

vector_fruitfly <- get_phylopic(fruitfly_uuid, format = 'vector')
original_fruitfly <- get_phylopic(fruitfly_uuid, format = 'source')
raster_fruitfly <- get_phylopic(fruitfly_uuid, format = 'raster')
raster2_fruitfly <- get_phylopic(fruitfly_uuid, format = 'raster', height = 1295)
# original_fruitfly <- get_phylopic(diptera[44], format = 'original')

# OK, now we've got the image we want... let's add it to a plot!
pdf('figures/figure1_silhouettes_test.pdf')

plot(NULL, xlim = c(0,3), ylim = c(0, 3), type = "n", ann = FALSE, axes = F)

# for (i in 1:length(imgs)){
#     ypos <- i %% 10
#     xpos <- floor(i / 10)
#     add_phylopic_base(img = imgs[[i]], x = xpos, y = ypos, height = 0.5)
#     text(xpos, ypos - 0.4, i)
# }

# https://images.phylopic.org/images/be776840-35e4-45a8-87be-b0a39d22bd56/vector.svg
# https://images.phylopic.org/images/be776840-35e4-45a8-87be-b0a39d22bd56/source.svg

add_phylopic_base(img = raster_fruitfly, x = 2, y = 1, height = 0.5)
add_phylopic_base(img = vector_fruitfly, x = 1, y = 1, height = 0.5)
add_phylopic_base(img = raster2_fruitfly, x = 2, y = 2, height = 0.5)
add_phylopic_base(img = original_fruitfly, x = 1, y = 2, height = 0.5)

dev.off()

# uuid <- diptera[44]
# fruit_links <- phy_GET(file.path("images", uuid))$`_links`
# url <- fruit_links$sourceFile$href
# img <- make_png(url, height)


plot_silhouettes <- TRUE
if ( plot_silhouettes ){
  require('rphylopic')
}

diptera <- get_uuid(name = "diptera", filter = 'by', n = 100)
imgs <- sapply(diptera, get_phylopic)

if ( plot_silhouettes ) {
  for 
}