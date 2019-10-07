#==========================================================================================================
#**********************************************************************************************************
# N-of-1 t-tests: Figures 1 and 2
# Run time: less than 1 minute
#           on a Dell Latitude 5590 with an 
#           Intel Core i7-8650U CPU @ 1.90 GHz 
#           with 32.0 GB of RAM
#           R version 3.5.3
#**********************************************************************************************************
#==========================================================================================================


#==========================================================================================================
#------------------------------------------------- Set-up -------------------------------------------------
#==========================================================================================================
# Change to correct folder   
setwd("<insert location where figures will be outputted>")

# Load library
library(MASS)

# Serially correlated errors generated with rho=.5
# epsilon1 and epsilon3 have m=5, epsilon2 has m=7
epsilon1 <- c(-1.1535625, -1.2990488, -0.4307985, -0.6368679, -0.3788991)
epsilon2 <- c(-0.5491445, -0.3542094, -0.6961421, 1.3830632, 0.1664940, -0.5992540, -0.3068639)
epsilon3 <- c(0.59285090, -0.17726845, 0.82539252, 0.37960591, -0.09169061)


#==========================================================================================================
#------------------------------------------------ Figures -------------------------------------------------
#==========================================================================================================

#------------------------------------------------ Figure 1 ------------------------------------------------
tiff("N-of-1 t-tests Figure 1.tiff", bg="white", width=900, height=600, compression="lzw")
layout(matrix(c(1,3,2,4),2,2))
resize <- 1
#       btm, lt, top, rt
par(mar=c(4.1,5,2.5,1.6), xpd=NA, cex=resize)

# (a) Paired level-change A ----
x <- 1:5
y0 <- c(1+(0:4/4)); y1 <- c(3+(0:4/4))
rand.ord <- c(1, -1, 1, 1, -1)*0.25
.x <- c(x+rand.ord, x-rand.ord); y <- c(y0,y1)
ydat <- y+c(epsilon1,epsilon3)/3
shapes <- c( rep(19,5), rep(1,5))
plot(.x, ydat, pch=shapes, xlim=c(0.5, 5.5), ylim=c(-1,5), axes=FALSE, 
     ylab='', xlab='', cex=1.5*resize)
segments(c(1:4)+.5, -1, c(1:4)+.5, 5, col='gray80')
lines(x,y0, lty=1, lwd=resize*2)
lines(x,y1, lty=3, lwd=resize*2)
text(3,y0[3]-.25,'A',cex=1.75*resize,pos=1)
text(3,y1[3],'B',cex=1.75*resize,pos=1)
mtext('y', side=2, line=1.5, cex=1.5*resize)
mtext('Block', side=1, line=2.5, cex=1.5*resize)
box()
axis(2, at=c(-1, 1, 3, 5), label=c(' ',' ',' ',' '))
axis(1, at=x, cex.axis=1.5*resize)
plotdim <- par("usr")
text(plotdim[2]-(plotdim[2]-plotdim[1])*1.15,(plotdim[4]-plotdim[3])*1.1+plotdim[3], "(a)",cex=1.25*resize)

# (b) Paired level-change B ----
d <- y1 - y0
ddat <- ydat[6:10] - ydat[1:5]
plot(x, ddat, pch=17, xlim=c(0.5, 5.5), ylim=c(0, 6), axes=FALSE, 
     ylab='', xlab='', cex=1.5*resize)
axis(2, at=c(0, 2, 4, 6), label=c(' ',' ',' ',' '), las=2)
axis(1, at=x, cex.axis=1.5*resize)
mtext(expression(paste("d = ", "y"["B"], " - ", "y"["A"])), side=2, line=1.5, cex=1.5*resize)
mtext('Block', side=1, line=2.5, cex=1.5*resize)
box()
lines(x,d, lty=4, lwd=resize*2)
text(3,d[3],'B - A',cex=1.75*resize,pos=1)
plotdim <- par("usr")
text(plotdim[2]-(plotdim[2]-plotdim[1])*1.15,(plotdim[4]-plotdim[3])*1.1+plotdim[3], "(b)",cex=1.25*resize)


# (c) Paired rate-change A ----
x <- 1:5
rand.ord <- c(-1, 1, 1, -1, 1)*0.25
.x <- c(x+rand.ord, x-rand.ord); y <- c(y0,y1)
y0 <- c(0:4)/4; y1 <- 0:4 
y <- c(y0,y1)
ydat <- y+c(epsilon1,epsilon3)/3
shapes <- c( rep(19,5), rep(1,5))
plot(.x, ydat, pch=shapes, xlim=c(0.5, 5.5), ylim=c(-1,5), axes=FALSE, 
     ylab='', xlab='', cex=1.5*resize)
segments(c(1:4)+.5, -1, c(1:4)+.5, 5, col='gray80')
lines(x,y0, lty=1, lwd=resize*2)
lines(x,y1, lty=3, lwd=resize*2)
mtext('y', side=2, line=1.5, cex=1.5*resize)
mtext('Block', side=1, line=2.5, cex=1.5*resize)
axis(1, at=x, cex.axis=1.5*resize)
axis(2, at=c(-1, 1, 3, 5), label=c(' ',' ',' ',' '))
text(3,y0[3]-.25,'A',cex=1.75*resize,pos=1)
text(3,y1[3],'B',cex=1.75*resize,pos=1)
box()
plotdim <- par("usr")
text(plotdim[2]-(plotdim[2]-plotdim[1])*1.15,(plotdim[4]-plotdim[3])*1.1+plotdim[3], "(c)",cex=1.25*resize)

# (d) Paired rate-change B ----
d <- y1 - y0
ddat <- ydat[6:10] - ydat[1:5]
plot(x, ddat, pch=17, xlim=c(0.5, 5.5), ylim=c(0,4), axes=FALSE, 
     ylab='', xlab='', cex=1.5*resize)
axis(2, at=c(0, 1.33, 2.67, 4), label=c(' ',' ',' ',' '), las=2)
axis(1, at=x, cex.axis=1.5*resize)
mtext(expression(paste("d = ", "y"["B"], " - ", "y"["A"])), side=2, line=1.5, cex=1.5*resize)
mtext('Block', side=1, line=2.5, cex=1.5*resize)
lines(x ,d, lty=4, lwd=resize*2)
box()
text(3,d[3]*.8,'B - A',cex=1.75*resize,pos=1)
plotdim <- par("usr")
text(plotdim[2]-(plotdim[2]-plotdim[1])*1.15,(plotdim[4]-plotdim[3])*1.1+plotdim[3], "(d)",cex=1.25*resize)
dev.off()

#------------------------------------------------ Figure 2 ------------------------------------------------
tiff("N-of-1 t-tests Figure 2.tiff", bg="white", width=900, height=300, compression="lzw")
layout(matrix(c(1,2),1,2))
resize <- 1
#       btm, lt, top, rt
par(mar=c(4.1,5,2.5,1.6), xpd=NA, cex=resize)
#par(mar=c(4.1,5,2.5,1.6), xpd=NA, cex=resize)

# (a) 2-sample level-change ----
x0 <- 1:5; x1 <- 7:13
y0 <- rep(0,5); y1 <- rep(1, 7)
x <- c(x0,x1); y <- c(y0,y1)
ydat <- y+c(epsilon1,epsilon2)/3
shapes <- c( rep(19,5), rep(1,7))
plot(x, ydat, pch=shapes, ylim=c(-1,2), axes=FALSE, xlab='', ylab='', cex=1.5*resize)
lines(x0,y0, lty=1, lwd=resize*2)
lines(x1,y1, lty=3, lwd=resize*2)
axis(1, at=c(1,3,5,7,9,11,13), label=c('1','3','5', '1','3','5','7'), cex.axis=1.5*resize)
axis(1, at=c(2,4,8,10,12), label=c('2','4', '2','4','6'), cex.axis=1.5*resize)
text(3,-.2,'A',cex=1.75*resize,pos=1)
text(10,1,'B',cex=1.75*resize,pos=1)
mtext('y', side=2, line=1.5, cex=1.5*resize)
mtext('Observation', side=1, line=2.5, cex=1.5*resize)
box()
axis(2, at=c(-1,0,1,2), label=c(' ',' ',' ',' '))
plotdim <- par("usr")
text(plotdim[2]-(plotdim[2]-plotdim[1])*1.15,(plotdim[4]-plotdim[3])*1.1+plotdim[3], "(a)",cex=1.25*resize)

# (b) 2-sample rate-change ----
x0 <- 1:5; x1 <- 7:13
y0 <- c(0:4)/4; y1 <- c(0:6)/1.5 
y <- c(y0,y1)
ydat <- y+c(epsilon1,epsilon2)/3
shapes <- c( rep(19,5), rep(1,7))
x <- c(x0,x1)
plot(x, ydat, pch=shapes, ylim=c(-1,5), axes=FALSE, ylab='', xlab='', cex=1.5*resize)
lines(x0,y0, lty=1, lwd=resize*2)
lines(x1,y1, lty=3, lwd=resize*2)
text(3,y0[3]*.5,'A',cex=1.75*resize,pos=1)
text(10,y1[4]*.9,'B',cex=1.75*resize,pos=1)
mtext('y', side=2, line=1.5, cex=1.5*resize)
mtext('Observation', side=1, line=2.5, cex=1.5*resize)
box()
axis(2, at=c(-1, 1, 3, 5), label=c(' ',' ',' ',' '))
axis(1, at=c(1,3,5,7,9,11,13), label=c('1','3','5', '1','3','5','7'), cex.axis=1.5*resize)
axis(1, at=c(2,4,8,10,12), label=c('2','4', '2','4','6'), cex.axis=1.5*resize)
plotdim <- par("usr")
text(plotdim[2]-(plotdim[2]-plotdim[1])*1.15,(plotdim[4]-plotdim[3])*1.1+plotdim[3], "(b)",cex=1.25*resize)
dev.off()





