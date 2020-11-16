# Jiaying's non-uniqueness of phat example
par(pty = "s")
z = c(-0.5, -1, 4, -3, 0.5)
v = c(-0.25, 0, 0, 0.4, 0.5)
x = c(-1,1)
plot(x,x, type = "n", xlab = "", ylab = "", axes = FALSE, bty = "n" )
verts = function(K, A, b){
    nv = nrow(K)
    B = matrix(0,nv,2)
    for(i in 1:nv)
	B[i,] = solve(A[K[i,],],b[K[i,]])
    B
}
M = 1000
A = cbind(1, -c(z,M,M,0,0))
b = c(v,2*M,-2*M,-1.2,1.2)
K = list(
    K1 = cbind(c(1,5,6),c(5,6,1)),
    K2 = cbind(c(1,2,7),c(2,7,1)),
    K3 = cbind(c(3,1,4,8),c(1,4,8,3)),
    K4 = cbind(c(4,5,2,6),c(5,2,6,4)),
    K5 = cbind(c(3,5,9),c(5,9,3)),
    K6 = cbind(c(2,3,4),c(3,4,2)))
for(i in 1:6){
    B = verts(K[[i]], A, b)
    polygon(B[,2],B[,1], col = "skyblue")
}
# Line Labels [paste/expression mysteries aargh!]
Labels = list( expression(L[1]),
	      expression(L[2]),
	      expression(L[3]),
	      expression(L[4]),
	      expression(L[5]))
# determined by lab = locator(5)
loc.Labels = list(x =  c(-1.004424, -1.01302343,  0.17934171, 
			 -0.08649861,  1.00436081),
		  y = c(0.3205637,  0.91636306,  1.00753127,
			1.02314035,  0.89690312))
for(i in 1:5)
    text(loc.Labels[[1]][i], loc.Labels[[2]][i], labels = Labels[[i]])

# Polygon Labels [paste/expression mysteries aargh!]
Labels = list( expression(P[1]),
	      expression(P[2]),
	      expression(P[3]),
	      expression(P[4]),
	      expression(P[6]),
	      expression(P[5]))
# determined by lab = locator(5)
loc.Labels = list(x = c(-1.00289155,  1.02002596,  0.06291857, 
			-0.53817104,  0.54508011,  0.05471212),
		  y = c(0.11818783, -0.88815825, -0.92915311,  
			0.94752513,  0.98451957,  0.07552926))
loc.counts = list(x = c(-0.88573936,  0.82733333,  0.11516691, 
			-0.40868114,  0.09882879,  0.34031147),
		  y = c(0.12749720, -0.73885324, -0.64723638,  
			0.74507448, -0.03584664,  0.81743703))
for(i in 1:6){
    text(loc.Labels[[1]][i], loc.Labels[[2]][i], labels = Labels[[i]])
    text(loc.counts[[1]][i], loc.counts[[2]][i], labels = "3", col = 2)
}
loc.counts = list(x = c(0.2923439, -0.2889679, 0.06056319),
		  y = c(-0.36331944,  0.04096497, 0.4360896))
for(i in 1:3)
    text(loc.counts[[1]][i], loc.counts[[2]][i], labels = "1", col = 2)
loc.counts = list(x = c(0.5650635, -0.5907933,  0.0510295, 
			-0.8496609, 0.5158996,  0.0536279, -0.09240575),
		  y = c(-0.7865217, -0.4373644, -0.2191364,  
			0.4922822,  0.1549019,  0.8226151, 0.2862798))
for(i in 1:7)
    text(loc.counts[[1]][i], loc.counts[[2]][i], labels = "2", col = 2)

# Arrows
x = c(-0.95, -0.6, 0.3, -0.2, 0.4)  
d = c(-0.07, 0.08, -.3, -0.3, 0.07)
y = v + z*x + d
for(i in 1:5) {
    abline(c(v[i], z[i]), col = 1)
    x0 = x[i]
    y0 = y[i]
    a = v[i]
    b = z[i]
    xhat = (y0*b - a*b + x0)/(1 + b^2)
    yhat = a + b*xhat
    arrows(xhat, yhat, x0, y0, length = .03, col = 1, lwd = 1.5)
}
