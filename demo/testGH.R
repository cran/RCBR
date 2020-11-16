# Groeneboom and Hendrikx profiling for Horowitz data
data(Horowitz93)
H <- Horowitz93
H <- data.frame(x = H$DOVTT, v = H$DCOST/100, y = H$DEPEND, k = H$CARS)
D <- H[H$k == 1, ]
f <- glm(y ~ x, offset = v,  data = D, family = binomial(link = "probit"))
g0 <- prcbr(y ~ v | x, data = D, b0 = 0)
g1 <- prcbr(y ~ v | x, data = D, b0 = -20:20/100)
h0 <- prcbr(y ~ v | x, data = D, b0 = .03, logL = FALSE, omethod = "L-BFGS-B", lo = 0.005, up = 0.1)
h1 <- prcbr(y ~ v | x, data = D, b0 = (pi/3 + 0:40)/400, logL = FALSE)
plot((pi/3 + 0:40)/400, h1$bopt, type = "l")
abline(v = h0$bopt$par)

