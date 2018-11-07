MASS
pracma
GPFDA
RandomFields
caret
rstan
Rcpp
ggplot2
viridis> package_name::method_name()> library(package_name)
> method_name()
library(ggplot2)
library(Rcpp) 
source("r_code/auxiliary_functions.R")
source("r_code/jupyter_stuff.R") # Only needed if using jupyter

options(repr.plot.width=5, repr.plot.height=4)

# Generate some toy data
reg_data <- lm_data(seed=0)
reg_data$type <- 'observed'

# Fit linear model
model1 <- lm(y ~ x, data=reg_data)

# Predictions
x_linspace <- seq(.5, 2.5 * pi, length.out=100) 
to_append <- data.frame(x=x_linspace, 
                        y=predict(model1, newdata=data.frame(x=x_linspace)),
                        type='order_1')
reg_data <- rbind(reg_data, to_append)

# Plot results
plt <- ggplot(data=subset(reg_data, type=='observed'), aes(x, y)) + 
        geom_point(col="steelblue", size=2) + 
        geom_line(data=subset(reg_data, type=='order_1'), col="red") +
        ggtitle(label = "Linear Regression Example") + 
        theme_minimal()
publish_gg(plt)

n_poly <- 4
poly_data <- reg_data
for (i in 2:n_poly) {
    new_dataset <- subset(reg_data, type=='observed')[,-c(3)]
    new_linspace <- data.frame(x=x_linspace)
    for (j in 2:i) {
        new_dataset[paste0("x", j)] <- new_dataset$x^j
        new_linspace[paste0("x", j)] <- new_linspace$x^j
    }
    m <- lm(y ~ ., data = new_dataset)
    to_append <- data.frame(x=x_linspace, 
                            y=predict(m, newdata=new_linspace), 
                            type=paste("order", i, sep="_"))
    poly_data <- rbind(poly_data, to_append)
}

# Plot results
plt <- ggplot(data=subset(poly_data, type=='observed'), aes(x, y)) + 
        geom_point(col="steelblue", size=2) + 
        geom_line(data=subset(poly_data, type!='observed'), aes(color=type)) +
        ggtitle(label = "Polynomial Regression Example") + 
        theme_minimal()
publish_gg(plt)

spatial_reg <- soil_data(n_peaks=2, n_data = 150, seed=0)
head(spatial_reg)

gg <- ggplot(spatial_reg, aes(lng, lat)) + 
        geom_point(aes(col=soiliness), size=2.5) +
        viridis::scale_color_viridis(option="plasma") +
        ggtitle("Clustered Spatial Data") +
        theme_minimal()
publish_gg(gg)

# Fit a model with a trend of order 1
surf_order1 <- lm(soiliness ~ lng * lat, data=spatial_reg)

# Fit a model with a trend of order 2
surf_order2 <- lm(soiliness ~ lng * lat + I(lng^2) + I(lat^2), data=spatial_reg)

# Grid for predictions
surf_grid <- as.data.frame(make_grid(size = 20))

surf_grid$surface_o1 <- predict(surf_order1, newdata=surf_grid)
surf_grid$surface_o2 <- predict(surf_order2, newdata=surf_grid)

plt <- ggplot(surf_grid, aes(lng, lat)) + 
        geom_raster(aes(fill=surface_o1)) +
        geom_contour(aes(z=surface_o1), col="white", linetype=1, alpha=.5) +
        geom_point(data=spatial_reg, aes(fill=soiliness), colour="white", pch=21, size=3) +
        viridis::scale_fill_viridis(option="plasma", na.value="darkblue") +
        ggtitle("Trend Surface: Order 1") +
        theme_minimal()
publish_gg(plt)

plt <- ggplot(surf_grid, aes(lng, lat)) + 
        geom_raster(aes(fill=surface_o2)) +
        geom_contour(aes(z=surface_o2), col="white", linetype=1, alpha=.5) +
        geom_point(data=spatial_reg, aes(fill=soiliness), colour="white", 
                   pch=21, size=3) +
        ggtitle("Trend Surface Order 2") +
        viridis::scale_fill_viridis(option="plasma", na.value="darkblue") +
        theme_minimal()
publish_gg(plt)

x_0 <- as.matrix(0, nrow=1, ncol=1)
x_linspace <- as.matrix(seq(-8, 8, length.out = 100), ncol=1)
kval <- c(as.numeric(rbf(x_0, x_linspace, sigma=1, lengthscale=1)),
       as.numeric(rqd(x_0, x_linspace, sigma=1, lengthscale=1, alpha=1)),
       as.numeric(pow(x_0, x_linspace, sigma=1, lengthscale=1, gamma=1)))
kname <- c(rep("exp_quad", 100), rep("rat_quad", 100), rep("pow_exp", 100))
covdf <- data.frame(x=rep(x_linspace, 3), k=kval, covariance= kname)


plt <- ggplot(covdf, aes(x, k)) + geom_line(aes(col=covariance)) +
        ggtitle("Different Covariance Functions") +
        theme_minimal()
publish_gg(plt)

xx <- random_points_ordered(n_points = 10, seed=0)

plt <- ggplot(as.data.frame(xx), aes(lng, lat)) + 
        geom_point(col="steelblue", size=2.5) + 
        theme_minimal()
publish_gg(plt)

K <- rqd(xx, xx, lengthscale=20, alpha=2)
kdf <- data.frame(k = as.vector(K), 
                  point_i = factor(rep(10:1, 10)), 
                  point_j = factor(sort(rep(1:10, 10), decreasing=TRUE)))

plt <- ggplot(kdf, aes(point_j, point_i)) +
        geom_raster(aes(fill=k)) +
        viridis::scale_fill_viridis(option="plasma") +
        ggtitle("Covariance Matrix") +
        scale_y_discrete(limits=rev(levels(kdf$point_i))) +
        theme_minimal()
publish_gg(plt)

# For the moment we will assume that the following quantities are known
rbf_sigma <- 20
rbf_lengthscale <- 25
epsilon_sigma <- 10

# These are the locations of the train set and predictions
X <- as.matrix(spatial_reg[, c("lng", "lat")])
X_star <- as.matrix(surf_grid[, c("lng", "lat")])

# Output value
Y <- as.matrix(spatial_reg$soiliness)

# Fit a Gaussian process
gp_fit <- gp_regression_rbf(X, Y, X_star, k_var=rbf_sigma^2, k_len=rbf_lengthscale, noise_var=epsilon_sigma^2)
surf_grid$surface_gp_mean <- gp_fit$gp_mean
surf_grid$surface_gp_var <- gp_fit$gp_var

plt <- ggplot(surf_grid, aes(lng, lat)) + 
        geom_raster(aes(fill=surface_gp_mean)) +
        geom_contour(aes(z=surface_gp_mean), col="white", linetype=1, alpha=.5) +
        #geom_point(data=spatial_reg, aes(fill=soiliness), colour="white", pch=21, size=3) +
        viridis::scale_fill_viridis(option="plasma", na.value="darkblue") +
        ggtitle("Gaussian Process - Predictive Mean") +
        theme_minimal()
publish_gg(plt)

plt <- ggplot(surf_grid, aes(lng, lat)) + 
      geom_raster(aes(fill=surface_gp_var)) +
      geom_contour(aes(z=surface_gp_var), col="white", linetype=1, alpha=.5) +
      geom_point(data=spatial_reg, colour="white", pch=21, size=3) +
      viridis::scale_fill_viridis(option="plasma", na.value="darkblue") +
      ggtitle("Gaussian Process - Predictive Variance") +
      theme_minimal()
publish_gg(plt)

x_0 <- as.matrix(0, nrow=1, ncol=1)
x_linspace <- as.matrix(seq(-8, 8, length.out = 100), ncol=1)

exp_quad <- c()
for (ell in c(.5, 1, 1.5, 2, 4)) {
    k <- rbf(x_0, x_linspace, sigma=1, lengthscale=ell)
    exp_quad <- rbind(exp_quad,
                      data.frame(x=x_linspace, 
                                 cov_function=as.numeric(k),
                                 lengthscale=paste0('lengthscale=',ell)))
}

plt <- ggplot(exp_quad, aes(x, cov_function)) + 
        geom_line(aes(color=lengthscale)) +
        ggtitle("Exponentiated Quadratic Covariance Function (sigma=1)") +
        theme_minimal()
publish_gg(plt)

x_0 <- as.matrix(0, nrow=1, ncol=1)
x_linspace <- as.matrix(seq(-8, 8, length.out = 100), ncol=1)

rat_quad <- c()
for (ell in c(.5, 1, 1.5, 2, 4)) {
    k <- rqd(x_0, x_linspace, sigma=1, lengthscale=ell, alpha=.1)
    rat_quad <- rbind(rat_quad,
                      data.frame(x=x_linspace, 
                                 cov_function=as.numeric(k),
                                 lengthscale=paste0('lengthscale=',ell)))
}

plt <- ggplot(rat_quad, aes(x, cov_function)) + 
        geom_line(aes(color=lengthscale)) +
        ggtitle("Rational Quadratic Covariance Function (sigma=1, alpha=.1)") +
        theme_minimal()
publish_gg(plt)

fit_hyper <- rstan::stan(file="stan_code/gp_gaussian_fit_hyper.stan",
                         data=list(n_data = nrow(spatial_reg),
                                   n_dim = 2,
                                   y_data = spatial_reg$soiliness,
                                   X_data = spatial_reg[, c("lng", "lat")],
                                   mu_data = rep(0, nrow(spatial_reg))),
                         pars = c("noise_var", "cov_var", "cov_length"),
                         warmup = 1000, iter = 5000, chains = 4,
                         verbose = FALSE, seed=123)

hyper_params <- rstan::extract(fit_hyper)
noise_var <- mean(cbind(hyper_params$noise_var))
cov_var <- mean(cbind(hyper_params$cov_var))
cov_length <- colMeans(cbind(hyper_params$cov_length))

# Fit a Gaussian process
gp_fit_stan <- gp_regression_rbf(X, Y, X_star, k_var=cov_var, k_len=mean(cov_length), noise_var=noise_var)
surf_grid$surface_gp_mean_stan <- gp_fit_stan$gp_mean
surf_grid$surface_gp_var_stan <- gp_fit_stan$gp_var

# Plot predictive mean
plt <- ggplot(surf_grid, aes(lng, lat)) + 
      geom_raster(aes(fill=surface_gp_mean_stan)) +
      geom_contour(aes(z=surface_gp_mean_stan), col="white", linetype=1, alpha=.5) +
      geom_point(data=spatial_reg, aes(fill=soiliness), colour="white", pch=21, size=3) +
      viridis::scale_fill_viridis(option="plasma", na.value="darkblue") +
      ggtitle("Gaussian Process - Predictive Mean") +
      theme_minimal()
publish_gg(plt)

# Fit a Gaussian process
gp_fit_stan <- gp_regression_rbf(X, Y, X_star, k_var=cov_var, k_len=mean(cov_length), noise_var=noise_var)
surf_grid$surface_gp_mean_stan <- gp_fit_stan$gp_mean
surf_grid$surface_gp_var_stan <- gp_fit_stan$gp_var

# Plot predictive mean
plt <- ggplot(surf_grid, aes(lng, lat)) + 
      geom_raster(aes(fill=surface_gp_mean_stan)) +
      geom_contour(aes(z=surface_gp_mean_stan), col="white", linetype=1, alpha=.5) +
      geom_point(data=spatial_reg, aes(fill=soiliness), colour="white", pch=21, size=3) +
      viridis::scale_fill_viridis(option="plasma", na.value="darkblue") +
      ggtitle("Gaussian Process - Predictive Mean") +
      theme_minimal()
publish_gg(plt)

# K-fold validation diagram
k = 5
label_ <- rep("Train", k^2)
for(i in 1:k){
    label_[k*(i-1) +(k+1-i)] <- "Validation"
}
cvplot <- data.frame(x=rep(1:k, k), y=sort(rep(1:k, k), decreasing = TRUE), label=label_)

plt <- ggplot(cvplot, aes(x, y)) + geom_tile(aes(fill=label)) + 
        geom_vline(xintercept = 1:k + .5, col='white') +
        geom_hline(yintercept = 1:k + .5, col="white", size=10) +
        scale_y_continuous(breaks = 1:k, labels = paste("fold", 1:k)) +
        scale_x_continuous(breaks = 1:k, labels = paste("subset", 1:k)) +
        ggtitle(paste0(k, "-Fold Validation Diagram")) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank()) +
        ylab('') + xlab('')
publish_gg(plt)

# RBF kernel
noise_var_1 <- 6
rbf_var <- 18
rbf_length <- 40

# Rational quadratic kernel
noise_var_2 <- 6
rqd_var <- 18
rqd_length <- 40
rqd_alpha <-1

surf_grid$model_1 <- gp_regression_rbf(X, Y, X_star, k_var=rbf_var, k_len=rbf_length, 
                                       noise_var=noise_var_1)$gp_mean
surf_grid$model_2 <- gp_regression_rqd(X, Y, X_star, k_var=rqd_var, k_len=rqd_length, 
                                       noise_var=noise_var_2, alpha=rqd_alpha)$gp_mean

# Change dimensions to allow side by side display
options(repr.plot.width=10, repr.plot.height=4)

# Plot predictive mean
plt1 <- ggplot(surf_grid, aes(lng, lat)) + 
  geom_raster(aes(fill=model_1)) +
  geom_contour(aes(z=model_1), col="white", linetype=1, alpha=.5) +
  viridis::scale_fill_viridis(option="plasma", na.value="darkblue", lim=c(5,55)) +
  ggtitle("Exponentiated Quadratic Kernel") +
  theme_minimal()

plt2 <- ggplot(surf_grid, aes(lng, lat)) + 
  geom_raster(aes(fill=model_2)) +
  geom_contour(aes(z=model_2), col="white", linetype=1, alpha=.5) +
  viridis::scale_fill_viridis(option="plasma", na.value="darkblue", lim=c(5,55)) +
  ggtitle("Rational Quadratic Kernel") +
  theme_minimal()

plt <- cowplot::plot_grid(plt1, plt2)
publish_gg(plt)

# Make an index with the k folds
ix = caret::createFolds(Y, k = 5)
sse1 = 0
sse2 = 0
for (i in 1:5){
    # Split data
    Xtrain <- X[-ix[[i]], ] 
    Ytrain <- Y[-ix[[i]]]
    Xval <- X[ix[[i]], ]
    Yval <- Y[ix[[i]]]
    # Train models
    m1 <- gp_regression_rbf(X=Xtrain, Y=Ytrain, X_pred=Xval, k_var=rbf_var, k_len=rbf_length, 
                            noise_var=noise_var_1)
    m2 <- gp_regression_rqd(X=Xtrain, Y=Ytrain, X_pred=Xval, k_var=rqd_var, k_len=rqd_length, 
                            noise_var=noise_var_1, alpha=rqd_alpha)
    # Predict hold-out data
    Y1_hat <- m1$gp_mean
    Y2_hat <- m2$gp_mean
    # Compute SSE
    sse1 <- sse1 + sum((Yval - Y1_hat)^2)
    sse2 <- sse2 + sum((Yval - Y2_hat)^2)
}
print(c(sse1, sse2))

K1 <- rbf(X, X, sigma=sqrt(rbf_var), lengthscale=rbf_length) + diag(x=noise_var_1, nrow=nrow(X))
K2 <- rqd(X, X, sigma=sqrt(rqd_var), lengthscale=rqd_length, alpha=rqd_alpha) + diag(x=noise_var_2, nrow=nrow(X))

W1 <- chol2inv(chol(K1))
W2 <- chol2inv(chol(K2))

s1 <- 1/diag(W1)
s2 <- 1/diag(W2)
mu1 <- Y - (W1 %*% Y) / diag(W1)
mu2 <- Y - (W2 %*% Y) / diag(W2)

log_pi1 <- .5 * sum(-log(s1) - (Y - mu1)^2/s1 - log(2*pi))
log_pi2 <- .5 * sum(-log(s2) - (Y - mu2)^2/s2 - log(2*pi))

print(c(log_pi1, log_pi2))

data(meuse, package="sp")
head(meuse)

options(repr.plot.width=5, repr.plot.height=4)
plt <- ggplot(meuse, aes(x, y)) + 
          geom_point(aes(fill=copper), colour="white", pch=21, size=3) +
          viridis::scale_fill_viridis(option="plasma", na.value="darkblue") +
          ggtitle("Meuse Data") +
          theme_minimal()
publish_gg(plt)
