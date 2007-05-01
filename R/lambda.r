lambda <- function(z, x, xi, h = NULL, kernel = c("Gaussian",
"Beta", "Triangle", "Cosinus", "Optcosinus"), g = 0){

    if (is.null(h)){cat("WARNING! BANDWIDTH MUST BE SPECIFIED!","\n")}
    kernel <- match.arg(kernel)
    inwindow <- (abs((xi-x)/h)<=1)
#   inwindow <- abs((xi-x)<=h)

    if (kernel == "Gaussian") {
        sum(z*kern.G(x,xi,h)*inwindow)/sum(kern.G(x,xi,h)*inwindow)}

    else if (kernel == "Beta") {
        sum(z*kern.B(x,xi,h,g)*inwindow)/sum(kern.B(x,xi,h,g)*inwindow)}

    else if (kernel == "Triangle") {
        sum(z*kern.T(x,xi,h)*inwindow)/sum(kern.T(x,xi,h)*inwindow)}

    else if (kernel == "Cosinus") {
        sum(z*kern.C(x,xi,h)*inwindow)/sum(kern.C(x,xi,h)*inwindow)}

    else if (kernel == "Optcosinus") {
        sum(z*kern.O(x,xi,h)*inwindow)/sum(kern.O(x,xi,h)*inwindow)}

}
