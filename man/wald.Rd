\name{wald.test}
\alias{wald.test}
\alias{print.wald.test}
\title{Wald Test for Model Coefficients}
\description{Computes a Wald \eqn{\chi^2}{chi-squared} test for 1 or more coefficients, given their variance-covariance matrix.}
\usage{
  wald.test(Sigma, b, Terms = NULL, L = NULL, H0 = NULL,  
            df = NULL, verbose = FALSE)
  \method{print}{wald.test}(x, digits = 2, ...)
  }

\arguments{
  \item{Sigma}{
    A var-cov matrix, usually extracted from one of the fitting functions (e.g., \code{lm}, \code{glm}, ...).
    }
  \item{b}{
    A vector of coefficients with var-cov matrix \code{Sigma}. These coefficients are usually extracted from 
    one of the fitting functions available in \R (e.g., \code{lm}, \code{glm},...).
    }
  \item{Terms}{
    An optional integer vector specifying which coefficients should be \emph{jointly} tested, using a Wald 
    \eqn{\chi^2}{chi-squared} or \eqn{F} test. Its elements correspond to the columns or rows of the var-cov 
    matrix given in \code{Sigma}. Default is \code{NULL}.
    }
  \item{L}{
    An optional matrix conformable to \code{b}, such as its product with \code{b} i.e., \code{L \%*\% b} 
    gives the linear combinations of the coefficients to be tested. Default is \code{NULL}.
    }
  \item{H0}{
    A numeric vector giving the null hypothesis for the test. It must be as long as \code{Terms} or 
    must have the same number of columns as \code{L}. Default to 0 for all the coefficients to be tested.
    }
  \item{df}{
    A numeric vector giving the degrees of freedom to be used in an \eqn{F} test, i.e. the degrees of freedom 
    of the residuals of the model from which \code{b} and \code{Sigma} were fitted. Default to NULL, for no 
    \eqn{F} test. See the section \bold{Details} for more information.
    }
  \item{verbose}{
    A logical scalar controlling the amount of output information. The default is \code{FALSE}, providing minimum output.
    }
  \item{x}{
    Object of class \dQuote{wald.test}
    }
  \item{digits}{
    Number of decimal places for displaying test results. Default to 2.
    }
  \item{...}{
    Additional arguments to \code{print}.
    }
  }

\details{
  The key assumption is that the coefficients asymptotically follow a (multivariate) normal distribution with mean = 
  model coefficients and variance = their var-cov matrix.\cr
  One (and only one) of \code{Terms} or \code{L} must be given. When \code{L} is given, it must have the same number of 
  columns as the length of \code{b}, and the same number of rows as the number of linear combinations of coefficients. 
  When \code{df} is given, the \eqn{\chi^2}{chi-squared} Wald statistic is divided by \code{m} = the number of 
  linear combinations of coefficients to be tested (i.e., \code{length(Terms)} or \code{nrow(L)}). Under the null 
  hypothesis \code{H0}, this new statistic follows an \eqn{F(m, df)} distribution.
  }

\value{
  An object of class \code{wald.test}, printed with \code{print.wald.test}.
  }

\references{
  Diggle, P.J., Liang, K.-Y., Zeger, S.L., 1994. Analysis of longitudinal data. Oxford, Clarendon Press, 253 p.\cr
  Draper, N.R., Smith, H., 1998. Applied Regression Analysis. New York, John Wiley & Sons, Inc., 706 p.
  }



\examples{
  data(DT_yatesoats)
DT <- DT_yatesoats

m3 <- mmer(fixed=Y ~ V + N + V:N-1,
           random = ~ B + B:MP,
           rcov=~units,
           data = DT)

wald.test(b = m3$Beta$Estimate, Sigma = m3$VarBeta, Terms = 2)

LL <- matrix(0,nrow=1, ncol=12)
LL[1,2] <- 1
LL[1,3] <- -1
LL

wald.test(b = m3$Beta$Estimate, Sigma = m3$VarBeta, L=LL)

  }

